#ifndef CORE_PARTICLE_CACHE_HPP
#define CORE_PARTICLE_CACHE_HPP

#include <algorithm>
#include <cassert>
#include <functional>
#include <iterator>
#include <map>
#include <memory>
#include <utility>
#include <vector>

#include <boost/container/flat_map.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/mpi/exception.hpp>
#include <boost/serialization/map.hpp>

#include "utils/NoOp.hpp"
#include "utils/mpi/gather_buffer.hpp"
#include "utils/parallel/Callback.hpp"
#include "utils/serialization/flat_map.hpp"

namespace detail {
class TakeSecond {
public:
  template <typename T, typename U> U &operator()(std::pair<T, U> &p) const {
    return p.second;
  }

  template <typename T, typename U>
  U const &operator()(std::pair<T, U> const &p) const {
    return p.second;
  }
};

/**
 * @brief Merge two ordered containers into a new one.
 *
 * This implementation has a different tradeoff than
 * flat_map::merge, here we use O(N) extra memory to
 * get O(N) time complexity, while the flat_map implementation
 * avoids extra memory usage, but will cause O(N^2) copies
 * on average.
 * Inspired by the example implementation in
 * http://en.cppreference.com/w/cpp/algorithm/merge.
 */
template <typename Container> class Merge {
public:
  Container operator()(Container const &a, Container const &b) const {
    using value_type = typename Container::value_type;
    Container ret;
    ret.reserve(a.size() + b.size());

    auto first1 = a.begin();
    auto last1 = a.end();
    auto first2 = b.begin();
    auto last2 = b.end();

    while (first1 != last1) {
      if (first2 == last2) {
        for (; first1 != last1; ++first1) {
          ret.emplace_hint(ret.end(), *first1);
        }
        break;
      }

      if (first2->first < first1->first) {
        ret.emplace_hint(ret.end(), *first2);
        ++first2;
      } else {
        ret.emplace_hint(ret.end(), *first1);
        ++first1;
      }
    }

    for (; first2 != last2; ++first2)
      ret.emplace_hint(ret.end(), *first2);

    return ret;
  }
};
}

template <typename Cells, typename UnaryOp = Utils::NoOp,
          typename Range = typename std::remove_reference<decltype(
              std::declval<Cells>().particles())>::type,
          typename Particle = typename std::iterator_traits<
              typename Range::iterator>::value_type>
class ParticleCache {
  using map_type = boost::container::flat_map<int, Particle>;

  map_type remote_parts;
  std::vector<int> bond_info;
  bool m_valid, m_valid_bonds;

  Utils::Parallel::Callback update_cb;
  Utils::Parallel::Callback update_bonds_cb;

  Cells &cells;
  UnaryOp m_op;

  void m_update_bonds() {
    std::vector<int> local_bonds;

    for (auto &p : cells.particles()) {
      std::copy(p.bl.begin(), p.bl.end(), std::back_inserter(local_bonds));
    }

    Utils::Mpi::gather_buffer(local_bonds,
                              Communication::mpiCallbacks().comm());
  }

  void m_update() {
    remote_parts.clear();
    for (auto const &p : cells.particles()) {
      typename map_type::iterator it;
      /* Add the particle to the map */
      std::tie(it, std::ignore) =
          remote_parts.emplace(std::make_pair(p.identity(), p));

      /* And run the op on it. */
      m_op(it->second);
    }

    /* Reduce data to the master by merging the flat_maps from the
     * nodes in a reduction tree. */
    boost::mpi::reduce(Communication::mpiCallbacks().comm(), remote_parts,
                       remote_parts, detail::Merge<map_type>(), 0);
  }

  void m_recv_bonds() {
    bond_info.clear();

    for (auto &p : cells.particles()) {
      std::copy(p.bl.begin(), p.bl.end(), std::back_inserter(bond_info));
    }

    Utils::Mpi::gather_buffer(bond_info, Communication::mpiCallbacks().comm());

    auto it = bond_info.begin();
    for (auto &e : remote_parts) {
      auto &p = e.second;
      p.bl.e = nullptr;
      p.bl.max = 0;
      p.bl.resize(p.bl.size());

      std::copy_n(it, p.bl.size(), p.bl.begin());
      it += p.bl.size();
    }
  }

public:
  using value_iterator =
      boost::transform_iterator<detail::TakeSecond,
                                typename map_type::const_iterator,
                                Particle const &, const Particle>;

  ParticleCache() = delete;
  ParticleCache(Cells &cells, UnaryOp &&op = UnaryOp{})
      : update_cb([this](int, int) { this->m_update(); }),
        update_bonds_cb([this](int, int) { this->m_update_bonds(); }),
        cells(cells), m_valid(false), m_valid_bonds(false),
        m_op(std::forward<UnaryOp>(op)) {}
  ParticleCache(ParticleCache const &) = delete;
  ParticleCache(ParticleCache &&) = delete;

  void clear() {
    remote_parts.clear();
    bond_info.clear();
  }

  value_iterator begin() {
    assert(Communication::mpiCallbacks().comm().rank() == 0);

    if (!m_valid)
      update();

    return value_iterator(remote_parts.begin());
  }

  value_iterator end() {
    assert(Communication::mpiCallbacks().comm().rank() == 0);

    if (!m_valid)
      update();

    return value_iterator(remote_parts.end());
  }

  bool valid() const { return m_valid; }
  bool valid_bonds() const { return m_valid_bonds; }

  void invalidate() {
    clear();
    /* Release memory */
    remote_parts.shrink_to_fit();
    bond_info.shrink_to_fit();
    /* Adjust state */
    m_valid = false;
    m_valid_bonds = false;
  }

  void update_bonds() {
    update();

    if (!m_valid_bonds) {
      update_bonds_cb.call();
      m_recv_bonds();
      m_valid_bonds = true;
    }
  }

  void update() {
    if (m_valid)
      return;

    update_cb.call();

    m_update();

    m_valid = true;
  }

  /** Number of particles in the config. */
  size_t size() {
    assert(Communication::mpiCallbacks().comm().rank() == 0);

    if (!m_valid)
      update();

    return remote_parts.size();
  }

  bool empty() {
    assert(Communication::mpiCallbacks().comm().rank() == 0);

    if (!m_valid)
      update();

    return remote_parts.empty();
  }

  /**
   * @brief Access particle by id.
   * If the particle config is not valid this will trigger
   * a global update.
   * Will throw std::out_of_range if the particle does
   * not exists.
   */
  Particle const &operator[](int id) {
    assert(Communication::mpiCallbacks().comm().rank() == 0);

    if (!m_valid)
      update();

    return remote_parts.at(id);
  }

  int min_id() const { return remote_parts.begin()->first; }
  int max_id() const { return remote_parts.rbegin()->first; }
};

#endif
