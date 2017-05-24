#ifndef CORE_PARTICLE_CACHE_HPP
#define CORE_PARTICLE_CACHE_HPP

#include <algorithm>
#include <cassert>
#include <functional>
#include <iterator>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

#include <boost/container/flat_set.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/mpi/collectives.hpp>

#include "utils/NoOp.hpp"
#include "utils/mpi/gather_buffer.hpp"
#include "utils/parallel/Callback.hpp"
#include "utils/serialization/flat_set.hpp"

namespace detail {
  /**
  * @brief Compare particles by id.
  */
class IdCompare {
public:
  template <typename Particle>
  bool operator()(Particle const &a, Particle const &b) const {
    return a.identity() < b.identity();
  }
};

/**
 * @brief Merge two ordered containers into a new one.
 *
 * This implementation has a different tradeoff than
 * flat_set::merge, here we use O(N) extra memory to
 * get O(N) time complexity, while the flat_map implementation
 * avoids extra memory usage, but will cause O(N^2) copies
 * on average.
 * Inspired by the example implementation in
 * http://en.cppreference.com/w/cpp/algorithm/merge.
 */
template <typename Container, typename Compare> class Merge {
  Compare m_comp;

public:
  Merge(Compare &&comp = Compare{}) : m_comp(comp) {}
  Container operator()(Container const &a, Container const &b) const {
    Container ret;
    ret.reserve(a.size() + b.size());

    auto first1 = a.begin();
    auto last1 = a.end();
    auto first2 = b.begin();
    auto last2 = b.end();

    while (first1 != last1) {
      if (first2 == last2) {
        /* The first range has no more elements, so we can
           just copy the rest of range 2. */
        for (; first1 != last1; ++first1) {
          ret.emplace_hint(ret.end(), *first1);
        }
        break;
      }

      if (m_comp(*first2, *first1)) {
        ret.emplace_hint(ret.end(), *first2);
        ++first2;
      } else {
        ret.emplace_hint(ret.end(), *first1);
        ++first1;
      }
    }

/* The first range has no more elements, so we can
   just copy the rest of range 2. */
    for (; first2 != last2; ++first2)
      ret.emplace_hint(ret.end(), *first2);

    return ret;
  }
};
}

/**
* @brief Particle cache on the master.
*
* This class implements cached access to all particles on the
* master node. This implementation fetches all particles to
* the master on first access. Updates of the particle data are
* triggered automatically on access. The data in the cache
* is invalidated automatically on_particle_change, and then
* updated on the next access.
* By default the particles do not have valid bond information
* on them. If bonds are needed, update_bonds() has to be called
* explicitly to update the bond cache.
*
* To update the cache particles are sorted by id on the nodes,
* and the sorted arrays a merged in a reduction tree, until the
* master node recives a complete and sorted particle array.
*
* This class can be customized by running a unary opration on
* the particles. This op is run on all the nodes. It can be used
* e.g. to fold or unfold the coordinates on the fly.
*/
template <typename Cells, typename UnaryOp = Utils::NoOp,
          typename Range = typename std::remove_reference<decltype(
              std::declval<Cells>().particles())>::type,
          typename Particle = typename std::iterator_traits<
              typename Range::iterator>::value_type>
class ParticleCache {
  using map_type = boost::container::flat_set<Particle, detail::IdCompare>;

  std::unordered_map<int, int> id_index;
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
      local_bonds.push_back(p.identity());
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
      std::tie(it, std::ignore) = remote_parts.emplace(p);

      /* And run the op on it. */
      m_op(*it);
    }

    /* Reduce data to the master by merging the flat_sets from the
     * nodes in a reduction tree. */
    boost::mpi::reduce(Communication::mpiCallbacks().comm(), remote_parts,
                       remote_parts,
                       detail::Merge<map_type, detail::IdCompare>(), 0);
  }

  void m_recv_bonds() {
    bond_info.clear();

    Utils::Mpi::gather_buffer(bond_info, Communication::mpiCallbacks().comm());

    for(auto it = bond_info.begin(); it != bond_info.end();)
    {
      auto &p = remote_parts.begin()[id_index[*it++]];
      p.bl.e = nullptr;
      p.bl.max = 0;
      p.bl.resize(p.bl.size());

      std::copy_n(it, p.bl.size(), p.bl.begin());
      it += p.bl.size();
    }
  }

  void m_update_index() {
    /* Try to avoid rehashing along the way */
    id_index.reserve(remote_parts.size()+1);

    int index = 0;
    for (auto const &p : remote_parts) {
      id_index.insert(std::make_pair(p.identity(), index++));
    }
  }

public:
  using value_iterator = typename map_type::const_iterator;

  ParticleCache() = delete;
  ParticleCache(Cells &cells, UnaryOp &&op = UnaryOp{})
      : update_cb([this](int, int) { this->m_update(); }),
        update_bonds_cb([this](int, int) { this->m_update_bonds(); }),
        cells(cells), m_valid(false), m_valid_bonds(false),
        m_op(std::forward<UnaryOp>(op)) {}
  ParticleCache(ParticleCache const &) = delete;
  ParticleCache(ParticleCache &&) = delete;
  ParticleCache operator=(ParticleCache const &) = delete;
  ParticleCache operator=(ParticleCache &&) = delete;

  void clear() {
    id_index.clear();
    remote_parts.clear();
    bond_info.clear();
  }

  /**
   * @brief Iterator pointing to the particle with the lowest id.
   *
   * Returns an random access iterator that traverses the particle
   * in order of ascending id. If the cache is not up-to-date,
   * an update is triggered. This iterator stays valid as long
   * as the cache is valid. Since the cache could be invalidated
   * and updated elsewhere, iterators into the cache should not be
   * stored.
   */
  value_iterator begin() {
    assert(Communication::mpiCallbacks().comm().rank() == 0);

    if (!m_valid)
      update();

    return value_iterator(remote_parts.begin());
  }

  /**
   * @brief Iterator pointing past the particle with the highest id.
   *
   * If the cache is not up-to-date,
   * an update is triggered.
   */
  value_iterator end() {
    assert(Communication::mpiCallbacks().comm().rank() == 0);

    if (!m_valid)
      update();

    return value_iterator(remote_parts.end());
  }

  /**
   * @brief Returns true if the cache is up-to-date.
   *
   * If false, particle access will trigger an update.
   */
  bool valid() const { return m_valid; }

  /**
   * @brief Returns true if the bond cache is up-to-date.
   *
   * If false, particle access will trigger an update.
   */
  bool valid_bonds() const { return m_valid_bonds; }

  /**
   * @brief Invalidate the cache and free memory.
   */
  void invalidate() {
    clear();
    /* Release memory */
    remote_parts.shrink_to_fit();
    bond_info.shrink_to_fit();
    /* Adjust state */
    m_valid = false;
    m_valid_bonds = false;
  }

  /**
   * @brief Update bond information.
   *
   * If the particle data is not valid,
   * it will be updated first.
   */
  void update_bonds() {
    update();

    if (!m_valid_bonds) {
      update_bonds_cb.call();
      m_recv_bonds();
      m_valid_bonds = true;
    }
  }

  /**
   * @brief Update particle information.
   *
   * This triggers a global update. All nodes
   * sort their particle by id, and send them
   * to the master.
   *
   * Complexity: N/P*log(N/P) for the sorting,
   * 2*N*(1 - 0.5^(log(p) + 1)) for the merging.
   */
  void update() {
    if (m_valid)
      return;

    update_cb.call();

    m_update();
    m_update_index();

    m_valid = true;
  }

  /** Number of particles in the config.
    *
     * Complexity: O(1)
  */
  size_t size() {
    assert(Communication::mpiCallbacks().comm().rank() == 0);

    if (!m_valid)
      update();

    return remote_parts.size();
  }

  /**
   * @brief size() == 0 ?
   *
   * Complexity: O(1)
   */
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
   *
   * Complexity: O(1)
   */
  Particle const &operator[](int id) {
    assert(Communication::mpiCallbacks().comm().rank() == 0);

    if (!m_valid)
      update();

    return remote_parts.begin()[id_index.at(id)];
  }
};

#endif
