#ifndef CORE_UTILS_PART_CFG_HPP
#define CORE_UTILS_PART_CFG_HPP

#include <algorithm>
#include <cassert>
#include <functional>
#include <iterator>
#include <map>
#include <memory>
#include <utility>
#include <vector>

#include <boost/iterator/transform_iterator.hpp>

#include "utils/parallel/Callback.hpp"

namespace detail {
class TakeSecond {
public:
  template <typename T, typename U> U &operator()(std::pair<T, U> &p) const {
    return p.second;
  }
};
}

template <typename Cells,
          typename Range = typename std::remove_reference<decltype(
              std::declval<Cells>().particles())>::type,
          typename Particle = typename std::iterator_traits<
              typename Range::iterator>::value_type>
class PartCfg {
  using map_type = std::map<int, std::reference_wrapper<Particle>>;

  map_type mapping;
  std::vector<Particle> remote_parts;
  std::vector<int> bond_info;
  bool m_valid, m_valid_bonds;

  Utils::Parallel::Callback update_cb;
  Utils::Parallel::Callback update_bonds_cb;

  Cells &cells;

  template <typename Container> void m_update_references(Container &range) {
    for (auto &p : range) {
      mapping.insert({p.identity(), std::ref(p)});
    }
  }

  void m_update_bonds() {
    std::vector<int> local_bonds;

    for (auto &p : cells.particles()) {
      std::copy(p.bl.begin(), p.bl.end(), std::back_inserter(local_bonds));
    }

    auto const size = local_bonds.size();
    MPI_Gather(&size, 1, MPI_INT, NULL, 1, MPI_INT, 0,
               Communication::mpiCallbacks().comm());
    MPI_Send(local_bonds.data(), local_bonds.size(), MPI_INT, 0, 42,
             Communication::mpiCallbacks().comm());
  }

  void m_update() {
    auto particles = cells.particles();
    auto const local_parts = particles.size();

    MPI_Gather(&local_parts, 1, MPI_INT, NULL, 1, MPI_INT, 0,
               Communication::mpiCallbacks().comm());

    remote_parts.clear();
    remote_parts.reserve(local_parts);
    std::copy(std::begin(particles), std::end(particles),
              std::back_inserter(remote_parts));

    MPI_Send(remote_parts.data(), remote_parts.size() * sizeof(Particle),
             MPI_BYTE, 0, 42, Communication::mpiCallbacks().comm());
  }

  void m_recv_bonds() {
    std::vector<int> sizes(Communication::mpiCallbacks().comm().size());
    int ign{0};

    MPI_Gather(&ign, 1, MPI_INT, sizes.data(), 1, MPI_INT, 0,
               Communication::mpiCallbacks().comm());
    auto const remote_size =
        std::accumulate(std::begin(sizes), std::end(sizes), 0);

    bond_info.resize(remote_size);

    int offset = 0;
    for (int i = 1; i < sizes.size(); i++) {
      MPI_Recv(&(bond_info[offset]), sizes[i], MPI_INT, i, 42,
               Communication::mpiCallbacks().comm(), MPI_STATUS_IGNORE);
      offset += sizes[i];
    }

    auto it = bond_info.begin();
    for (auto &p : remote_parts) {
      p.bl.e = nullptr;
      p.bl.max = 0;
      p.bl.resize(p.bl.size());

      std::copy_n(it, p.bl.size(), p.bl.begin());
      it += p.bl.size();
    }
  }

public:
  using value_iterator =
      boost::transform_iterator<detail::TakeSecond, typename map_type::iterator,
                                Particle &, Particle>;

  PartCfg() = delete;
  PartCfg(Cells &cells)
      : update_cb([this](int, int) { this->m_update(); }),
        update_bonds_cb([this](int, int) { this->m_update_bonds(); }),
        cells(cells), m_valid(false), m_valid_bonds(false) {}
  PartCfg(PartCfg const &) = delete;
  PartCfg(PartCfg &&) = delete;

  void clear() {
    mapping.clear();
    remote_parts.clear();
    bond_info.clear();
  }

  value_iterator begin() {
    assert(Communication::mpiCallbacks().comm().rank() == 0);

    if (!m_valid)
      update();

    return value_iterator(mapping.begin());
  }

  value_iterator end() {
    assert(Communication::mpiCallbacks().comm().rank() == 0);

    if (!m_valid)
      update();

    return value_iterator(mapping.end());
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

    std::vector<int> sizes(Communication::mpiCallbacks().comm().size());
    int ign{0};
    MPI_Gather(&ign, 1, MPI_INT, sizes.data(), 1, MPI_INT, 0,
               Communication::mpiCallbacks().comm());
    auto const remote_size =
        std::accumulate(std::begin(sizes), std::end(sizes), 0);

    remote_parts.resize(remote_size);

    std::vector<MPI_Request> reqs(sizes.size() - 1);

    int offset = 0;
    for (int i = 1; i < sizes.size(); i++) {
      MPI_Irecv(&(remote_parts[offset]), sizes[i] * sizeof(Particle), MPI_BYTE,
                i, 42, Communication::mpiCallbacks().comm(), &reqs[i - 1]);
      offset += sizes[i];
    }

    m_update_references(cells.particles());

    MPI_Waitall(reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);

    m_update_references(remote_parts);

    m_valid = true;
  }

  /** Number of particles in the config. */
  size_t size() {
    assert(Communication::mpiCallbacks().comm().rank() == 0);

    if (!m_valid)
      update();

    return mapping.size();
  }
  
  bool empty() {
    assert(Communication::mpiCallbacks().comm().rank() == 0);

    if (!m_valid)
      update();

    return mapping.empty();
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

    return mapping.at(id);
  }
};

#endif
