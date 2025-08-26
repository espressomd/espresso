/*
 * Copyright (C) 2010-2022 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "cell_system/ParticleDecomposition.hpp"

#include "BoxGeometry.hpp"
#include "LocalBox.hpp"
#include "Particle.hpp"
#include "ParticleList.hpp"
#include "ParticleRange.hpp"
#include "algorithm/link_cell.hpp"
#include "bond_error.hpp"
#include "cell_system/Cell.hpp"
#include "cell_system/CellStructureType.hpp"
#include "config/config.hpp"
#include "ghosts.hpp"
#include "system/Leaf.hpp"

#include <utils/Vector.hpp>

#include <boost/container/static_vector.hpp>
#include <boost/iterator/indirect_iterator.hpp>
#include <boost/range/algorithm/transform.hpp>

#include <algorithm>
#include <cassert>
#include <concepts>
#include <functional>
#include <iterator>
#include <memory>
#include <optional>
#include <set>
#include <span>
#include <stdexcept>
#include <unordered_set>
#include <utility>
#include <vector>

#ifdef ESPRESSO_CALIPER
#include <caliper/cali.h>
#endif

// forward declarations
#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
namespace Kokkos {
template <class DataType, class... Properties> class View;
class HostSpace;
struct LayoutRight;
template <unsigned T> struct MemoryTraits;
} // namespace Kokkos
namespace Cabana {
class HalfNeighborTag;
struct VerletLayout2D;
class TeamVectorOpTag;
template <typename... Types> struct MemberTypes;
template <class DataType, class MemorySpace, int, class MemoryTraits>
class AoSoA;
} // namespace Cabana
namespace Communication {
struct KokkosHandle;
} // namespace Communication
template <class MemorySpace, class ListAlgorithm, class Layout, class BuildTag>
class CustomVerletList;
#endif // ESPRESSO_SHARED_MEMORY_PARALLELISM

template <typename Callable>
concept ParticleCallback = requires(Callable c, Particle &p) {
  { c(p) } -> std::same_as<void>;
};

using ParticleUnaryOp = std::function<void(Particle &)>;

namespace Cells {
enum Resort : unsigned {
  RESORT_NONE = 0u,
  RESORT_LOCAL = 1u,
  RESORT_GLOBAL = 2u
};

/**
 * @brief Flags to select particle parts for communication.
 */
enum DataPart : unsigned {
  DATA_PART_NONE = 0u,       /**< Nothing */
  DATA_PART_PROPERTIES = 1u, /**< Particle::p */
  DATA_PART_POSITION = 2u,   /**< Particle::r */
  DATA_PART_MOMENTUM = 8u,   /**< Particle::m */
  DATA_PART_FORCE = 16u,     /**< Particle::f */
#ifdef ESPRESSO_BOND_CONSTRAINT
  DATA_PART_RATTLE = 32u, /**< Particle::rattle */
#endif
  DATA_PART_BONDS = 64u /**< Particle::bonds */
};
} // namespace Cells

/**
 * @brief Map the data parts flags from cells to those
 *        used internally by the ghost communication.
 *
 * @param data_parts data parts flags
 * @return ghost communication flags
 */
unsigned map_data_parts(unsigned data_parts);

namespace Cells {
inline ParticleRange particles(std::span<Cell *const> cells) {
  /* Find first non-empty cell */
  auto first_non_empty = std::ranges::find_if(
      cells, [](const Cell *c) { return not c->particles().empty(); });

  return {CellParticleIterator(first_non_empty, cells.end()),
          CellParticleIterator(cells.end())};
}
} // namespace Cells

/**
 * @brief Distance vector and length handed to pair kernels.
 */
struct Distance {
  explicit Distance(Utils::Vector3d const &vec21)
      : vec21(vec21), dist2(vec21.norm2()) {}

  Utils::Vector3d vec21;
  double dist2;
};

namespace detail {
// NOLINTNEXTLINE(bugprone-exception-escape)
struct MinimalImageDistance {
  BoxGeometry const box;

  Distance operator()(Particle const &p1, Particle const &p2) const {
    return Distance(box.get_mi_vector(p1.pos(), p2.pos()));
  }
};

struct EuclidianDistance {
  Distance operator()(Particle const &p1, Particle const &p2) const {
    return Distance(p1.pos() - p2.pos());
  }
};
} // namespace detail

/** Describes a cell structure / cell system. Contains information
 *  about the communication of cell contents (particles, ghosts, ...)
 *  between different nodes and the relation between particle
 *  positions and the cell system. All other properties of the cell
 *  system which are not common between different cell systems have to
 *  be stored in separate structures.
 */
class CellStructure : public System::Leaf<CellStructure> {
#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
public:
  static constexpr auto vector_length = 1;
  struct AoSoA_pack;
  using ForceType = Kokkos::View<double **[3], Kokkos::LayoutRight>;
  using VirialType = Kokkos::View<double *[3], Kokkos::LayoutRight>;
  using data_types = Cabana::MemberTypes<double[3], double, int, int>;
  using memory_space = Kokkos::HostSpace;
  using AoSoAType = Cabana::AoSoA<data_types, memory_space, vector_length,
                                  Kokkos::MemoryTraits<0>>;
  using ListAlgorithm = Cabana::HalfNeighborTag;
  using ListType =
      CustomVerletList<Kokkos::HostSpace, ListAlgorithm, Cabana::VerletLayout2D,
                       Cabana::TeamVectorOpTag>;
#endif // ESPRESSO_SHARED_MEMORY_PARALLELISM

private:
  /** The local id-to-particle index */
  std::vector<Particle *> m_particle_index;
  /** Implementation of the primary particle decomposition */
  std::unique_ptr<ParticleDecomposition> m_decomposition;
  /** Active type in m_decomposition */
  CellStructureType m_type = CellStructureType::NSQUARE;
  /** One of @ref Cells::Resort, announces the level of resort needed.
   */
  unsigned m_resort_particles = Cells::RESORT_NONE;
  bool m_verlet_skin_set = false;
  bool m_rebuild_verlet_list = true;
  bool m_rebuild_verlet_list_cabana = true;
  std::vector<std::pair<Particle *, Particle *>> m_verlet_list;
  double m_le_pos_offset_at_last_resort = 0.;
  /** @brief Verlet list skin. */
  double m_verlet_skin = 0.;
  double m_verlet_reuse = 0.;
#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
  int m_cached_max_local_particle_id = 0;
  int m_max_id = 0;
  std::unique_ptr<ForceType> m_local_force;
#ifdef ESPRESSO_ROTATION
  std::unique_ptr<ForceType> m_local_torque;
#endif
#ifdef ESPRESSO_NPT
  std::unique_ptr<VirialType> m_local_virial;
#endif
  std::unique_ptr<ListType> m_verlet_list_cabana;
  std::unique_ptr<AoSoAType> m_particle_storage;
  /** particle properties for Cabana */
  std::unique_ptr<AoSoA_pack> m_aosoa;
  /** The local id-to-index for aosoa data */
  std::vector<Particle *> m_unique_particles;
  std::shared_ptr<Communication::KokkosHandle> m_kokkos_handle;
#endif // ESPRESSO_SHARED_MEMORY_PARALLELISM

public:
  CellStructure(BoxGeometry const &box);
  virtual ~CellStructure();

  bool use_verlet_list = true;

  /**
   * @brief Update local particle index.
   *
   * Update the entry for a particle in the local particle
   * index.
   *
   * @param id Entry to update.
   * @param p Pointer to the particle.
   */
  void update_particle_index(int id, Particle *p) {
    assert(id >= 0);
    // cppcheck-suppress assertWithSideEffect
    assert(not p or p->id() == id);

    if (static_cast<unsigned int>(id) >= m_particle_index.size())
      m_particle_index.resize(static_cast<unsigned int>(id + 1));

    m_particle_index[static_cast<unsigned int>(id)] = p;
  }

  /**
   * @brief Update local particle index.
   *
   * Update the entry for a particle in the local particle
   * index.
   *
   * @param p Pointer to the particle.
   */
  void update_particle_index(Particle &p) {
    update_particle_index(p.id(), std::addressof(p));
  }

  /**
   * @brief Update local particle index.
   *
   * @param pl List of particles whose index entries should be updated.
   */
  void update_particle_index(ParticleList &pl) {
    for (auto &p : pl) {
      update_particle_index(p.id(), std::addressof(p));
    }
  }

  /**
   * @brief Clear the particles index.
   */
  void clear_particle_index() { m_particle_index.clear(); }

private:
  /**
   * @brief Append a particle to a list and update this
   *        particle index accordingly.
   * @param pl List to add the particle to.
   * @param p Particle to add.
   */
  Particle &append_indexed_particle(ParticleList &pl, Particle &&p) {
    /* Check if cell may reallocate, in which case the index
     * entries for all particles in this cell have to be
     * updated. */
    auto const may_reallocate = pl.size() >= pl.capacity();
    auto &new_part = pl.insert(std::move(p));

    if (may_reallocate)
      update_particle_index(pl);
    else {
      update_particle_index(new_part);
    }

    return new_part;
  }

public:
  /**
   * @brief Get a local particle by id.
   *
   * @param id Particle to get.
   * @return Pointer to particle if it is local,
   *         nullptr otherwise.
   */
  Particle *get_local_particle(int id) {
    assert(id >= 0);

    if (static_cast<unsigned int>(id) >= m_particle_index.size())
      return nullptr;

    return m_particle_index[static_cast<unsigned int>(id)];
  }

  /** @overload */
  const Particle *get_local_particle(int id) const {
    assert(id >= 0);

    if (static_cast<unsigned int>(id) >= m_particle_index.size())
      return nullptr;

    return m_particle_index[static_cast<unsigned int>(id)];
  }

  template <class InputRange, class OutputIterator>
  void get_local_particles(InputRange ids, OutputIterator out) {
    std::ranges::transform(ids, out,
                           [this](int id) { return get_local_particle(id); });
  }

  CellStructureType decomposition_type() const { return m_type; }

  /** Maximal cutoff supported by current cell system. */
  Utils::Vector3d max_cutoff() const { return decomposition().max_cutoff(); }

  /** Maximal pair range supported by current cell system. */
  Utils::Vector3d max_range() const { return decomposition().max_range(); }

  ParticleRange local_particles() const {
    return Cells::particles(decomposition().local_cells());
  }

  ParticleRange ghost_particles() const {
    return Cells::particles(decomposition().ghost_cells());
  }

  std::size_t count_local_particles() const {
    std::size_t count = 0;
    for (auto const &cell : m_decomposition->local_cells()) {
      count += cell->particles().size();
    }
    return count;
  }

  /** @brief whether to use parallel version of @ref for_each_local_particle */
  bool use_parallel_for_each_local_particle() const {
#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
    return true;
#else
    return false;
#endif
  }

  /**
   * @brief Run a kernel on all local particles.
   * The kernel is assumed to be thread-safe.
   */
  void for_each_local_particle(ParticleUnaryOp &&f) const {
#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
    if (use_parallel_for_each_local_particle()) {
      parallel_for_each_particle_impl(decomposition().local_cells(), f);
      return;
    }
#endif
    for (auto &p : local_particles()) {
      f(p);
    }
  }

  /**
   * @brief Run a kernel on all ghost particles.
   * The kernel is assumed to be thread-safe.
   */
  void for_each_ghost_particle(ParticleUnaryOp &&f) const {
    for (auto &p : ghost_particles()) {
      f(p);
    }
  }

private:
  /** Cell system dependent function to find the right cell for a
   *  particle.
   *  \param  p Particle.
   *  \return pointer to cell where to put the particle, nullptr
   *          if the particle does not belong on this node.
   */
  Cell *particle_to_cell(const Particle &p) {
    return decomposition().particle_to_cell(p);
  }
  Cell const *particle_to_cell(const Particle &p) const {
    return decomposition().particle_to_cell(p);
  }

#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
  void parallel_for_each_particle_impl(std::span<Cell *const> cells,
                                       ParticleUnaryOp &f) const;
#endif

public:
  /**
   * @brief Add a particle.
   *
   * Moves a particle into the cell system. This adds
   * a particle to the local node, irrespective of where
   * it belongs.
   *
   * @param p Particle to add.
   * @return Pointer to the particle in the cell
   *         system.
   */
  Particle *add_particle(Particle &&p);

  /**
   * @brief Add a particle.
   *
   * Moves a particle into the cell system, if it
   * belongs to this node. Otherwise this does not
   * have an effect and the particle is discarded.
   * This can be used to add a particle without
   * knowledge where it should be placed by calling
   * the function on all nodes, it will then add
   * the particle in exactly one place.
   *
   * @param p Particle to add.
   * @return Pointer to particle if it is local, null
   *         otherwise.
   */
  Particle *add_local_particle(Particle &&p);

  /**
   * @brief Remove a particle.
   *
   * Removes a particle and all bonds pointing
   * to it. This is a collective call.
   *
   * @param id Id of particle to remove.
   */
  void remove_particle(int id);

  /**
   * @brief Get the maximal particle id on this node.
   *
   * This returns the highest particle id on
   * this node, or -1 if there are no particles on this node.
   */
  int get_max_local_particle_id() const;
#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
  int get_cached_max_local_particle_id() const {
    return m_cached_max_local_particle_id;
  }
#endif

  /**
   * @brief Remove all particles from the cell system.
   *
   * This allows linear time removal of all particles from
   * the system, removing each particle individually would
   * be quadratic.
   */
  void remove_all_particles();

  /**
   * @brief Get the underlying particle decomposition.
   *
   * Should be used solely for informative purposes.
   *
   * @return The active particle decomposition.
   */
  ParticleDecomposition const &decomposition() const {
    return assert(m_decomposition), *m_decomposition;
  }

private:
  ParticleDecomposition &decomposition() {
    return assert(m_decomposition), *m_decomposition;
  }

public:
  /**
   * @brief Increase the local resort level at least to @p level.
   */
  void set_resort_particles(Cells::Resort level) {
    m_resort_particles |= level;
    assert(m_resort_particles >= level);
  }

  /**
   * @brief Get the currently scheduled resort level.
   */
  unsigned get_resort_particles() const { return m_resort_particles; }

  /**
   * @brief Set the resort level to sorted.
   */
  void clear_resort_particles() { m_resort_particles = Cells::RESORT_NONE; }

  /**
   * @brief Check whether a particle has moved further than half the skin
   * since the last Verlet list update, thus requiring a resort.
   * @param additional_offset   Offset which is added to the distance the
   *                            particle has travelled when comparing to half
   *                            the Verlet skin (e.g., for Lees-Edwards BC).
   * @return Whether a resort is needed.
   */
  bool
  check_resort_required(Utils::Vector3d const &additional_offset = {}) const;

  auto get_le_pos_offset_at_last_resort() const {
    return m_le_pos_offset_at_last_resort;
  }

  /**
   * @brief Synchronize number of ghosts.
   */
  void ghosts_count();

  /**
   * @brief Update ghost particles.
   *
   * Update ghost particles with data from the real particles.
   *
   * @param data_parts Particle parts to update, combination of @ref
   * Cells::DataPart
   */
  void ghosts_update(unsigned data_parts);

  /**
   * @brief Update ghost particles, with particle resort if needed.
   *
   * Update ghost particles with data from the real particles.
   * Resort particles if a resort is due.
   *
   * @param data_parts Particle parts to update, combination of @ref
   * Cells::DataPart
   */
  void update_ghosts_and_resort_particle(unsigned data_parts);

  /**
   * @brief Add forces from ghost particles to real particles.
   */
  void ghosts_reduce_forces();

#ifdef ESPRESSO_BOND_CONSTRAINT
  /**
   * @brief Add rattle corrections from ghost particles to real particles.
   */
  void ghosts_reduce_rattle_correction();
#endif

  /**
   * @brief Resort particles.
   */
  void resort_particles(bool global_flag);

  /** @brief Whether the Verlet skin is set. */
  auto is_verlet_skin_set() const { return m_verlet_skin_set; }

  /** @brief Get the Verlet skin. */
  auto get_verlet_skin() const { return m_verlet_skin; }

  /** @brief Set the Verlet skin. */
  void set_verlet_skin(double value);

  /** @brief Set the Verlet skin using a heuristic. */
  void set_verlet_skin_heuristic();

  void update_verlet_stats(int n_steps, int n_verlet_updates) {
    if (n_verlet_updates > 0) {
      m_verlet_reuse = n_steps / static_cast<double>(n_verlet_updates);
    } else {
      m_verlet_reuse = 0.;
    }
  }

  /** @brief Average number of integration steps the Verlet list was re-used */
  auto get_verlet_reuse() const { return m_verlet_reuse; }

  /**
   * @brief Resolve ids to particles.
   *
   * @throws BondResolutionError if one of the ids
   *         was not found.
   *
   * @param partner_ids Ids to resolve.
   * @return Vector of Particle pointers.
   */
  auto resolve_bond_partners(std::span<const int> partner_ids) {
    boost::container::static_vector<Particle *, 4> partners;
    get_local_particles(partner_ids, std::back_inserter(partners));

    /* Check if id resolution failed for any partner */
    if (std::ranges::find(partners, nullptr) != partners.end()) {
      throw BondResolutionError{};
    }

    return partners;
  }

private:
  /**
   * @brief Execute kernel for every bond on particle.
   * @tparam Handler Callable, which can be invoked with
   *                 (Particle, int, std::span<Particle *>),
   *                 returning a bool.
   * @param p Particles for whom the bonds are evaluated.
   * @param handler is called for every bond, and handed
   *                p, the bond id and a span with the bond
   *                partners as arguments. Its return value
   *                should indicate if the bond was broken.
   */
  template <class Handler>
  void execute_bond_handler(Particle &p, Handler handler) {
    for (const BondView bond : p.bonds()) {
      auto const partner_ids = bond.partner_ids();

      try {
        auto partners = resolve_bond_partners(partner_ids);
        auto const partners_span = std::span(partners.data(), partners.size());
        auto const bond_broken = handler(p, bond.bond_id(), partners_span);
        if (bond_broken) {
          bond_broken_error(p.id(), partner_ids);
        }
      } catch (const BondResolutionError &) {
        bond_broken_error(p.id(), partner_ids);
      }
    }
  }

  /**
   * @brief Go through ghost cells and remove the ghost entries from the
   * local particle index.
   */
  void invalidate_ghosts() {
    for (auto const &p : ghost_particles()) {
      if (get_local_particle(p.id()) == &p) {
        update_particle_index(p.id(), nullptr);
      }
    }
  }

  /** @brief Set the particle decomposition, keeping the particles. */
  void set_particle_decomposition(
      std::unique_ptr<ParticleDecomposition> &&decomposition) {
    clear_particle_index();

    /* Swap in new cell system */
    std::swap(m_decomposition, decomposition);

    /* Add particles to new system */
    for (auto &p : Cells::particles(decomposition->local_cells())) {
      add_particle(std::move(p));
    }
  }

public:
  /**
   * @brief Set the particle decomposition to @ref AtomDecomposition.
   */
  void set_atom_decomposition();

  /**
   * @brief Set the particle decomposition to @ref RegularDecomposition.
   *
   * @param range Interaction range.
   * @param fully_connected_boundary neighbor cell directions for Lees-Edwards.
   */
  void set_regular_decomposition(
      double range,
      std::optional<std::pair<int, int>> fully_connected_boundary);

  /**
   * @brief Set the particle decomposition to @ref HybridDecomposition.
   *
   * @param cutoff_regular Interaction cutoff_regular.
   * @param n_square_types Particle types to put into n_square decomposition.
   */
  void set_hybrid_decomposition(double cutoff_regular,
                                std::set<int> n_square_types);

private:
  /**
   * @brief Run link_cell algorithm for local cells.
   *
   * @tparam Kernel Needs to be callable with (Particle, Particle, Distance).
   * @param kernel Pair kernel functor.
   */
  template <class Kernel> void link_cell(Kernel kernel) {
    auto const maybe_box = decomposition().minimum_image_distance();
    auto const local_cells_span = decomposition().local_cells();
    auto const first = boost::make_indirect_iterator(local_cells_span.begin());
    auto const last = boost::make_indirect_iterator(local_cells_span.end());

    if (maybe_box) {
      Algorithm::link_cell(
          first, last,
          [&kernel, df = detail::MinimalImageDistance{decomposition().box()}](
              Particle &p1, Particle &p2) { kernel(p1, p2, df(p1, p2)); });
    } else {
      if (decomposition().box().type() != BoxType::CUBOID) {
        throw std::runtime_error("Non-cuboid box type is not compatible with a "
                                 "particle decomposition that relies on "
                                 "EuclideanDistance for distance calculation.");
      }
      Algorithm::link_cell(
          first, last,
          [&kernel, df = detail::EuclidianDistance{}](
              Particle &p1, Particle &p2) { kernel(p1, p2, df(p1, p2)); });
    }
  }

#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
public:
  auto get_max_id() const { return m_max_id; }

  void set_kokkos_handle(std::shared_ptr<Communication::KokkosHandle> handle);
  void rebuild_local_properties(double pair_cutoff);
  void reset_local_properties();
  void reset_local_force();

  auto &get_local_force() { return *m_local_force; }
#ifdef ESPRESSO_ROTATION
  auto &get_local_torque() { return *m_local_torque; }
#endif
#ifdef ESPRESSO_NPT
  auto &get_local_virial() { return *m_local_virial; }
#endif
  auto &get_aosoa() { return *m_aosoa; }
  auto const &get_unique_particles() const { return m_unique_particles; }
  auto const &get_verlet_list_cabana() const { return *m_verlet_list_cabana; }
  void clear_local_properties();

  [[nodiscard]] auto is_verlet_list_cabana_rebuild_needed() const {
    return m_rebuild_verlet_list_cabana or (not use_verlet_list);
  }

  /**
   * @brief Reset local properties of the Verlet list.
   * @param cutoff    Pair interaction cutoff.
   * @return True if a rebuild is needed.
   */
  [[nodiscard]] auto prepare_verlet_list_cabana(double cutoff) {
    auto const rebuild = is_verlet_list_cabana_rebuild_needed();
    if (rebuild) {
      // If we have to rebuild, we need to count the particles
      set_index_map(); // parallelized index_map
      // Create essential variables for MD
      rebuild_local_properties(cutoff);
    } else {
      // If we do not rebuild we can use the saved map
      reset_local_properties();
    }
    return rebuild;
  }

  void rebuild_verlet_list_cabana(auto &&kernel) {
    assert(is_verlet_list_cabana_rebuild_needed());
    kernel(m_decomposition->local_cells(), m_decomposition->box(),
           *m_verlet_list_cabana);
    m_rebuild_verlet_list_cabana = false;
  }

  void set_index_map();
#endif

private:
  /** Non-bonded pair loop with verlet lists.
   *
   * @param pair_kernel Kernel to apply
   * @param verlet_criterion Filter for verlet lists.
   */
  template <class PairKernel, class VerletCriterion>
  void verlet_list_loop(PairKernel pair_kernel,
                        const VerletCriterion &verlet_criterion) {
    /* In this case the verlet list update is attached to
     * the pair kernel, and the verlet list is rebuilt as
     * we go. */
    if (m_rebuild_verlet_list) {
      m_verlet_list.clear();

      link_cell([&](Particle &p1, Particle &p2, Distance const &d) {
        if (verlet_criterion(p1, p2, d)) {
          m_verlet_list.emplace_back(&p1, &p2);
          pair_kernel(p1, p2, d);
        }
      });

      m_rebuild_verlet_list = false;
      m_rebuild_verlet_list_cabana = true;
    } else {
      auto const maybe_box = decomposition().minimum_image_distance();
      /* In this case the pair kernel is just run over the verlet list. */
      if (maybe_box) {
        auto const distance_function =
            detail::MinimalImageDistance{decomposition().box()};
        for (auto const &[p1, p2] : m_verlet_list) {
          pair_kernel(*p1, *p2, distance_function(*p1, *p2));
        }
      } else {
        auto const distance_function = detail::EuclidianDistance{};
        for (auto const &[p1, p2] : m_verlet_list) {
          pair_kernel(*p1, *p2, distance_function(*p1, *p2));
        }
      }
    }
  }

public:
  /** Bonded pair loop.
   * @param bond_kernel Kernel to apply
   */
  template <class BondKernel> void bond_loop(BondKernel const &bond_kernel) {
    for (auto &p : local_particles()) {
      execute_bond_handler(p, bond_kernel);
    }
  }

  /** Non-bonded pair loop.
   * @param pair_kernel Kernel to apply
   */
  template <class PairKernel> void non_bonded_loop(PairKernel pair_kernel) {
    link_cell(pair_kernel);
  }

  /** Non-bonded pair loop with potential use
   * of verlet lists.
   * @param pair_kernel Kernel to apply
   * @param verlet_criterion Filter for verlet lists.
   */
  template <class PairKernel, class VerletCriterion>
  void non_bonded_loop(PairKernel pair_kernel,
                       const VerletCriterion &verlet_criterion) {
    if (use_verlet_list) {
      verlet_list_loop(pair_kernel, verlet_criterion);
    } else {
      /* No verlet lists, just run the kernel with pairs from the cells. */
      link_cell(pair_kernel);
    }
  }

  /**
   * @brief Check that particle index is commensurate with particles.
   *
   * For each local particles is checked that has a correct entry
   * in the particles index, and that there are no excess (non-existing)
   * particles in the index.
   */
  void check_particle_index() const;

  /**
   * @brief Check that particles are in the correct cell.
   *
   * This checks for all local particles that the result
   * of particles_to_cell is the cell the particles is
   * actually in, e.g. that the particles are sorted according
   * to particles_to_cell.
   */
  void check_particle_sorting() const;

public:
  /**
   * @brief Find cell a particle is stored in.
   *
   * For local particles, this returns the cell they
   * are stored in, otherwise nullptr is returned.
   *
   * @param p Particle to find cell for
   * @return Cell for particle or nullptr.
   */
  Cell *find_current_cell(const Particle &p) {
    assert(not get_resort_particles());

    if (p.is_ghost()) {
      return nullptr;
    }

    return particle_to_cell(p);
  }

  /**
   * @brief Run kernel on all particles inside local cell and its neighbors.
   *
   * @param p      Particle to find cell for
   * @param kernel Function with signature <tt>double(Particle const&,
   *               Particle const&, Utils::Vector3d const&)</tt>
   * @return false if cell is not found, otherwise true
   */
  template <class Kernel>
  bool run_on_particle_short_range_neighbors(Particle const &p,
                                             Kernel &kernel) {
    auto const cell = find_current_cell(p);

    if (cell == nullptr) {
      return false;
    }

    auto const maybe_box = decomposition().minimum_image_distance();

    if (maybe_box) {
      auto const distance_function =
          detail::MinimalImageDistance{decomposition().box()};
      short_range_neighbor_loop(p, cell, kernel, distance_function);
    } else {
      auto const distance_function = detail::EuclidianDistance{};
      short_range_neighbor_loop(p, cell, kernel, distance_function);
    }
    return true;
  }

private:
  template <class Kernel, class DistanceFunc>
  void short_range_neighbor_loop(Particle const &p1, Cell *const cell,
                                 Kernel &kernel, DistanceFunc const &df) {
    /* Iterate over particles inside cell */
    for (auto const &p2 : cell->particles()) {
      if (p1.id() != p2.id()) {
        auto const vec = df(p1, p2).vec21;
        kernel(p1, p2, vec);
      }
    }
    /* Iterate over all neighbors */
    for (auto const neighbor : cell->neighbors().all()) {
      /* Iterate over particles in neighbors */
      if (neighbor != cell) {
        for (auto const &p2 : neighbor->particles()) {
          auto const vec = df(p1, p2).vec21;
          kernel(p1, p2, vec);
        }
      }
    }
  }
};
