/*
 * Copyright (C) 2025-2026 The ESPResSo project
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

#include <config/config.hpp>

#include "cell_system/CellStructure.hpp"

#include "aosoa_pack.hpp"
#include "bond_forces_kokkos.hpp"
#include "custom_verlet_list.hpp"
#include "forces_cabana.hpp"
#include "kokkos_helpers.hpp"

#include <Cabana_Core.hpp>
#include <Cabana_NeighborList.hpp>

#ifdef ESPRESSO_CALIPER
#include "caliper_utils.hpp"
#endif

#include <functional>
#include <iterator>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

ESPRESSO_ATTR_ALWAYS_INLINE inline void
commit_particle(Particle const &p, auto const index,
                CellStructure::AoSoA_pack &aosoa, bool const rebuild) {
  // position/image/director/velocity/id/mass alias the authoritative
  // ParticleStore columns and are read by *store row* via row(i); they are
  // never copied here.
  // charge/dipm are refreshed per step (when a solver is active) in
  // refresh_pack_charges / refresh_pack_dipm, not here.
  // `type` is PACK-OWNED and written here on rebuild (read pack-indexed by the
  // hot pair kernels). The ParticleStore type column stays authoritative; a
  // mid-run type change forces a rebuild (System::on_particle_change ->
  // set_resort_particles), so this cache is refreshed before it is next read.
  // The has-exclusion flag (like `type`) is written ONLY on rebuild. Exclusions
  // cannot change between rebuilds: every add/delete goes through
  // on_particle_change -> set_resort_particles(RESORT_LOCAL), which forces the
  // next update_cabana_state onto the full-commit branch. The pack-owned
  // `flags` column persists across partial-update steps (per-step force reset
  // does not touch it; clear_local_properties zeroes it but also forces a
  // resort/rebuild). Reading the exclusions sidecar per particle per step was
  // pure waste -- keeping it out of the partial-update hot loop matters.
  if (rebuild) {
    aosoa.type(index) = p.type();
#ifdef ESPRESSO_EXCLUSIONS
    bool const has_exclusion = not p.exclusions().empty();
    aosoa.set_has_exclusion(index, has_exclusion);
    // Record the any-exclusion aggregate on the host. This is intentionally NOT
    // done inside the device-qualified set_has_exclusion: this commit sweep
    // runs on the host execution space, and the aggregate is a host
    // std::atomic. It is cleared before the sweep in update_cabana_state.
    if (has_exclusion) {
      aosoa.mark_any_exclusion();
    }
#else
    aosoa.flags(index) = 0;
#endif
  }
}

ESPRESSO_ATTR_ALWAYS_INLINE inline void
link_cell_kokkos(std::span<Cell *const> cells, BoxGeometry const &box_geo,
                 auto const &verlet_criterion,
                 Kokkos::View<int *, Kokkos::HostSpace> const &id_to_index,
                 int const max_id, auto const &intra_operator,
                 auto const &inter_operator) {

  // implementation detail: max_id refers to the max local particle id,
  // but ghost particles from other ranks may have larger particle ids;
  // -1 is used as a sentinel value for particle ids from other threads

  // Component-major layout only: gather the position column into an
  // INTERLEAVED scratch buffer once per Verlet build. The build examines
  // O(N * candidates) positions by row; under LayoutLeft each such read
  // touches three far-apart component streams (three cache lines per
  // candidate), which perf measured at ~2x the build cost at 4000
  // particles/core. The O(N) sequential gather below is prefetch-friendly
  // and runs at rebuild cadence only. The kernels consume base + runtime
  // strides, so they read the scratch ({3, 1}) and the native column
  // (view strides) through the same code path. Under particle-major the
  // column is already interleaved and the scratch is skipped entirely.
  std::vector<double> interleaved_positions;
  if constexpr (std::is_same_v<ParticleStore::StateVectorLayout,
                               Kokkos::LayoutLeft>) {
    if (not cells.empty()) {
      auto &store = cells.front()->store();
      auto const &position_view = store.position_view();
      auto const total = store.number_of_particles();
      interleaved_positions.resize(3u * total);
      for (std::size_t row = 0u; row < total; ++row) {
        interleaved_positions[3u * row + 0u] = position_view(row, 0);
        interleaved_positions[3u * row + 1u] = position_view(row, 1);
        interleaved_positions[3u * row + 2u] = position_view(row, 2);
      }
    }
  }

  // Hoist the cuboid minimum-image parameters by value: the fold runs once
  // per candidate pair and must not chase the BoxGeometry reference for its
  // box lengths every time. Lees-Edwards boxes take the full BoxGeometry
  // path (shear offset handling).
  bool const has_lees_edwards = box_geo.type() == BoxType::LEES_EDWARDS;
  auto const cuboid_minimum_image = box_geo.cuboid_minimum_image();
  auto const minimum_image_dist2 =
      [&box_geo, has_lees_edwards, cuboid_minimum_image](
          Utils::Vector3d const &a, Utils::Vector3d const &b) {
        return has_lees_edwards ? box_geo.get_mi_dist2(a, b)
                                : cuboid_minimum_image.dist2(a, b);
      };

  // Iterate the cells' store-ROW bags directly and REBIND two cached views
  // (p1 + partner) per work item via Particle::attach_to_store, instead of
  // driving RowParticleRange iterators (each embeds a Particle by value, so
  // std::next(it) and per-neighbour range begin()/end() would build fresh
  // Particles). One reused view per role per work item (one cell per Kokkos
  // work item) is thread-safe. Iteration ORDER is unchanged.
  auto intra_kernel = [&cells, minimum_image_dist2, &verlet_criterion,
                       &id_to_index, &intra_operator, &interleaved_positions,
                       max_id](const int i) {
    auto &store = cells[i]->store();
    // Contiguous store-row range; clean store, so index directly.
    auto const offset = cells[i]->offset();
    auto const n = cells[i]->count();
    // Hoist raw id/position column pointers once per work item -- the
    // per-candidate id and position reads are the hottest loads of the Verlet
    // build, and the Particle accessors cost an attached-branch + view deref
    // per candidate. LayoutRight vector columns are row-contiguous (element
    // (row, j) at data() + 3 * row + j). The Particle views stay attached for
    // the verlet_criterion (it may read types/flags); iteration order and
    // arithmetic are unchanged.
    auto const *const id_column = store.id_view().data();
    // Layout-agnostic hoisted position access: base pointer plus the view's
    // run-time strides (row stride, component stride). Particle-major
    // (LayoutRight) gives {3, 1}; component-major (LayoutLeft) gives
    // {1, padded extent}. See the StateVectorLayout toggle in ParticleStore.
    auto const &position_view = store.position_view();
    // Prefer the interleaved rebuild-cadence scratch when populated
    // (component-major layout); otherwise read the native column. Same
    // base-plus-strides consumption either way.
    bool const use_scratch = not interleaved_positions.empty();
    auto const *const position_column =
        use_scratch ? interleaved_positions.data() : position_view.data();
    auto const pos_row_stride =
        use_scratch ? std::size_t{3u}
                    : static_cast<std::size_t>(position_view.stride(0));
    auto const pos_comp_stride =
        use_scratch ? std::size_t{1u}
                    : static_cast<std::size_t>(position_view.stride(1));
    Particle p1, p2;
    for (std::size_t a = 0u; a < n; ++a) {
      auto const row_a = offset + a;
      p1.attach_to_store(store, static_cast<int>(row_a));
      if (id_column[row_a] <= max_id) {
        auto const ii = id_to_index(id_column[row_a]);
        if (ii >= 0) {
          // Hoist p1's position out of the inner loop (one read per outer
          // particle instead of per pair candidate).
          auto const *const p1_base = position_column + row_a * pos_row_stride;
          auto const p1_pos =
              Utils::Vector3d{p1_base[0u], p1_base[pos_comp_stride],
                              p1_base[2u * pos_comp_stride]};
          // pairs in this cell (j > i), same order as before
          for (std::size_t b = a + 1u; b < n; ++b) {
            auto const row_b = offset + b;
            if (id_column[row_b] <= max_id) {
              p2.attach_to_store(store, static_cast<int>(row_b));
              auto const *const p2_base =
                  position_column + row_b * pos_row_stride;
              auto const p2_pos =
                  Utils::Vector3d{p2_base[0u], p2_base[pos_comp_stride],
                                  p2_base[2u * pos_comp_stride]};
              if (verlet_criterion(p1, p2,
                                   minimum_image_dist2(p1_pos, p2_pos))) {
                auto const jj = id_to_index(id_column[row_b]);
                if (jj >= 0) {
                  intra_operator(ii, jj);
                }
              }
            }
          }
        }
      }
    }
  };

  auto inter_kernel = [&cells, minimum_image_dist2, &verlet_criterion,
                       &id_to_index, &inter_operator, &interleaved_positions,
                       max_id](const int i) {
    auto &store = cells[i]->store();
    // Contiguous store-row range; clean store, so index directly.
    auto const offset = cells[i]->offset();
    auto const n = cells[i]->count();
    // Hoisted raw column pointers: see intra_kernel (incl. the layout-agnostic
    // stride note). All cells share the one active ParticleStore, so the
    // pointers are valid for neighbor cells too.
    auto const *const id_column = store.id_view().data();
    auto const &position_view = store.position_view();
    // Prefer the interleaved rebuild-cadence scratch when populated
    // (component-major layout); otherwise read the native column. Same
    // base-plus-strides consumption either way.
    bool const use_scratch = not interleaved_positions.empty();
    auto const *const position_column =
        use_scratch ? interleaved_positions.data() : position_view.data();
    auto const pos_row_stride =
        use_scratch ? std::size_t{3u}
                    : static_cast<std::size_t>(position_view.stride(0));
    auto const pos_comp_stride =
        use_scratch ? std::size_t{1u}
                    : static_cast<std::size_t>(position_view.stride(1));
    Particle p1, p2;
    for (std::size_t a = 0u; a < n; ++a) {
      auto const row_a = offset + a;
      p1.attach_to_store(store, static_cast<int>(row_a));
      if (id_column[row_a] <= max_id) {
        auto const ii = id_to_index(id_column[row_a]);
        if (ii >= 0) {
          // Hoist p1's position out of the inner loops (see intra_kernel).
          auto const *const p1_base = position_column + row_a * pos_row_stride;
          auto const p1_pos =
              Utils::Vector3d{p1_base[0u], p1_base[pos_comp_stride],
                              p1_base[2u * pos_comp_stride]};
          // pairs with neighboring cells, same order as before
          for (auto *neighbor : cells[i]->neighbors().red()) {
            auto &nb_store = neighbor->store();
            auto const nb_offset = neighbor->offset();
            auto const nb_n = neighbor->count();
            for (std::size_t k = 0u; k < nb_n; ++k) {
              auto const row_k = nb_offset + k;
              if (id_column[row_k] <= max_id) {
                p2.attach_to_store(nb_store, static_cast<int>(row_k));
                auto const *const p2_base =
                    position_column + row_k * pos_row_stride;
                auto const p2_pos =
                    Utils::Vector3d{p2_base[0u], p2_base[pos_comp_stride],
                                    p2_base[2u * pos_comp_stride]};
                if (verlet_criterion(p1, p2,
                                     minimum_image_dist2(p1_pos, p2_pos))) {
                  auto const jj = id_to_index(id_column[row_k]);
                  if (jj >= 0) {
                    inter_operator(ii, jj);
                  }
                }
              }
            }
          }
        }
      }
    }
  };

  kokkos_parallel_range_for<Kokkos::DefaultHostExecutionSpace>(
      "intra", std::size_t{0}, cells.size(), intra_kernel);
  Kokkos::fence();

  kokkos_parallel_range_for<Kokkos::DefaultHostExecutionSpace>(
      "inter", std::size_t{0}, cells.size(), inter_kernel);
  Kokkos::fence();
}

// @p make_verlet_criterion is a nullary factory: constructing the criterion
// fills an O(n_types^2) cutoff table, so it is only invoked when the Verlet
// list is actually rebuilt, not on every force call.
template <class execution_space = Kokkos::DefaultHostExecutionSpace>
ESPRESSO_ATTR_ALWAYS_INLINE inline void
update_cabana_state(CellStructure &cell_structure,
                    auto const &make_verlet_criterion, double const pair_cutoff,
                    auto const integ_switch) {
#ifdef ESPRESSO_CALIPER
  ESPRESSO_CALI_MARK_FUNCTION;
#endif
  auto const rebuild = cell_structure.prepare_verlet_list_cabana(pair_cutoff);
  auto const &unique_particles = cell_structure.get_unique_particles();
  auto const n_part = unique_particles.size();
  auto const max_id = cell_structure.get_cached_max_local_particle_id();
  // Recompute the store-side derived director, then point the pack's
  // store-aliased views at the current store columns + translation/director
  // views.
#if defined(ESPRESSO_GAY_BERNE) or defined(ESPRESSO_DIPOLES)
  cell_structure.update_director_view();
#endif
  cell_structure.bind_pack_store_views();
  auto &aosoa = cell_structure.get_aosoa();

  if (rebuild) {
    auto &id_to_index = cell_structure.get_id_to_index();

    // ===================================================
    // Fill particle storage (full commit)
    // ===================================================
#ifdef ESPRESSO_CALIPER
    ESPRESSO_CALI_MARK_BEGIN("AoSoA commit full");
#endif
    int pair_count = 0;
    int angle_count = 0;
    int dihedral_count = 0;
#ifdef ESPRESSO_EXCLUSIONS
    // commit_particle accumulates the any-exclusion aggregate below; clear it
    // before the sweep repopulates it (read O(1) at the dispatch gate).
    aosoa.reset_any_exclusion();
#endif
    kokkos_parallel_range_for<execution_space>(
        "AoSoA write", std::size_t{0}, n_part,
        [&unique_particles, &aosoa, &id_to_index, &cell_structure, &pair_count,
         &angle_count, &dihedral_count](int const index) {
          auto const &p = *unique_particles.at(index);
          commit_particle(p, index, aosoa, true);
          id_to_index(p.id()) = index;
          if (not p.is_ghost()) {
            cell_structure.update_bond_storage(pair_count, angle_count,
                                               dihedral_count, p);
          }
        });
    Kokkos::fence();
    using host_space = Kokkos::DefaultHostExecutionSpace;
    auto &bs = cell_structure.bond_state();
    if (pair_count) {
      auto &pair_bond_list = bs.pair_list;
      kokkos_parallel_range_for<host_space>(
          "resolve_pair_bond_indices", std::size_t{0}, pair_count,
          [&pair_bond_list, &id_to_index](int idx) {
            for (int col = 0; col < 2; ++col) {
              pair_bond_list(idx, col) = id_to_index(pair_bond_list(idx, col));
            }
          });
    }
    if (angle_count) {
      auto &angle_bond_list = bs.angle_list;
      kokkos_parallel_range_for<host_space>(
          "resolve_angle_bond_indices", std::size_t{0}, angle_count,
          [&angle_bond_list, &id_to_index](int idx) {
            for (int col = 0; col < 3; ++col) {
              angle_bond_list(idx, col) =
                  id_to_index(angle_bond_list(idx, col));
            }
          });
    }
    if (dihedral_count) {
      auto &dihedral_bond_list = bs.dihedral_list;
      kokkos_parallel_range_for<host_space>(
          "resolve_dihedral_bond_indices", std::size_t{0}, dihedral_count,
          [&dihedral_bond_list, &id_to_index](int idx) {
            for (int col = 0; col < 4; ++col) {
              dihedral_bond_list(idx, col) =
                  id_to_index(dihedral_bond_list(idx, col));
            }
          });
    }
    if (pair_count != 0 or angle_count != 0 or dihedral_count != 0) {
      Kokkos::fence();
    }
#ifdef ESPRESSO_CALIPER
    ESPRESSO_CALI_MARK_END("AoSoA commit full");
#endif

    // ===================================================
    // Get Verlet pairs and fill Verlet list
    // ===================================================
    bool rebuild_vl = (integ_switch != INTEG_METHOD_STEEPEST_DESCENT and
                       cell_structure.use_verlet_list);
#ifdef ESPRESSO_CALIPER
    ESPRESSO_CALI_MARK_BEGIN("Verlet list creation");
#endif
    cell_structure.rebuild_verlet_list_cabana(
        [&](std::span<Cell *const> cells, BoxGeometry const &box,
            CellStructure::ListType &verlet_list) {
          auto const verlet_criterion = make_verlet_criterion();
          link_cell_kokkos(
              std::move(cells), box, verlet_criterion, id_to_index, max_id,
              [&](const int i, const int j) {
                // intra cell loop
                verlet_list.addNeighborLB(i, j);
              },
              [&](const int i, const int j) {
                // inter cell loop
                verlet_list.addNeighbor(i, j);
              });

          if (verlet_list.hasOverflow()) {
            cell_structure.use_verlet_list = false;
            runtimeWarningMsg()
                << "Verlet list overflow detected: neighbor count exceeded "
                   "max_counts. Falling back to the link cell algorithm. "
                   "Configured max is "
                << Cabana::NeighborList<CellStructure::ListType>::maxNeighbor(
                       verlet_list);
          }
        },
        rebuild_vl);
#ifdef ESPRESSO_CALIPER
    ESPRESSO_CALI_MARK_END("Verlet list creation");
#endif
  }
  // Partial-update (no-rebuild) branch: NOTHING to commit. The pack does not
  // own per-step particle data (position/velocity/etc. alias the ParticleStore
  // columns, read by store row), and the only pack-owned commit writes --
  // `type` and the has-exclusion `flags` -- are rebuild-cadence data written
  // on the full-commit branch above (both are refreshed by a forced rebuild
  // whenever they can change). The store-view rebind (bind_pack_store_views
  // above) already repointed the aliased views for this step.
}

#ifdef ESPRESSO_ELECTROSTATICS
// Refresh the PACK-OWNED contiguous charge column from the authoritative
// ParticleStore q column. This runs once per step (O(N)) ONLY when a coulomb
// actor is active; the hot pair loop and the P3M gather/spread loops then read
// `aosoa.charge(pack_index)` contiguously, which is O(pairs) >> O(N).
// Pure-LJ runs never call this and pay zero cost.
ESPRESSO_ATTR_ALWAYS_INLINE inline void
refresh_pack_charges(CellStructure &cell_structure) {
  using execution_space = Kokkos::DefaultHostExecutionSpace;
  auto const &unique_particles = cell_structure.get_unique_particles();
  auto const n_part = unique_particles.size();
  auto &aosoa = cell_structure.get_aosoa();
  kokkos_parallel_range_for<execution_space>(
      "refresh pack charges", std::size_t{0}, n_part,
      [&unique_particles, &aosoa](std::size_t const index) {
        aosoa.pair_charge(index) = unique_particles.at(index)->q();
      });
}
#endif

#ifdef ESPRESSO_DIPOLES
// Pack-owned dipm column refresh, guarded by an active dipolar actor
// (see refresh_pack_charges for the rationale).
ESPRESSO_ATTR_ALWAYS_INLINE inline void
refresh_pack_dipm(CellStructure &cell_structure) {
  using execution_space = Kokkos::DefaultHostExecutionSpace;
  auto const &unique_particles = cell_structure.get_unique_particles();
  auto const n_part = unique_particles.size();
  auto &aosoa = cell_structure.get_aosoa();
  kokkos_parallel_range_for<execution_space>(
      "refresh pack dipm", std::size_t{0}, n_part,
      [&unique_particles, &aosoa](std::size_t const index) {
        aosoa.pair_dipm(index) = unique_particles.at(index)->dipm();
      });
}
#endif

// Optional replacement for the Verlet-list pair loop. When set, it is invoked
// with the current Verlet list and the pack particle count instead of the
// generic Cabana::neighbor_parallel_for over @p nonbonded_kernel. forces.cpp
// installs the compile-time-specialized pair kernel this way when the active
// feature set allows it; energy/pressure leave it empty and keep the generic
// loop.
using ShortRangeVerletPairLoop =
    std::function<void(CellStructure::ListType const &, std::size_t)>;

#ifdef ESPRESSO_CUDA
namespace System {
class System;
} // namespace System

// Opt-in DEVICE (GPU) short-range pair-force path. Defined in
// forces_lj_device.cu (compiled by the CUDA compiler) so the heavy
// Kokkos/Cabana/System headers it needs never reach forces.cpp. Lennard-Jones
// is the first potential implemented; see forces_lj_device.cu for the gate.
bool gpu_core_enabled();
ShortRangeVerletPairLoop
create_device_short_range_pair_loop(System::System &system, bool has_coulomb,
                                    bool has_dipoles, bool has_elc,
                                    bool has_dpd, bool has_npt_virial);
#endif // ESPRESSO_CUDA

// @p make_verlet_criterion is a nullary factory: constructing the criterion
// fills an O(n_types^2) cutoff table, so it is only invoked on the link-cell
// fallback path, which is the only consumer here.
template <class execution_space = Kokkos::DefaultHostExecutionSpace>
void cabana_short_range(auto const &pair_bonds_kernel,
                        auto const &angle_bonds_kernel,
                        auto const &dihedral_bonds_kernel,
                        auto const &nonbonded_kernel,
                        CellStructure &cell_structure, double pair_cutoff,
                        double bond_cutoff, auto const &make_verlet_criterion,
                        auto const integ_switch,
                        ShortRangeVerletPairLoop const &verlet_pair_loop = {}) {
  assert(cell_structure.get_resort_particles() == Cells::RESORT_NONE);

  if (bond_cutoff >= 0.) {
#ifdef ESPRESSO_CALIPER
    ESPRESSO_CALI_MARK_BEGIN("cabana_bond_loop");
#endif
    using host_space = Kokkos::DefaultHostExecutionSpace;
    auto const n_pair_bonds = cell_structure.get_local_pair_bond_numbers();
    auto const n_angle_bonds = cell_structure.get_local_angle_bond_numbers();
    auto const n_dihedral_bonds =
        cell_structure.get_local_dihedral_bond_numbers();
    if (n_pair_bonds > 0) {
      kokkos_parallel_range_for<host_space>("for_each_local_pair_bonds",
                                            std::size_t{0}, n_pair_bonds,
                                            pair_bonds_kernel);
    }
    if (n_angle_bonds > 0) {
      kokkos_parallel_range_for<host_space>("for_each_local_angle_bonds",
                                            std::size_t{0}, n_angle_bonds,
                                            angle_bonds_kernel);
    }
    if (n_dihedral_bonds > 0) {
      kokkos_parallel_range_for<host_space>("for_each_local_dihedral_bonds",
                                            std::size_t{0}, n_dihedral_bonds,
                                            dihedral_bonds_kernel);
    }
    if (n_pair_bonds != 0 or n_angle_bonds != 0 or n_dihedral_bonds != 0) {
      Kokkos::fence();
    }
#ifdef ESPRESSO_CALIPER
    ESPRESSO_CALI_MARK_END("cabana_bond_loop");
#endif
  }

  // Cabana short range loop
  if (pair_cutoff > 0.) {
#ifdef ESPRESSO_CALIPER
    ESPRESSO_CALI_MARK_BEGIN("cabana_pair_loop");
#endif
    if (integ_switch != INTEG_METHOD_STEEPEST_DESCENT and
        cell_structure.use_verlet_list) {
      auto const &verlet_list = cell_structure.get_verlet_list_cabana();
      auto const n_particles = cell_structure.get_unique_particles().size();
      if (verlet_pair_loop) {
        verlet_pair_loop(verlet_list, n_particles);
      } else if (Kokkos::num_threads() == 1) {
        using NL = Cabana::NeighborList<CellStructure::ListType>;
        for (std::size_t i = 0; i < n_particles; ++i) {
          auto const nn = NL::numNeighbor(verlet_list, i);
          for (std::size_t n = 0; n < nn; ++n) {
            nonbonded_kernel(i, NL::getNeighbor(verlet_list, i, n));
          }
        }
      } else {
        Kokkos::RangePolicy<execution_space> policy(std::size_t{0},
                                                    n_particles);
        Cabana::neighbor_parallel_for(policy, nonbonded_kernel, verlet_list,
                                      Cabana::FirstNeighborsTag(),
                                      Cabana::SerialOpTag());
      }
    } else {
      cell_structure.cell_list_loop(
          [&](std::span<Cell *const> cells, BoxGeometry const &box) {
            auto const verlet_criterion = make_verlet_criterion();
            link_cell_kokkos(
                std::move(cells), box, verlet_criterion,
                cell_structure.get_id_to_index(),
                cell_structure.get_cached_max_local_particle_id(),
                [&](const int i, const int j) {
                  // intra cell loop
                  nonbonded_kernel(i, j);
                },
                [&](const int i, const int j) {
                  // inter cell loop
                  nonbonded_kernel(i, j);
                });
          });
    }
    Kokkos::fence();
#ifdef ESPRESSO_CALIPER
    ESPRESSO_CALI_MARK_END("cabana_pair_loop");
#endif
  }
}
