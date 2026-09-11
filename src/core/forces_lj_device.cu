/*
 * Copyright (C) 2026 The ESPResSo project
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

/** @file
 *  Opt-in DEVICE (GPU / Kokkos CUDA) short-range pair-force path.
 *
 *  Implemented on the device short-range path so far: Lennard-Jones (nonbonded)
 *  and P3M real-space Coulomb. @ref create_device_short_range_pair_loop gates
 *  on that feature subset (LJ-only among nonbonded pair potentials; P3M or no
 *  Coulomb; no dipoles/ELC/DPD/Thole/Gay-Berne/exclusions/NPT; cuboid box).
 *  Extending it to the full short-range kernel -- more nonbonded pair
 *  potentials, exclusions / Thole -- means generalizing @ref LJParamsDevice and
 *  the flat param table it feeds, and relaxing the gate below accordingly. The
 *  file and the factory are named generically ("short_range"); the LJ
 *  specifics (@ref LJParamsDevice, @ref lj_force_factor_device, the LJ-only
 *  @c active_pair_mask scan) are clearly the current subset.
 *
 *  This translation unit is compiled by the CUDA compiler (see
 *  src/core/CMakeLists.txt) so the heavy Kokkos/Cabana/System headers it needs
 *  never reach forces.cpp. forces.cpp calls the factory through the
 *  declaration in short_range_cabana.hpp; with the opt-in flag off (the
 *  default) the factory returns an empty @ref ShortRangeVerletPairLoop and the
 *  host path runs unchanged, so the CPU identity is preserved bit-for-bit.
 *
 *  The kernel runs on @c Kokkos::DefaultExecutionSpace (the device under
 *  @c Kokkos_ENABLE_CUDA, OpenMP on a host-only build). It is device-clean:
 *  it captures only Kokkos device Views and POD scalars by value, calls no
 *  host-only code, and throws nothing.
 */

#include <config/config.hpp>

#include "BoxGeometry.hpp"
#include "cell_system/CellStructure.hpp"
#include "electrostatics/coulomb.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "short_range_cabana.hpp"
#include "system/System.hpp"
#ifdef ESPRESSO_P3M
#include "electrostatics/p3m.hpp"
#endif

#include <utils/Vector.hpp>
#include <utils/math/AS_erfc_part.hpp>
#include <utils/math/int_pow.hpp>
#include <utils/math/sqr.hpp>

#include <Kokkos_Core.hpp>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <memory>
#include <numbers>
#include <variant>

/**
 * @brief Persistent device buffers for the GPU short-range pair loop, reused
 * across steps so no device memory is (re)allocated per force call. The
 * neighbor list / pack maps are refreshed only on a Verlet rebuild (keyed on
 * the Cabana generation); the force accumulator is re-zeroed each step.
 * Owned via a shared_ptr on CellStructure so it is destroyed before
 * Kokkos::finalize.
 */
struct DeviceShortRangeBuffers {
  using Exec = Kokkos::DefaultExecutionSpace;
  using Mem = Exec::memory_space;
  Kokkos::View<int *, Mem> counts;
  Kokkos::View<int **, Kokkos::LayoutRight, Mem> neighbors;
  Kokkos::View<int *, Mem> type;
  Kokkos::View<int *, Mem> row_map;
  Kokkos::View<double *[3], Kokkos::LayoutRight, Exec> force;
  // Distinct from any real generation, so the first call always copies.
  std::uint64_t cached_generation = std::numeric_limits<std::uint64_t>::max();

  void ensure_list(std::size_t n_counts, std::size_t n_nbr_rows,
                   std::size_t n_nbr_cols, std::size_t n_type,
                   std::size_t n_row_map) {
    if (counts.extent(0) != n_counts) {
      counts = decltype(counts)(
          Kokkos::view_alloc(Kokkos::WithoutInitializing, "sr_counts"),
          n_counts);
    }
    if (neighbors.extent(0) != n_nbr_rows or
        neighbors.extent(1) != n_nbr_cols) {
      neighbors = decltype(neighbors)(
          Kokkos::view_alloc(Kokkos::WithoutInitializing, "sr_neighbors"),
          n_nbr_rows, n_nbr_cols);
    }
    if (type.extent(0) != n_type) {
      type = decltype(type)(
          Kokkos::view_alloc(Kokkos::WithoutInitializing, "sr_type"), n_type);
    }
    if (row_map.extent(0) != n_row_map) {
      row_map = decltype(row_map)(
          Kokkos::view_alloc(Kokkos::WithoutInitializing, "sr_row_map"),
          n_row_map);
    }
  }
  void ensure_force(std::size_t n) {
    if (force.extent(0) != n) {
      force = decltype(force)(
          Kokkos::view_alloc(Kokkos::WithoutInitializing, "sr_force"), n);
    }
  }
};

/**
 * @brief Read the opt-in GPU-core flag once (cached).
 *
 * Enabled iff the environment variable @c ESPRESSO_GPU_CORE is exactly
 * @c "1". Read once at first call and cached for the process lifetime, so the
 * value is stable across steps. Defaults to @c false when the variable is
 * unset or has any other value, keeping the host path the default.
 */
bool gpu_core_enabled() {
  static bool const enabled = []() {
    char const *const value = std::getenv("ESPRESSO_GPU_CORE");
    return value != nullptr and std::strcmp(value, "1") == 0;
  }();
  return enabled;
}

/**
 * @brief POD Lennard-Jones parameters for one type pair, device-capturable.
 *
 * Mirrors @ref LJ_Parameters reduced to exactly what the force factor needs.
 * @c max_cut / @c min_cut are the *shifted* cutoffs
 * (@ref LJ_Parameters::max_cutoff / @ref LJ_Parameters::min_cutoff).
 */
struct LJParamsDevice {
  double eps;
  double sig;
  double offset;
  double max_cut;
  double min_cut;
};

/**
 * @brief Device Lennard-Jones force factor, matching @ref lj_pair_force_factor.
 *
 * Returns the scalar factor @c f such that the force on particle @c i is
 * <tt>f * d</tt> with <tt>d = pos_i - pos_j</tt> (minimum image), i.e. the
 * exact double-precision formula of the host path.
 */
KOKKOS_INLINE_FUNCTION double lj_force_factor_device(LJParamsDevice const &lj,
                                                     double dist) {
  if (dist < lj.max_cut and dist > lj.min_cut) {
    auto const r_off = dist - lj.offset;
    // int_pow<6> matches the host lj_pair_force_factor bit-for-bit (same
    // multiply tree); it is DEVICE_QUALIFIER constexpr, so device-callable.
    auto const frac6 = Utils::int_pow<6>(lj.sig / r_off);
    return 48.0 * lj.eps * frac6 * (frac6 - 0.5) / (r_off * dist);
  }
  return 0.0;
}

/**
 * @brief P3M real-space Coulomb parameters, device-capturable POD.
 * @c active is false when no P3M Coulomb solver is present (the coulomb term is
 * then skipped and the kernel is pure LJ).
 */
struct CoulombP3MParamsDevice {
  double prefactor;
  double alpha;
  double r_cut;
  bool active;
};

/**
 * @brief Device P3M real-space Coulomb force factor, matching
 * @ref CoulombP3M::pair_force with @c USE_ERFC_APPROXIMATION == 1 (the active
 * branch; see p3m/math.hpp). Returns @c f such that the force on @c i is
 * <tt>f * d</tt> (@c d = pos_i - pos_j, minimum image). @c AS_erfc_part is
 * @c constexpr, hence device-callable, so this is bit-for-bit the host formula.
 */
KOKKOS_INLINE_FUNCTION double
coulomb_p3m_force_factor(double q1q2, double dist,
                         CoulombP3MParamsDevice const &c) {
  if (q1q2 == 0. or dist >= c.r_cut or dist <= 0.) {
    return 0.;
  }
  auto const adist = c.alpha * dist;
  auto const exp_adist_sq = std::exp(-adist * adist);
  auto const dist_sq = dist * dist;
  auto const two_a_sqrt_pi_i = 2. * c.alpha * std::numbers::inv_sqrtpi;
  auto const erfc_part_ri = Utils::AS_erfc_part(adist) / dist;
  auto const fac = exp_adist_sq * (erfc_part_ri + two_a_sqrt_pi_i) / dist_sq;
  return fac * c.prefactor * q1q2;
}

/**
 * @brief Minimum-image coordinate difference (a - b), device-local copy of
 * detail::get_mi_coord_masked. Kept here so the shared, widely-included
 * BoxGeometry.hpp stays free of Kokkos/device annotations (which perturbed the
 * host FP codegen of the Coulomb short-range path). @p length_inv_masked is 0
 * for non-periodic directions, disabling the fold there.
 */
KOKKOS_INLINE_FUNCTION double device_mi_coord(double a, double b, double length,
                                              double length_inv_masked) {
  auto const dx = a - b;
  return dx - std::rint(dx * length_inv_masked) * length;
}

/**
 * @brief Build the flat device LJ parameter table.
 *
 * Fills a host mirror of a @c n_types*n_types table (row-major, index
 * <tt>ti*n_types + tj</tt>) from @p nonbonded_ias, then deep-copies it to the
 * device. The table is symmetric because @ref
 * InteractionsNonBonded::get_ia_param is symmetric.
 */
static Kokkos::View<LJParamsDevice *, Kokkos::DefaultExecutionSpace>
build_lj_device_param_table(InteractionsNonBonded const &nonbonded_ias,
                            int n_types) {
  auto const n =
      static_cast<std::size_t>(n_types) * static_cast<std::size_t>(n_types);
  Kokkos::View<LJParamsDevice *, Kokkos::DefaultExecutionSpace> table_device(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "lj_device_param_table"),
      n);
  auto table_host = Kokkos::create_mirror_view(table_device);
  for (int ti = 0; ti < n_types; ++ti) {
    for (int tj = 0; tj < n_types; ++tj) {
#ifdef ESPRESSO_LENNARD_JONES
      auto const &lj = nonbonded_ias.get_ia_param(ti, tj).lj;
      table_host(static_cast<std::size_t>(ti) *
                     static_cast<std::size_t>(n_types) +
                 static_cast<std::size_t>(tj)) = LJParamsDevice{
          lj.eps, lj.sig, lj.offset, lj.max_cutoff(), lj.min_cutoff()};
#else
      static_cast<void>(nonbonded_ias);
      table_host(
          static_cast<std::size_t>(ti) * static_cast<std::size_t>(n_types) +
          static_cast<std::size_t>(tj)) = LJParamsDevice{0., 0., 0., -1., 0.};
#endif
    }
  }
  Kokkos::deep_copy(table_device, table_host);
  return table_device;
}

/**
 * @brief Build the opt-in device short-range (pure-LJ) Verlet pair loop, or an
 * empty loop.
 *
 * Returns an empty @ref ShortRangeVerletPairLoop (leaving forces.cpp on the
 * host path) unless ALL of the following hold:
 *   - @ref gpu_core_enabled is @c true (env @c ESPRESSO_GPU_CORE == "1");
 *   - the box is a CUBOID;
 *   - no NPT virial is being accumulated (@p has_npt_virial is @c false);
 *   - no dipolar, ELC or DPD kernel is active
 *     (@p has_dipoles / @p has_elc / @p has_dpd are @c false);
 *   - no Thole or Gay-Berne pair potential is active;
 *   - no particle carries an exclusion;
 *   - NO Coulomb solver is active (@p has_coulomb is @c false);
 *   - every type pair's @c active_pair_mask is exactly the Lennard-Jones bit
 *     (no other pair potential contributes).
 *
 * The gate mirrors @ref create_specialized_verlet_pair_loop, tightened for the
 * pure-LJ device kernel. The caller passes the same feature predicates it uses
 * for that factory as concrete bools (no templates), so this factory can live
 * in a compiled translation unit and be declared in a header.
 */
ShortRangeVerletPairLoop create_device_short_range_pair_loop(
    System::System &system, [[maybe_unused]] bool has_coulomb,
    [[maybe_unused]] bool has_dipoles, [[maybe_unused]] bool has_elc,
    [[maybe_unused]] bool has_dpd, [[maybe_unused]] bool has_npt_virial) {
  if (not gpu_core_enabled())
    return {};
#ifndef ESPRESSO_LENNARD_JONES
  return {};
#else
  if (system.box_geo->type() != BoxType::CUBOID)
    return {};
#ifdef ESPRESSO_NPT
  if (has_npt_virial)
    return {};
#endif
#ifdef ESPRESSO_DIPOLES
  if (has_dipoles)
    return {};
#endif
  // P3M real-space Coulomb is supported on the device; any other Coulomb
  // method (or ELC) falls back to the host path.
  CoulombP3MParamsDevice coulomb_params{0., 0., 0., false};
#ifdef ESPRESSO_ELECTROSTATICS
  if (has_elc)
    return {};
  if (has_coulomb) {
#ifdef ESPRESSO_P3M
    CoulombP3M const *p3m = nullptr;
    if (auto &solver = system.coulomb.impl->solver; solver.has_value()) {
      if (std::holds_alternative<std::shared_ptr<CoulombP3M>>(*solver)) {
        p3m = std::get<std::shared_ptr<CoulombP3M>>(*solver).get();
      }
    }
    if (p3m == nullptr)
      return {}; // non-P3M Coulomb solver: host fallback
    coulomb_params = CoulombP3MParamsDevice{
        p3m->prefactor, p3m->p3m_params.alpha, p3m->p3m_params.r_cut, true};
#else
    return {}; // Coulomb active but P3M not compiled in: host fallback
#endif
  }
#else
  static_cast<void>(has_coulomb);
#endif
#ifdef ESPRESSO_DPD
  if (has_dpd)
    return {};
#endif

  auto const &nonbonded_ias = *system.nonbonded_ias;
  auto const max_type = nonbonded_ias.get_max_seen_particle_type();
  auto const n_types = max_type + 1;

  // Every active type pair must be *exactly* Lennard-Jones: the LJ bit set and
  // no other pair-potential bit. This also rejects Thole and Gay-Berne.
  constexpr auto lj_bit = pair_potential_bit(PairPotential::LennardJones);
  for (int type_i = 0; type_i <= max_type; ++type_i) {
    for (int type_j = type_i; type_j <= max_type; ++type_j) {
      auto const mask =
          nonbonded_ias.get_ia_param(type_i, type_j).active_pair_mask;
      if (mask != 0u and mask != lj_bit)
        return {};
    }
  }

  auto &cell_structure = *system.cell_structure;
#ifdef ESPRESSO_EXCLUSIONS
  {
    auto const &aosoa = cell_structure.get_aosoa();
    auto const n_pack = cell_structure.get_unique_particles().size();
    for (std::size_t i = 0u; i < n_pack; ++i) {
      if (aosoa.has_exclusion(i))
        return {};
    }
  }
#endif

  // POD cuboid fold parameters, computed host-side and captured by value into
  // the device kernel; box_l_inv is 0 for non-periodic directions (fold off).
  auto const box_l = system.box_geo->length();
  Utils::Vector3d box_l_inv{};
  for (auto d = 0u; d < 3u; ++d) {
    box_l_inv[d] = system.box_geo->periodic(d) ? 1.0 / box_l[d] : 0.0;
  }
  auto const max_cutoff_sq = Utils::sqr(system.maximal_cutoff());

  // Build the device param table once per loop construction (i.e. per step).
  auto const lj_table = build_lj_device_param_table(nonbonded_ias, n_types);

  return [&cell_structure, lj_table, n_types, box_l, box_l_inv, max_cutoff_sq,
          coulomb_params](CellStructure::ListType const &verlet_list,
                          std::size_t const n) {
    using ExecSpace = Kokkos::DefaultExecutionSpace;

    if (n == 0)
      return;

    // The kernel reads pack-indexed data: counts/neighbors (neighbors(i,k) is a
    // pack index j), type (pack-indexed), row_map (pack index -> store row).
    // These change only on a Verlet rebuild, so we copy them to the device once
    // per rebuild (keyed on the Cabana generation); positions change every step
    // and are synced on their own below.
    auto const &counts_host = verlet_list.counts;
    auto const &neighbors_host = verlet_list.neighbors;
    auto const &aosoa = cell_structure.get_aosoa();

    // Persistent, reused device buffers: no per-step device (re)allocation.
    auto &handle = cell_structure.device_sr_buffers();
    if (not handle) {
      handle = std::make_shared<DeviceShortRangeBuffers>();
    }
    auto &buf = *handle;
    buf.ensure_force(n);

    auto const generation = cell_structure.verlet_list_cabana_generation();
    if (buf.cached_generation != generation) {
      buf.ensure_list(counts_host.extent(0), neighbors_host.extent(0),
                      neighbors_host.extent(1), aosoa.type.extent(0),
                      aosoa.row_map.extent(0));
      Kokkos::deep_copy(buf.counts, counts_host);
      Kokkos::deep_copy(buf.neighbors, neighbors_host);
      Kokkos::deep_copy(buf.type, aosoa.type);
      Kokkos::deep_copy(buf.row_map, aosoa.row_map);
      buf.cached_generation = generation;
    }

    // Sync only the per-step columns the kernel reads: the authoritative
    // positions every step, charges only when Coulomb is active.
    auto &store = cell_structure.particle_store();
    Kokkos::deep_copy(store.position_view_device(), store.position_view());
    auto const position_device = store.position_view_device();
#ifdef ESPRESSO_ELECTROSTATICS
    if (coulomb_params.active) {
      Kokkos::deep_copy(store.q_view_device(), store.q_view());
    }
    auto const q_device = store.q_view_device();
#endif

    auto const counts_device = buf.counts;
    auto const neighbors_device = buf.neighbors;
    auto const type_device = buf.type;
    auto const row_map_device = buf.row_map;
    // Persistent device force accumulator (pack-indexed): re-zero each step.
    auto const force_device = buf.force;
    Kokkos::deep_copy(force_device, 0.0);

    auto const lj_local = lj_table;
    auto const n_types_local = n_types;
    auto const cutoff_sq = max_cutoff_sq;
    auto const coulomb = coulomb_params;

    Kokkos::parallel_for(
        "device_lj_nonbonded_pairs",
        Kokkos::RangePolicy<ExecSpace>(std::size_t{0}, n),
        KOKKOS_LAMBDA(std::size_t const i) {
          auto const n_neighbors = counts_device(i);
          if (n_neighbors == 0)
            return;

          auto const row_i = row_map_device(i);
          double const xi = position_device(row_i, 0);
          double const yi = position_device(row_i, 1);
          double const zi = position_device(row_i, 2);
          auto const type_i = type_device(i);
#ifdef ESPRESSO_ELECTROSTATICS
          double const qi = coulomb.active ? q_device(row_i) : 0.0;
#endif

          for (int k = 0; k < n_neighbors; ++k) {
            auto const j = neighbors_device(i, k);
            auto const row_j = row_map_device(j);

            // Minimum-image vector i - j, matching the host pair kernel fold.
            double const dx = device_mi_coord(xi, position_device(row_j, 0),
                                              box_l[0], box_l_inv[0]);
            double const dy = device_mi_coord(yi, position_device(row_j, 1),
                                              box_l[1], box_l_inv[1]);
            double const dz = device_mi_coord(zi, position_device(row_j, 2),
                                              box_l[2], box_l_inv[2]);
            double const dsq = dx * dx + dy * dy + dz * dz;
            if (dsq > cutoff_sq)
              continue;

            auto const type_j = type_device(j);
            auto const &lj =
                lj_local(static_cast<std::size_t>(type_i) *
                             static_cast<std::size_t>(n_types_local) +
                         static_cast<std::size_t>(type_j));
            double const dist = Kokkos::sqrt(dsq);
            // LJ (0 outside its shifted cutoff) plus, if active, the P3M
            // real-space Coulomb term (0 outside r_cut). Both share the
            // neighbor cutoff gate above; each self-gates on its own range, so
            // no early-out on the LJ cutoff (that would drop Coulomb pairs
            // beyond it).
            double factor = lj_force_factor_device(lj, dist);
#ifdef ESPRESSO_ELECTROSTATICS
            if (coulomb.active) {
              factor +=
                  coulomb_p3m_force_factor(qi * q_device(row_j), dist, coulomb);
            }
#endif
            if (factor == 0.0)
              continue;

            double const fx = factor * dx;
            double const fy = factor * dy;
            double const fz = factor * dz;
            // Half list: update both i and j (atomic; another thread owns j).
            Kokkos::atomic_add(&force_device(i, 0), fx);
            Kokkos::atomic_add(&force_device(i, 1), fy);
            Kokkos::atomic_add(&force_device(i, 2), fz);
            Kokkos::atomic_add(&force_device(j, 0), -fx);
            Kokkos::atomic_add(&force_device(j, 1), -fy);
            Kokkos::atomic_add(&force_device(j, 2), -fz);
          }
        });
    Kokkos::fence();

    // Write the device LJ forces into the host local_force (pack-indexed). It
    // is ZERO on entry (init_forces_and_thermostat), and bonds live in
    // scatter_force, so the reduce in forces.cpp sums the two correctly.
    Kokkos::deep_copy(cell_structure.get_local_force(), force_device);
  };
#endif // ESPRESSO_LENNARD_JONES
}
