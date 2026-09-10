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

#include "aosoa_pack.hpp"
#include "forces_inline.hpp"
#include "short_range_cabana_helpers.hpp"

#include <utils/Vector.hpp>

#include <Cabana_Core.hpp>
#include <Kokkos_ScatterView.hpp>

#include <cstddef>
#include <memory>
#include <optional>
#include <variant>
#include <vector>

#ifdef ESPRESSO_P3M
/** @brief The active P3M solver, or nullptr. Deliberately checks only the
 *  top level of the solver variant (unlike @ref get_actor_by_type, which
 *  recurses into layer corrections): under ELC the pair force must go through
 *  the generic coulomb kernel plus the ELC corrections, not the bare P3M
 *  fast path.
 */
inline CoulombP3M const *
get_toplevel_p3m_solver(Coulomb::Solver const &coulomb) {
  if (auto const &solver = coulomb.impl->solver; solver.has_value()) {
    if (std::holds_alternative<std::shared_ptr<CoulombP3M>>(*solver)) {
      return std::get<std::shared_ptr<CoulombP3M>>(*solver).get();
    }
  }
  return nullptr;
}
#endif

#ifdef ESPRESSO_ELECTROSTATICS
/** @brief Real-space charge-charge pair force: the P3M fast path when the
 *  solver is (top-level) P3M, the type-erased coulomb kernel otherwise.
 *  Shared by @ref ForcesKernel and @ref SpecializedForcesKernel so the
 *  dispatch cannot drift between them.
 */
ESPRESSO_ATTR_ALWAYS_INLINE inline Utils::Vector3d coulomb_pair_force(
    double q1q2, Utils::Vector3d const &d, double dist,
    Coulomb::ShortRangeForceKernel::kernel_type const *coulomb_kernel
#ifdef ESPRESSO_P3M
    ,
    CoulombP3M const *p3m
#endif
) {
#ifdef ESPRESSO_P3M
  if (p3m) [[likely]] {
    return p3m->pair_force(q1q2, d, dist);
  }
#endif
  return (*coulomb_kernel)(q1q2, d, dist);
}
#endif

struct ForcesKernel {
  BondedInteractionsMap const &bonded_ias;
  InteractionsNonBonded const &nonbonded_ias;
  Coulomb::ShortRangeForceKernel::kernel_type const *const coulomb_kernel;
  Dipoles::ShortRangeForceKernel::kernel_type const *const dipoles_kernel;
  Coulomb::ShortRangeForceCorrectionsKernel::kernel_type const *elc_kernel;
  Coulomb::ShortRangeEnergyKernel::kernel_type const *const coulomb_u_kernel;
  Thermostat::Thermostat const &thermostat;
  BoxGeometry const &box_geo;
  std::vector<Particle *> const &unique_particles;
  CellStructure::ScatterForce local_force;
#ifdef ESPRESSO_ROTATION
  CellStructure::ScatterForce local_torque;
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  CellStructure::ScatterForce local_dip_fld;
#endif
#ifdef ESPRESSO_NPT
  Utils::Vector3d *const global_virial;
  CellStructure::ScatterVirial local_virial;
#endif
  CellStructure::AoSoA_pack const &aosoa;
#ifdef ESPRESSO_P3M
  CoulombP3M const *p3m;
#endif
  double system_max_cutoff_sq;

  ForcesKernel(
      BondedInteractionsMap const &bonded_ias_,
      InteractionsNonBonded const &nonbonded_ias_,
      Coulomb::ShortRangeForceKernel::kernel_type const *coulomb_kernel_,
      Dipoles::ShortRangeForceKernel::kernel_type const *dipoles_kernel_,
      Coulomb::ShortRangeForceCorrectionsKernel::kernel_type const *elc_kernel_,
      Coulomb::ShortRangeEnergyKernel::kernel_type const *coulomb_u_kernel_,
      Coulomb::Solver const &coulomb_,
      Thermostat::Thermostat const &thermostat_, BoxGeometry const &box_geo_,
      std::vector<Particle *> const &unique_particles_,
      CellStructure::ScatterForce local_force_,
#ifdef ESPRESSO_ROTATION
      CellStructure::ScatterForce local_torque_,
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
      CellStructure::ScatterForce local_dip_fld_,
#endif
#ifdef ESPRESSO_NPT
      Utils::Vector3d *const global_virial_,
      CellStructure::ScatterVirial local_virial_,
#endif
      CellStructure::AoSoA_pack const &aosoa_, double system_max_cutoff_)
      : bonded_ias(bonded_ias_), nonbonded_ias(nonbonded_ias_),
        coulomb_kernel(coulomb_kernel_), dipoles_kernel(dipoles_kernel_),
        elc_kernel(elc_kernel_), coulomb_u_kernel(coulomb_u_kernel_),
        thermostat(thermostat_), box_geo(box_geo_),
        unique_particles(unique_particles_),
        local_force(std::move(local_force_)),
#ifdef ESPRESSO_ROTATION
        local_torque(std::move(local_torque_)),
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
        local_dip_fld(std::move(local_dip_fld_)),
#endif
#ifdef ESPRESSO_NPT
        global_virial(global_virial_), local_virial(std::move(local_virial_)),
#endif
        aosoa(aosoa_), system_max_cutoff_sq(Utils::sqr(system_max_cutoff_)) {
#ifdef ESPRESSO_P3M
    p3m = get_toplevel_p3m_solver(coulomb_);
#endif
  }

#ifdef ESPRESSO_NPT
  ESPRESSO_ATTR_ALWAYS_INLINE KOKKOS_INLINE_FUNCTION bool npt_active() const {
    return global_virial != nullptr;
  }
#endif

  ESPRESSO_ATTR_ALWAYS_INLINE inline void operator()(std::size_t i,
                                                     std::size_t j) const {

    // calc distance (component-wise, avoids constructing pos1/pos2 Vector3d
    // on the hot early-exit path; pos1/pos2 are built lazily below only
    // where kernels actually require them)
    auto const d = box_geo.get_mi_vector(
        aosoa.position(i, 0), aosoa.position(i, 1), aosoa.position(i, 2),
        aosoa.position(j, 0), aosoa.position(j, 1), aosoa.position(j, 2));
    auto const dist_sq = d.norm2();

    // Early exit if distance > maximal global cutoff
    if (dist_sq > system_max_cutoff_sq)
      return;
    auto const dist = std::sqrt(dist_sq);
    auto const &ia_params =
        nonbonded_ias.get_ia_param(aosoa.type(i), aosoa.type(j));

    ParticleForce pf{};

    // Determine which data needs to be loaded based on active algorithms
#if defined(ESPRESSO_EXCLUSIONS) or defined(ESPRESSO_THOLE)
    bool need_particle_pointers = false;
#ifdef ESPRESSO_EXCLUSIONS
    need_particle_pointers |= aosoa.has_exclusion(i) or aosoa.has_exclusion(j);
#endif
#ifdef ESPRESSO_THOLE
    need_particle_pointers |=
        thole_active(ia_params, coulomb_kernel != nullptr);
#endif

    Particle const *p1_ptr = nullptr;
    Particle const *p2_ptr = nullptr;
    if (need_particle_pointers) {
      p1_ptr = unique_particles.at(i);
      p2_ptr = unique_particles.at(j);
    }
#endif

    /***********************************************/
    /* non-bonded pair potentials                  */
    /***********************************************/

    if (dist <= ia_params.max_cut) {
#ifdef ESPRESSO_EXCLUSIONS
      bool skip_non_bonded = false;
      if (aosoa.has_exclusion(i) or aosoa.has_exclusion(j)) {
        skip_non_bonded = not do_nonbonded(*p1_ptr, *p2_ptr);
      }
#else
      constexpr bool skip_non_bonded = false;
#endif
      if (not skip_non_bonded) {
        pf.f += calc_central_radial_force(ia_params, d, dist);

        // Only call Thole force kernel if active
#ifdef ESPRESSO_THOLE
        if (thole_active(ia_params, coulomb_kernel != nullptr)) {
          pf.f += thole_pair_force(*p1_ptr, *p2_ptr, ia_params, d, dist,
                                   bonded_ias, coulomb_kernel);
        }
#endif
        // Only call Gay-Berne force kernel if active
#ifdef ESPRESSO_GAY_BERNE
        if (gay_berne_active(dist, ia_params)) {
          auto const dir1 = aosoa.get_vector_at(aosoa.director, i);
          auto const dir2 = aosoa.get_vector_at(aosoa.director, j);
          pf += gb_pair_force(dir1, dir2, ia_params, d, dist);
        }
#endif
#ifdef ESPRESSO_DPD
        if (dpd_active(ia_params, thermostat.thermo_switch)) {
          auto const pos1 = aosoa.get_vector_at(aosoa.position, i);
          auto const pos2 = aosoa.get_vector_at(aosoa.position, j);
          auto const vel1 = aosoa.get_vector_at(aosoa.velocity, i);
          auto const vel2 = aosoa.get_vector_at(aosoa.velocity, j);
          auto const force = dpd_pair_force(
              pos1, vel1, aosoa.id(i), pos2, vel2, aosoa.id(j), *thermostat.dpd,
              box_geo, ia_params, d, dist, dist_sq);
          pf += force;
        }
#endif // ESPRESSO_DPD

      } // not skip_non_bonded
    } // not dist > ia_params.max_cut

    /*********************************************************************/
    /* everything before this contributes to the virial pressure in NpT  */
    /* via d (x) pf.f; electrostatic and dipolar real-space contributions */
    /* are added in explicitly below instead: Coulomb reuses the pair    */
    /* energy as a virial proxy, dipoles compute d . F directly (see     */
    /* rationale below)                                                 */
    /*********************************************************************/
#ifdef ESPRESSO_NPT
    Utils::Vector3d virial{};
    if (npt_active()) {
      virial = hadamard_product(pf.f, d);
    }
#endif // ESPRESSO_NPT

#ifdef ESPRESSO_ELECTROSTATICS
    Utils::Vector3d f1_asym{};
    Utils::Vector3d f2_asym{};
    // real-space electrostatic charge-charge interaction
    if (coulomb_kernel != nullptr) {
      if ((aosoa.charge(i) != 0.) and (aosoa.charge(j) != 0.)) {
        auto const q1q2 = aosoa.charge(i) * aosoa.charge(j);
        pf.f += coulomb_pair_force(q1q2, d, dist, coulomb_kernel
#ifdef ESPRESSO_P3M
                                   ,
                                   p3m
#endif
        );
        if (elc_kernel) {
          auto const pos1 = aosoa.get_vector_at(aosoa.position, i);
          auto const pos2 = aosoa.get_vector_at(aosoa.position, j);
          (*elc_kernel)(pos1, pos2, f1_asym, f2_asym, q1q2);
        }
#ifdef ESPRESSO_NPT
        if (npt_active()) {
          auto const pos1 = aosoa.get_vector_at(aosoa.position, i);
          auto const pos2 = aosoa.get_vector_at(aosoa.position, j);
          virial[0] += (*coulomb_u_kernel)(pos1, pos2, q1q2, d, dist);
        }
#endif // ESPRESSO_NPT
      }
    }
#endif // ESPRESSO_ELECTROSTATICS

#ifdef ESPRESSO_DIPOLES
    if (dipoles_kernel != nullptr) {
      auto const d1d2 = aosoa.dipm(i) * aosoa.dipm(j);
      if (d1d2 != 0.) {
        auto const dir1 = aosoa.get_vector_at(aosoa.director, i);
        auto const dir2 = aosoa.get_vector_at(aosoa.director, j);
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
        Utils::Vector3d dip_fld_i{};
        Utils::Vector3d dip_fld_j{};
#endif
        auto const dip_pf =
            (*dipoles_kernel)(d1d2, aosoa.dipm(i) * dir1, aosoa.dipm(j) * dir2,
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
                              dip_fld_i, dip_fld_j,
#endif
                              d, dist, dist_sq);
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
        auto access_dip_fld = local_dip_fld.access();
        access_dip_fld(i, 0) += dip_fld_i[0];
        access_dip_fld(i, 1) += dip_fld_i[1];
        access_dip_fld(i, 2) += dip_fld_i[2];
        access_dip_fld(j, 0) += dip_fld_j[0];
        access_dip_fld(j, 1) += dip_fld_j[1];
        access_dip_fld(j, 2) += dip_fld_j[2];
#endif
#ifdef ESPRESSO_NPT
        if (npt_active()) {
          // d . F = -n * U for a homogeneous potential of degree n
          // (Euler's theorem, independent of centrality); n=-3 here vs
          // n=-1 for Coulomb. Ewald screening makes that only
          // approximate, and for dipoles the approximation measurably
          // fails NpT pressure consistency (see test_pressure_with_dp3m),
          // so d . F is computed explicitly here instead.
          virial[0] += d * dip_pf.f;
        }
#endif // ESPRESSO_NPT
        pf += dip_pf;
      }
    }
#endif // ESPRESSO_DIPOLES

    auto opf = calc_opposing_force(pf, d);
#ifdef ESPRESSO_ELECTROSTATICS
    pf.f += f1_asym;
    opf.f += f2_asym;
#endif // ESPRESSO_ELECTROSTATICS

    auto access_force = local_force.access();

    access_force(i, 0) += pf.f[0];
    access_force(i, 1) += pf.f[1];
    access_force(i, 2) += pf.f[2];
#ifdef ESPRESSO_ROTATION
    auto access_torque = local_torque.access();
    access_torque(i, 0) += pf.torque[0];
    access_torque(i, 1) += pf.torque[1];
    access_torque(i, 2) += pf.torque[2];
#endif

    access_force(j, 0) += opf.f[0];
    access_force(j, 1) += opf.f[1];
    access_force(j, 2) += opf.f[2];
#ifdef ESPRESSO_ROTATION
    access_torque(j, 0) += opf.torque[0];
    access_torque(j, 1) += opf.torque[1];
    access_torque(j, 2) += opf.torque[2];
#endif
#ifdef ESPRESSO_NPT
    if (npt_active()) {
      auto access_virial = local_virial.access();
      access_virial(0) += virial[0];
      access_virial(1) += virial[1];
      access_virial(2) += virial[2];
    }
#endif
  }
};

/** @brief Pair potentials fully handled by @ref SpecializedForcesKernel,
 *  i.e. exactly the central-radial family dispatched through
 *  @ref calc_central_radial_force. Potentials with their own branch in
 *  @ref ForcesKernel (Gay-Berne, DPD) are deliberately absent, and so is any
 *  newly added potential until it is proven compatible -- the dispatch gate
 *  in forces.cpp rejects any pair mask with a bit outside this allowlist, so
 *  new potentials fall back to the generic kernel by default.
 */
constexpr unsigned specialized_kernel_pair_mask =
    pair_potential_bit(PairPotential::LennardJones) |
    pair_potential_bit(PairPotential::WCA) |
    pair_potential_bit(PairPotential::LennardJonesGeneric) |
    pair_potential_bit(PairPotential::SmoothStep) |
    pair_potential_bit(PairPotential::Hertzian) |
    pair_potential_bit(PairPotential::Gaussian) |
    pair_potential_bit(PairPotential::BMHTF) |
    pair_potential_bit(PairPotential::Buckingham) |
    pair_potential_bit(PairPotential::Morse) |
    pair_potential_bit(PairPotential::SoftSphere) |
    pair_potential_bit(PairPotential::Hat) |
    pair_potential_bit(PairPotential::LJCos) |
    pair_potential_bit(PairPotential::LJCos2) |
    pair_potential_bit(PairPotential::Tabulated);

/**
 * @brief Own-the-loop specialization of the non-bonded pair kernel.
 *
 * Handles the common case of a cuboid box with only central radial pair
 * forces active, optionally with real-space electrostatics (@p HasCoulomb).
 * The dispatch in forces.cpp selects it only when no torque-producing,
 * asymmetric, thermostat or exclusion interaction is active (no dipoles, ELC,
 * DPD, NPT virial, exclusions, Thole or Gay-Berne, no Lees-Edwards), and
 * falls back to @ref ForcesKernel otherwise.
 *
 * Unlike @ref ForcesKernel (invoked once per pair by
 * Cabana::neighbor_parallel_for), this owns the Verlet-list loop and runs once
 * per particle. That lets it hoist the per-particle position, type and charge,
 * capture the cuboid box parameters by value, and obtain the ScatterView
 * accessor once per particle instead of once per pair -- the per-pair
 * `access()` is an omp_get_thread_num call that dominates the LJ pair cost. It
 * also accumulates the i-side force in a local register over all of particle
 * i's neighbors and updates the ScatterView once per particle instead of once
 * per pair; the j-side write stays per-pair (each neighbor j is distinct). This
 * removes the per-pair i-side accessor update at the cost of changing the
 * i-side summation order (a single register sum added once, rather than
 * incremental scatter adds), so the result matches @ref ForcesKernel to
 * floating-point round-off rather than bitwise. `if constexpr (HasCoulomb)`
 * compiles the electrostatics path in or out entirely.
 *
 * The neighbors are processed in fixed-size tiles, each in three passes: a
 * scalar gather of the neighbor positions into SoA scratch, a vectorized
 * minimum-image pass (@ref CuboidMinimumImage::batch_vector_dist2) that folds
 * the whole tile at once and squares the distances, and a scalar pass that
 * runs the short-range force only for the pairs whose squared distance passes
 * the cutoff gate. The vectorized pass yields the same per-pair fold vector and
 * squared distance as the scalar path, so identity is preserved.
 */
template <bool HasCoulomb> struct SpecializedForcesKernel {
  static constexpr int tile_size = 64;

  InteractionsNonBonded const &nonbonded_ias;
  CellStructure::AoSoA_pack const &aosoa;
  CellStructure::ScatterForce local_force;
  Kokkos::View<int const *, Kokkos::HostSpace> counts;
  Kokkos::View<int const **, Kokkos::LayoutRight, Kokkos::HostSpace> neighbors;
  CuboidMinimumImage minimum_image;
  double system_max_cutoff_sq;
#ifdef ESPRESSO_ELECTROSTATICS
  Coulomb::ShortRangeForceKernel::kernel_type const *coulomb_kernel;
#ifdef ESPRESSO_P3M
  CoulombP3M const *p3m;
#endif
#endif

  ESPRESSO_ATTR_ALWAYS_INLINE inline void
  operator()(std::size_t const i) const {
    auto const n_neighbors = counts(i);
    if (n_neighbors == 0)
      return;

    auto const x_i = aosoa.position(i, 0);
    auto const y_i = aosoa.position(i, 1);
    auto const z_i = aosoa.position(i, 2);
    auto const type_i = aosoa.type(i);
#ifdef ESPRESSO_ELECTROSTATICS
    double charge_i = 0.;
    if constexpr (HasCoulomb) {
      charge_i = aosoa.charge(i);
    }
#endif

    // One ScatterView accessor for all of this particle's pairs.
    auto access_force = local_force.access();

    // Accumulate the i-side force over all neighbors in a local register and
    // write it to the ScatterView once at the end (see class docs).
    Utils::Vector3d f_i{};

    // Per-tile SoA scratch (thread-local, on the stack).
    int js[tile_size];
    double sx[tile_size], sy[tile_size], sz[tile_size];
    double dx0[tile_size], dx1[tile_size], dx2[tile_size], dsq[tile_size];

    for (int base = 0; base < n_neighbors; base += tile_size) {
      // ``+tile_size`` creates a prvalue and avoids an ODR-use of a host-space
      // variable from device code (Kokkos::min() takes arguments by const &T)
      auto const m = Kokkos::min(+tile_size, n_neighbors - base);

      // Pass 1: scalar gather of the tile's neighbor positions.
      for (int t = 0; t < m; ++t) {
        auto const j = neighbors(i, base + t);
        js[t] = j;
        auto const row_j = static_cast<std::size_t>(j);
        sx[t] = aosoa.position(row_j, 0);
        sy[t] = aosoa.position(row_j, 1);
        sz[t] = aosoa.position(row_j, 2);
      }

      // Pass 2: vectorized minimum-image fold + squared distance.
      minimum_image.batch_vector_dist2(x_i, y_i, z_i, m, sx, sy, sz, dx0, dx1,
                                       dx2, dsq);

      // Pass 3: scalar short-range force for the pairs that pass the gate,
      // in neighbor order (i before j) to preserve the accumulation order.
      for (int t = 0; t < m; ++t) {
        if (dsq[t] > system_max_cutoff_sq)
          continue;
        auto const j = static_cast<std::size_t>(js[t]);
        Utils::Vector3d const d{dx0[t], dx1[t], dx2[t]};
        auto const dist = std::sqrt(dsq[t]);
        auto const &ia_params =
            nonbonded_ias.get_ia_param(type_i, aosoa.type(j));

        Utils::Vector3d f{};
        if (dist <= ia_params.max_cut) {
          f += calc_central_radial_force(ia_params, d, dist);
        }
#ifdef ESPRESSO_ELECTROSTATICS
        if constexpr (HasCoulomb) {
          auto const charge_j = aosoa.charge(j);
          if (charge_i != 0. and charge_j != 0.) {
            auto const q1q2 = charge_i * charge_j;
            f += coulomb_pair_force(q1q2, d, dist, coulomb_kernel
#ifdef ESPRESSO_P3M
                                    ,
                                    p3m
#endif
            );
          }
        }
#endif

        f_i += f;
        access_force(j, 0) -= f[0];
        access_force(j, 1) -= f[1];
        access_force(j, 2) -= f[2];
      }
    }

    // Single i-side ScatterView update for the whole neighbor loop.
    access_force(i, 0) += f_i[0];
    access_force(i, 1) += f_i[1];
    access_force(i, 2) += f_i[2];
  }
};
