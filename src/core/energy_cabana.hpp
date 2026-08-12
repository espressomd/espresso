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

#pragma once

#include <config/config.hpp>

#include "aosoa_pack.hpp"
#include "energy_inline.hpp"
#include "short_range_cabana_helpers.hpp"

#include <utils/Vector.hpp>

#include <Cabana_Core.hpp>

#include <omp.h>

#include <cstddef>
#include <memory>
#include <optional>
#include <variant>
#include <vector>

struct EnergyBinLayout {
  std::size_t n_bonded;       // initialized by bonded_ias->get_next_key()
  std::size_t n_types;        // max_seen_particle_type
  std::size_t off_bonded = 0; // [0, n_bonded)
  std::size_t off_nb_inter;
  std::size_t off_nb_intra;
  std::size_t off_coulomb;
  std::size_t off_dipolar;
  std::size_t total;

  EnergyBinLayout(std::size_t n_bonded_, std::size_t n_types_)
      : n_bonded(n_bonded_), n_types(n_types_) {
    auto const n_nb = n_types * (n_types + 1) / 2;
    off_nb_inter = off_bonded + n_bonded;
    off_nb_intra = off_nb_inter + n_nb;
    off_coulomb = off_nb_intra + n_nb;
    off_dipolar = off_coulomb + 1;
    total = off_dipolar + 1;
  }

  ESPRESSO_ATTR_ALWAYS_INLINE KOKKOS_INLINE_FUNCTION std::size_t
  nb_inter_idx(int t1, int t2) const {
    return off_nb_inter +
           Utils::lower_triangular(std::max(t1, t2), std::min(t1, t2));
  }

  ESPRESSO_ATTR_ALWAYS_INLINE KOKKOS_INLINE_FUNCTION std::size_t
  nb_intra_idx(int t1, int t2) const {
    return off_nb_intra +
           Utils::lower_triangular(std::max(t1, t2), std::min(t1, t2));
  }

  KOKKOS_INLINE_FUNCTION std::size_t dipolar_idx() const { return off_dipolar; }
  KOKKOS_INLINE_FUNCTION std::size_t coulomb_idx() const { return off_coulomb; }
  KOKKOS_INLINE_FUNCTION std::size_t bonded_idx(int b) const {
    return off_bonded + b;
  }
};

struct EnergyKernel {
  using memory_space = Kokkos::HostSpace;
  BondedInteractionsMap const &bonded_ias;
  InteractionsNonBonded const &nonbonded_ias;
  Coulomb::Solver const &coulomb;
  Coulomb::ShortRangeEnergyKernel::kernel_type const *coulomb_u_kernel;
  Dipoles::ShortRangeEnergyKernel::kernel_type const *dipoles_u_kernel;
  BoxGeometry const &box_geo;
  std::vector<Particle *> const &unique_particles;
  Kokkos::View<double **, Kokkos::LayoutRight, memory_space> local_energy;
  EnergyBinLayout layout;
  CellStructure::AoSoA_pack const &aosoa;
  Kokkos::View<int *, Kokkos::LayoutRight, memory_space> mol_id_view;
  double system_max_cutoff_sq;

  EnergyKernel(
      BondedInteractionsMap const &bonded_ias_,
      InteractionsNonBonded const &nonbonded_ias_,
      Coulomb::Solver const &coulomb_,
      Coulomb::ShortRangeEnergyKernel::kernel_type const *coulomb_u_kernel_,
      Dipoles::ShortRangeEnergyKernel::kernel_type const *dipoles_u_kernel_,
      BoxGeometry const &box_geo_,
      std::vector<Particle *> const &unique_particles_,
      Kokkos::View<double **, Kokkos::LayoutRight, memory_space> const
          &local_energy_,
      EnergyBinLayout layout_, CellStructure::AoSoA_pack const &aosoa_,
      Kokkos::View<int *, memory_space> mol_id_view_, double system_max_cutoff_)
      : bonded_ias(bonded_ias_), nonbonded_ias(nonbonded_ias_),
        coulomb(coulomb_), coulomb_u_kernel(coulomb_u_kernel_),
        dipoles_u_kernel(dipoles_u_kernel_), box_geo(box_geo_),
        unique_particles(unique_particles_), local_energy(local_energy_),
        layout(layout_), aosoa(aosoa_), mol_id_view(std::move(mol_id_view_)),
        system_max_cutoff_sq(system_max_cutoff_ * system_max_cutoff_) {}

  ESPRESSO_ATTR_ALWAYS_INLINE inline void operator()(std::size_t i,
                                                     std::size_t j) const {
    auto const d = box_geo.get_mi_vector(
        aosoa.position(i, 0), aosoa.position(i, 1), aosoa.position(i, 2),
        aosoa.position(j, 0), aosoa.position(j, 1), aosoa.position(j, 2));
    auto const dist_sq = d.norm2();
    if (dist_sq > system_max_cutoff_sq)
      return;
    auto const dist = std::sqrt(dist_sq);

    auto const t1 = aosoa.type(i);
    auto const t2 = aosoa.type(j);
    auto const &ia_params = nonbonded_ias.get_ia_param(t1, t2);

    // Determine which data needs to be loaded based on active algorithms
#if defined(ESPRESSO_EXCLUSIONS) or defined(ESPRESSO_THOLE)
    bool need_particle_pointers = false;
#ifdef ESPRESSO_EXCLUSIONS
    need_particle_pointers |= aosoa.has_exclusion(i) or aosoa.has_exclusion(j);
#endif
#ifdef ESPRESSO_THOLE
    need_particle_pointers |=
        thole_active(ia_params, coulomb_u_kernel != nullptr);
#endif

    Particle const *p1_ptr = nullptr;
    Particle const *p2_ptr = nullptr;
    if (need_particle_pointers) {
      p1_ptr = unique_particles.at(i);
      p2_ptr = unique_particles.at(j);
    }
#endif

    auto const tid = omp_get_thread_num();
    double e_nb = 0.0;

    if (dist <= ia_params.max_cut) {
#ifdef ESPRESSO_EXCLUSIONS
      bool skip = false;
      if (aosoa.has_exclusion(i) or aosoa.has_exclusion(j))
        skip = not do_nonbonded(*p1_ptr, *p2_ptr);
      if (not skip)
#endif
      {
        e_nb += calc_central_radial_energy(ia_params, dist);

        // Only call Thole energy kernel if active
#ifdef ESPRESSO_THOLE
        if (thole_active(ia_params, coulomb_u_kernel != nullptr)) {
          e_nb += thole_pair_energy(*p1_ptr, *p2_ptr, ia_params, d, dist,
                                    bonded_ias, coulomb, coulomb_u_kernel);
        }
#endif
        // Only call Gay-Berne energy kernel if active
#ifdef ESPRESSO_GAY_BERNE
        if (gay_berne_active(dist, ia_params)) {
          auto const dir1 = aosoa.get_vector_at(aosoa.director, i);
          auto const dir2 = aosoa.get_vector_at(aosoa.director, j);
          e_nb += gb_pair_energy(dir1, dir2, ia_params, d, dist);
        }
#endif
      }
    }
    // pick inter vs intra bin like Observable_stat::non_bonded_contribution
    // does
    auto const bin = (mol_id_view(i) == mol_id_view(j))
                         ? layout.nb_intra_idx(t1, t2)
                         : layout.nb_inter_idx(t1, t2);
    local_energy(tid, bin) += e_nb;

#ifdef ESPRESSO_ELECTROSTATICS
    if (coulomb_u_kernel != nullptr) {
      auto const q1 = aosoa.charge(i), q2 = aosoa.charge(j);
      if (q1 != 0. and q2 != 0.) {
        auto const pos1 = aosoa.get_vector_at(aosoa.position, i);
        auto const pos2 = aosoa.get_vector_at(aosoa.position, j);
        double const e_c = (*coulomb_u_kernel)(pos1, pos2, q1 * q2, d, dist);
        local_energy(tid, layout.coulomb_idx()) += e_c;
      }
    }
#endif

#ifdef ESPRESSO_DIPOLES
    if (dipoles_u_kernel != nullptr) {
      if (aosoa.dipm(i) != 0. and aosoa.dipm(j) != 0.) {
        auto const dir1 = aosoa.get_vector_at(aosoa.director, i);
        auto const dir2 = aosoa.get_vector_at(aosoa.director, j);
        double const e_d = (*dipoles_u_kernel)(
            aosoa.dipm(i) * dir1, aosoa.dipm(j) * dir2, d, dist, dist_sq);
        local_energy(tid, layout.dipolar_idx()) += e_d;
      }
    }
#endif
  }
};

static void reduce_cabana_energy(
    Kokkos::View<double **, Kokkos::LayoutRight, Kokkos::HostSpace> const
        &local_energy,
    EnergyBinLayout const &layout, Observable_stat &obs,
    BondedInteractionsMap const &bonded_ias, int n_types) {
  auto const nthreads = local_energy.extent(0);
  auto host =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, local_energy);

  auto sum_bin = [&](std::size_t bin) {
    double s = 0.;
    for (std::size_t t = 0; t < nthreads; ++t)
      s += host(t, bin);
    return s;
  };

  for (int b = 0; b < int(layout.n_bonded); ++b)
    obs.bonded_contribution(b)[0] += sum_bin(layout.bonded_idx(b));

  for (int t1 = 0; t1 < n_types; ++t1)
    for (int t2 = 0; t2 <= t1; ++t2)
      obs.non_bonded_inter_contribution(t1, t2)[0] +=
          sum_bin(layout.nb_inter_idx(t1, t2));

  for (int t1 = 0; t1 < n_types; ++t1)
    for (int t2 = 0; t2 <= t1; ++t2)
      obs.non_bonded_intra_contribution(t1, t2)[0] +=
          sum_bin(layout.nb_intra_idx(t1, t2));

  obs.coulomb[0] += sum_bin(layout.coulomb_idx());
  obs.dipolar[0] += sum_bin(layout.dipolar_idx());
}
