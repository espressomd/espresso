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

#include "BoxGeometry.hpp"
#include "Observable_stat.hpp"
#include "Particle.hpp"
#include "aosoa_pack.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "electrostatics/coulomb.hpp"
#include "exclusions.hpp"
#include "forces_inline.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "pressure_inline.hpp"
#include "short_range_cabana_helpers.hpp"

#ifdef ESPRESSO_DPD
#include "dpd.hpp"
#endif

#include <utils/Vector.hpp>
#include <utils/math/tensor_product.hpp>

#include <Cabana_Core.hpp>

#include <omp.h>

#include <cstddef>
#include <span>
#include <vector>

struct PressureBinLayout {
  std::size_t n_bonded;
  std::size_t n_types;
  std::size_t off_bonded = 0;
  std::size_t off_nb_inter;
  std::size_t off_nb_intra;
  std::size_t off_coulomb;
  std::size_t off_dipolar;
  std::size_t off_dpd;
  std::size_t total;

  PressureBinLayout(std::size_t n_bonded_, std::size_t n_types_)
      : n_bonded(n_bonded_), n_types(n_types_) {
    auto const n_nb = n_types * (n_types + 1) / 2;
    off_nb_inter = off_bonded + n_bonded;
    off_nb_intra = off_nb_inter + n_nb;
    off_coulomb = off_nb_intra + n_nb;
    off_dipolar = off_coulomb + 1;
    off_dpd = off_dipolar + 1;
    total = off_dpd + 1;
  }

  KOKKOS_INLINE_FUNCTION
  std::size_t tensor_offset(std::size_t bin, std::size_t k) const {
    return bin * 9 + k;
  }

  KOKKOS_INLINE_FUNCTION
  std::size_t nb_inter_idx(int t1, int t2) const {
    return off_nb_inter +
           Utils::lower_triangular(std::max(t1, t2), std::min(t1, t2));
  }

  KOKKOS_INLINE_FUNCTION
  std::size_t nb_intra_idx(int t1, int t2) const {
    return off_nb_intra +
           Utils::lower_triangular(std::max(t1, t2), std::min(t1, t2));
  }

  KOKKOS_INLINE_FUNCTION std::size_t coulomb_idx() const { return off_coulomb; }
  KOKKOS_INLINE_FUNCTION std::size_t dipolar_idx() const { return off_dipolar; }
  KOKKOS_INLINE_FUNCTION std::size_t dpd_idx() const { return off_dpd; }
  KOKKOS_INLINE_FUNCTION std::size_t bonded_idx(int b) const {
    return off_bonded + static_cast<std::size_t>(b);
  }
};

struct PressureKernel {
  BondedInteractionsMap const &bonded_ias;
  InteractionsNonBonded const &nonbonded_ias;
  Coulomb::Solver const &coulomb;
  Coulomb::ShortRangeForceKernel::kernel_type const *coulomb_f_kernel;
  Coulomb::ShortRangePressureKernel::kernel_type const *coulomb_p_kernel;
  BoxGeometry const &box_geo;
#ifdef ESPRESSO_DPD
  DPDThermostat const *dpd;
#endif
  std::vector<Particle *> const &unique_particles;
  Kokkos::View<double **, Kokkos::LayoutRight> local_pressure;
  PressureBinLayout layout;
  CellStructure::AoSoA_pack const &aosoa;
  Kokkos::View<int *> mol_id_view;
  double system_max_cutoff;
  int thermo_switch;

  PressureKernel(
      BondedInteractionsMap const &bonded_ias_,
      InteractionsNonBonded const &nonbonded_ias_,
      Coulomb::Solver const &coulomb_,
      Coulomb::ShortRangeForceKernel::kernel_type const *coulomb_f_kernel_,
      Coulomb::ShortRangePressureKernel::kernel_type const *coulomb_p_kernel_,
      BoxGeometry const &box_geo_,
#ifdef ESPRESSO_DPD
      DPDThermostat const *dpd_,
#endif
      std::vector<Particle *> const &unique_particles_,
      Kokkos::View<double **, Kokkos::LayoutRight> const &local_pressure_,
      PressureBinLayout layout_, CellStructure::AoSoA_pack const &aosoa_,
      Kokkos::View<int *> mol_id_view_, double system_max_cutoff_,
      int thermo_switch_)
      : bonded_ias(bonded_ias_), nonbonded_ias(nonbonded_ias_),
        coulomb(coulomb_), coulomb_f_kernel(coulomb_f_kernel_),
        coulomb_p_kernel(coulomb_p_kernel_), box_geo(box_geo_),
#ifdef ESPRESSO_DPD
        dpd(dpd_),
#endif
        unique_particles(unique_particles_), local_pressure(local_pressure_),
        layout(layout_), aosoa(aosoa_), mol_id_view(std::move(mol_id_view_)),
        system_max_cutoff(system_max_cutoff_), thermo_switch(thermo_switch_) {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(std::size_t i, std::size_t j) const {
    auto const pos1 = aosoa.get_vector_at(aosoa.position, i);
    auto const pos2 = aosoa.get_vector_at(aosoa.position, j);
    auto const d = box_geo.get_mi_vector(pos1, pos2);
    auto const dist = d.norm();
    if (dist > system_max_cutoff)
      return;

    auto const t1 = aosoa.type(i);
    auto const t2 = aosoa.type(j);
    auto const &ia_params = nonbonded_ias.get_ia_param(t1, t2);
    auto const tid = omp_get_thread_num();

    // --- Non-bonded virial ---
    if (dist <= ia_params.max_cut) {
      bool skip = false;
#ifdef ESPRESSO_EXCLUSIONS
      if (aosoa.has_exclusion(i) or aosoa.has_exclusion(j)) {
        skip =
            not do_nonbonded(*unique_particles.at(i), *unique_particles.at(j));
      }
#endif
      if (not skip) {
        Utils::Vector3d f = calc_central_radial_force(ia_params, d, dist);

#ifdef ESPRESSO_GAY_BERNE
        if (gay_berne_active(dist, ia_params)) {
          auto const dir1 = aosoa.get_vector_at(aosoa.director, i);
          auto const dir2 = aosoa.get_vector_at(aosoa.director, j);
          f += calc_non_central_force(dir1, dir2, ia_params, d, dist);
        }
#endif
#ifdef ESPRESSO_THOLE
        if (thole_active(ia_params, coulomb_f_kernel != nullptr)) {
          f += thole_pair_force(*unique_particles.at(i),
                                *unique_particles.at(j), ia_params, d, dist,
                                bonded_ias, coulomb_f_kernel);
        }
#endif

        auto const stress = Utils::flatten(Utils::tensor_product(d, f));
        auto const bin = (mol_id_view(i) == mol_id_view(j))
                             ? layout.nb_intra_idx(t1, t2)
                             : layout.nb_inter_idx(t1, t2);
        for (std::size_t k = 0; k < 9; ++k)
          local_pressure(tid, layout.tensor_offset(bin, k)) += stress[k];

#ifdef ESPRESSO_DPD
        if (dpd_active(ia_params, thermo_switch)) {
          auto const pid1 = aosoa.id(i);
          auto const pid2 = aosoa.id(j);
          auto const noise_vec = (ia_params.dpd.radial.pref > 0.0 ||
                                  ia_params.dpd.trans.pref > 0.0)
                                     ? dpd_noise(*dpd, pid1, pid2)
                                     : Utils::Vector3d{};
          auto const vel1 = aosoa.get_vector_at(aosoa.velocity, i);
          auto const vel2 = aosoa.get_vector_at(aosoa.velocity, j);
          auto const v21 = box_geo.velocity_difference(pos1, pos2, vel1, vel2);
          auto const dist2 = d.norm2();
          // f_r/f_t: dissipative force from radial/transverse DPD channel
          auto const f_r =
              dpd_pair_force(ia_params.dpd.radial, v21, dist, noise_vec);
          auto const f_t =
              dpd_pair_force(ia_params.dpd.trans, v21, dist, noise_vec);
          auto const P = Utils::tensor_product(d / dist2, d);
          auto const f_d = P * (f_r - f_t) + f_t;
          auto const s = Utils::flatten(Utils::tensor_product(d, f_d));
          for (std::size_t k = 0; k < 9; ++k)
            local_pressure(tid, layout.tensor_offset(layout.dpd_idx(), k)) +=
                s[k];
        }
#endif
      }
    }

#ifdef ESPRESSO_ELECTROSTATICS
    if (coulomb_p_kernel != nullptr) {
      auto const q1 = aosoa.charge(i), q2 = aosoa.charge(j);
      if (q1 != 0. and q2 != 0.) {
        auto const p_c = Utils::flatten((*coulomb_p_kernel)(q1 * q2, d, dist));
        for (std::size_t k = 0; k < 9; ++k)
          local_pressure(tid, layout.tensor_offset(layout.coulomb_idx(), k)) +=
              p_c[k];
      }
    }
#endif
  }
};

static void reduce_cabana_pressure(
    Kokkos::View<double **, Kokkos::LayoutRight> const &local_pressure,
    PressureBinLayout const &layout, Observable_stat &obs,
    [[maybe_unused]] BondedInteractionsMap const &bonded_ias, int n_types) {
  auto const nthreads = local_pressure.extent(0);
  auto host =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, local_pressure);

  auto sum_tensor = [&](std::size_t bin, std::span<double> dest) {
    for (std::size_t k = 0; k < 9; ++k) {
      double s = 0.;
      for (std::size_t t = 0; t < nthreads; ++t)
        s += host(t, layout.tensor_offset(bin, k));
      dest[k] += s;
    }
  };

  for (int b = 0; b < int(layout.n_bonded); ++b)
    sum_tensor(layout.bonded_idx(b), obs.bonded_contribution(b));

  for (int t1 = 0; t1 < n_types; ++t1)
    for (int t2 = 0; t2 <= t1; ++t2)
      sum_tensor(layout.nb_inter_idx(t1, t2),
                 obs.non_bonded_inter_contribution(t1, t2));

  for (int t1 = 0; t1 < n_types; ++t1)
    for (int t2 = 0; t2 <= t1; ++t2)
      sum_tensor(layout.nb_intra_idx(t1, t2),
                 obs.non_bonded_intra_contribution(t1, t2));

  if (!obs.coulomb.empty())
    sum_tensor(layout.coulomb_idx(), obs.coulomb.subspan(0, 9));
  if (!obs.dipolar.empty())
    sum_tensor(layout.dipolar_idx(), obs.dipolar.subspan(0, 9));
  if (!obs.dpd.empty())
    sum_tensor(layout.dpd_idx(), obs.dpd);
}
