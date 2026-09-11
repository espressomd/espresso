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
#include "bond_error.hpp"
#include "cell_system/LocalBondState.hpp"
#include "energy_cabana.hpp"
#include "energy_inline.hpp"

#include <utils/Vector.hpp>

#include <Kokkos_Core.hpp>

#include <omp.h>

#include <cstddef>
#include <optional>
#include <variant>

struct BondsEnergyKernelData {
  using execution_space = Kokkos::HostSpace;
  BondedInteractionsMap const &bonded_ias;
  BoxGeometry const &box_geo;
  Kokkos::View<double **, Kokkos::LayoutRight, execution_space> local_energy;
  EnergyBinLayout layout;
  CellStructure::AoSoA_pack const &aosoa;
};

struct PairBondsEnergyKernel {
  BondsEnergyKernelData data;
  LocalBondState::PairBondlistType bond_list;
  LocalBondState::PairBondIDType bond_ids;
  Coulomb::ShortRangeEnergyKernel::kernel_type const *const coulomb_u_kernel;

  PairBondsEnergyKernel(
      BondsEnergyKernelData data_, LocalBondState::PairBondlistType bond_list_,
      LocalBondState::PairBondIDType bond_ids_,
      Coulomb::ShortRangeEnergyKernel::kernel_type const *coulomb_u_kernel_)
      : data(std::move(data_)), bond_list(std::move(bond_list_)),
        bond_ids(std::move(bond_ids_)), coulomb_u_kernel(coulomb_u_kernel_) {}

  ESPRESSO_ATTR_ALWAYS_INLINE inline void operator()(std::size_t idx) const {
    auto const &bonded_ias = data.bonded_ias;
    auto const &box_geo = data.box_geo;
    auto &local_energy = data.local_energy;
    auto const &layout = data.layout;
    auto const &aosoa = data.aosoa;
    auto const bond_id = bond_ids(idx);

    // TODO: omp_get_thread_num() is only available for the OpenMP backend.
    // This should be updated when using other Kokkos backends.
    auto const thread_id = omp_get_thread_num();

    auto const i = bond_list(idx, 0);
    auto const j = bond_list(idx, 1);
    // Translate pack indices to ParticleStore rows once.
    auto const row_i = aosoa.row(i);
    auto const row_j = aosoa.row(j);
    auto const &iaparams = *bonded_ias.at(bond_id);

    auto const pos1 = aosoa.get_vector_at(aosoa.position, row_i);
    auto const pos2 = aosoa.get_vector_at(aosoa.position, row_j);
    auto const dx = box_geo.get_mi_vector(pos1, pos2);

    std::optional<double> energy = calc_pair_bonded_energy(
        iaparams, dx, pos1, pos2,
#ifdef ESPRESSO_ELECTROSTATICS
        // charge aliases the store column; read by *store row*.
        aosoa.charge(row_i) * aosoa.charge(row_j), coulomb_u_kernel
#else
        0.0, nullptr
#endif
    );

    if (energy) {
      local_energy(thread_id, layout.bonded_idx(bond_id)) += energy.value();
    } else {
      auto partner_id = aosoa.id(row_j);
      bond_broken_error(aosoa.id(row_i), {&partner_id, 1});
    }
  }
};

struct AngleBondsEnergyKernel {
  BondsEnergyKernelData data;
  LocalBondState::AngleBondlistType bond_list;
  LocalBondState::AngleBondIDType bond_ids;

  AngleBondsEnergyKernel(BondsEnergyKernelData data_,
                         LocalBondState::AngleBondlistType bond_list_,
                         LocalBondState::AngleBondIDType bond_ids_)
      : data(std::move(data_)), bond_list(std::move(bond_list_)),
        bond_ids(std::move(bond_ids_)) {}

  ESPRESSO_ATTR_ALWAYS_INLINE inline void operator()(std::size_t idx) const {
    auto const &bonded_ias = data.bonded_ias;
    auto const &box_geo = data.box_geo;
    auto &local_energy = data.local_energy;
    auto const &layout = data.layout;
    auto const &aosoa = data.aosoa;
    auto const bond_id = bond_ids(idx);

    // TODO: omp_get_thread_num() is only available for the OpenMP backend.
    // This should be updated when using other Kokkos backends.
    auto const thread_id = omp_get_thread_num();

    auto const i = bond_list(idx, 0);
    auto const j = bond_list(idx, 1);
    auto const k = bond_list(idx, 2);
    // Translate pack indices to ParticleStore rows once.
    auto const row_i = aosoa.row(i);
    auto const row_j = aosoa.row(j);
    auto const row_k = aosoa.row(k);
    auto const &iaparams = *bonded_ias.at(bond_id);

    auto const pos1 = aosoa.get_vector_at(aosoa.position, row_i);
    auto const pos2 = aosoa.get_vector_at(aosoa.position, row_j);
    auto const pos3 = aosoa.get_vector_at(aosoa.position, row_k);
    auto const vec1 = box_geo.get_mi_vector(pos2, pos1);
    auto const vec2 = box_geo.get_mi_vector(pos3, pos1);

    std::optional<double> energy =
        calc_angle_bonded_energy(iaparams, vec1, vec2);

    if (energy) {
      local_energy(thread_id, layout.bonded_idx(bond_id)) += energy.value();
    } else {
      std::array<int, 2> pids = {aosoa.id(row_j), aosoa.id(row_k)};
      bond_broken_error(aosoa.id(row_i), {pids.data(), 2});
    }
  }
};

struct DihedralBondsEnergyKernel {
  BondsEnergyKernelData data;
  LocalBondState::DihedralBondlistType bond_list;
  LocalBondState::DihedralBondIDType bond_ids;

  DihedralBondsEnergyKernel(BondsEnergyKernelData data_,
                            LocalBondState::DihedralBondlistType bond_list_,
                            LocalBondState::DihedralBondIDType bond_ids_)
      : data(std::move(data_)), bond_list(std::move(bond_list_)),
        bond_ids(std::move(bond_ids_)) {}

  ESPRESSO_ATTR_ALWAYS_INLINE inline void operator()(std::size_t idx) const {
    auto const &bonded_ias = data.bonded_ias;
    auto const &box_geo = data.box_geo;
    auto &local_energy = data.local_energy;
    auto const &layout = data.layout;
    auto const &aosoa = data.aosoa;
    auto const bond_id = bond_ids(idx);

    // TODO: omp_get_thread_num() is only available for the OpenMP backend.
    // This should be updated when using other Kokkos backends.
    auto const thread_id = omp_get_thread_num();

    auto const i = bond_list(idx, 0);
    auto const j = bond_list(idx, 1);
    auto const k = bond_list(idx, 2);
    auto const m = bond_list(idx, 3);
    // Translate pack indices to ParticleStore rows once.
    auto const row_i = aosoa.row(i);
    auto const row_j = aosoa.row(j);
    auto const row_k = aosoa.row(k);
    auto const row_m = aosoa.row(m);
    auto const &iaparams = *bonded_ias.at(bond_id);

    auto const pos1 = aosoa.get_vector_at(aosoa.position, row_i);
    auto const pos2 = aosoa.get_vector_at(aosoa.position, row_j);
    auto const pos3 = aosoa.get_vector_at(aosoa.position, row_k);
    auto const pos4 = aosoa.get_vector_at(aosoa.position, row_m);
    auto const v12 = box_geo.get_mi_vector(pos1, pos2);
    auto const v23 = box_geo.get_mi_vector(pos3, pos1);
    auto const v34 = box_geo.get_mi_vector(pos4, pos3);

    std::optional<double> energy =
        calc_dihedral_bonded_energy(iaparams, v12, v23, v34);

    if (energy) {
      local_energy(thread_id, layout.bonded_idx(bond_id)) += energy.value();
    } else {
      std::array<int, 3> pids = {aosoa.id(row_j), aosoa.id(row_k),
                                 aosoa.id(row_m)};
      bond_broken_error(aosoa.id(row_i), {pids.data(), 3});
    }
  }
};
