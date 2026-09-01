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

// Dispatched particle-parallel like the force kernels (see
// bond_forces_kokkos.hpp), for consistency and so that a single per-particle
// gather structure (LocalBondState::pp_*) drives every bond kernel. Energy
// never needed a ScatterView (it already accumulated into a per-thread
// private (thread_id, bin) array), so this is purely a dispatch-axis change:
// each row is only a candidate contribution for the participant holding its
// bond's primary entry -- mirror rows are skipped -- which reproduces
// exactly the old per-bond kernel's single evaluation with "self" playing
// the role of the old kernel's owner particle "i".
struct PairBondsEnergyKernel {
  BondsEnergyKernelData data;
  LocalBondState::PPPairDegreeType pp_pair_degree;
  LocalBondState::PPPairSlotType pp_pair_slots;
  Coulomb::ShortRangeEnergyKernel::kernel_type const *const coulomb_u_kernel;

  PairBondsEnergyKernel(
      BondsEnergyKernelData data_,
      LocalBondState::PPPairDegreeType pp_pair_degree_,
      LocalBondState::PPPairSlotType pp_pair_slots_,
      Coulomb::ShortRangeEnergyKernel::kernel_type const *coulomb_u_kernel_)
      : data(std::move(data_)), pp_pair_degree(std::move(pp_pair_degree_)),
        pp_pair_slots(std::move(pp_pair_slots_)),
        coulomb_u_kernel(coulomb_u_kernel_) {}

  ESPRESSO_ATTR_ALWAYS_INLINE inline void operator()(std::size_t self) const {
    auto const &bonded_ias = data.bonded_ias;
    auto const &box_geo = data.box_geo;
    auto &local_energy = data.local_energy;
    auto const &layout = data.layout;
    auto const &aosoa = data.aosoa;
    // TODO: omp_get_thread_num() is only available for the OpenMP backend.
    // This should be updated when using other Kokkos backends.
    auto const thread_id = omp_get_thread_num();
    auto const degree = pp_pair_degree(self);

    for (int slot = 0; slot < degree; ++slot) {
      auto const j = pp_pair_slots(self, slot, 0);
      // mirror, or an unresolvable row (-1, see
      // CellStructure::update_bond_storage) -- skip either way.
      if (pp_pair_slots(self, slot, 2) == 0 or j < 0) {
        continue;
      }
      auto const i = self;
      auto const bond_id = pp_pair_slots(self, slot, 1);
      auto const &iaparams = *bonded_ias.at(bond_id);

      auto const pos1 = aosoa.get_vector_at(aosoa.position, i);
      auto const pos2 = aosoa.get_vector_at(aosoa.position, j);
      auto const dx = box_geo.get_mi_vector(pos1, pos2);

      std::optional<double> energy = calc_pair_bonded_energy(
          iaparams, dx, pos1, pos2,
#ifdef ESPRESSO_ELECTROSTATICS
          aosoa.charge(i) * aosoa.charge(j), coulomb_u_kernel
#else
          0.0, nullptr
#endif
      );

      if (energy) {
        local_energy(thread_id, layout.bonded_idx(bond_id)) += energy.value();
      } else {
        auto partner_id = aosoa.id(j);
        bond_broken_error(aosoa.id(i), {&partner_id, 1});
      }
    }
  }
};

struct AngleBondsEnergyKernel {
  BondsEnergyKernelData data;
  LocalBondState::PPAngleDegreeType pp_angle_degree;
  LocalBondState::PPAngleSlotType pp_angle_slots;

  AngleBondsEnergyKernel(BondsEnergyKernelData data_,
                         LocalBondState::PPAngleDegreeType pp_angle_degree_,
                         LocalBondState::PPAngleSlotType pp_angle_slots_)
      : data(std::move(data_)), pp_angle_degree(std::move(pp_angle_degree_)),
        pp_angle_slots(std::move(pp_angle_slots_)) {}

  ESPRESSO_ATTR_ALWAYS_INLINE inline void operator()(std::size_t self) const {
    auto const &bonded_ias = data.bonded_ias;
    auto const &box_geo = data.box_geo;
    auto &local_energy = data.local_energy;
    auto const &layout = data.layout;
    auto const &aosoa = data.aosoa;
    auto const thread_id = omp_get_thread_num();
    auto const degree = pp_angle_degree(self);

    for (int slot = 0; slot < degree; ++slot) {
      // Only self_slot==0 (self == vertex == column 0) contributes, to
      // avoid triple-counting; -1 (unresolvable row) also fails this.
      if (pp_angle_slots(self, slot, 4) != 0) {
        continue;
      }
      auto const i = self;
      auto const j = pp_angle_slots(self, slot, 1);
      auto const k = pp_angle_slots(self, slot, 2);
      auto const bond_id = pp_angle_slots(self, slot, 3);
      auto const &iaparams = *bonded_ias.at(bond_id);

      auto const pos1 = aosoa.get_vector_at(aosoa.position, i);
      auto const pos2 = aosoa.get_vector_at(aosoa.position, j);
      auto const pos3 = aosoa.get_vector_at(aosoa.position, k);
      auto const vec1 = box_geo.get_mi_vector(pos2, pos1);
      auto const vec2 = box_geo.get_mi_vector(pos3, pos1);

      std::optional<double> energy =
          calc_angle_bonded_energy(iaparams, vec1, vec2);

      if (energy) {
        local_energy(thread_id, layout.bonded_idx(bond_id)) += energy.value();
      } else {
        std::array<int, 2> pids = {aosoa.id(j), aosoa.id(k)};
        bond_broken_error(aosoa.id(i), {pids.data(), 2});
      }
    }
  }
};

struct DihedralBondsEnergyKernel {
  BondsEnergyKernelData data;
  LocalBondState::PPDihedralDegreeType pp_dihedral_degree;
  LocalBondState::PPDihedralSlotType pp_dihedral_slots;

  DihedralBondsEnergyKernel(
      BondsEnergyKernelData data_,
      LocalBondState::PPDihedralDegreeType pp_dihedral_degree_,
      LocalBondState::PPDihedralSlotType pp_dihedral_slots_)
      : data(std::move(data_)),
        pp_dihedral_degree(std::move(pp_dihedral_degree_)),
        pp_dihedral_slots(std::move(pp_dihedral_slots_)) {}

  ESPRESSO_ATTR_ALWAYS_INLINE inline void operator()(std::size_t self) const {
    auto const &bonded_ias = data.bonded_ias;
    auto const &box_geo = data.box_geo;
    auto &local_energy = data.local_energy;
    auto const &layout = data.layout;
    auto const &aosoa = data.aosoa;
    auto const thread_id = omp_get_thread_num();
    auto const degree = pp_dihedral_degree(self);

    for (int slot = 0; slot < degree; ++slot) {
      if (pp_dihedral_slots(self, slot, 5) != 0) { // not chain position 0
        continue;
      }
      auto const i = self;
      auto const j = pp_dihedral_slots(self, slot, 1);
      auto const k = pp_dihedral_slots(self, slot, 2);
      auto const m = pp_dihedral_slots(self, slot, 3);
      auto const bond_id = pp_dihedral_slots(self, slot, 4);
      auto const &iaparams = *bonded_ias.at(bond_id);

      auto const pos1 = aosoa.get_vector_at(aosoa.position, i);
      auto const pos2 = aosoa.get_vector_at(aosoa.position, j);
      auto const pos3 = aosoa.get_vector_at(aosoa.position, k);
      auto const pos4 = aosoa.get_vector_at(aosoa.position, m);
      auto const v12 = box_geo.get_mi_vector(pos1, pos2);
      auto const v23 = box_geo.get_mi_vector(pos3, pos1);
      auto const v34 = box_geo.get_mi_vector(pos4, pos3);

      std::optional<double> energy =
          calc_dihedral_bonded_energy(iaparams, v12, v23, v34);

      if (energy) {
        local_energy(thread_id, layout.bonded_idx(bond_id)) += energy.value();
      } else {
        std::array<int, 3> pids = {aosoa.id(j), aosoa.id(k), aosoa.id(m)};
        bond_broken_error(aosoa.id(i), {pids.data(), 3});
      }
    }
  }
};
