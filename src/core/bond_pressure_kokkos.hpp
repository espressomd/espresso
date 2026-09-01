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
#include "pressure_cabana.hpp"
#include "pressure_inline.hpp"

#include <utils/Vector.hpp>
#include <utils/math/tensor_product.hpp>

#include <Kokkos_Core.hpp>

#include <omp.h>

#include <array>
#include <cstddef>

struct BondsPressureKernelData {
  using execution_space = Kokkos::HostSpace;
  BondedInteractionsMap const &bonded_ias;
  BoxGeometry const &box_geo;
  Kokkos::View<double **, Kokkos::LayoutRight, execution_space> local_pressure;
  PressureBinLayout layout;
  CellStructure::AoSoA_pack const &aosoa;
};

// Dispatched particle-parallel, gated to primary rows -- see
// bond_energy_kokkos.hpp's PairBondsEnergyKernel doc comment for the
// rationale; pressure never used a ScatterView either, so this is again
// purely a dispatch-axis change with a mirror-skipping gate.
struct PairBondsPressureKernel {
  BondsPressureKernelData data;
  LocalBondState::PPPairDegreeType pp_pair_degree;
  LocalBondState::PPPairSlotType pp_pair_slots;
  Coulomb::ShortRangeForceKernel::kernel_type const *const coulomb_f_kernel;

  PairBondsPressureKernel(
      BondsPressureKernelData data_,
      LocalBondState::PPPairDegreeType pp_pair_degree_,
      LocalBondState::PPPairSlotType pp_pair_slots_,
      Coulomb::ShortRangeForceKernel::kernel_type const *coulomb_f_kernel_)
      : data(std::move(data_)), pp_pair_degree(std::move(pp_pair_degree_)),
        pp_pair_slots(std::move(pp_pair_slots_)),
        coulomb_f_kernel(coulomb_f_kernel_) {}

  ESPRESSO_ATTR_ALWAYS_INLINE inline void operator()(std::size_t self) const {
    auto const &bonded_ias = data.bonded_ias;
    auto const &box_geo = data.box_geo;
    auto &local_pressure = data.local_pressure;
    auto const &layout = data.layout;
    auto const &aosoa = data.aosoa;
    auto const thread_id = omp_get_thread_num();
    auto const degree = pp_pair_degree(self);

    for (int slot = 0; slot < degree; ++slot) {
      auto const j = pp_pair_slots(self, slot, 0);
      // mirror, or an unresolvable row (-1) -- skip either way.
      if (pp_pair_slots(self, slot, 2) == 0 or j < 0) {
        continue;
      }
      auto const i = self;
      auto const bond_id = pp_pair_slots(self, slot, 1);
      auto const &iaparams = *bonded_ias.at(bond_id);

      auto const pos1 = aosoa.get_vector_at(aosoa.position, i);
      auto const pos2 = aosoa.get_vector_at(aosoa.position, j);

      std::optional<Utils::Matrix<double, 3, 3>> pressure =
          calc_bonded_virial_pressure_tensor(iaparams, pos1, pos2, box_geo,
                                             coulomb_f_kernel,
#ifdef ESPRESSO_ELECTROSTATICS
                                             aosoa.charge(i) * aosoa.charge(j)
#else
                                             0.0
#endif
          );

      if (pressure) {
        auto const flat = Utils::flatten(*pressure);
        for (std::size_t k = 0; k < 9; ++k)
          local_pressure(thread_id,
                         layout.tensor_offset(layout.bonded_idx(bond_id), k)) +=
              flat[k];
      } else {
        auto partner_id = aosoa.id(j);
        bond_broken_error(aosoa.id(i), {&partner_id, 1});
      }
    }
  }
};

struct AngleBondsPressureKernel {
  BondsPressureKernelData data;
  LocalBondState::PPAngleDegreeType pp_angle_degree;
  LocalBondState::PPAngleSlotType pp_angle_slots;

  AngleBondsPressureKernel(BondsPressureKernelData data_,
                           LocalBondState::PPAngleDegreeType pp_angle_degree_,
                           LocalBondState::PPAngleSlotType pp_angle_slots_)
      : data(std::move(data_)), pp_angle_degree(std::move(pp_angle_degree_)),
        pp_angle_slots(std::move(pp_angle_slots_)) {}

  ESPRESSO_ATTR_ALWAYS_INLINE inline void operator()(std::size_t self) const {
    auto const &bonded_ias = data.bonded_ias;
    auto const &box_geo = data.box_geo;
    auto &local_pressure = data.local_pressure;
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

      std::optional<Utils::Matrix<double, 3, 3>> pressure =
          calc_bonded_three_body_pressure_tensor(iaparams, pos1, pos2, pos3,
                                                 box_geo);

      if (pressure) {
        auto const flat = Utils::flatten(*pressure);
        for (std::size_t k2 = 0; k2 < 9; ++k2)
          local_pressure(
              thread_id,
              layout.tensor_offset(layout.bonded_idx(bond_id), k2)) += flat[k2];
      } else {
        std::array<int, 2> pids = {aosoa.id(j), aosoa.id(k)};
        bond_broken_error(aosoa.id(i), {pids.data(), 2});
      }
    }
  }
};

struct DihedralBondsPressureKernel {
  BondsPressureKernelData data;
  LocalBondState::PPDihedralDegreeType pp_dihedral_degree;
  LocalBondState::PPDihedralSlotType pp_dihedral_slots;

  DihedralBondsPressureKernel(
      BondsPressureKernelData data_,
      LocalBondState::PPDihedralDegreeType pp_dihedral_degree_,
      LocalBondState::PPDihedralSlotType pp_dihedral_slots_)
      : data(std::move(data_)),
        pp_dihedral_degree(std::move(pp_dihedral_degree_)),
        pp_dihedral_slots(std::move(pp_dihedral_slots_)) {}

  ESPRESSO_ATTR_ALWAYS_INLINE inline void operator()(std::size_t self) const {
    auto const &bonded_ias = data.bonded_ias;
    auto const &box_geo = data.box_geo;
    auto &local_pressure = data.local_pressure;
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

      std::optional<Utils::Matrix<double, 3, 3>> pressure =
          calc_bonded_four_body_pressure_tensor(iaparams, pos1, pos2, pos3,
                                                pos4, box_geo);

      if (pressure) {
        auto const flat = Utils::flatten(*pressure);
        for (std::size_t k3 = 0; k3 < 9; ++k3)
          local_pressure(
              thread_id,
              layout.tensor_offset(layout.bonded_idx(bond_id), k3)) += flat[k3];
      } else {
        std::array<int, 3> pids = {aosoa.id(j), aosoa.id(k), aosoa.id(m)};
        bond_broken_error(aosoa.id(i), {pids.data(), 3});
      }
    }
  }
};
