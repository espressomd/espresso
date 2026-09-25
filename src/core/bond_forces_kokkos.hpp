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
#include "cell_system/LocalBondState.hpp"
#include "forces_inline.hpp"

#include <utils/Vector.hpp>

#include <Kokkos_Core.hpp>
#include <Kokkos_ScatterView.hpp>

#include <cstddef>
#include <optional>
#include <variant>

struct BondsKernelData {
  BondedInteractionsMap const &bonded_ias;
  BondBreakage::BondBreakage &bond_breakage;
  BoxGeometry const &box_geo;
  CellStructure::ScatterForce local_force;
#ifdef ESPRESSO_NPT
  CellStructure::ScatterVirial local_virial;
#endif
  CellStructure::AoSoA_pack const &aosoa;
  bool const has_breakage_specs;
};

struct PairBondsKernel {
  BondsKernelData data;
  LocalBondState::PairBondlistType bond_list;
  LocalBondState::PairBondIDType bond_ids;
  Coulomb::ShortRangeForceKernel::kernel_type const *const coulomb_kernel;

  PairBondsKernel(
      BondsKernelData data_, LocalBondState::PairBondlistType bond_list_,
      LocalBondState::PairBondIDType bond_ids_,
      Coulomb::ShortRangeForceKernel::kernel_type const *coulomb_kernel_)
      : data(std::move(data_)), bond_list(std::move(bond_list_)),
        bond_ids(std::move(bond_ids_)), coulomb_kernel(coulomb_kernel_) {}

  ESPRESSO_ATTR_ALWAYS_INLINE inline void operator()(std::size_t idx) const {
    auto const &bonded_ias = data.bonded_ias;
    auto const &box_geo = data.box_geo;
    auto local_force = data.local_force.access();
    auto const &aosoa = data.aosoa;
    auto &bond_breakage = data.bond_breakage;
#ifdef ESPRESSO_NPT
    auto local_virial = data.local_virial.access();
#endif
    auto const has_breakage_specs = data.has_breakage_specs;
    auto const bond_id = bond_ids(idx);

    auto const i = bond_list(idx, 0);
    auto const j = bond_list(idx, 1);
    auto const &iaparams = *bonded_ias.at(bond_id);

    auto const dx =
        box_geo.get_mi_vector(aosoa.get_vector_at(aosoa.position, i),
                              aosoa.get_vector_at(aosoa.position, j));
    //  Consider for bond breakage
    if (has_breakage_specs &&
        bond_breakage.check_and_handle_breakage(
            aosoa.id(i), {{aosoa.id(j), std::nullopt}}, bond_id, dx.norm())) {
      return;
    }

    if (auto const *iap = std::get_if<ThermalizedBond>(&iaparams)) {
      auto const result = iap->forces(
#ifdef ESPRESSO_MASS
          aosoa.mass(i), aosoa.mass(j),
#else
          1.0, 1.0,
#endif
          aosoa.get_vector_at(aosoa.velocity, i),
          aosoa.get_vector_at(aosoa.velocity, j), aosoa.id(i), aosoa.id(j), dx);
      if (result) {
        auto const &forces = result.value();

        local_force(i, 0) += std::get<0>(forces)[0];
        local_force(i, 1) += std::get<0>(forces)[1];
        local_force(i, 2) += std::get<0>(forces)[2];
        local_force(j, 0) += std::get<1>(forces)[0];
        local_force(j, 1) += std::get<1>(forces)[1];
        local_force(j, 2) += std::get<1>(forces)[2];
      } else {
        auto partner_id = aosoa.id(j);
        bond_broken_error(aosoa.id(i), {&partner_id, 1});
      }
      return;
    }

    auto const result =
        calc_bond_pair_force(iaparams, dx,
#ifdef ESPRESSO_ELECTROSTATICS
                             aosoa.charge(i) * aosoa.charge(j), coulomb_kernel
#else
                             0.0, nullptr
#endif
        );

    if (result) {
      auto const f = result.value();
      local_force(i, 0) += f[0];
      local_force(i, 1) += f[1];
      local_force(i, 2) += f[2];
      local_force(j, 0) -= f[0];
      local_force(j, 1) -= f[1];
      local_force(j, 2) -= f[2];
#ifdef ESPRESSO_NPT
      auto const virial = hadamard_product(f, dx);
      local_virial(0) += virial[0];
      local_virial(1) += virial[1];
      local_virial(2) += virial[2];
#endif
    } else {
      auto partner_id = aosoa.id(j);
      bond_broken_error(aosoa.id(i), {&partner_id, 1});
    }
  }
};

struct AngleBondsKernel {
  BondsKernelData data;
  LocalBondState::AngleBondlistType bond_list;
  LocalBondState::AngleBondIDType bond_ids;

  AngleBondsKernel(BondsKernelData data_,
                   LocalBondState::AngleBondlistType bond_list_,
                   LocalBondState::AngleBondIDType bond_ids_)
      : data(std::move(data_)), bond_list(std::move(bond_list_)),
        bond_ids(std::move(bond_ids_)) {}

  ESPRESSO_ATTR_ALWAYS_INLINE inline void operator()(std::size_t idx) const {
    auto const &bonded_ias = data.bonded_ias;
    auto const &box_geo = data.box_geo;
    auto local_force = data.local_force.access();
    auto const &aosoa = data.aosoa;
    auto &bond_breakage = data.bond_breakage;
    auto const has_breakage_specs = data.has_breakage_specs;
    auto const bond_id = bond_ids(idx);

    auto const i = bond_list(idx, 0);
    auto const j = bond_list(idx, 1);
    auto const k = bond_list(idx, 2);
    auto const &iaparams = *bonded_ias.at(bond_id);

    auto const pos1 = aosoa.get_vector_at(aosoa.position, i);
    auto const pos2 = aosoa.get_vector_at(aosoa.position, j);
    auto const pos3 = aosoa.get_vector_at(aosoa.position, k);
    auto const vec1 = box_geo.get_mi_vector(pos2, pos1);
    auto const vec2 = box_geo.get_mi_vector(pos3, pos1);

    //  Consider for bond breakage
    if (has_breakage_specs &&
        bond_breakage.check_and_handle_breakage(
            aosoa.id(i), {{aosoa.id(j), aosoa.id(k)}}, bond_id,
            box_geo.get_mi_vector(pos2, pos3).norm())) {
      return;
    }
    if (std::get_if<OifGlobalForcesBond>(&iaparams)) {
      return;
    }

    auto const result = calc_bonded_three_body_force(iaparams, vec1, vec2);

    if (result) {
      auto const &forces = result.value();

      local_force(i, 0) += std::get<0>(forces)[0];
      local_force(i, 1) += std::get<0>(forces)[1];
      local_force(i, 2) += std::get<0>(forces)[2];
      local_force(j, 0) += std::get<1>(forces)[0];
      local_force(j, 1) += std::get<1>(forces)[1];
      local_force(j, 2) += std::get<1>(forces)[2];
      local_force(k, 0) += std::get<2>(forces)[0];
      local_force(k, 1) += std::get<2>(forces)[1];
      local_force(k, 2) += std::get<2>(forces)[2];
    } else {
      std::array<int, 2> pids = {aosoa.id(j), aosoa.id(k)};
      bond_broken_error(aosoa.id(i), {pids.data(), 2});
    }
  }
};

struct DihedralBondsKernel {
  BondsKernelData data;
  LocalBondState::DihedralBondlistType bond_list;
  LocalBondState::DihedralBondIDType bond_ids;

  DihedralBondsKernel(BondsKernelData data_,
                      LocalBondState::DihedralBondlistType bond_list_,
                      LocalBondState::DihedralBondIDType bond_ids_)
      : data(std::move(data_)), bond_list(std::move(bond_list_)),
        bond_ids(std::move(bond_ids_)) {}

  ESPRESSO_ATTR_ALWAYS_INLINE inline void operator()(std::size_t idx) const {
    auto const &bonded_ias = data.bonded_ias;
    auto const &box_geo = data.box_geo;
    auto local_force = data.local_force.access();
    auto const &aosoa = data.aosoa;
    auto const bond_id = bond_ids(idx);

    auto const i = bond_list(idx, 0);
    auto const j = bond_list(idx, 1);
    auto const k = bond_list(idx, 2);
    auto const m = bond_list(idx, 3);
    auto const &iaparams = *bonded_ias.at(bond_id);

    auto const pos1 = aosoa.get_vector_at(aosoa.position, i);
    auto const pos2 = aosoa.get_vector_at(aosoa.position, j);
    auto const pos3 = aosoa.get_vector_at(aosoa.position, k);
    auto const pos4 = aosoa.get_vector_at(aosoa.position, m);
    auto const vel1 = aosoa.get_vector_at(aosoa.velocity, i);
    auto const vel3 = aosoa.get_vector_at(aosoa.velocity, k);
    auto const image1 = aosoa.get_vector_at(aosoa.image, i);

    auto const result = calc_bonded_four_body_force(
        iaparams, box_geo, pos1, pos2, pos3, pos4, vel1, vel3, image1);

    if (result) {
      auto const &forces = result.value();

      local_force(i, 0) += std::get<0>(forces)[0];
      local_force(i, 1) += std::get<0>(forces)[1];
      local_force(i, 2) += std::get<0>(forces)[2];
      local_force(j, 0) += std::get<1>(forces)[0];
      local_force(j, 1) += std::get<1>(forces)[1];
      local_force(j, 2) += std::get<1>(forces)[2];
      local_force(k, 0) += std::get<2>(forces)[0];
      local_force(k, 1) += std::get<2>(forces)[1];
      local_force(k, 2) += std::get<2>(forces)[2];
      local_force(m, 0) += std::get<3>(forces)[0];
      local_force(m, 1) += std::get<3>(forces)[1];
      local_force(m, 2) += std::get<3>(forces)[2];
    } else {
      std::array<int, 3> pids = {aosoa.id(j), aosoa.id(k), aosoa.id(m)};
      bond_broken_error(aosoa.id(i), {pids.data(), 3});
    }
  }
};
