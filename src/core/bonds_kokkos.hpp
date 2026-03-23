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

#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM

#include "aosoa_pack.hpp"
#include "forces_inline.hpp"

#include <utils/Vector.hpp>

#include <Kokkos_Core.hpp>

#include <omp.h>

#include <cstddef>
#include <memory>
#include <optional>
#include <variant>
#include <vector>

struct BondsKernelData {
  BondedInteractionsMap const &bonded_ias;
  BondBreakage::BondBreakage &bond_breakage;
  BoxGeometry const &box_geo;
  CellStructure::ForceType &local_force;
#ifdef ESPRESSO_NPT
  CellStructure::VirialType &local_virial;
#endif
  CellStructure::AoSoA_pack const &aosoa;
  bool const has_breakage_specs;
};

struct PairBondsKernel {
  BondsKernelData const &data;
  CellStructure::PairBondlistType const &bond_list;
  CellStructure::PairBondIDType const &bond_ids;
  Coulomb::ShortRangeForceKernel::kernel_type const *const coulomb_kernel;

  PairBondsKernel(
      BondsKernelData const &data_,
      CellStructure::PairBondlistType const &bond_list_,
      CellStructure::PairBondIDType const &bond_ids_,
      Coulomb::ShortRangeForceKernel::kernel_type const *coulomb_kernel_)
      : data(data_), bond_list(bond_list_), bond_ids(bond_ids_),
        coulomb_kernel(coulomb_kernel_) {}

  ESPRESSO_ATTR_ALWAYS_INLINE KOKKOS_INLINE_FUNCTION void
  operator()(std::size_t idx) const {
    auto const &bonded_ias = data.bonded_ias;
    auto const &box_geo = data.box_geo;
    auto &local_force = data.local_force;
    auto const &aosoa = data.aosoa;
    auto &bond_breakage = data.bond_breakage;
#ifdef ESPRESSO_NPT
    auto &local_virial = data.local_virial;
#endif
    auto const has_breakage_specs = data.has_breakage_specs;

    auto const &partners = Kokkos::subview(bond_list, idx, Kokkos::ALL);
    auto const &bond_id = bond_ids(idx);

    auto const i = partners(0);

    auto const &iaparams = *bonded_ias.at(bond_id);
    // TODO: omp_get_thread_num() is only available for the OpenMP backend.
    // This should be updated when using other Kokkos backends.
    auto const thread_id = omp_get_thread_num();

    auto const j = partners(1);
    auto const dx =
        box_geo.get_mi_vector(aosoa.get_vector_at(aosoa.position, i),
                              aosoa.get_vector_at(aosoa.position, j));
    std::optional<Utils::Vector3d> result;
    // Consider for bond breakage
    if (has_breakage_specs &&
        bond_breakage.check_and_handle_breakage(
            aosoa.id(i), {{aosoa.id(j), std::nullopt}}, bond_id, dx.norm())) {
      return;
    }

    if (auto const *iap = std::get_if<ThermalizedBond>(&iaparams)) {
      auto const res = iap->forces(
#ifdef ESPRESSO_MASS
          aosoa.mass(i), aosoa.mass(j),
#else
	  1.0, 1.0,
#endif
	  aosoa.get_vector_at(aosoa.velocity, i),
          aosoa.get_vector_at(aosoa.velocity, j), aosoa.id(i), aosoa.id(j), dx);
      if (res) {
        auto const &forces = res.value();

        local_force(i, thread_id, 0) += std::get<0>(forces)[0];
        local_force(i, thread_id, 1) += std::get<0>(forces)[1];
        local_force(i, thread_id, 2) += std::get<0>(forces)[2];
        local_force(j, thread_id, 0) += std::get<1>(forces)[0];
        local_force(j, thread_id, 1) += std::get<1>(forces)[1];
        local_force(j, thread_id, 2) += std::get<1>(forces)[2];
      } else {
        auto partner_id = aosoa.id(j);
        bond_broken_error(aosoa.id(i), {&partner_id, 1});
      }
      return;
    }

    result =
        calc_bond_pair_force(iaparams, dx
#ifdef ESPRESSO_ELECTROSTATICS
                             ,
                             aosoa.charge(i) * aosoa.charge(j), coulomb_kernel
#endif
        );

    if (result) {
      auto const f = result.value();
      local_force(i, thread_id, 0) += f[0];
      local_force(i, thread_id, 1) += f[1];
      local_force(i, thread_id, 2) += f[2];
      local_force(j, thread_id, 0) -= f[0];
      local_force(j, thread_id, 1) -= f[1];
      local_force(j, thread_id, 2) -= f[2];
#ifdef ESPRESSO_NPT
      auto virial = hadamard_product(result.value(), dx);
      local_virial(thread_id, 0) += virial[0];
      local_virial(thread_id, 1) += virial[1];
      local_virial(thread_id, 2) += virial[2];
#endif
    } else {
      auto partner_id = aosoa.id(j);
      bond_broken_error(aosoa.id(i), {&partner_id, 1});
    }
  }
};

struct AngleBondsKernel {
  BondsKernelData const &data;
  CellStructure::AngleBondlistType const &bond_list;
  CellStructure::AngleBondIDType const &bond_ids;

  AngleBondsKernel(BondsKernelData const &data_,
                   CellStructure::AngleBondlistType const &bond_list_,
                   CellStructure::AngleBondIDType const &bond_ids_)
      : data(data_), bond_list(bond_list_), bond_ids(bond_ids_) {}

  ESPRESSO_ATTR_ALWAYS_INLINE KOKKOS_INLINE_FUNCTION void
  operator()(std::size_t idx) const {
    auto const &bonded_ias = data.bonded_ias;
    auto const &box_geo = data.box_geo;
    auto const &local_force = data.local_force;
    auto const &aosoa = data.aosoa;
    auto &bond_breakage = data.bond_breakage;
    auto const has_breakage_specs = data.has_breakage_specs;

    auto const &partners = Kokkos::subview(bond_list, idx, Kokkos::ALL);
    auto const &bond_id = bond_ids(idx);

    auto const i = partners(0);

    auto const &iaparams = *bonded_ias.at(bond_id);
    // TODO: omp_get_thread_num() is only available for the OpenMP backend.
    // This should be updated when using other Kokkos backends.
    auto const thread_id = omp_get_thread_num();

    auto const j = partners(1);
    auto const k = partners(2);
    auto const pos1 = aosoa.get_vector_at(aosoa.position, i);
    auto const pos2 = aosoa.get_vector_at(aosoa.position, j);
    auto const pos3 = aosoa.get_vector_at(aosoa.position, k);
    auto const vec1 = box_geo.get_mi_vector(pos2, pos1);
    auto const vec2 = box_geo.get_mi_vector(pos3, pos1);
    std::optional<std::tuple<Utils::Vector3d, Utils::Vector3d, Utils::Vector3d>>
        result;
    // Consider for bond breakage
    if (has_breakage_specs &&
        bond_breakage.check_and_handle_breakage(
            aosoa.id(i), {{aosoa.id(j), aosoa.id(k)}}, bond_id,
            box_geo.get_mi_vector(pos2, pos3).norm())) {
      return;
    }
    if (std::get_if<OifGlobalForcesBond>(&iaparams)) {
      return;
    }

    result = calc_bonded_three_body_force(iaparams, vec1, vec2);

    if (result) {
      auto const &forces = result.value();

      local_force(i, thread_id, 0) += std::get<0>(forces)[0];
      local_force(i, thread_id, 1) += std::get<0>(forces)[1];
      local_force(i, thread_id, 2) += std::get<0>(forces)[2];
      local_force(j, thread_id, 0) += std::get<1>(forces)[0];
      local_force(j, thread_id, 1) += std::get<1>(forces)[1];
      local_force(j, thread_id, 2) += std::get<1>(forces)[2];
      local_force(k, thread_id, 0) += std::get<2>(forces)[0];
      local_force(k, thread_id, 1) += std::get<2>(forces)[1];
      local_force(k, thread_id, 2) += std::get<2>(forces)[2];
    } else {
      std::array<int, 2> pids = {aosoa.id(j), aosoa.id(k)};
      bond_broken_error(aosoa.id(i), {pids.data(), 2});
    }
  }
};

struct DihedralBondsKernel {
  BondsKernelData const &data;
  CellStructure::DihedralBondlistType const &bond_list;
  CellStructure::DihedralBondIDType const &bond_ids;

  DihedralBondsKernel(BondsKernelData const &data_,
                      CellStructure::DihedralBondlistType const &bond_list_,
                      CellStructure::DihedralBondIDType const &bond_ids_)
      : data(data_), bond_list(bond_list_), bond_ids(bond_ids_) {}

  ESPRESSO_ATTR_ALWAYS_INLINE KOKKOS_INLINE_FUNCTION void
  operator()(std::size_t idx) const {
    auto const &bonded_ias = data.bonded_ias;
    auto const &box_geo = data.box_geo;
    auto const &local_force = data.local_force;
    auto const &aosoa = data.aosoa;

    auto const &partners = Kokkos::subview(bond_list, idx, Kokkos::ALL);
    auto const &bond_id = bond_ids(idx);

    auto const i = partners(0);

    auto const &iaparams = *bonded_ias.at(bond_id);
    // TODO: omp_get_thread_num() is only available for the OpenMP backend.
    // This should be updated when using other Kokkos backends.
    auto const thread_id = omp_get_thread_num();

    auto const j = partners(1);
    auto const k = partners(2);
    auto const m = partners(3);
    auto const pos1 = aosoa.get_vector_at(aosoa.position, i);
    auto const pos2 = aosoa.get_vector_at(aosoa.position, j);
    auto const pos3 = aosoa.get_vector_at(aosoa.position, k);
    auto const pos4 = aosoa.get_vector_at(aosoa.position, m);
    auto const vel1 = aosoa.get_vector_at(aosoa.velocity, i);
    auto const vel3 = aosoa.get_vector_at(aosoa.velocity, k);
    auto const image1 = aosoa.get_vector_at(aosoa.image, i);

    std::optional<std::tuple<Utils::Vector3d, Utils::Vector3d, Utils::Vector3d,
                             Utils::Vector3d>>
        result = calc_bonded_four_body_force(iaparams, box_geo, pos1, pos2,
                                             pos3, pos4, vel1, vel3, image1);

    if (result) {
      auto const &forces = result.value();

      local_force(i, thread_id, 0) += std::get<0>(forces)[0];
      local_force(i, thread_id, 1) += std::get<0>(forces)[1];
      local_force(i, thread_id, 2) += std::get<0>(forces)[2];
      local_force(j, thread_id, 0) += std::get<1>(forces)[0];
      local_force(j, thread_id, 1) += std::get<1>(forces)[1];
      local_force(j, thread_id, 2) += std::get<1>(forces)[2];
      local_force(k, thread_id, 0) += std::get<2>(forces)[0];
      local_force(k, thread_id, 1) += std::get<2>(forces)[1];
      local_force(k, thread_id, 2) += std::get<2>(forces)[2];
      local_force(m, thread_id, 0) += std::get<3>(forces)[0];
      local_force(m, thread_id, 1) += std::get<3>(forces)[1];
      local_force(m, thread_id, 2) += std::get<3>(forces)[2];
    } else {
      std::array<int, 3> pids = {aosoa.id(j), aosoa.id(k), aosoa.id(m)};
      bond_broken_error(aosoa.id(i), {pids.data(), 3});
    }
  }
};

#endif // ESPRESSO_SHARED_MEMORY_PARALLELISM
