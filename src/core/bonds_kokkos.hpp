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

#if defined(__GNUG__) or defined(__clang__)
#define ESPRESSO_ATTR_ALWAYS_INLINE [[gnu::always_inline]]
#else
#define ESPRESSO_ATTR_ALWAYS_INLINE
#endif

struct BondsKernel {
  BondedInteractionsMap const &bonded_ias;
  BondBreakage::BondBreakage &bond_breakage;
  Coulomb::ShortRangeForceKernel::kernel_type const *const coulomb_kernel;
  BoxGeometry const &box_geo;
  Kokkos::View<int *> const &id_to_index;
  CellStructure::ForceType const &local_force;
#ifdef ESPRESSO_NPT
  CellStructure::VirialType const &local_virial;
#endif
  CellStructure::BondlistType const &bond_list;
  CellStructure::BondIDType const &bond_ids;
  CellStructure::AoSoA_pack const &aosoa;

  BondsKernel(
      BondedInteractionsMap const &bonded_ias_,
      BondBreakage::BondBreakage &bond_breakage_,
      Coulomb::ShortRangeForceKernel::kernel_type const *coulomb_kernel_,
      BoxGeometry const &box_geo_,
      Kokkos::View<int *> const &id_to_index_,
      CellStructure::ForceType const &local_force_,
#ifdef ESPRESSO_NPT
      CellStructure::VirialType const &local_virial_,
#endif
      CellStructure::BondlistType const &bond_list_,
      CellStructure::BondIDType const &bond_ids_,
      CellStructure::AoSoA_pack const &aosoa_)
      : bonded_ias(bonded_ias_), bond_breakage(bond_breakage_),
        coulomb_kernel(coulomb_kernel_), box_geo(box_geo_),
	id_to_index(id_to_index_),
        local_force(local_force_),
#ifdef ESPRESSO_NPT
        local_virial(local_virial_),
#endif
        bond_list(bond_list_), bond_ids(bond_ids_), aosoa(aosoa_) {
  }

  ESPRESSO_ATTR_ALWAYS_INLINE KOKKOS_INLINE_FUNCTION bool
  check_breakage(Kokkos::View<int *> const &partners, int const bond_id) const {
    auto const i = id_to_index(partners(0));
    // Consider for bond breakage
    if (partners(2) == -1) { // pair bonds
      auto const j = id_to_index(partners(1));
      auto const pos1 = aosoa.get_vector_at(aosoa.position, i);
      auto const pos2 = aosoa.get_vector_at(aosoa.position, j);
      auto d = box_geo.get_mi_vector(pos1, pos2).norm();
      if (bond_breakage.check_and_handle_breakage(
              aosoa.id(i), {{aosoa.id(j), std::nullopt}}, bond_id, d)) {
        return true;
      }
    } else if (partners(3) == -1) { // angle bond
      auto const j = id_to_index(partners(1));
      auto const k = id_to_index(partners(2));
      auto const pos2 = aosoa.get_vector_at(aosoa.position, j);
      auto const pos3 = aosoa.get_vector_at(aosoa.position, k);
      auto d = box_geo.get_mi_vector(pos2, pos3).norm();
      if (bond_breakage.check_and_handle_breakage(
              aosoa.id(i), {{aosoa.id(j), aosoa.id(k)}}, bond_id, d)) {
        return true;
      }
    }
    return false;
  }

  ESPRESSO_ATTR_ALWAYS_INLINE KOKKOS_INLINE_FUNCTION bool
  calculate_bond_forces(Kokkos::View<int *> const &partners,
                        int const bond_id) const {
    auto const &iaparams = *bonded_ias.at(bond_id);
    auto const thread_id = omp_get_thread_num();
    auto const i = id_to_index(partners(0));

    switch (number_of_partners(iaparams)) {

    case 1: {
      auto const j = id_to_index(partners(1));
      auto const pos1 = aosoa.get_vector_at(aosoa.position, i);
      auto const pos2 = aosoa.get_vector_at(aosoa.position, j);
      auto const vel1 = aosoa.get_vector_at(aosoa.velocity, i);
      auto const vel2 = aosoa.get_vector_at(aosoa.velocity, j);
      auto const dx = box_geo.get_mi_vector(pos1, pos2);

      if (auto const *iap = std::get_if<ThermalizedBond>(&iaparams)) {
        auto result = iap->forces(aosoa.mass(i), aosoa.mass(j), vel1, vel2,
                                  aosoa.id(i), aosoa.id(j), dx);
        if (result) {
          auto const &forces = result.value();

          local_force(i, thread_id, 0) += std::get<0>(forces)[0];
          local_force(i, thread_id, 1) += std::get<0>(forces)[1];
          local_force(i, thread_id, 2) += std::get<0>(forces)[2];
          local_force(j, thread_id, 0) += std::get<1>(forces)[0];
          local_force(j, thread_id, 1) += std::get<1>(forces)[1];
          local_force(j, thread_id, 2) += std::get<1>(forces)[2];

          return false;
        }
      } else {
        auto result = calc_bond_pair_force(iaparams, dx
#ifdef ESPRESSO_ELECTROSTATICS
                                           ,
                                           aosoa.charge(i), aosoa.charge(j),
                                           coulomb_kernel
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
          return false;
        }
      }
      return true;
    }

    case 2: {
      auto const j = id_to_index(partners(1));
      auto const k = id_to_index(partners(2));

      if (std::get_if<OifGlobalForcesBond>(&iaparams)) {
        return false;
      }
      auto const pos1 = aosoa.get_vector_at(aosoa.position, i);
      auto const pos2 = aosoa.get_vector_at(aosoa.position, j);
      auto const pos3 = aosoa.get_vector_at(aosoa.position, k);
      auto const vec1 = box_geo.get_mi_vector(pos2, pos1);
      auto const vec2 = box_geo.get_mi_vector(pos3, pos1);
      auto const result = calc_bonded_three_body_force(iaparams, vec1, vec2);
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

        return false;
      }
      return true;
    }
    case 3: {
      auto const j = id_to_index(partners(1));
      auto const k = id_to_index(partners(2));
      auto const m = id_to_index(partners(3));
      auto const pos1 = aosoa.get_vector_at(aosoa.position, i);
      auto const pos2 = aosoa.get_vector_at(aosoa.position, j);
      auto const pos3 = aosoa.get_vector_at(aosoa.position, k);
      auto const pos4 = aosoa.get_vector_at(aosoa.position, m);

      std::optional<std::tuple<Utils::Vector3d, Utils::Vector3d,
                               Utils::Vector3d, Utils::Vector3d>> result;
      if (auto const *iap = std::get_if<OifLocalForcesBond>(&iaparams)) {
	auto const fp2 = box_geo.unfolded_position(pos1, aosoa.get_vector_at(aosoa.image, i));
	auto const fp1 = fp2 + box_geo.get_mi_vector(pos2, fp2);
	auto const fp3 = fp2 + box_geo.get_mi_vector(pos3, fp2);
	auto const fp4 = fp2 + box_geo.get_mi_vector(pos4, fp2);
	auto const vel2 = aosoa.get_vector_at(aosoa.velocity, i);
	auto const vel3 = aosoa.get_vector_at(aosoa.velocity, k);

	result = iap->calc_forces(fp2, fp1, fp3, fp4, vel2, vel3);
      }
      if (auto const *iap = std::get_if<IBMTribend>(&iaparams)) {
	//result = iap->calc_forces(box_geo, p1, p2, p3, p4);
	result = iap->calc_forces(box_geo, pos1, pos2, pos3, pos4);
      }
      // note: particles in a dihedral bond are ordered as p2-p1-p3-p4
      auto const v12 = box_geo.get_mi_vector(pos1, pos2);
      auto const v23 = box_geo.get_mi_vector(pos3, pos1);
      auto const v34 = box_geo.get_mi_vector(pos4, pos3);
      if (auto const *iap = std::get_if<DihedralBond>(&iaparams)) {
	result = iap->forces(v12, v23, v34);
      }
#ifdef ESPRESSO_TABULATED
      if (auto const *iap = std::get_if<TabulatedDihedralBond>(&iaparams)) {
	result = iap->forces(v12, v23, v34);
      }
#endif
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

        return false;
      }
      //throw BondUnknownTypeError();
      return true;
    }
    default:
      return false;
    }
  }

  ESPRESSO_ATTR_ALWAYS_INLINE KOKKOS_INLINE_FUNCTION void
  operator()(std::size_t idx) const {
    auto const &partners = Kokkos::subview(bond_list, idx, Kokkos::ALL);
    auto const &bond_id = bond_ids(idx);
    auto bond_broken = false;

    auto breakage = check_breakage(partners, bond_id);
    if (not breakage) {
      bond_broken = calculate_bond_forces(partners, bond_id);
    }

    if (bond_broken) {
      std::span<int> s(partners.data(), partners.extent(0));
      bond_broken_error(s);
    }
  }
};

#endif // ESPRESSO_SHARED_MEMORY_PARALLELISM
