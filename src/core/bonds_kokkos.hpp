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

struct BondsKernel {
  BondedInteractionsMap const &bonded_ias;
  BondBreakage::BondBreakage &bond_breakage;
  Coulomb::ShortRangeForceKernel::kernel_type const *const coulomb_kernel;
  BoxGeometry const &box_geo;
  CellStructure::ForceType const &local_force;
#ifdef ESPRESSO_NPT
  CellStructure::VirialType const &local_virial;
#endif
  CellStructure::BondlistType const &bond_list;
  CellStructure::BondIDType const &bond_ids;
  CellStructure::AoSoA_pack const &aosoa;
  bool const has_breakage_specs;

  BondsKernel(
      BondedInteractionsMap const &bonded_ias_,
      BondBreakage::BondBreakage &bond_breakage_,
      Coulomb::ShortRangeForceKernel::kernel_type const *coulomb_kernel_,
      BoxGeometry const &box_geo_, CellStructure::ForceType const &local_force_,
#ifdef ESPRESSO_NPT
      CellStructure::VirialType const &local_virial_,
#endif
      CellStructure::BondlistType const &bond_list_,
      CellStructure::BondIDType const &bond_ids_,
      CellStructure::AoSoA_pack const &aosoa_)
      : bonded_ias(bonded_ias_), bond_breakage(bond_breakage_),
        coulomb_kernel(coulomb_kernel_), box_geo(box_geo_),
        local_force(local_force_),
#ifdef ESPRESSO_NPT
        local_virial(local_virial_),
#endif
        bond_list(bond_list_), bond_ids(bond_ids_), aosoa(aosoa_),
        has_breakage_specs(!bond_breakage.breakage_specs.empty()) {
  }

  ESPRESSO_ATTR_ALWAYS_INLINE KOKKOS_INLINE_FUNCTION void
  operator()(std::size_t idx) const {
    auto const &partners = Kokkos::subview(bond_list, idx, Kokkos::ALL);
    auto const &bond_id = bond_ids(idx);

    auto const i = partners(0);

    auto const &iaparams = *bonded_ias.at(bond_id);
    // TODO: omp_get_thread_num() is only available for the OpenMP backend.
    // This should be updated when using other Kokkos backends.
    auto const thread_id = omp_get_thread_num();

    switch (number_of_partners(iaparams)) {
    // case 0: zero-partner bonds are implicitly skipped
    case 1: {
      auto const j = partners(1);
      auto const dx =
          box_geo.get_mi_vector(aosoa.get_vector_at(aosoa.position, i),
                                aosoa.get_vector_at(aosoa.position, j));
      std::optional<Utils::Vector3d> result;
      // Consider for bond breakage
      if (has_breakage_specs &&
          bond_breakage.check_and_handle_breakage(
              aosoa.id(i), {{aosoa.id(j), std::nullopt}}, bond_id, dx.norm())) {
        break;
      }
      if (auto const *iap = std::get_if<ThermalizedBond>(&iaparams)) {
        auto const res = iap->forces(aosoa.mass(i), aosoa.mass(j),
                                     aosoa.get_vector_at(aosoa.velocity, i),
                                     aosoa.get_vector_at(aosoa.velocity, j),
                                     aosoa.id(i), aosoa.id(j), dx);
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
        break;
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
      break;
    }

    case 2: {
      auto const j = partners(1);
      auto const k = partners(2);
      auto const pos1 = aosoa.get_vector_at(aosoa.position, i);
      auto const pos2 = aosoa.get_vector_at(aosoa.position, j);
      auto const pos3 = aosoa.get_vector_at(aosoa.position, k);
      auto const vec1 = box_geo.get_mi_vector(pos2, pos1);
      auto const vec2 = box_geo.get_mi_vector(pos3, pos1);
      std::optional<
          std::tuple<Utils::Vector3d, Utils::Vector3d, Utils::Vector3d>>
          result;
      // Consider for bond breakage
      // if (bond_breakage.check_and_handle_breakage(
      if (has_breakage_specs &&
          bond_breakage.check_and_handle_breakage(
              aosoa.id(i), {{aosoa.id(j), aosoa.id(k)}}, bond_id,
              box_geo.get_mi_vector(pos2, pos3).norm())) {
        break;
      }
      if (std::get_if<OifGlobalForcesBond>(&iaparams)) {
        break;
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
        // std::span<int> s(partners.data(), partners.extent(0));
        // bond_broken_error(s);
        std::array<int, 2> pids = {aosoa.id(j), aosoa.id(k)};
        bond_broken_error(aosoa.id(i), {pids.data(), 2});
      }
      break;
    }
    case 3: {
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

      std::optional<std::tuple<Utils::Vector3d, Utils::Vector3d,
                               Utils::Vector3d, Utils::Vector3d>>
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
        // std::span<int> s(partners.data(), partners.extent(0));
        // bond_broken_error(s);
        std::array<int, 3> pids = {aosoa.id(j), aosoa.id(k), aosoa.id(m)};
        bond_broken_error(aosoa.id(i), {pids.data(), 3});
      }
      break;
    }
      // no default: bond_list construction only includes 1-, 2-, and 3-partner
      // bonds
    }
  }
};

#endif // ESPRESSO_SHARED_MEMORY_PARALLELISM
