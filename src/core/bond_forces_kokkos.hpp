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

// Pair bonds are dispatched particle-parallel (see PairBondsKernel below):
// each work-item writes only to its own row of a plain force View, so no
// ScatterView is needed there. The NPT virial fold-in is still a genuine
// cross-particle reduction (many particles' primary pair bonds add into the
// same 3 global components), so it keeps using a ScatterView.
struct PairBondsKernelData {
  BondedInteractionsMap const &bonded_ias;
  BondBreakage::BondBreakage &bond_breakage;
  BoxGeometry const &box_geo;
  CellStructure::ForceType local_force;
#ifdef ESPRESSO_NPT
  CellStructure::ScatterVirial local_virial;
#endif
  CellStructure::AoSoA_pack const &aosoa;
  bool const has_breakage_specs;
};

// Particle-parallel gather kernel: one work-item per LOCAL particle
// ("self"), which loops over its own row of the per-particle pair-bond
// structure (built by CellStructure::update_bond_storage(), one slot per
// BondList entry -- primary or mirror -- the particle holds) and adds only
// its own force contribution to its own row of a plain force View. Pair
// force is antisymmetric (force(-dx) == -force(dx)), so a mirror-holding
// participant gets the correct sign automatically from its own
// self-referenced @c dx, with no role bookkeeping needed for the force
// itself. Bond breakage and the NPT virial fold-in are bond-level (not
// per-participant) quantities and are only evaluated by the primary-holding
// participant, to avoid counting them twice.
struct PairBondsKernel {
  PairBondsKernelData data;
  LocalBondState::PPPairDegreeType pp_pair_degree;
  LocalBondState::PPPairSlotType pp_pair_slots;
  Coulomb::ShortRangeForceKernel::kernel_type const *const coulomb_kernel;

  PairBondsKernel(
      PairBondsKernelData data_,
      LocalBondState::PPPairDegreeType pp_pair_degree_,
      LocalBondState::PPPairSlotType pp_pair_slots_,
      Coulomb::ShortRangeForceKernel::kernel_type const *coulomb_kernel_)
      : data(std::move(data_)), pp_pair_degree(std::move(pp_pair_degree_)),
        pp_pair_slots(std::move(pp_pair_slots_)),
        coulomb_kernel(coulomb_kernel_) {}

  ESPRESSO_ATTR_ALWAYS_INLINE inline void operator()(std::size_t self) const {
    auto const &bonded_ias = data.bonded_ias;
    auto const &box_geo = data.box_geo;
    auto &local_force = data.local_force;
    auto const &aosoa = data.aosoa;
    auto &bond_breakage = data.bond_breakage;
#ifdef ESPRESSO_NPT
    auto local_virial = data.local_virial.access();
#endif
    auto const has_breakage_specs = data.has_breakage_specs;
    auto const degree = pp_pair_degree(self);
    Utils::Vector3d total_force{};

    for (int slot = 0; slot < degree; ++slot) {
      auto const other = pp_pair_slots(self, slot, 0);
      // -1 marks a row whose partner failed to resolve on this rank (see
      // CellStructure::update_bond_storage) -- this participant's
      // authoritative home rank always resolves it and applies the force
      // there, so skipping here loses nothing.
      if (other < 0) {
        continue;
      }

      auto const bond_id = pp_pair_slots(self, slot, 1);
      auto const is_primary = pp_pair_slots(self, slot, 2) != 0;
      auto const &iaparams = *bonded_ias.at(bond_id);

      auto const dx =
          box_geo.get_mi_vector(aosoa.get_vector_at(aosoa.position, self),
                                aosoa.get_vector_at(aosoa.position, other));

      // Bond breakage is checked from every participant's own point of
      // view (own id first, own self-referenced distance), not gated on
      // is_primary: the old per-bond kernel gave both participants zero
      // force in the same evaluation where the bond broke, and a mirror
      // that skipped this check would otherwise keep adding a (possibly
      // unbounded, e.g. HarmonicBond) force for one extra evaluation.
      // Both sides may enqueue a breakage record for the same bond; the
      // second ::remove_bond() call in process_queue_impl() is then just
      // a harmless no-op (bond_breakage.cpp, BondBreakage::execute).
      if (has_breakage_specs &&
          bond_breakage.check_and_handle_breakage(
              aosoa.id(self), {{aosoa.id(other), std::nullopt}}, bond_id,
              dx.norm())) {
        continue;
      }

      if (auto const *iap = std::get_if<ThermalizedBond>(&iaparams)) {
        // ThermalizedBond's noise draw is keyed by (id1, id2) through a
        // Philox counter, which is NOT symmetric under swapping the two
        // ids. Every participant must therefore evaluate it with the same
        // canonical (primary-holder first) argument order regardless of
        // which one is "self", then pick out its own slot of the returned
        // (force_on_id1, force_on_id2) pair -- this reproduces exactly the
        // single evaluation the old per-bond kernel used to make.
        auto const canonical_dx = is_primary ? dx : -dx;
        auto const result =
            is_primary ? iap->forces(
#ifdef ESPRESSO_MASS
                             aosoa.mass(self), aosoa.mass(other),
#else
                             1.0, 1.0,
#endif
                             aosoa.get_vector_at(aosoa.velocity, self),
                             aosoa.get_vector_at(aosoa.velocity, other),
                             aosoa.id(self), aosoa.id(other), canonical_dx)
                       : iap->forces(
#ifdef ESPRESSO_MASS
                             aosoa.mass(other), aosoa.mass(self),
#else
                             1.0, 1.0,
#endif
                             aosoa.get_vector_at(aosoa.velocity, other),
                             aosoa.get_vector_at(aosoa.velocity, self),
                             aosoa.id(other), aosoa.id(self), canonical_dx);
        if (result) {
          total_force += is_primary ? std::get<0>(result.value())
                                    : std::get<1>(result.value());
        } else if (is_primary) {
          auto partner_id = aosoa.id(other);
          bond_broken_error(aosoa.id(self), {&partner_id, 1});
        }
        continue;
      }

      auto const result = calc_bond_pair_force(
          iaparams, dx,
#ifdef ESPRESSO_ELECTROSTATICS
          aosoa.charge(self) * aosoa.charge(other), coulomb_kernel
#else
          0.0, nullptr
#endif
      );

      if (result) {
        auto const f = result.value();
        total_force += f;
#ifdef ESPRESSO_NPT
        if (is_primary) {
          auto const virial = hadamard_product(f, dx);
          local_virial(0) += virial[0];
          local_virial(1) += virial[1];
          local_virial(2) += virial[2];
        }
#endif
      } else if (is_primary) {
        auto partner_id = aosoa.id(other);
        bond_broken_error(aosoa.id(self), {&partner_id, 1});
      }
    }
    local_force(self, 0) += total_force[0];
    local_force(self, 1) += total_force[1];
    local_force(self, 2) += total_force[2];
  }
};

// Angle bonds are dispatched particle-parallel, same rationale as
// PairBondsKernel above. No NPT virial fold-in exists for angle bonds today
// (BondsKernelData's local_virial was already unused by the old
// bond-indexed AngleBondsKernel), so this kernel's data has no virial field.
struct AngleBondsKernelData {
  BondedInteractionsMap const &bonded_ias;
  BondBreakage::BondBreakage &bond_breakage;
  BoxGeometry const &box_geo;
  CellStructure::ForceType local_force;
  CellStructure::AoSoA_pack const &aosoa;
  bool const has_breakage_specs;
  CellStructure::ScatterForce scatter_force;
  int pp_num_particles;
};

// Evaluates every rank-locally-owned angle bond's force exactly once
// (dispatched over angle_count, i.e. angle_list/angle_ids -- the same
// compact per-bond list energy/pressure calculation already uses), instead
// of once per participant like AngleBondsKernel below used to do
// unconditionally. This is the compute side of the split: it does the
// actual (trig-heavy) geometry + force evaluation and scatter-adds the 3
// resulting force vectors directly into the (real, non-ghost) participants'
// own local_force rows via a ScatterView, exactly like the pre-particle-
// parallel-gather bond-indexed AngleBondsKernel used to -- so
// AngleBondsKernel's per-participant gather pass can just skip a row
// instead of recomputing it. A ghost participant's row is skipped here (its
// share has nowhere atomics-safe to go without re-introducing the MPI
// ghost-force-reduction this architecture avoids); it is instead applied,
// independently and redundantly, by its own real-owner rank's
// AngleBondsKernel fallback path.
struct AngleBondsForceComputeKernel {
  AngleBondsKernelData data;
  LocalBondState::AngleBondlistType bond_list;
  LocalBondState::AngleBondIDType bond_ids;

  AngleBondsForceComputeKernel(AngleBondsKernelData data_,
                               LocalBondState::AngleBondlistType bond_list_,
                               LocalBondState::AngleBondIDType bond_ids_)
      : data(std::move(data_)), bond_list(std::move(bond_list_)),
        bond_ids(std::move(bond_ids_)) {}

  ESPRESSO_ATTR_ALWAYS_INLINE inline void operator()(std::size_t idx) const {
    auto const &bonded_ias = data.bonded_ias;
    auto const &box_geo = data.box_geo;
    auto const &aosoa = data.aosoa;
    auto &bond_breakage = data.bond_breakage;
    auto const has_breakage_specs = data.has_breakage_specs;
    auto scatter_force = data.scatter_force.access();
    auto const pp_num_particles = data.pp_num_particles;
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

    // Consider for bond breakage. Evaluated exactly once per bond here
    // (unlike AngleBondsKernel's fallback path, which every participant
    // checks independently), matching the pre-gather kernel's behavior.
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
      auto const &f0 = std::get<0>(forces);
      auto const &f1 = std::get<1>(forces);
      auto const &f2 = std::get<2>(forces);
      // i (the vertex) is always real/local by construction (angle_list
      // only ever gets a primary entry from a non-ghost owner -- see
      // CellStructure::update_bond_storage); j/k may be ghosts here.
      scatter_force(i, 0) += f0[0];
      scatter_force(i, 1) += f0[1];
      scatter_force(i, 2) += f0[2];
      if (j < pp_num_particles) {
        scatter_force(j, 0) += f1[0];
        scatter_force(j, 1) += f1[1];
        scatter_force(j, 2) += f1[2];
      }
      if (k < pp_num_particles) {
        scatter_force(k, 0) += f2[0];
        scatter_force(k, 1) += f2[1];
        scatter_force(k, 2) += f2[2];
      }
    } else {
      std::array<int, 2> pids = {aosoa.id(j), aosoa.id(k)};
      bond_broken_error(aosoa.id(i), {pids.data(), 2});
    }
  }
};

// Particle-parallel gather kernel for angle bonds: one work-item per LOCAL
// particle, looping its own row of the per-particle angle-bond structure.
// Every row stores the full original (vertex, arm1, arm2) plus which
// position "self" occupies ("self_slot", derived at build time -- see
// PPAngleSlotType's doc comment): not every 3-body force function reachable
// through calc_bonded_three_body_force shares angle_generic_force's
// f_left(a, b) == f_right(b, a) symmetry (IBMTriel's does not), so unlike an
// earlier version of this kernel, an arm mirror cannot just treat itself as
// the "left" participant -- the geometry is always set up vertex-referenced,
// exactly as the old bond-indexed kernel did, and only the choice of which
// tuple element is "self"'s force depends on self_slot.
//
// Each row's "bond_index" column (see PPAngleSlotType's doc comment) picks
// one of two paths: if AngleBondsForceComputeKernel already evaluated this
// bond (the common case), it already scatter-added this participant's share
// directly into local_force, so this row has nothing left to do -- skip it.
// Otherwise (this bond's primary entry lives on another rank) fall back to
// evaluating it directly from this row, exactly as this kernel always used
// to.
struct AngleBondsKernel {
  AngleBondsKernelData data;
  LocalBondState::PPAngleDegreeType pp_angle_degree;
  LocalBondState::PPAngleSlotType pp_angle_slots;

  AngleBondsKernel(AngleBondsKernelData data_,
                   LocalBondState::PPAngleDegreeType pp_angle_degree_,
                   LocalBondState::PPAngleSlotType pp_angle_slots_)
      : data(std::move(data_)), pp_angle_degree(std::move(pp_angle_degree_)),
        pp_angle_slots(std::move(pp_angle_slots_)) {}

  ESPRESSO_ATTR_ALWAYS_INLINE inline void operator()(std::size_t self) const {
    auto const &bonded_ias = data.bonded_ias;
    auto const &box_geo = data.box_geo;
    auto &local_force = data.local_force;
    auto const &aosoa = data.aosoa;
    auto &bond_breakage = data.bond_breakage;
    auto const has_breakage_specs = data.has_breakage_specs;
    auto const degree = pp_angle_degree(self);
    Utils::Vector3d total_force{};

    for (int slot = 0; slot < degree; ++slot) {
      auto const self_slot = pp_angle_slots(self, slot, 4);
      // -1 marks an unresolvable row -- see PairBondsKernel's comment.
      if (self_slot < 0) {
        continue;
      }

      auto const bond_index = pp_angle_slots(self, slot, 5);
      if (bond_index >= 0) {
        // Already applied directly by AngleBondsForceComputeKernel's
        // scatter-add -- nothing left to do for this row.
        continue;
      }

      // Fallback: this bond's primary entry isn't resolvable on this
      // rank (it straddles a rank boundary and the owner lives
      // elsewhere), so AngleBondsForceComputeKernel never evaluated it --
      // evaluate it directly.
      auto const vertex = pp_angle_slots(self, slot, 0);
      auto const arm1 = pp_angle_slots(self, slot, 1);
      auto const arm2 = pp_angle_slots(self, slot, 2);
      auto const bond_id = pp_angle_slots(self, slot, 3);
      auto const &iaparams = *bonded_ias.at(bond_id);

      if (std::get_if<OifGlobalForcesBond>(&iaparams)) {
        continue;
      }

      auto const pos_vertex = aosoa.get_vector_at(aosoa.position, vertex);
      auto const pos_arm1 = aosoa.get_vector_at(aosoa.position, arm1);
      auto const pos_arm2 = aosoa.get_vector_at(aosoa.position, arm2);

      auto const vec1 = box_geo.get_mi_vector(pos_arm1, pos_vertex);
      auto const vec2 = box_geo.get_mi_vector(pos_arm2, pos_vertex);

      // Bond breakage is checked from every participant's row (all 3
      // independently reconstruct the same vertex/arm1/arm2 geometry and
      // distance, so this is deterministically the same check either
      // way), always reported with the vertex as the primary participant
      // -- REVERT_BIND_AT_POINT_OF_COLLISION in bond_breakage.cpp assumes
      // the reported particle_id is the vertex/virtual-site. Up to 3
      // duplicate queue entries for the same bond are a harmless no-op
      // (see PairBondsKernel's comment on the analogous pair-bond case).
      if (has_breakage_specs &&
          bond_breakage.check_and_handle_breakage(
              aosoa.id(vertex), {{aosoa.id(arm1), aosoa.id(arm2)}}, bond_id,
              box_geo.get_mi_vector(pos_arm1, pos_arm2).norm())) {
        continue;
      }

      auto const result = calc_bonded_three_body_force(iaparams, vec1, vec2);

      if (result) {
        auto const &forces = result.value();
        switch (self_slot) {
        case 0:
          total_force += std::get<0>(forces);
          break;
        case 1:
          total_force += std::get<1>(forces);
          break;
        default:
          total_force += std::get<2>(forces);
          break;
        }
      } else if (self_slot == 0) {
        std::array<int, 2> pids = {aosoa.id(arm1), aosoa.id(arm2)};
        bond_broken_error(aosoa.id(vertex), {pids.data(), 2});
      }
    }
    local_force(self, 0) += total_force[0];
    local_force(self, 1) += total_force[1];
    local_force(self, 2) += total_force[2];
  }
};

// Dihedral bonds are dispatched particle-parallel, same rationale as
// PairBondsKernel/AngleBondsKernel above. Dihedral bonds don't support
// bond_breakage (the old bond-indexed kernel never checked it either, and
// bond_breakage.cpp's bond_handler() only handles 1- and 2-partner bonds),
// so this kernel's data carries no bond_breakage field.
struct DihedralBondsKernelData {
  BondedInteractionsMap const &bonded_ias;
  BoxGeometry const &box_geo;
  CellStructure::ForceType local_force;
  CellStructure::AoSoA_pack const &aosoa;
  CellStructure::ScatterForce scatter_force;
  int pp_num_particles;
};

// Evaluates every rank-locally-owned dihedral bond's force exactly once
// (dispatched over dihedral_count), scatter-adding the 4 resulting force
// vectors directly into the (real, non-ghost) participants' own local_force
// rows -- see AngleBondsForceComputeKernel's doc comment for the rationale;
// same split, applied to the 4-body case.
struct DihedralBondsForceComputeKernel {
  DihedralBondsKernelData data;
  LocalBondState::DihedralBondlistType bond_list;
  LocalBondState::DihedralBondIDType bond_ids;

  DihedralBondsForceComputeKernel(
      DihedralBondsKernelData data_,
      LocalBondState::DihedralBondlistType bond_list_,
      LocalBondState::DihedralBondIDType bond_ids_)
      : data(std::move(data_)), bond_list(std::move(bond_list_)),
        bond_ids(std::move(bond_ids_)) {}

  ESPRESSO_ATTR_ALWAYS_INLINE inline void operator()(std::size_t idx) const {
    auto const &bonded_ias = data.bonded_ias;
    auto const &box_geo = data.box_geo;
    auto const &aosoa = data.aosoa;
    auto scatter_force = data.scatter_force.access();
    auto const pp_num_particles = data.pp_num_particles;
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
      auto const &f0 = std::get<0>(forces);
      auto const &f1 = std::get<1>(forces);
      auto const &f2 = std::get<2>(forces);
      auto const &f3 = std::get<3>(forces);
      // i is always real/local by construction (see
      // AngleBondsForceComputeKernel's comment); j/k/m may be ghosts here.
      scatter_force(i, 0) += f0[0];
      scatter_force(i, 1) += f0[1];
      scatter_force(i, 2) += f0[2];
      if (j < pp_num_particles) {
        scatter_force(j, 0) += f1[0];
        scatter_force(j, 1) += f1[1];
        scatter_force(j, 2) += f1[2];
      }
      if (k < pp_num_particles) {
        scatter_force(k, 0) += f2[0];
        scatter_force(k, 1) += f2[1];
        scatter_force(k, 2) += f2[2];
      }
      if (m < pp_num_particles) {
        scatter_force(m, 0) += f3[0];
        scatter_force(m, 1) += f3[1];
        scatter_force(m, 2) += f3[2];
      }
    } else {
      std::array<int, 3> pids = {aosoa.id(j), aosoa.id(k), aosoa.id(m)};
      bond_broken_error(aosoa.id(i), {pids.data(), 3});
    }
  }
};

// Particle-parallel gather kernel for dihedral bonds: one work-item per
// LOCAL particle, looping its own row of the per-particle dihedral-bond
// structure. Every row stores the full chain [i, j, k, m] in the bond's
// original creation order plus which of the 4 positions this row's own
// particle occupies ("chain_slot", derived at build time -- see
// PPDihedralSlotType's doc comment, since unlike pair/angle bonds the
// dihedral force formula has no self-referenced symmetry to exploit).
// vel1/vel3/image1 always refer to chain positions 0/2/0 respectively,
// regardless of which position "self" is -- calc_bonded_four_body_force
// needs all 4 chain positions' geometry either way, so this is no more
// redundant per participant than the position lookups themselves.
//
// See AngleBondsKernel's doc comment for the bond_index fast-path/fallback
// split: in the common case, DihedralBondsForceComputeKernel already
// scatter-added this participant's share directly into local_force, so this
// row has nothing left to do.
struct DihedralBondsKernel {
  DihedralBondsKernelData data;
  LocalBondState::PPDihedralDegreeType pp_dihedral_degree;
  LocalBondState::PPDihedralSlotType pp_dihedral_slots;

  DihedralBondsKernel(DihedralBondsKernelData data_,
                      LocalBondState::PPDihedralDegreeType pp_dihedral_degree_,
                      LocalBondState::PPDihedralSlotType pp_dihedral_slots_)
      : data(std::move(data_)),
        pp_dihedral_degree(std::move(pp_dihedral_degree_)),
        pp_dihedral_slots(std::move(pp_dihedral_slots_)) {}

  ESPRESSO_ATTR_ALWAYS_INLINE inline void operator()(std::size_t self) const {
    auto const &bonded_ias = data.bonded_ias;
    auto const &box_geo = data.box_geo;
    auto &local_force = data.local_force;
    auto const &aosoa = data.aosoa;
    auto const degree = pp_dihedral_degree(self);
    Utils::Vector3d total_force{};

    for (int slot = 0; slot < degree; ++slot) {
      auto const chain_slot = pp_dihedral_slots(self, slot, 5);
      // -1 marks an unresolvable row -- see PairBondsKernel's comment.
      if (chain_slot < 0) {
        continue;
      }

      auto const bond_index = pp_dihedral_slots(self, slot, 6);
      if (bond_index >= 0) {
        // Already applied directly by DihedralBondsForceComputeKernel's
        // scatter-add -- nothing left to do for this row.
        continue;
      }

      // Fallback: this bond's primary entry isn't resolvable on this
      // rank -- see AngleBondsKernel's fallback comment.
      auto const i = pp_dihedral_slots(self, slot, 0);
      auto const j = pp_dihedral_slots(self, slot, 1);
      auto const k = pp_dihedral_slots(self, slot, 2);
      auto const m = pp_dihedral_slots(self, slot, 3);
      auto const bond_id = pp_dihedral_slots(self, slot, 4);
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
        switch (chain_slot) {
        case 0:
          total_force += std::get<0>(forces);
          break;
        case 1:
          total_force += std::get<1>(forces);
          break;
        case 2:
          total_force += std::get<2>(forces);
          break;
        default:
          total_force += std::get<3>(forces);
          break;
        }
      } else if (chain_slot == 0) {
        std::array<int, 3> pids = {aosoa.id(j), aosoa.id(k), aosoa.id(m)};
        bond_broken_error(aosoa.id(i), {pids.data(), 3});
      }
    }
    local_force(self, 0) += total_force[0];
    local_force(self, 1) += total_force[1];
    local_force(self, 2) += total_force[2];
  }
};
