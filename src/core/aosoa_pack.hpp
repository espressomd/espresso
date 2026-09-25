/*
 * Copyright (C) 2025-2026 The ESPResSo project
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

#include "cell_system/CellStructure.hpp"

#include <utils/attributes.hpp>
#include <utils/device_qualifier.hpp>

#include <Kokkos_Core.hpp>

#include <omp.h>

#include <atomic>
#include <cstdint>
#include <span>

struct CellStructure::AoSoA_pack {
  // Particle-major (@ref ParticleStore::StateVectorLayout, LayoutRight) to
  // match the ParticleStore columns these views alias / are derived from.
  //
  // position/image/director/velocity/id/mass alias the authoritative
  // ParticleStore host columns and are indexed by *store row*, obtained by
  // translating a pack index i through row(i). For i < n_local the translation
  // is the identity (both are built in cell-traversal order); only the deduped
  // ghost tail is remapped.
  //
  // `type` is a PACK-OWNED contiguous array, written at pack-rebuild time in
  // commit_particle (`type(index)=p.type()`) and read pack-indexed by the hot
  // pair kernels (forces/energy/pressure_cabana). The ParticleStore type column
  // REMAINS authoritative; this pack copy is a derived cache refreshed on
  // rebuild (a mid-run type change goes through System::on_particle_change,
  // which calls set_resort_particles).
  //
  // Likewise `charge`/`dipm` are pack-owned contiguous arrays, refreshed per
  // step ONLY when a coulomb (resp. dipolar) actor is active (see
  // refresh_pack_charges / refresh_pack_dipm). The hot pair kernels and the P3M
  // gather/spread loops read them PACK-INDEXED (contiguous). Pure-LJ runs never
  // touch them, paying zero cost. The store q/dipm columns stay authoritative.
  //
  // The layout of the store-aliased vector views (position/image/director/
  // velocity) MUST match the store columns' layout.
  using PositionViewType =
      Kokkos::View<double *[3], ParticleStore::StateVectorLayout,
                   Kokkos::HostSpace>;
  using VelocityViewType =
      Kokkos::View<double *[3], ParticleStore::StateVectorLayout,
                   Kokkos::HostSpace>;
  using DirectorViewType =
      Kokkos::View<double *[3], ParticleStore::StateVectorLayout,
                   Kokkos::HostSpace>;
  using ImageViewType = Kokkos::View<int *[3], ParticleStore::StateVectorLayout,
                                     Kokkos::HostSpace>;
  using RowMapViewType = Kokkos::View<int const *, Kokkos::HostSpace>;
  using ChargeViewType = Kokkos::View<double *, Kokkos::HostSpace>;
  using DipmViewType = Kokkos::View<double *, Kokkos::HostSpace>;
  using IdViewType = Kokkos::View<int *, Kokkos::HostSpace>;
  using TypeViewType = Kokkos::View<int *, Kokkos::HostSpace>;
  using MassViewType = Kokkos::View<double *, Kokkos::HostSpace>;
  using FlagsViewType = Kokkos::View<uint8_t *, Kokkos::HostSpace>;

  // Store-aliased / store-derived columns (indexed by store row).
  PositionViewType position;
  ImageViewType image;
  DirectorViewType director;
  VelocityViewType velocity;
  IdViewType id;
  MassViewType mass;
  // `charge` aliases the authoritative ParticleStore q column and is read by
  // *store row* on the cold BOND path (BondedCoulomb needs charges even with NO
  // coulomb solver, so it must always be valid — hence it aliases the store
  // rather than the guarded pack-owned pair_charge below).
  ChargeViewType charge;
  // Pack-index -> store-row translation. Identity on the local prefix.
  RowMapViewType row_map;

  // Pack-owned columns (indexed by pack index).
  // `type` is written on rebuild (read by the hot pair kernels). `pair_charge`/
  // `pair_dipm` are the contiguous hot-path charge/dipm columns, refreshed per
  // step ONLY when the respective solver is active (see refresh_pack_charges /
  // refresh_pack_dipm); the hot pair kernels and P3M gather/spread read them
  // pack-indexed. `flags` is written every commit and is the allocation
  // sentinel (type/pair_charge/pair_dipm/flags are resized together).
  TypeViewType type;
  ChargeViewType pair_charge;
  DipmViewType pair_dipm;
  FlagsViewType flags;

  AoSoA_pack() = default;

  AoSoA_pack(std::size_t num_particles) { resize(num_particles); }

  /** @brief Translate a pack index to its ParticleStore row. */
  ESPRESSO_ATTR_ALWAYS_INLINE KOKKOS_INLINE_FUNCTION int
  row(std::size_t i) const {
    return row_map(i);
  }

  HOST_ONLY_QUALIFIER
  void resize(std::size_t num_particles) {
    // id/mass/charge are store-aliased (bound in bind_pack_store_views), not
    // allocated here. `type`/`pair_charge`/`pair_dipm`/`flags` are pack-owned
    // contiguous columns allocated here; `flags` serves as the allocation
    // sentinel (all four are resized together).
    if (flags.extent(0) == 0) {
      // First allocation
      type =
          TypeViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing, "type"),
                       num_particles);
      pair_charge = ChargeViewType(
          Kokkos::view_alloc(Kokkos::WithoutInitializing, "pair_charge"),
          num_particles);
      pair_dipm = DipmViewType(
          Kokkos::view_alloc(Kokkos::WithoutInitializing, "pair_dipm"),
          num_particles);
      flags = FlagsViewType("flags", num_particles);
    } else {
      // Reallocation
      Kokkos::realloc(Kokkos::WithoutInitializing, type, num_particles);
      Kokkos::realloc(Kokkos::WithoutInitializing, pair_charge, num_particles);
      Kokkos::realloc(Kokkos::WithoutInitializing, pair_dipm, num_particles);
      Kokkos::realloc(flags, num_particles);
    }
  }

  template <typename array_layout, typename T, std::size_t N>
  DEVICE_QUALIFIER std::span<T, N>
  get_span_at(Kokkos::View<T *[N], array_layout, Kokkos::HostSpace> const &view,
              std::size_t i) const {
    return std::span<T, N>(const_cast<T *>(&view(i, 0)), N);
  }

  template <typename array_layout, typename T, std::size_t N>
  DEVICE_QUALIFIER Utils::Vector<T, N> get_vector_at(
      Kokkos::View<T *[N], array_layout, Kokkos::HostSpace> const &view,
      std::size_t i) const {
    Utils::Vector<T, N> result;
    auto const data = result.data();
#if !defined(__NVCOMPILER) && !defined(__CUDACC__)
#if defined(__clang__)
#pragma omp unroll
#elif defined(__GNUC__) or defined(__GNUG__)
#pragma GCC unroll 8
#endif
#endif
    for (std::size_t j = 0ul; j < N; j += 1ul) {
      data[j] = view(i, j);
    }
    return result;
  }

  template <typename array_layout, typename T, std::size_t N>
  DEVICE_QUALIFIER void
  set_vector_at(Kokkos::View<T *[N], array_layout, Kokkos::HostSpace> &view,
                std::size_t i, Utils::Vector<T, N> const &value) {
#if !defined(__NVCOMPILER) && !defined(__CUDACC__)
#if defined(__clang__)
#pragma omp unroll
#elif defined(__GNUC__) or defined(__GNUG__)
#pragma GCC unroll 8
#endif
#endif
    for (std::size_t j = 0ul; j < N; j += 1ul) {
      view(i, j) = value[j];
    }
  }

  // Aggregate of the per-particle exclusion flag over the whole pack (local +
  // ghost). It lets the specialized-kernel dispatch decide in O(1) whether any
  // packed particle carries an exclusion, instead of sweeping every flag on
  // every force call. It is only ever transitioned false->true during a commit
  // sweep (which runs under Kokkos::parallel_for on the host execution space),
  // so a relaxed atomic store suffices; reads happen after the sweep's
  // Kokkos::fence() and are therefore race-free. The atomic member makes
  // AoSoA_pack non-copyable and non-movable, which is fine: it is only ever
  // default-constructed via std::make_unique and accessed by reference (see
  // CellStructure::m_aosoa).
  std::atomic<bool> any_exclusion{false};

  void reset_any_exclusion() {
    any_exclusion.store(false, std::memory_order_relaxed);
  }

  bool has_any_exclusion() const {
    return any_exclusion.load(std::memory_order_relaxed);
  }

  // Host-only: record that some particle carries an exclusion. Kept out of the
  // device-qualified set_has_exclusion so the atomic never appears in device
  // code; commit_particle calls this from the host commit sweep. Test before
  // setting: after the first write the flag's cache line stays shared, so
  // exclusion-heavy sweeps don't ping-pong it across threads.
  void mark_any_exclusion() {
    if (not any_exclusion.load(std::memory_order_relaxed)) {
      any_exclusion.store(true, std::memory_order_relaxed);
    }
  }

  DEVICE_QUALIFIER void set_has_exclusion(std::size_t i, bool value) {
    flags(i) = value ? uint8_t{1} : uint8_t{0};
  }

  DEVICE_QUALIFIER bool has_exclusion(std::size_t i) const {
    return flags(i) == uint8_t{1};
  }

#ifdef ESPRESSO_ELECTROSTATICS
  /**
   * @brief Debug freshness check for the pack-owned @ref pair_charge column.
   *
   * The hot Coulomb pair kernels and the P3M real-space charge gather read
   * @ref pair_charge PACK-INDEXED, trusting that it was resynced from the
   * authoritative ParticleStore q column this step (via refresh_pack_charges).
   * The store-aliased @ref charge view (indexed by store row) IS the current
   * authoritative charge, so a fresh @ref pair_charge must equal it
   * element-for-element. O(n_part) and debug-only (ESPRESSO_ADDITIONAL_CHECKS),
   * i.e. negligible vs the O(mesh)/O(pairs) P3M work. @p n_part is the local
   * particle count the caller is about to gather over.
   */
  void assert_pair_charge_fresh([[maybe_unused]] std::size_t n_part) const {
#ifdef ESPRESSO_ADDITIONAL_CHECKS
    for (std::size_t i = 0ul; i < n_part; ++i) {
      assert(pair_charge(i) == charge(row(i)) &&
             "P3M/Coulomb gather read a STALE pack charge column: "
             "refresh_pack_charges was not called after the store q column "
             "last changed");
    }
#endif
  }
#endif
};
