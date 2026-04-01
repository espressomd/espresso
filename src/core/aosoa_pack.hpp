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

#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM

#include "cell_system/CellStructure.hpp"

#include <Kokkos_Core.hpp>

#include <omp.h>

#include <cstdint>
#include <span>

struct CellStructure::AoSoA_pack {
  using PositionViewType =
      Kokkos::View<double *[3], Kokkos::LayoutRight, Kokkos::HostSpace>;
  using VelocityViewType =
      Kokkos::View<double *[3], Kokkos::LayoutRight, Kokkos::HostSpace>;
  using DirectorViewType =
      Kokkos::View<double *[3], Kokkos::LayoutRight, Kokkos::HostSpace>;
  using ImageViewType =
      Kokkos::View<int *[3], Kokkos::LayoutRight, Kokkos::HostSpace>;
  using ChargeViewType = Kokkos::View<double *, Kokkos::HostSpace>;
  using DipmViewType = Kokkos::View<double *, Kokkos::HostSpace>;
  using IdViewType = Kokkos::View<int *, Kokkos::HostSpace>;
  using TypeViewType = Kokkos::View<int *, Kokkos::HostSpace>;
  using MassViewType = Kokkos::View<double *, Kokkos::HostSpace>;
  using FlagsViewType = Kokkos::View<uint8_t *, Kokkos::HostSpace>;

  PositionViewType position;
  VelocityViewType velocity;
  DirectorViewType director;
  ImageViewType image;
  ChargeViewType charge;
  DipmViewType dipm;
  IdViewType id;
  TypeViewType type;
  MassViewType mass;
  FlagsViewType flags;

  AoSoA_pack() = default;

  AoSoA_pack(std::size_t num_particles) { resize(num_particles); }

  void resize(std::size_t num_particles) {
    if (position.extent(0) == 0) {
      // First allocation
      position = PositionViewType("position", num_particles);
      image = ImageViewType("image", num_particles);
#ifdef ESPRESSO_ELECTROSTATICS
      charge = ChargeViewType("charge", num_particles);
#endif
      id = IdViewType("id", num_particles);
      type = TypeViewType("type", num_particles);
#ifdef ESPRESSO_MASS
      mass = MassViewType("mass", num_particles);
#endif
      flags = FlagsViewType("flags", num_particles);
      velocity = PositionViewType("velocity", num_particles);
#if defined(ESPRESSO_GAY_BERNE) or defined(ESPRESSO_DIPOLES)
      director = DirectorViewType("director", num_particles);
#endif
#ifdef ESPRESSO_DIPOLES
      dipm = DipmViewType("dipm", num_particles);
#endif
    } else {
      // Reallocation
      Kokkos::realloc(position, num_particles);
      Kokkos::realloc(image, num_particles);
#ifdef ESPRESSO_ELECTROSTATICS
      Kokkos::realloc(charge, num_particles);
#endif
      Kokkos::realloc(id, num_particles);
      Kokkos::realloc(type, num_particles);
#ifdef ESPRESSO_MASS
      Kokkos::realloc(mass, num_particles);
#endif
      Kokkos::realloc(flags, num_particles);
      Kokkos::realloc(velocity, num_particles);
#if defined(ESPRESSO_GAY_BERNE) or defined(ESPRESSO_DIPOLES)
      Kokkos::realloc(director, num_particles);
#endif
#ifdef ESPRESSO_DIPOLES
      Kokkos::realloc(dipm, num_particles);
#endif
    }
  }

  template <typename array_layout, typename T, std::size_t N>
  std::span<T, N>
  get_span_at(Kokkos::View<T *[N], array_layout, Kokkos::HostSpace> const &view,
              std::size_t i) const {
    return std::span<T, N>(const_cast<T *>(&view(i, 0)), N);
  }

  template <typename array_layout, typename T, std::size_t N>
  Utils::Vector<T, N> get_vector_at(
      Kokkos::View<T *[N], array_layout, Kokkos::HostSpace> const &view,
      std::size_t i) const {
    Utils::Vector<T, N> result;
    auto const data = result.data();
#if (defined(__GNUC__) or defined(__GNUG__)) && !defined(__clang__)
#pragma GCC unroll 8
#else
#pragma omp unroll
#endif
    for (std::size_t j = 0ul; j < N; j += 1ul) {
      data[j] = view(i, j);
    }
    return result;
  }

  template <typename array_layout, typename T, std::size_t N>
  void
  set_vector_at(Kokkos::View<T *[N], array_layout, Kokkos::HostSpace> &view,
                std::size_t i, Utils::Vector<T, N> const &value) {
#if (defined(__GNUC__) or defined(__GNUG__)) && !defined(__clang__)
#pragma GCC unroll 8
#else
#pragma omp unroll
#endif
    for (std::size_t j = 0ul; j < N; j += 1ul) {
      view(i, j) = value[j];
    }
  }

  void set_has_exclusion(std::size_t i, bool value) {
    flags(i) = value ? uint8_t{1} : uint8_t{0};
  }

  bool has_exclusion(std::size_t i) const { return flags(i) == uint8_t{1}; }
};

#endif // ESPRESSO_SHARED_MEMORY_PARALLELISM
