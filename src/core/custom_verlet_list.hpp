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

#include <Cabana_VerletList.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>

template <class MemorySpace, class ListAlgorithm, class Layout, class BuildTag>
class CustomVerletList;
namespace Cabana {
template <class MemorySpace, class AlgorithmTag, class BuildTag>
class NeighborList<
    CustomVerletList<MemorySpace, AlgorithmTag, VerletLayout2D, BuildTag>>;
}

// ONLY FOR 2D LAYOUT, OTHERWISE NEIGHBOR LIST INTERFACE IMPLEMENTATION WILL
// CAUSE PROBLEMS (NOT IMPLEMENTED)
template <class MemorySpace, class AlgorithmTag, class LayoutTag,
          class BuildTag = Cabana::TeamVectorOpTag>
class CustomVerletList : public Cabana::VerletList<MemorySpace, AlgorithmTag,
                                                   LayoutTag, BuildTag> {
public:
  CustomVerletList() = default;
  CustomVerletList(std::size_t const begin, std::size_t const end,
                   std::size_t const max_neigh) {
    initializeData(end - begin, max_neigh);
  }

  Kokkos::View<int *, MemorySpace> counts;
  Kokkos::View<int **, Kokkos::LayoutRight, MemorySpace> neighbors;

  // Note: writing to 'overflow' from multiple threads by 'setOverflow()'
  // without synchronization can lead to a data race (unspecified behavior).
  // Since the same value is written from multiple threads concurrently,
  // this should not affect program behavior.
  // https://www.openmp.org/spec-html/5.0/openmpsu9.html
  bool overflow = false;

  // Method to initialize _data without filling neighbors
  inline void initializeData(std::size_t const num_particles,
                             std::size_t const max_neigh) {
    counts = Kokkos::View<int *, MemorySpace>("num_neighbors", num_particles);
    neighbors = Kokkos::View<int **, Kokkos::LayoutRight, MemorySpace>(
        Kokkos::view_alloc(MemorySpace{}, Kokkos::WithoutInitializing,
                           "neighbors"),
        num_particles, max_neigh);
  }

  // Method to realloc _data
  inline void reallocData(std::size_t const num_particles,
                          std::size_t const max_neigh) {
    Kokkos::realloc(counts, num_particles);
    Kokkos::realloc(Kokkos::WithoutInitializing, neighbors, num_particles,
                    max_neigh);
  }

  // Method to add a neighbor
  KOKKOS_INLINE_FUNCTION
  void addNeighborAtomicLB(int pid, int nid) {
    auto count = counts(pid);
    auto count_n = counts(nid);

    if (count > count_n) {
      std::swap(pid, nid);
    }
    count = Kokkos::atomic_fetch_add(&counts(pid), 1);
    auto overflow = count >= neighbors.extent(1);
    if (overflow) {
      setOverflow();
    } else {
      neighbors(pid, count) = nid;
    }
  }

  // Thread-safe but non-atomic method to add a neighbor
  KOKKOS_INLINE_FUNCTION
  void addNeighbor(int pid, int nid) {
    auto const count = counts(pid);

    auto overflow = count >= neighbors.extent(1);
    if (overflow) {
      setOverflow();
    } else {
      neighbors(pid, count) = nid;
      counts(pid) += 1;
    }
  }

  // Non-atomic and load-balanced method to add a neighbor
  KOKKOS_INLINE_FUNCTION
  void addNeighborLB(int pid, int nid) {
    auto count = counts(pid);
    auto count_n = counts(nid);

    if (count > count_n) {
      std::swap(pid, nid);
      count = counts(pid);
    }
    auto overflow = count >= neighbors.extent(1);
    if (overflow) {
      setOverflow();
    } else {
      neighbors(pid, count) = nid;
      counts(pid) += 1;
    }
  }

  // Sorting a neighbor
  KOKKOS_INLINE_FUNCTION
  void sortNeighbors() {
    Kokkos::parallel_for("custom_verlet_list::sort_neighbors",
                         Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(
                             std::size_t{0}, counts.size()),
                         [&](std::size_t const i) {
                           auto const count = counts(i);
                           auto *ptr = &neighbors(i, 0);
                           std::sort(ptr, ptr + count);
                         });
    Kokkos::fence();
  }

  // Find max counts
  KOKKOS_INLINE_FUNCTION
  auto get_variance_max_counts(auto &ostream) {
    auto max_counts = 0l;
    auto ave_counts = 0l;
    auto ave_sq_counts = 0l;
    for (int pid = 0; pid < counts.extent(0); ++pid) {
      auto const count = static_cast<long>(counts(pid));
      if (max_counts < count)
        max_counts = count;
      ave_counts += count;
      ave_sq_counts += count * count;
    }
    if (counts.extent(0) != 0) {
      ave_counts /= static_cast<long>(counts.extent(0));
      ave_sq_counts /= static_cast<long>(counts.extent(0));
      ave_sq_counts -= ave_counts * ave_counts;
    }
    ostream << "max:" << max_counts << " ave:" << ave_counts
            << " var:" << ave_sq_counts << std::endl;
    return static_cast<int>(max_counts);
  }

  KOKKOS_INLINE_FUNCTION
  auto get_max_counts() {
    int max_counts;
    Kokkos::Max<int> max_reduce(max_counts);
    Kokkos::parallel_reduce(
        "custom_verlet_list::reduce_max",
        Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(std::size_t{0},
                                                               counts.size()),
        [&](std::size_t const i, int &value) {
          if (counts(i) > value)
            value = counts(i);
        },
        max_reduce);
    Kokkos::fence();
    return max_counts;
  }

  KOKKOS_INLINE_FUNCTION bool hasOverflow() const { return overflow; }

private:
  KOKKOS_INLINE_FUNCTION void setOverflow() { overflow = true; }
};

template <class MemorySpace, class AlgorithmTag, class BuildTag>
class Cabana::NeighborList<CustomVerletList<MemorySpace, AlgorithmTag,
                                            Cabana::VerletLayout2D, BuildTag>> {
public:
  //! Kokkos memory space.
  using memory_space = MemorySpace;
  //! Neighbor list type.
  using list_type = CustomVerletList<MemorySpace, AlgorithmTag,
                                     Cabana::VerletLayout2D, BuildTag>;

  //! Get the total number of neighbors across all particles.
  KOKKOS_INLINE_FUNCTION
  static std::size_t totalNeighbor(list_type const &list) {
    std::size_t const num_p = list.counts.size();
    std::size_t total_n = 0;
    for (std::size_t i = 0; i < num_p; ++i)
      total_n += list.counts(i);
    return total_n;
  }

  //! Get the maximum number of neighbors per particle.
  KOKKOS_INLINE_FUNCTION
  static std::size_t maxNeighbor(list_type const &list) {
    return list.neighbors.extent(1);
  }

  //! Get the number of neighbors for a given particle index.
  KOKKOS_INLINE_FUNCTION
  static std::size_t numNeighbor(list_type const &list,
                                 std::size_t const particle_index) {
    return list.counts(particle_index);
  }

  //! Get the id for a neighbor for a given particle index and the index of
  //! the neighbor relative to the particle.
  KOKKOS_INLINE_FUNCTION
  static std::size_t getNeighbor(list_type const &list,
                                 std::size_t const particle_index,
                                 std::size_t const count) {
    return list.neighbors(particle_index, count);
  }
};
