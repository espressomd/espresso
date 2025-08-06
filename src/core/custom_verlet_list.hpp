/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#ifdef SHARED_MEMORY_PARALLELISM

#include <Cabana_VerletList.hpp>
#include <algorithm>

namespace Cabana {
// ONLY FOR 2D LAYOUT, OTHERWISE NEIGHBOR LIST INTERFACE IMPLEMENTATION WILL
// CAUSE PROBLEMS (NOT IMPLEMENTED)
template <class MemorySpace, class AlgorithmTag, class LayoutTag,
          class BuildTag = TeamVectorOpTag>
class CustomVerletList
    : public VerletList<MemorySpace, AlgorithmTag, LayoutTag, BuildTag> {
public:
  using Base = VerletList<MemorySpace, AlgorithmTag, LayoutTag, BuildTag>;

  // Default constructor
  CustomVerletList() : Base() {}

  // Custom constructor
  CustomVerletList(const std::size_t begin, const std::size_t end,
                   const std::size_t max_neigh) {
    initializeData(end - begin, max_neigh);
  }
  virtual ~CustomVerletList() {};

public:
  Kokkos::View<int *, MemorySpace> counts;
  Kokkos::View<int **, Kokkos::LayoutRight, MemorySpace> neighbors;

  // Method to initialize _data without filling neighbors
  KOKKOS_INLINE_FUNCTION
  void initializeData(const std::size_t num_particles,
                      const std::size_t max_neigh) {
    counts = Kokkos::View<int *, MemorySpace>("num_neighbors", num_particles);
    neighbors = Kokkos::View<int **, Kokkos::LayoutRight, MemorySpace>(
        Kokkos::ViewAllocateWithoutInitializing("neighbors"), num_particles,
        max_neigh);
  }

  // Method to add a neighbor
  KOKKOS_INLINE_FUNCTION
  void addNeighborAtomicLB(int pid, int nid) {
    std::size_t count = counts(pid);
    std::size_t count_n = counts(nid);

    if (count > count_n) {
      int tmp = pid;
      pid = nid;
      nid = tmp;
    }
    count = Kokkos::atomic_fetch_add(&counts(pid), 1);
#ifndef NDEBUG
    if (count >= neighbors.extent(1)) {
      throw std::runtime_error(
          "Number of count is larger than VerletList size.");
    }
#endif
    neighbors(pid, count) = nid;
  }

  // Thread safe but non atomic method to add a neighbor
  KOKKOS_INLINE_FUNCTION
  void addNeighbor(int pid, int nid) {
    std::size_t count = counts(pid);

#ifndef NDEBUG
    if (count >= neighbors.extent(1)) {
      throw std::runtime_error(
          "Number of count is larger than VerletList size.");
    }
#endif
    neighbors(pid, count) = nid;
    counts(pid) += 1;
  }

  // Non atomic and load balancing method to add a neighbor
  KOKKOS_INLINE_FUNCTION
  void addNeighborLB(int pid, int nid) {
    std::size_t count = counts(pid);
    std::size_t count_n = counts(nid);

    if (count > count_n) {
      int tmp = pid;
      pid = nid;
      nid = tmp;
      count = counts(pid);
    }
#ifndef NDEBUG
    if (count >= neighbors.extent(1)) {
      throw std::runtime_error(
          "Number of count is larger than VerletList size.");
    }
#endif
    neighbors(pid, count) = nid;
    counts(pid) += 1;
  }

  // Sorting a neighbor
  KOKKOS_INLINE_FUNCTION
  void sortNeighbors() {
    Kokkos::parallel_for(
        "custom_velet_list::sort_neighbors",
        Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, counts.size()),
        [&](const int i) {
          const int count = counts(i);
          int *ptr = &neighbors(i, 0);
          std::sort(ptr, ptr + count);
        });
    Kokkos::fence();
  }

  // Find max counts
  KOKKOS_INLINE_FUNCTION
  std::size_t get_variance_max_counts() {
    std::size_t max_counts = 0;
    std::size_t ave_counts = 0;
    std::size_t ave_sq_counts = 0;
    for (int pid = 0; pid < counts.extent(0); ++pid) {
      std::size_t count = counts(pid);
      if (max_counts < count)
        max_counts = count;
      ave_counts += count;
      ave_sq_counts += count * count;
    }
    if (counts.extent(0) != 0) {
      ave_counts /= counts.extent(0);
      ave_sq_counts /= counts.extent(0);
      ave_sq_counts -= ave_counts * ave_counts;
      std::cout << "max:" << max_counts << " ave:" << ave_counts
                << " var:" << ave_sq_counts << std::endl;
    }
    return max_counts;
  }

  KOKKOS_INLINE_FUNCTION
  std::size_t get_max_counts() {
    int max;
    Kokkos::Max<int> max_reduce(max);
    Kokkos::parallel_reduce(
        "custom_velet_list::reduce_max",
        Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, counts.size()),
        [&](const int i, int &value) {
          if (counts(i) > value)
            value = counts(i);
        },
        max_reduce);
    Kokkos::fence();
    return static_cast<std::size_t>(max);
  }
};

template <class MemorySpace, class AlgorithmTag, class BuildTag>
class NeighborList<
    CustomVerletList<MemorySpace, AlgorithmTag, VerletLayout2D, BuildTag>> {
public:
  //! Kokkos memory space.
  using memory_space = MemorySpace;
  //! Neighbor list type.
  using list_type =
      CustomVerletList<MemorySpace, AlgorithmTag, VerletLayout2D, BuildTag>;

  //! Get the total number of neighbors across all particles.
  KOKKOS_INLINE_FUNCTION
  static std::size_t totalNeighbor(const list_type &list) {
    std::size_t total_n = 0;
    std::size_t num_p = list.counts.size();
    for (std::size_t i = 0; i < num_p; ++i)
      total_n += list.counts(i);
    return total_n;
  }

  //! Get the maximum number of neighbors per particle.
  KOKKOS_INLINE_FUNCTION
  static std::size_t maxNeighbor(const list_type &list) {
    // Stored during neighbor search.
    return list.max_n;
  }

  //! Get the number of neighbors for a given particle index.
  KOKKOS_INLINE_FUNCTION
  static std::size_t numNeighbor(const list_type &list,
                                 const std::size_t particle_index) {
    return list.counts(particle_index);
  }

  //! Get the id for a neighbor for a given particle index and the index of
  //! the neighbor relative to the particle.
  KOKKOS_INLINE_FUNCTION
  static std::size_t getNeighbor(const list_type &list,
                                 const std::size_t particle_index,
                                 const std::size_t count) {
    return list.neighbors(particle_index, count);
  }
};

} // namespace Cabana

#endif
