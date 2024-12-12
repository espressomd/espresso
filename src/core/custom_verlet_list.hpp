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

#ifdef CABANA 

#include <Cabana_VerletList.hpp>

namespace Cabana
{
// ONLY FOR 2D LAYOUT, OTHERWISE NEIGHBOR LIST INTERFACE IMPLEMENTATION WILL CAUSE PROBLEMS (NOT IMPLEMENTED)
template <class MemorySpace, class AlgorithmTag, class LayoutTag, class BuildTag = TeamVectorOpTag>
class CustomVerletList : public VerletList<MemorySpace, AlgorithmTag, LayoutTag, BuildTag>
{
  public:
    using Base = VerletList<MemorySpace, AlgorithmTag, LayoutTag, BuildTag>;

    // Default constructor
    CustomVerletList() : Base() {}

    // Custom constructor
    template <class PositionSlice>
    CustomVerletList(PositionSlice x, const std::size_t begin, const std::size_t end,
                     const std::size_t max_neigh)
    {
        initializeData(x.size(), max_neigh);
    }


public:
    Kokkos::View<int*, MemorySpace> counts;
    Kokkos::View<int**, MemorySpace> neighbors;

    // Method to initialize _data without filling neighbors
    KOKKOS_INLINE_FUNCTION
    void initializeData(const std::size_t num_particles, const std::size_t max_neigh)
    {
        counts = Kokkos::View<int*, MemorySpace>("num_neighbors", num_particles);
        neighbors = Kokkos::View<int**, MemorySpace>(
            Kokkos::ViewAllocateWithoutInitializing("neighbors"),
            num_particles, max_neigh);
    }

    // Method to dynamically expand the size of max_neighbors
    KOKKOS_INLINE_FUNCTION
    void expandMaxNeighbors(const std::size_t new_max_neigh)
    {
        // Create a new view with the larger size
        Kokkos::View<int**, MemorySpace> new_neighbors(
            Kokkos::ViewAllocateWithoutInitializing("neighbors"),
            neighbors.extent(0), new_max_neigh);

        // Copy existing data to the new view
        Kokkos::parallel_for("copy_neighbors", neighbors.extent(0), KOKKOS_LAMBDA(const int i) {
            for (std::size_t j = 0; j < counts(i); ++j) {
                new_neighbors(i, j) = neighbors(i, j);
            }
        });

        // Replace the old view with the new view
        neighbors = new_neighbors;
    }

    // Method to add a neighbor
    KOKKOS_INLINE_FUNCTION
    void addNeighbor(const int pid, const int nid)
    {
        std::size_t count = Kokkos::atomic_fetch_add(&counts(pid), 1);
        if (count >= neighbors.extent(1)) {
            expandMaxNeighbors(neighbors.extent(1) * 2);
        }
        neighbors(pid, count) = nid;
    }
};

template <class MemorySpace, class AlgorithmTag, class BuildTag>
class NeighborList<
    CustomVerletList<MemorySpace, AlgorithmTag, VerletLayout2D, BuildTag>>
{
  public:
    //! Kokkos memory space.
    using memory_space = MemorySpace;
    //! Neighbor list type.
    using list_type =
        CustomVerletList<MemorySpace, AlgorithmTag, VerletLayout2D, BuildTag>;

    //! Get the total number of neighbors across all particles.
    KOKKOS_INLINE_FUNCTION
    static std::size_t totalNeighbor( const list_type& list )
    {
        std::size_t num_p = list._data.counts.size();
        for ( std::size_t i = 0; i < num_p; ++i )
            num_p += list.counts( i );
        return num_p;
    }

    //! Get the maximum number of neighbors per particle.
    KOKKOS_INLINE_FUNCTION
    static std::size_t maxNeighbor( const list_type& list )
    {
        // Stored during neighbor search.
        return list.max_n;
    }

    //! Get the number of neighbors for a given particle index.
    KOKKOS_INLINE_FUNCTION
    static std::size_t numNeighbor( const list_type& list,
                                    const std::size_t particle_index )
    {
        return list.counts( particle_index );
    }

    //! Get the id for a neighbor for a given particle index and the index of
    //! the neighbor relative to the particle.
    KOKKOS_INLINE_FUNCTION
    static std::size_t getNeighbor( const list_type& list,
                                    const std::size_t particle_index,
                                    const std::size_t count)
    {
        return list.neighbors( particle_index, count );
    }
};

}

#endif