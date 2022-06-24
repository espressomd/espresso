/*
 * Copyright (C) 2010-2022 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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
#ifndef CORE_GRID_BASED_ALGORITHMS_HALO_HPP
#define CORE_GRID_BASED_ALGORITHMS_HALO_HPP
/** \file
 *
 * Halo scheme for parallelization of lattice algorithms.
 * Header file for \ref halo.cpp.
 *
 */

#include "grid_based_algorithms/lattice.hpp"

#include <utils/Vector.hpp>

#include <mpi.h>

#include <memory>
#include <vector>

/** \name Types of halo communications */
/**@{*/
#define HALO_LOCL                                                              \
  0 /**< Tag for local exchange of halo regions on the same processor */
#define HALO_SENDRECV                                                          \
  1                 /**< Tag for halo exchange between different processors */
#define HALO_SEND 2 /**< Tag for halo send only */
#define HALO_RECV 3 /**< Tag for halo receive only */
#define HALO_OPEN 4 /**< Tag for halo open boundary */
/**@}*/

/** \name Tags for halo communications */
/**@{*/
#define REQ_HALO_SPREAD 501 /**< Tag for halo update */
#define REQ_HALO_CHECK 599  /**< Tag for consistency check of halo regions */
/**@}*/

/** Layout of the lattice data.
 *  The description is similar to MPI datatypes but a bit more compact.
 */
struct FieldType {
  FieldType(int new_extent)
      : count(0), disps({}), lengths({}), extent(new_extent), vblocks(0),
        vstride(0), vskip(0), vflag(false), subtype(nullptr) {}
  FieldType(int new_vblocks, int new_vstride, int new_vskip, bool new_vflag,
            std::shared_ptr<FieldType> oldtype)
      : count(oldtype->count), disps(oldtype->disps), lengths(oldtype->lengths),
        extent(0), vblocks(new_vblocks), vstride(new_vstride), vskip(new_vskip),
        vflag(new_vflag), subtype(oldtype) {
    if (vflag) {
      extent = oldtype->extent * ((vblocks - 1) * vskip + vstride);
    } else {
      extent = oldtype->extent * vstride + (vblocks - 1) * vskip;
    }
  }
  int count;                /**< number of subtypes in fieldtype */
  std::vector<int> disps;   /**< displacements of the subtypes */
  std::vector<int> lengths; /**< lengths of the subtypes */
  int extent;  /**< extent of the complete fieldtype including gaps */
  int vblocks; /**< number of blocks in field vectors */
  int vstride; /**< size of strides in field vectors */
  int vskip;   /**< displacement between strides in field vectors */
  bool vflag;
  std::shared_ptr<FieldType> subtype;
};

/** Structure describing a Halo region */
struct HaloInfo {

  int type; /**< type of halo communication */

  int source_node; /**< index of processor which sends halo data */
  int dest_node;   /**< index of processor receiving halo data */

  unsigned long s_offset; /**< offset for send buffer */
  unsigned long r_offset; /**< offset for receive buffer */

  std::shared_ptr<FieldType>
      fieldtype;         /**< type layout of the data being exchanged */
  MPI_Datatype datatype; /**< MPI datatype of data being communicated */
};

/** Structure holding a set of \ref HaloInfo which comprise a certain
 *  parallelization scheme */
class HaloCommunicator {
public:
  HaloCommunicator(int num) : num(num) {}

  int num; /**< number of halo communications in the scheme */

  std::vector<HaloInfo> halo_info; /**< set of halo communications */
};

/** Preparation of the halo parallelization scheme. Sets up the
 *  necessary data structures for \ref halo_communication
 *  @param[in,out] hc       halo communicator being created
 *  @param[in]     lattice  lattice the communication is created for
 *  @param datatype         MPI datatype for the lattice data
 *  @param local_node_grid  Number of nodes in each spatial dimension
 */
void prepare_halo_communication(HaloCommunicator &hc, const Lattice &lattice,
                                MPI_Datatype datatype,
                                const Utils::Vector3i &local_node_grid);

/** Frees data structures associated with a halo communicator
 *  @param[in,out] hc  halo communicator to be released
 */
void release_halo_communication(HaloCommunicator &hc);

/** Perform communication according to the parallelization scheme
 *  described by the halo communicator
 *  @param[in]  hc    halo communicator describing the parallelization scheme
 *  @param[in]  base  base plane of local node
 */
void halo_communication(const HaloCommunicator &hc, char *base);

#endif /* CORE_GRID_BASED_ALGORITHMS_HALO_HPP */
