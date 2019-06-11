/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/** \file
 *
 * Halo scheme for parallelization of lattice algorithms.
 * Header file for \ref halo.cpp.
 *
 */

#ifndef _HALO_HPP
#define _HALO_HPP

#include "grid_based_algorithms/lattice.hpp"
#include <mpi.h>

/** \name Types of halo communications */
/*@{*/
#define HALO_LOCL                                                              \
  0 /**< Tag for local exchange of halo regions on the same processor */
#define HALO_SENDRECV                                                          \
  1                 /**< Tag for halo exchange between different processors */
#define HALO_SEND 2 /**< Tag for halo send only */
#define HALO_RECV 3 /**< Tag for halo receive only */
#define HALO_OPEN 4 /**< Tag for halo open boundary */
/*@}*/

/** \name Tags for halo communications */
/*@{*/
#define REQ_HALO_SPREAD 501 /**< Tag for halo update */
#define REQ_HALO_CHECK 599  /**< Tag for consistency check of halo regions */
/*@}*/

/** Layout of the lattice data.
 *  The description is similar to MPI datatypes but a bit more compact.
 *  See \ref halo_create_fieldtype, \ref halo_create_field_vector and \ref
 *  halo_dtcopy to understand how it works.
 */
typedef struct _Fieldtype *Fieldtype;
struct _Fieldtype {
  int count;    /**< number of subtypes in fieldtype */
  int *disps;   /**< displacements of the subtypes */
  int *lengths; /**< lengths of the subtypes */
  int extent;   /**< extent of the complete fieldtype including gaps */
  int vblocks;  /**< number of blocks in field vectors */
  int vstride;  /**< size of strides in field vectors */
  int vskip;    /**< displacement between strides in field vectors */
  int vflag;
  Fieldtype subtype;
};

/** Predefined fieldtypes */
extern struct _Fieldtype fieldtype_double;
#define FIELDTYPE_DOUBLE (&fieldtype_double)

/** Structure describing a Halo region */
typedef struct {

  int type; /**< type of halo communication */

  int source_node; /**< index of processor which sends halo data */
  int dest_node;   /**< index of processor receiving halo data */

  // void *send_buffer;      /**< pointer to data being sent */
  // void *recv_buffer;      /**< pointer to data being received */

  unsigned long s_offset; /**< offset for send buffer */
  unsigned long r_offset; /**< offset for receive buffer */

  Fieldtype fieldtype;   /**< type layout of the data being exchanged */
  MPI_Datatype datatype; /**< MPI datatype of data being communicated */

} HaloInfo;

/** Structure holding a set of \ref HaloInfo which comprise a certain
 *  parallelization scheme */
typedef struct {

  int num; /**< number of halo communications in the scheme */

  HaloInfo *halo_info; /**< set of halo communications */

} HaloCommunicator;

/** Creates a fieldtype describing the data layout
 *  @param      count   number of subtypes
 *  @param[in]  lengths array of lengths of the subtypes
 *  @param[in]  disps   array of displacements the subtypes
 *  @param      extent  extent of the whole new fieldtype
 *  @param[out] newtype newly created fieldtype
 */
void halo_create_fieldtype(int count, int const *lengths, int const *disps,
                           int extent, Fieldtype *newtype);

/** Creates a field vector layout
 *  @param vblocks       number of vector blocks
 *  @param vstride       size of strides in field vector
 *  @param vskip         displacements of strides in field vector
 *  @param oldtype       fieldtype the vector is composed of
 *  @param[out] newtype  newly created fieldtype
 */
void halo_create_field_vector(int vblocks, int vstride, int vskip,
                              Fieldtype oldtype, Fieldtype *newtype);
void halo_create_field_hvector(int vblocks, int vstride, int vskip,
                               Fieldtype oldtype, Fieldtype *newtype);

/** Frees a fieldtype
 *  @param ftype pointer to the type to be freed
 */
void halo_free_fieldtype(Fieldtype *ftype);

/** Preparation of the halo parallelization scheme. Sets up the
 *  necessary datastructures for \ref halo_communication
 *  @param[in,out] hc       halo communicator being created
 *  @param[in]     lattice  lattice the communication is created for
 *  @param fieldtype        field layout of the lattice data
 *  @param datatype         MPI datatype for the lattice data
 *  @param local_node_grid  Number of nodes in each spatial dimension
 */
void prepare_halo_communication(HaloCommunicator *hc, Lattice const *lattice,
                                Fieldtype fieldtype, MPI_Datatype datatype,
                                const Utils::Vector3i &local_node_grid);

/** Frees datastructures associated with a halo communicator
 *  @param[in,out] hc  halo communicator to be released
 */
void release_halo_communication(HaloCommunicator *hc);

/** Perform communication according to the parallelization scheme
 *  described by the halo communicator
 *  @param[in]  hc    halo communicator describing the parallelization scheme
 *  @param[in]  base  base plane of local node
 */
void halo_communication(HaloCommunicator const *hc, char *base);

#endif /* HALO_H */
