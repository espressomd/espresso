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
 * Implementation of file \ref halo.hpp.
 *
 */

#include "config.hpp"
#include <cstring>

#include "communication.hpp"
#include "debug.hpp"
#include "grid.hpp"
#include "grid_based_algorithms/lattice.hpp"
#include "halo.hpp"

/** Primitive fieldtypes and their initializers */
struct _Fieldtype fieldtype_double = {0, nullptr, nullptr, sizeof(double), 0,
                                      0, 0,       0,       nullptr};

void halo_create_field_vector(int vblocks, int vstride, int vskip,
                              Fieldtype oldtype, Fieldtype *const newtype) {
  int i;

  Fieldtype ntype = *newtype = (Fieldtype)Utils::malloc(sizeof(*ntype));

  ntype->subtype = oldtype;
  ntype->vflag = 1;

  ntype->vblocks = vblocks;
  ntype->vstride = vstride;
  ntype->vskip = vskip;

  ntype->extent = oldtype->extent * ((vblocks - 1) * vskip + vstride);

  int count = ntype->count = oldtype->count;
  ntype->lengths = (int *)Utils::malloc(count * 2 * sizeof(int));
  ntype->disps = (int *)((char *)ntype->lengths + count * sizeof(int));

  for (i = 0; i < count; i++) {
    ntype->disps[i] = oldtype->disps[i];
    ntype->lengths[i] = oldtype->lengths[i];
  }
}

void halo_create_field_hvector(int vblocks, int vstride, int vskip,
                               Fieldtype oldtype, Fieldtype *const newtype) {
  int i;

  Fieldtype ntype = *newtype = (Fieldtype)Utils::malloc(sizeof(*ntype));

  ntype->subtype = oldtype;
  ntype->vflag = 0;

  ntype->vblocks = vblocks;
  ntype->vstride = vstride;
  ntype->vskip = vskip;

  ntype->extent = oldtype->extent * vstride + (vblocks - 1) * vskip;

  int count = ntype->count = oldtype->count;
  ntype->lengths = (int *)Utils::malloc(count * 2 * sizeof(int));
  ntype->disps = (int *)((char *)ntype->lengths + count * sizeof(int));

  for (i = 0; i < count; i++) {
    ntype->disps[i] = oldtype->disps[i];
    ntype->lengths[i] = oldtype->lengths[i];
  }
}

void halo_free_fieldtype(Fieldtype *const ftype) {
  if ((*ftype)->count > 0) {
    free((*ftype)->lengths);
    (*ftype)->lengths = nullptr;
  }
  free(*ftype);
}

/** Set halo region to a given value
 * @param[out] dest pointer to the halo buffer
 * @param value integer value to write into the halo buffer
 * @param type halo field layout description
 */
void halo_dtset(char *dest, int value, Fieldtype type) {
  int vblocks = type->vblocks;
  int vstride = type->vstride;
  int vskip = type->vskip;
  int count = type->count;
  int *lens = type->lengths;
  int *disps = type->disps;
  int extent = type->extent;

  for (int i = 0; i < vblocks; i++) {
    for (int j = 0; j < vstride; j++) {
      for (int k = 0; k < count; k++)
        memset(dest + disps[k], value, lens[k]);
    }
    dest += vskip * extent;
  }
}

void halo_dtcopy(char *r_buffer, char *s_buffer, int count, Fieldtype type);

void halo_copy_vector(char *r_buffer, char *s_buffer, int count, Fieldtype type,
                      int vflag) {
  int i, j;
  char *dest, *src;

  int vblocks = type->vblocks;
  int vstride = type->vstride;
  int vskip = type->vskip;
  int extent = type->extent;

  HALO_TRACE(fprintf(stderr,
                     "%d: halo_copy_vector %p %p vblocks=%d vstride=%d "
                     "vskip=%d extent=%d subtype_extent=%d\n",
                     this_node, static_cast<void *>(r_buffer),
                     static_cast<void *>(s_buffer), vblocks, vstride, vskip,
                     extent, type->subtype->extent));

  if (vflag) {
    vskip *= type->subtype->extent;
  }

  for (i = 0; i < count; i++, s_buffer += extent, r_buffer += extent) {
    for (j = 0, dest = r_buffer, src = s_buffer; j < vblocks;
         j++, dest += vskip, src += vskip) {
      halo_dtcopy(dest, src, vstride, type->subtype);
    }
  }
}

/** Copy lattice data with layout described by fieldtype.
 * @param r_buffer data destination
 * @param s_buffer data source
 * @param count    amount of data to copy
 * @param type     field layout type
 */
void halo_dtcopy(char *r_buffer, char *s_buffer, int count, Fieldtype type) {

  HALO_TRACE(fprintf(
      stderr,
      "%d: halo_dtcopy r_buffer=%p s_buffer=%p blocks=%d stride=%d skip=%d\n",
      this_node, static_cast<void *>(r_buffer), static_cast<void *>(s_buffer),
      type->vblocks, type->vstride, type->vskip));

  if (type->subtype) {
    halo_copy_vector(r_buffer, s_buffer, count, type, type->vflag);
  } else {

    for (int i = 0; i < count;
         i++, s_buffer += type->extent, r_buffer += type->extent) {
      if (!type->count) {
        memmove(r_buffer, s_buffer, type->extent);
      } else {

        for (int j = 0; j < type->count; j++) {
          memmove(r_buffer + type->disps[j], s_buffer + type->disps[j],
                  type->lengths[j]);
        }
      }
    }
  }
}

void prepare_halo_communication(HaloCommunicator *const hc,
                                Lattice const *const lattice,
                                Fieldtype fieldtype, MPI_Datatype datatype,
                                const Utils::Vector3i &local_node_grid) {
  int k, n, dir, lr, cnt, num = 0;
  const auto grid = lattice->grid;
  const auto period = lattice->halo_grid;

  for (n = 0; n < hc->num; n++) {
    MPI_Type_free(&(hc->halo_info[n].datatype));
  }

  num = 2 * 3; /* two communications in each space direction */

  hc->num = num;
  hc->halo_info = Utils::realloc(hc->halo_info, num * sizeof(HaloInfo));

  int extent = fieldtype->extent;

  auto const node_neighbors = calc_node_neighbors(comm_cart);

  cnt = 0;
  for (dir = 0; dir < 3; dir++) {
    for (lr = 0; lr < 2; lr++) {

      HaloInfo *hinfo = &(hc->halo_info[cnt]);

      int nblocks = 1;
      for (k = dir + 1; k < 3; k++) {
        nblocks *= period[k];
      }
      int stride = 1;
      for (k = 0; k < dir; k++) {
        stride *= period[k];
      }
      int skip = 1;
      for (k = 0; k < dir + 1 && k < 2; k++) {
        skip *= period[k];
      }

      if (lr == 0) {
        /* send to left, recv from right */
        hinfo->s_offset = extent * stride * 1;
        hinfo->r_offset = extent * stride * (grid[dir] + 1);
      } else {
        /* send to right, recv from left */
        hinfo->s_offset = extent * stride * grid[dir];
        hinfo->r_offset = extent * stride * 0;
      }

      hinfo->source_node = node_neighbors[2 * dir + 1 - lr];
      hinfo->dest_node = node_neighbors[2 * dir + lr];

      halo_create_field_vector(nblocks, stride, skip, fieldtype,
                               &hinfo->fieldtype);

      MPI_Type_vector(nblocks, stride, skip, datatype, &hinfo->datatype);
      MPI_Type_commit(&hinfo->datatype);

      if (!box_geo.periodic(dir) &&
          (local_geo.boundary()[2 * dir + lr] != 0 ||
           local_geo.boundary()[2 * dir + 1 - lr] != 0)) {
        if (local_node_grid[dir] == 1) {
          hinfo->type = HALO_OPEN;
        } else if (lr == 0) {
          if (local_geo.boundary()[2 * dir + lr] == 1) {
            hinfo->type = HALO_RECV;
          } else {
            hinfo->type = HALO_SEND;
          }
        } else {
          if (local_geo.boundary()[2 * dir + lr] == -1) {
            hinfo->type = HALO_RECV;
          } else {
            hinfo->type = HALO_SEND;
          }
        }
      } else {
        if (local_node_grid[dir] == 1) {
          hc->halo_info[cnt].type = HALO_LOCL;
        } else {
          hc->halo_info[cnt].type = HALO_SENDRECV;
        }
      }

      HALO_TRACE(
          fprintf(stderr,
                  "%d: prepare_halo_communication dir=%d lr=%d s_offset=%ld "
                  "r_offset=%ld s_node=%d d_node=%d type=%d\n",
                  this_node, dir, lr, hinfo->s_offset, hinfo->r_offset,
                  hinfo->source_node, hinfo->dest_node, hinfo->type));

      cnt++;
    }
  }
}

void release_halo_communication(HaloCommunicator *const hc) {
  int n;

  for (n = 0; n < hc->num; n++) {
    MPI_Type_free(&(hc->halo_info[n].datatype));
  }

  free(hc->halo_info);
}

void halo_communication(HaloCommunicator const *const hc, char *const base) {
  int s_node, r_node;

  Fieldtype fieldtype;
  MPI_Datatype datatype;
  MPI_Request request;
  MPI_Status status;

  HALO_TRACE(fprintf(stderr, "%d: halo_comm base=%p num=%d\n", this_node,
                     static_cast<void *>(base), hc->num));

  for (int n = 0; n < hc->num; n++) {

    HALO_TRACE(fprintf(stderr, "%d: halo_comm round %d\n", this_node, n));

    int comm_type = hc->halo_info[n].type;
    char *s_buffer = (char *)base + hc->halo_info[n].s_offset;
    char *r_buffer = (char *)base + hc->halo_info[n].r_offset;

    switch (comm_type) {

    case HALO_LOCL:
      fieldtype = hc->halo_info[n].fieldtype;
      halo_dtcopy(r_buffer, s_buffer, 1, fieldtype);
      break;

    case HALO_SENDRECV:
      datatype = hc->halo_info[n].datatype;
      s_node = hc->halo_info[n].source_node;
      r_node = hc->halo_info[n].dest_node;

      HALO_TRACE(fprintf(stderr, "%d: halo_comm sendrecv %d to %d (%d) (%p)\n",
                         this_node, s_node, r_node, REQ_HALO_SPREAD,
                         (void *)&datatype));

      MPI_Sendrecv(s_buffer, 1, datatype, r_node, REQ_HALO_SPREAD, r_buffer, 1,
                   datatype, s_node, REQ_HALO_SPREAD, comm_cart, &status);
      break;

    case HALO_SEND:
      datatype = hc->halo_info[n].datatype;
      fieldtype = hc->halo_info[n].fieldtype;
      r_node = hc->halo_info[n].dest_node;

      HALO_TRACE(
          fprintf(stderr, "%d: halo_comm send to %d.\n", this_node, r_node));

      MPI_Isend(s_buffer, 1, datatype, r_node, REQ_HALO_SPREAD, comm_cart,
                &request);
      halo_dtset(r_buffer, 0, fieldtype);
      MPI_Wait(&request, &status);
      break;

    case HALO_RECV:
      datatype = hc->halo_info[n].datatype;
      s_node = hc->halo_info[n].source_node;

      HALO_TRACE(
          fprintf(stderr, "%d: halo_comm recv from %d.\n", this_node, s_node));

      MPI_Irecv(r_buffer, 1, datatype, s_node, REQ_HALO_SPREAD, comm_cart,
                &request);
      MPI_Wait(&request, &status);
      break;

    case HALO_OPEN:
      fieldtype = hc->halo_info[n].fieldtype;

      HALO_TRACE(fprintf(stderr, "%d: halo_comm open boundaries\n", this_node));

      /* \todo this does not work for the n_i - <n_i> */
      halo_dtset(r_buffer, 0, fieldtype);
      break;
    }
  }
}
