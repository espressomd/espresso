/* $Id$
 *
 * This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
 * It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
 * and by which you are legally bound while utilizing this file in any form or way.
 * There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * You should have received a copy of that license along with this program;
 * if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
 * write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
 * Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
 */

/** \file halo.c
 * 
 * Halo scheme for parallelization of lattice algorithms. 
 * Implementation of file \ref halo.h.
 *
 */

#include "utils.h"
#include "grid.h"
#include "lattice.h"
#include "halo.h"

#ifdef LATTICE

/** Creates a fieldtype describing the data layout 
 *  @param count   number of subtypes (Input)
 *  @param lengths array of lenghts of the subtytpes (Input)
 *  @param disps   array of displacements the subtypes (Input)
 *  @param extent  extent of the whole new fieldtype (Input)
 *  @param newtype newly created fieldtype (Input/Output)
 */
void halo_create_fieldtype(int count, int* lengths, int *disps, int extent, Fieldtype *newtype) {
  int i;

  Fieldtype ntype = *newtype = malloc(sizeof(*ntype));

  ntype->vblocks = 1;
  ntype->vstride = 1;
  ntype->vskip = 1;

  ntype->count = count;
  ntype->extent = extent;

  if (count>0) {

      ntype->lengths = malloc(count*2*sizeof(int));
      ntype->disps = (int *)((char *)ntype->lengths + count*sizeof(int));

      for (i=0;i<count;i++) {
	  ntype->disps[i] = disps[i];
	  ntype->lengths[i] = lengths[i];

      }

  }

}

/** Creates a field vector layout
 *  @param vblocks number of vector blocks(Input) 
 *  @param vstride size of strides in field vector (Input)
 *  @param vskip   displacements of strides in field vector (Input)
 *  @param oldtype fieldtype the vector is composed of (Input)
 *  @param newtype newly created fieldtype (Input/Output)
 */
void halo_create_field_vector(int vblocks, int vstride, int vskip, Fieldtype oldtype, Fieldtype *newtype) {
  int i;

  Fieldtype ntype = *newtype = malloc(sizeof(*ntype));

  ntype->vblocks = vblocks;
  ntype->vstride = vstride;
  ntype->vskip = vskip;

  ntype->extent = oldtype->extent;

  int count = ntype->count = oldtype->count;
  ntype->lengths = malloc(count*2*sizeof(int));
  ntype->disps = (int *)((char *)ntype->lengths + count*sizeof(int));

  for (i=0;i<count;i++) {
      ntype->disps[i] = oldtype->disps[i];
      ntype->lengths[i] = oldtype->lengths[i];
  }

}

/** Frees a fieldtype
 * @param ftype pointer to the type to be freed (Input)
 */
void halo_free_fieldtype(Fieldtype *ftype) {
  if ((*ftype)->count>0) {
    free((*ftype)->lengths);
    (*ftype)->lengths = NULL;
  }
  free(*ftype);
}

/** Copy lattice data with layout described by fieldtype.
 * @param r_buffer data destination
 * @param s_buffer data source
 * @param type     field layout type
 */
MDINLINE void halo_dtcopy(void *r_buffer, void *s_buffer, Fieldtype type) { 
    int i, j, k;
    void *dest, *src;

    int vblocks = type->vblocks;
    int vstride = type->vstride;
    int vskip   = type->vskip;
    int count   = type->count;
    int *lens   = type->lengths;
    int *disps  = type->disps;
    int extent  = type->extent;

    HALO_TRACE(fprintf(stderr, "%d: halo comm local copy r_buffer=%p s_buffer=%p\n",this_node,r_buffer,s_buffer));

    for (i=0; i<vblocks; i++, r_buffer+=vskip*extent, s_buffer+=vskip*extent) {
	dest=r_buffer; src=s_buffer;
	for (j=0; j<vstride; j++, dest+=extent, src+=extent) {
	    for (k=0; k<count; k++) {
		memcpy(dest+disps[k],src+disps[k],lens[k]);
	    }
	}
    }

}

/** Preparation of the halo parallelization scheme. Sets up the
 *  necessary datastructures for \ref halo_communication
 * @param hc         halo communicator beeing created (Input/Output)
 * @param lattice    lattice the communcation is created for (Input)
 * @param fieldtype  field layout of the lattice data (Input)
 * @param datatype   MPI datatype for the lattice data (Input)
 */
void prepare_halo_communication(HaloCommunicator *hc, Lattice *lattice, Fieldtype fieldtype, MPI_Datatype datatype) {
  int k, n, dir, lr, cnt, num = 0 ;
  int *grid  = lattice->grid ;
  int *period = lattice->halo_grid ;
  void *data = lattice->data;

  for (n=0; n<hc->num; n++) {
    MPI_Type_free(&(hc->halo_info[n].datatype));
  }

  for (dir=0; dir<3; dir++) {
    for (lr=0; lr<2; lr++) {

#ifdef PARTIAL_PERIODIC
      if ( PERIODIC(dir) || (boundary[2*dir+lr] == 0) )
#endif
	{
	    num += 1 ;
	}

    }
  }

  hc->num = num ;
  hc->halo_info = realloc(hc->halo_info,num*sizeof(HaloInfo)) ;

  int extent = fieldtype->extent;

  cnt = 0 ;
  for (dir=0; dir<3; dir++) {
    for (lr=0; lr<2; lr++) {

#ifdef PARTIAL_PERIODIC
	if ( PERIODIC(dir) || (boundary[2*dir+lr] == 0) )
#endif
	  {
	      HaloInfo *hinfo = &(hc->halo_info[cnt]) ;

	      int nblocks = 1 ;
	      for (k=dir+1;k<3;k++) {
		  nblocks *= period[k] ;
	      }
	      int stride = 1 ;
	      for (k=0;k<dir;k++) {
		  stride *= period[k] ;
	      }
	      int skip = 1 ;
	      for (k=0;k<dir+1 && k<2;k++) {
		  skip *= period[k] ;
	      }	      

	      if (lr==0) {
		/* send to left, recv from right */
		hinfo->send_buffer = &(((char *)data)[extent * stride * 1]);
		hinfo->recv_buffer = &(((char *)data)[extent * stride * (grid[dir]+1)]);
	      } else {
		/* send to right, recv from left */
		hinfo->send_buffer = &(((char *)data)[extent * stride * grid[dir]]);
		hinfo->recv_buffer   = &(((char *)data)[extent * stride * 0]);
	      } 

	      hinfo->source_node = node_neighbors[2*dir+1-lr];
	      hinfo->dest_node = node_neighbors[2*dir+lr];

	      halo_create_field_vector(nblocks, stride, skip, fieldtype, &hinfo->fieldtype);
	      
	      MPI_Type_vector(nblocks, stride, skip, datatype, &hinfo->datatype);
	      MPI_Type_commit(&hinfo->datatype);
			       
	      if (node_grid[dir] == 1) {
		  hc->halo_info[cnt].type = HALO_LOCL ;
	      }
	      else {
		  hc->halo_info[cnt].type = HALO_SENDRECV ;
	      }

	      HALO_TRACE(fprintf(stderr,"%d: prepare_halo_communication dir=%d lr=%d s_buffer=%p r_buffer=%p, s_node=%d d_node=%d extent=%d\n",this_node,dir,lr,hinfo->send_buffer,hinfo->recv_buffer,hinfo->source_node,hinfo->dest_node,(int)extent)) ;

	  }
	
	cnt++ ;

    }
  }

}

/** Frees datastrutures associated with a halo communicator 
 * @param hc halo communicator to be released
 */
void release_halo_communication(HaloCommunicator *hc) {
  int n;

  for (n=0; n<hc->num; n++) {
    HALO_TRACE(fprintf(stderr,"%d: freeing %p\n",this_node,&(hc->halo_info[n].datatype)));
    MPI_Type_free(&(hc->halo_info[n].datatype));
  }

  free(hc->halo_info);

}

/** Perform communication according to the parallelization scheme
 *  described by the halo communicator
 * @param hc halo communicator describing the parallelization scheme
 */
void halo_communication(HaloCommunicator *hc) {
  int n, comm_type, s_node, r_node;
  void *s_buffer, *r_buffer ;

  Fieldtype fieldtype;
  MPI_Datatype datatype;
  MPI_Status status[2] ;

    HALO_TRACE(fprintf(stderr, "%d: halo_comm %p (num=%d)\n", this_node, hc, hc->num)) ;

    for (n = 0; n < hc->num; n++) {

	HALO_TRACE(fprintf(stderr, "%d: halo_comm round %d\n", this_node, n)) ;

	comm_type = hc->halo_info[n].type ;
	s_buffer = hc->halo_info[n].send_buffer ;
	r_buffer = hc->halo_info[n].recv_buffer ;

	switch (comm_type) {

	    case HALO_LOCL:
	      fieldtype = hc->halo_info[n].fieldtype;
	      halo_dtcopy(r_buffer,s_buffer,fieldtype);
	      break ;

	    case HALO_SENDRECV:
	      datatype = hc->halo_info[n].datatype;
	      s_node = hc->halo_info[n].source_node ;
	      r_node = hc->halo_info[n].dest_node ;
	      
	      HALO_TRACE(fprintf(stderr,"%d: halo_comm sendrecv %d to %d (%d) (%p)\n",this_node,s_node,r_node,REQ_HALO_SPREAD,&datatype));

	      MPI_Sendrecv(s_buffer, 1, datatype, r_node, REQ_HALO_SPREAD,
			   r_buffer, 1, datatype, s_node, REQ_HALO_SPREAD,
			   MPI_COMM_WORLD, &status[0]);
	      break ;

	    case HALO_SEND:
	      datatype = hc->halo_info[n].datatype;
	      fieldtype = hc->halo_info[n].fieldtype;
	      s_node = hc->halo_info[n].source_node ;
	      r_node = hc->halo_info[n].dest_node ;
	      
	      HALO_TRACE(fprintf(stderr,"%d: halo_comm send to %d.\n",this_node,r_node));

	      MPI_Isend(s_buffer, 1, datatype, r_node, REQ_HALO_SPREAD, MPI_COMM_WORLD, &request);
	      halo_dtset(r_buffer,0,fieldtype);
	      MPI_Wait(&request,&status);
	      break;

	    case HALO_RECV:
	      datatype = hc->halo_info[n].datatype;
	      s_node = hc->halo_info[n].source_node ;
	      r_node = hc->halo_info[n].dest_node ;

	      HALO_TRACE(fprintf(stderr,"%d: halo_comm recv from %d.\n",this_node,s_node));

	      MPI_Irecv(r_buffer, 1, datatype, s_node, REQ_HALO_SPREAD, MPI_COMM_WORLD, &request);
	      MPI_Wait(&request,&status);
	      break;

	    case HALO_OPEN:
	      fieldtype = hc->halo_info[n].fieldtype;

	      HALO_TRACE(fprintf(stderr,"%d: halo_comm open boundaries\n",this_node));

	      //halo_dtset(r_buffer,0,fieldtype);
	      break;
	      
	}

    }

    HALO_TRACE(fprintf(stderr, "%d: halo_comm %p finished\n", this_node, hc));

}

#endif /* LATTICE */
