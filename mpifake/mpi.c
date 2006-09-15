// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
/** \file mpi.c
 *
 *  For more information about MPIFake, see \ref mpifake.h "mpifake.h".
 */
#include "mpi.h"

struct mpifake_dtype 
  mpifake_dtype_lb = { 0, 0, 0, 0, 0, 0, 0, NULL, NULL, NULL, NULL },
  mpifake_dtype_ub = { 0, 0, 0, 0, 0, 0, 0, NULL, NULL, NULL, NULL },
  mpifake_dtype_byte   = { 0, 0, sizeof(char), sizeof(char), 1, 1, sizeof(char), NULL, NULL, NULL, NULL },
  mpifake_dtype_char   = { 0, 0, sizeof(char), sizeof(char), 1, 1, sizeof(char), NULL, NULL, NULL, NULL },
  mpifake_dtype_int    = { 0, 0, sizeof(int), sizeof(int), 1, 1, sizeof(int), NULL, NULL, NULL, NULL },
  mpifake_dtype_long   = { 0, 0, sizeof(long), sizeof(long), 1, 1, sizeof(long), NULL, NULL, NULL, NULL },
  mpifake_dtype_double = { 0, 0, sizeof(double), sizeof(double), 1, 1, sizeof(double), NULL, NULL, NULL, NULL };

static void mpifake_dtblock(MPI_Datatype newtype, MPI_Datatype oldtype, int count, int disp);

int MPI_Type_struct(int count, int *lengths, MPI_Aint *disps, MPI_Datatype *oldtypes, MPI_Datatype *newtype) {

  int i;
  struct mpifake_dtype *ntype;

  ntype = *newtype = malloc(sizeof(struct mpifake_dtype));

  ntype->format = LAM_DTSTRUCT;

  /* initialize defaults */
  ntype->upper = 0;
  ntype->lower = 0;
  ntype->size = 0;

  ntype->count = count;

  if (count > 0) {
    ntype->dtypes = malloc(count * (sizeof(MPI_Datatype)+sizeof(int)+sizeof(int)));
    ntype->disps = (int *)((char *)ntype->dtypes + count*sizeof(MPI_Datatype));
    ntype->lengths = (int *)((char *)ntype->disps + count*sizeof(int));
  } else {
    ntype->size = 0;
  }

  for (i=0;i<count;i++) {

    ntype->disps[i] = disps[i];
    ntype->lengths[i] = lengths[i];
    ntype->dtypes[i] = oldtypes[i];

    mpifake_dtblock(ntype, oldtypes[i], lengths[i], disps[i]);

  }

  return MPI_SUCCESS;

}

int MPI_Type_contiguous(int count, MPI_Datatype oldtype, MPI_Datatype *newtype)
{
  struct mpifake_dtype *ntype;

  ntype = *newtype = malloc(sizeof(struct mpifake_dtype));

  ntype->format = LAM_DTCONTIG;

  /* initialize defaults */
  ntype->upper = 0;
  ntype->lower = 0;
  ntype->size = 0;

  ntype->count  = count;
  ntype->dtype  = oldtype;

  mpifake_dtblock(ntype, oldtype, count, 0);

  return MPI_SUCCESS;
}

int MPI_Type_vector(int count, int length, int stride,
		    MPI_Datatype oldtype, MPI_Datatype *newtype)
{
  struct mpifake_dtype *ntype;

  ntype = *newtype = malloc(sizeof(struct mpifake_dtype));

  ntype->format = LAM_DTVECTOR;

  /* initialize defaults */
  ntype->upper = 0;
  ntype->lower = 0;
  ntype->size = 0;

  ntype->count  = count;
  ntype->length = length;
  ntype->stride = stride;
  ntype->dtype  = oldtype;

  mpifake_dtblock(ntype, oldtype, length, 0);

  ntype->size *= count;

  stride *= (oldtype->upper - oldtype->lower) * (count - 1);

  if (stride > 0) {
    ntype->upper += stride;
  } else {
    ntype->lower += stride;
  }

  return MPI_SUCCESS;

}

int MPI_Type_hvector(int count, int length, int stride,
		     MPI_Datatype oldtype, MPI_Datatype *newtype)
{
  struct mpifake_dtype *ntype;

  ntype = *newtype = malloc(sizeof(struct mpifake_dtype));

  ntype->format = LAM_DTHVECTOR;

  /* initialize defaults */
  ntype->upper = 0;
  ntype->lower = 0;
  ntype->size = 0;

  ntype->count  = count;
  ntype->length = length;
  ntype->stride = stride;
  ntype->dtype  = oldtype;

  mpifake_dtblock(ntype, oldtype, length, 0);

  ntype->size *= count;

  stride *= count - 1;

  if (stride >0) {
    ntype->upper += stride;
  } else {
    ntype->lower += stride;
  }

  return MPI_SUCCESS;

}

static void mpifake_dtblock(MPI_Datatype newtype, MPI_Datatype oldtype, int count, int disp) {

  int extent, upper, lower;

  if (count > 0) {

    extent = (oldtype->upper - oldtype->lower) * (count - 1);
    if (extent > 0) {
      upper = oldtype->upper + extent + disp;
      lower = oldtype->lower + disp;
    } else {
      upper = oldtype->upper + disp;
      lower = oldtype->lower + extent + disp;
    }

  } else {

    upper = 0;
    lower = 0;

  }

  if (upper > newtype->upper) newtype->upper = upper;
  if (lower < newtype->lower) newtype->lower = lower;

  newtype->size += count * oldtype->size;

}

static void mpifake_cpy_hvector(void *dest, void *src, int num, MPI_Datatype dtype, int vflag) {

  int i, j;
  int extent, stride;
  MPI_Datatype subtype = dtype->dtype;

  extent = dtype->upper - dtype->lower ;

  stride = dtype->stride;
  if (vflag) {
    stride *= subtype->upper - subtype->lower;
  }

  for (i=0; i<num; i++, src+=extent, dest+=extent) {
    for (j=0; j<dtype->count; j++) {
      mpifake_copy(src+j*stride,dest+j*stride,&dtype->length,&subtype);
    }
  }

}

static void mpifake_cpy_struct(void *dest, void *src, int count, MPI_Datatype dtype){

  int i, j;
  int extent, *len, *disp;
  MPI_Datatype *type;

  extent = dtype->upper - dtype->lower;

  for (i=0; i<count; i++, src+=extent, dest+=extent) {
    len = dtype->lengths;
    disp = dtype->disps;
    type = dtype->dtypes;
    for (j=0; j<dtype->count; j++, len++, disp++, type++) {
      mpifake_copy(src+*disp,dest+*disp,len,type);
    }
  }

}

void mpifake_copy(void *src, void *dest, int *count, MPI_Datatype *dtype) {

  int num;

  switch ((*dtype)->format) {
    
  case LAM_DTSTRUCT:
    mpifake_cpy_struct(dest, src, *count, *dtype);
    break;

  case LAM_DTCONTIG:
    num = *count*(*dtype)->count;
    mpifake_copy(src, dest, &num, &((*dtype)->dtype));
    break;

  case LAM_DTVECTOR:
    mpifake_cpy_hvector(dest, src, *count, *dtype, 1);
    break;

  case LAM_DTHVECTOR:
    mpifake_cpy_hvector(dest, src, *count, *dtype, 0);
    break;

  default:
    memcpy((char *)dest, (char *)src, *count*(*dtype)->size);

  }

}

int mpifake_checked_copy(void *s, int scount, MPI_Datatype sdtype,
			 void *r, int rcount, MPI_Datatype rdtype)
{
  if (sdtype != rdtype || scount != rcount) {
    fprintf(stderr, "MPI_Gather: send type != recv type || scount != rcount\n");
    errexit();
  }
  mpifake_copy(s, r, &rcount, &rdtype);
  return MPI_SUCCESS;
}

