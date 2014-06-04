/*
  Copyright (C) 2010,2012,2013 The ESPResSo project
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
/** \file mpi.cpp
 *
 *  For more information about MPIFake, see \ref mpi.h "mpi.h".
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
static void mpifake_pack(void *dest, void *src, int num, MPI_Datatype dtype);
static int mpifake_unpack(void *dest, void *src, int num, MPI_Datatype dtype);

int MPI_Type_struct(int count, int *lengths, MPI_Aint *disps, MPI_Datatype *oldtypes, MPI_Datatype *newtype) {

  int i;
  struct mpifake_dtype *ntype;

  ntype = *newtype = static_cast<MPI_Datatype>(malloc(sizeof(struct mpifake_dtype)));

  ntype->format = LAM_DTSTRUCT;

  /* initialize defaults */
  ntype->upper = 0;
  ntype->lower = 0;
  ntype->size = 0;

  ntype->count = count;

  if (count > 0) {
    ntype->dtypes = 
      static_cast<mpifake_dtype**>
      (malloc(count * (sizeof(MPI_Datatype)+sizeof(int)+sizeof(int))));
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

  ntype = *newtype = static_cast<MPI_Datatype>(malloc(sizeof(struct mpifake_dtype)));

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

  ntype = *newtype = static_cast<MPI_Datatype>(malloc(sizeof(struct mpifake_dtype)));

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
		     MPI_Datatype oldtype, MPI_Datatype *newtype) {
  return MPI_Type_create_hvector(count, length, stride, oldtype, newtype);
}

int MPI_Type_create_hvector(int count, int length, int stride,
		     MPI_Datatype oldtype, MPI_Datatype *newtype)
{
  struct mpifake_dtype *ntype;

  ntype = *newtype = static_cast<MPI_Datatype>(malloc(sizeof(struct mpifake_dtype)));

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

static void mpifake_pack_hvector(void *_dest, void *_src, int num, MPI_Datatype dtype, int vflag) {

  char* dest = static_cast<char*>(_dest);
  char* src = static_cast<char*>(_src);

  MPI_Datatype subtype = dtype->dtype;
  char *s;
  int count   = dtype->count;
  int stride  = dtype->stride;
  int extent  = dtype->upper - dtype->lower;
  int blksize = dtype->length * subtype->size; 
  int i, j;

  if (vflag) {
    stride *= subtype->upper - subtype->lower;
  }

  for (i=0; i<num; i++, src += extent) {
    s = src;
    for (j=0; j<count; j++) {
      mpifake_pack(dest, s, dtype->length, dtype->dtype);
      dest += blksize;
      s += stride;
    }
  }

}

static void mpifake_pack_struct(void *_dest, void *_src, int num, MPI_Datatype dtype) {

  char *s;
  int *len;
  int *disp;
  MPI_Datatype *type;
  int extent = dtype->upper - dtype->lower;
  int blksize;
  int i, j;

  char* dest = static_cast<char*>(_dest);
  char* src = static_cast<char*>(_src);

  for (i=0; i<num; i++, src += extent) {
    s = src;
    len = dtype->lengths;
    disp = dtype->disps;
    type = dtype->dtypes;
    for (j=0; j<dtype->count; j++, len++, disp++, type++) {
      blksize = *len * (*type)->size;
      if (blksize >0) {
	mpifake_pack(dest, s + *disp, *len, *type);
	dest += blksize;
      }
    }
  }

}

static void mpifake_pack(void *dest, void *src, int count, MPI_Datatype dtype) {

  switch (dtype->format) {

  case LAM_DTSTRUCT:
    mpifake_pack_struct(dest, src, count, dtype);
    break;

  case LAM_DTCONTIG:
    mpifake_pack(dest, src, count * dtype->count, dtype->dtype);
    break;

  case LAM_DTVECTOR:
    mpifake_pack_hvector(dest, src, count, dtype, 1);
    break;

  case LAM_DTHVECTOR:
    mpifake_pack_hvector(dest, src, count, dtype, 0);
    break;

  default:
    memcpy((char *)dest, (char *)src, count *dtype->size);

  }

}   

static int mpifake_unpack_hvector(void *_dest, void *_src, int num, 
                                  MPI_Datatype dtype, int vflag) {

  MPI_Datatype subtype = dtype->dtype;
  char *d;
  char *src = static_cast<char*>(_src);
  char *dest = static_cast<char*>(_dest);
  char *start = src;
  int count   = dtype->count;
  int stride  = dtype->stride;
  int extent  = dtype->upper - dtype->lower;
  int blksize = dtype->length * subtype->size;
  int size;
  int i, j;

  if (vflag) {
    stride *= subtype->upper - subtype->lower;
  }

  for (i=0; i<num; i++, dest += extent) {
    d = dest;
    for (j=0; j<count; j++) {
      size = mpifake_unpack(d, src, dtype->length, subtype);
      src += size;

      if (size != blksize) {
	fprintf(stderr,"mpifake_unpack_hvector: size != blksize\n");
	errexit();
      }

      d += stride;
    }
  }

  return (char *)src - start;
}

static int mpifake_unpack_struct(void *_dest, void *_src, int num, MPI_Datatype dtype) {

  char *d;
  char *src = static_cast<char*>(_src);
  char *dest = static_cast<char*>(_dest);
  char *start = src;
  int *len;
  int *disp;
  MPI_Datatype *type;
  int extent = dtype->upper - dtype->lower;
  int size;
  int blksize;
  int i, j;

  for (i=0; i<num; i++, dest+=extent) {
    d =dest;
    len = dtype->lengths;
    disp = dtype->disps;
    type = dtype->dtypes;
    for (j=0; j<dtype->count; j++, len++, disp++, type++) {
      blksize = *len * (*type)->size;
      if (blksize >0) {
	size = mpifake_unpack(d + *disp, src, *len, *type);
	src += size;

	if (size != blksize) {
	  fprintf(stderr,"mpifake_unpack_struct: size != blksize\n");
	  errexit();
	}

      }
    }
  }
  
  return (char *)src - start;
}

static int mpifake_unpack(void *_dest, void *_src, int count, MPI_Datatype dtype) {
  char *src = static_cast<char*>(_src);
  char *dest = static_cast<char*>(_dest);
  int size;

  switch (dtype->format) {

  case LAM_DTSTRUCT:
    return mpifake_unpack_struct(dest, src, count, dtype);
    
  case LAM_DTCONTIG:
    return mpifake_unpack(dest, src, count * dtype->count, dtype->dtype);

  case LAM_DTVECTOR:
    return mpifake_unpack_hvector(dest, src, count, dtype, 1);

  case LAM_DTHVECTOR:
    return mpifake_unpack_hvector(dest, src, count, dtype, 0);

  default:
    size = count * dtype->size;
    memcpy((char *)dest, (char *)src, size);
    return size;

  }

}

static void mpifake_cpy_hvector(void *_dest, void *_src, int num, MPI_Datatype dtype, int vflag) {

  int i, j;
  char *src = static_cast<char*>(_src);
  char *dest = static_cast<char*>(_dest);
  int extent, stride;
  MPI_Datatype subtype = dtype->dtype;

  extent = dtype->upper - dtype->lower ;

  stride = dtype->stride;
  if (vflag) {
    stride *= subtype->upper - subtype->lower;
  }

  for (i=0; i<num; i++, src += extent, dest += extent) {
    for (j=0; j<dtype->count; j++) {
      mpifake_copy(src+j*stride,dest+j*stride,&dtype->length,&subtype);
    }
  }

}

static void mpifake_cpy_struct(void *_dest, void *_src, int count, MPI_Datatype dtype){

  char *src = static_cast<char*>(_src);
  char *dest = static_cast<char*>(_dest);
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

int mpifake_sendrecv(void *_s, int scount, MPI_Datatype sdtype,
                     void *_r, int rcount, MPI_Datatype rdtype)
{
  char *s = static_cast<char*>(_s);
  char *r = static_cast<char*>(_r);

  char *packbuf;

  if (sdtype == rdtype) {
    if (scount > rcount) {
      fprintf(stderr, "MPI_Gather: scount > rcount\n");
      errexit();
    } else {
      mpifake_copy(s, r, &rcount, &rdtype);
    }
  } else {
    packbuf = static_cast<char*>(malloc(scount * sdtype->size));
    mpifake_pack(packbuf, s, scount, sdtype);
    mpifake_unpack(r, packbuf, rcount, rdtype);
    free(packbuf);
  }

  return MPI_SUCCESS;
}

