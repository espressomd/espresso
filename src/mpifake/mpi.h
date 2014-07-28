/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
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
/** \file mpi.h
    This is the MPIfake implementation. This is NOT a real MPI implementation, but rather implements
    a subset of the MPI commands such that Espresso is able to run with a single processor. For this to work,
    you do not have to have any MPI implementation like LAM or MPICH installed.
*/

#ifndef _MPI_H
#define _MPI_H

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

//#include "utils.h"

void errexit();


/********************************** REMARK **********************/
/* This is the fake MPI header of Espresso, and has nothing to  */
/* with Espresso's MPI handling using other mpi variants        */
/********************************** REMARK **********************/

typedef struct mpifake_dtype *MPI_Datatype;

struct mpifake_dtype {
  int		format;
#define LAM_DTBASIC		0		/* basic datatype */
#define LAM_DTCONTIG		1		/* contiguous */
#define LAM_DTVECTOR		2		/* vector */
#define LAM_DTHVECTOR		3		/* hvector */
#define LAM_DTINDEXED		4		/* indexed */
#define LAM_DTHINDEXED		5		/* hindexed */
#define LAM_DTSTRUCT		6		/* struct */
#define LAM_DTHVECTORCREAT	7		/* extended vector */
#define	LAM_DTHINDEXEDCREAT	8		/* extended indexed */
#define	LAM_DTSTRUCTCREAT 	9		/* extended struct */
#define LAM_DTINDEXEDBLK	10		/* indexed block */
#define LAM_DTSUBARRAY		11		/* local array */
#define LAM_DTDARRAY		12		/* distributed array */
  int		lower;		/* lower extent */
  int		upper;		/* upper extent */
  int		size;		/* basic size */
  int		count;		/* count */
  int		length;		/* vector length */
  int	        stride;		/* vector stride */
  MPI_Datatype	dtype;		/* c/v/i datatype */
  int		*lengths;	/* i/s lengths */
  int	        *disps;		/* i/s displacements */
  MPI_Datatype	*dtypes;	/* struct datatypes */
};

typedef void *MPI_Status;
typedef void *MPI_Comm;
typedef void *MPI_Errhandler;
typedef void *MPI_Request;
typedef long MPI_Aint;

typedef void (MPI_User_function)(void *, void *, int *, MPI_Datatype *);
typedef void (MPI_Handler_function)(MPI_Comm *, int *, ...);

typedef MPI_User_function *MPI_Op;

void mpifake_copy(void *from, void *to, int *count, MPI_Datatype *dtype);
int mpifake_sendrecv(void *s, int scount, MPI_Datatype sdtype,
			 void *r, int rcount, MPI_Datatype rdtype);

#define MPI_LOR mpifake_copy
#define MPI_SUM mpifake_copy
#define MPI_MAX mpifake_copy
#define MPI_COPY mpifake_copy

#define MPI_STATUS_IGNORE NULL
#define MPI_SUCCESS 1

#define MPI_COMM_WORLD NULL

#define MPI_REQUEST_NULL NULL

#define MPI_ANY_SOURCE 0

#define MPI_IN_PLACE (void*)0x1

extern struct mpifake_dtype mpifake_dtype_int;
extern struct mpifake_dtype mpifake_dtype_double;
extern struct mpifake_dtype mpifake_dtype_byte;
extern struct mpifake_dtype mpifake_dtype_long;
extern struct mpifake_dtype mpifake_dtype_char;
extern struct mpifake_dtype mpifake_dtype_ub;
extern struct mpifake_dtype mpifake_dtype_lb;

#define MPI_INT    (&mpifake_dtype_int)
#define MPI_DOUBLE (&mpifake_dtype_double)
#define MPI_BYTE   (&mpifake_dtype_byte)
#define MPI_LONG   (&mpifake_dtype_long)
#define MPI_CHAR   (&mpifake_dtype_char)
#define MPI_LB     (&mpifake_dtype_lb)
#define MPI_UB     (&mpifake_dtype_ub)
#define MPI_DATATYPE_NULL NULL

int MPI_Type_struct(int count, int *lengths, MPI_Aint *disps, MPI_Datatype *oldtypes, MPI_Datatype *newtype);
int MPI_Type_contiguous(int count, MPI_Datatype oldtype, MPI_Datatype *newtype);
int MPI_Type_vector(int count, int length, int stride, MPI_Datatype oldtype, MPI_Datatype *newtype);
int MPI_Type_hvector(int count, int length, int stride, MPI_Datatype oldtype, MPI_Datatype *newtype);
int MPI_Type_create_hvector(int count, int length, int stride, MPI_Datatype oldtype, MPI_Datatype *newtype);

inline int MPI_Init(int *a, char ***b) { return MPI_SUCCESS; }
inline int MPI_Finalize(void) { return MPI_SUCCESS; }
inline int MPI_Comm_size(MPI_Comm comm, int *psize) { *psize = 1; return MPI_SUCCESS; }
inline int MPI_Comm_rank(MPI_Comm comm, int *rank) { *rank = 0; return MPI_SUCCESS; }
inline int MPI_Cart_rank(MPI_Comm comm, int *coords, int *rank) { *rank = 0; return MPI_SUCCESS; }
inline int MPI_Cart_coords(MPI_Comm comm, int rank, int maxdims, int *coords) { coords[0] = coords[1] = coords[2] = 0; return MPI_SUCCESS; }
inline int MPI_Cart_shift(MPI_Comm comm, int direction, int disp, \
			    int *rank_source, int *rank_dest){ *rank_source = *rank_dest = 0; return MPI_SUCCESS; }
inline int MPI_Cart_create(MPI_Comm comm_old, int ndims, int *dims, int *periods, int reorder, MPI_Comm *comm_cart)
{ *comm_cart = NULL; return MPI_SUCCESS; }
inline int MPI_Dims_create(int nnodes, int ndims, int *dims){ dims[0] = dims[1] = dims[2] = 1; return MPI_SUCCESS; }
inline int MPI_Comm_split(MPI_Comm comm, int colour, int key, MPI_Comm *newcomm) { return MPI_SUCCESS; }
inline int MPI_Comm_free(MPI_Comm *comm) { return MPI_SUCCESS; }
inline int MPI_Type_commit(MPI_Datatype *dtype) { return MPI_SUCCESS; }
inline int MPI_Type_free(MPI_Datatype *dtype) { free(*dtype); *dtype = NULL; return MPI_SUCCESS; }
inline int MPI_Type_extent(MPI_Datatype dtype, MPI_Aint *pextent) { *pextent = dtype->upper - dtype->lower; return MPI_SUCCESS; }
inline int MPI_Type_get_extent(MPI_Datatype dtype, MPI_Aint *lower, MPI_Aint *pextent) { 
  *lower = dtype->lower;
  *pextent = dtype->upper - dtype->lower; 
  return MPI_SUCCESS; 
}
inline int MPI_Barrier(MPI_Comm comm) { return MPI_SUCCESS; }
inline int MPI_Waitall(int count, MPI_Request *reqs, MPI_Status *stats) { return MPI_SUCCESS; }
inline int MPI_Wait(MPI_Request *reqs, MPI_Status *stats) { return MPI_SUCCESS; }
inline int MPI_Comm_create_errhandler(MPI_Handler_function *errfunc, MPI_Errhandler *errhdl) { return MPI_SUCCESS; }
inline int MPI_Comm_set_errhandler(MPI_Comm comm, MPI_Errhandler errhdl) { return MPI_SUCCESS; }
inline int MPI_Bcast(void *buff, int count, MPI_Datatype datatype, int root, MPI_Comm comm) { return MPI_SUCCESS; }

#ifndef GNU_MPIFAKE_DEBUG

inline int MPI_Recv(void *buf, int count, MPI_Datatype dtype, int src, int tag, MPI_Comm comm, MPI_Status *stat) {
  fprintf(stderr, "MPI_Recv on a single node\n"); errexit(); return MPI_SUCCESS; }
inline int MPI_Irecv(void *buf, int count, MPI_Datatype dtype, int src, int tag, MPI_Comm comm, MPI_Request *req) {
  fprintf(stderr, "MPI_Recv on a single node\n"); errexit(); return MPI_SUCCESS; }
inline int MPI_Send(void *buf, int count, MPI_Datatype dtype, int dst, int tag, MPI_Comm comm) {
  fprintf(stderr, "MPI_Recv on a single node\n"); errexit(); return MPI_SUCCESS; }
inline int MPI_Sendrecv(void *sbuf, int scount, MPI_Datatype stype, int dst, int stag, void *rbuf, int rcount, MPI_Datatype rtype, int src, int rtag, MPI_Comm comm, MPI_Status *stat) {
  fprintf(stderr, "MPI_Sendrecv on a single node\n"); errexit(); return MPI_SUCCESS;
}
inline int MPI_Isend(void *buf, int count, MPI_Datatype dtype, int dst, int tag, MPI_Comm comm, MPI_Request *req) {
  fprintf(stderr, "MPI_Recv on a single node\n"); errexit(); return MPI_SUCCESS; }

#else

inline int __MPI_ERR(char *func, char *file, int line) {
  fprintf(stderr, "%s on a single node at %s:%d\n", func, file, line);
  errexit(); return MPI_ERR_RANK;
}

#define MPI_Recv(buf, count, dtype, src, tag, comm, stat)  __MPI_ERR("MPI_Recv", __FILE__, __LINE__)
#define MPI_Irecv(buf, count, dtype, src, tag, comm, req) __MPI_ERR("MPI_IRecv", __FILE__, __LINE__)
#define MPI_Send(buf, count, dtype, dst, tag, comm) __MPI_ERR("MPI_Send", __FILE__, __LINE__)
#define MPI_Isend(buf, count, dtype, dst, tag, comm, req) __MPI_ERR("MPI_Isend", __FILE__, __LINE__)
#define MPI_Sendrecv(sbuf, scount, stype, dst, stag, rbuf, rcount, rtype, src, rtag, comm, stat) \
  __MPI_ERR("MPI_Sendrecv", __FILE__, __LINE__)

#endif

inline int MPI_Gather(void *sbuf, int scount, MPI_Datatype sdtype,
			void *rbuf, int rcount, MPI_Datatype rdtype,
			int root, MPI_Comm comm)
{ return mpifake_sendrecv(sbuf, scount, sdtype, rbuf, rcount, rdtype); }
inline int MPI_Allgather(void *sbuf, int scount, MPI_Datatype sdtype,
			   void *rbuf, int rcount, MPI_Datatype rdtype,
			   MPI_Comm comm)
{ return mpifake_sendrecv(sbuf, scount, sdtype, rbuf, rcount, rdtype); }
inline int MPI_Scatter(void *sbuf, int scount, MPI_Datatype sdtype,
			 void *rbuf, int rcount, MPI_Datatype rdtype,
			 int root, MPI_Comm comm)
{ return mpifake_sendrecv(sbuf, scount, sdtype, rbuf, rcount, rdtype); }
inline int MPI_Op_create(MPI_User_function func, int commute, MPI_Op *pop) { *pop = func; return MPI_SUCCESS; }
inline int MPI_Reduce(void *sbuf, void* rbuf, int count, MPI_Datatype dtype, MPI_Op op, int root, MPI_Comm comm)
{ op(sbuf, rbuf, &count, &dtype); return MPI_SUCCESS; }
inline int MPI_Allreduce(void *sbuf, void *rbuf, int count, MPI_Datatype dtype, MPI_Op op, MPI_Comm comm)
{ if(sbuf == MPI_IN_PLACE)
    return MPI_SUCCESS; 
  op(sbuf, rbuf, &count, &dtype); return MPI_SUCCESS; }
inline int MPI_Error_string(int errcode, char *string, int *len) { *string = 0; *len = 0; return MPI_SUCCESS; }

#endif
