// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
/** \file mpi.h
    This is the MPIfake implementation. This is NOT a real MPI implementation, but rather implements
    a subset of the MPI commands such that Espresso is able to run with a single processor. For this to work,
    you do not have to have any MPI implementation like LAM or MPICH installed.
*/

#ifndef MPI_H
#define MPI_H
#include <string.h>
#include "utils.h"

/********************************** REMARK **********************/
/* This is the fake MPI header of Espresso, and has nothing to  */
/* with Espresso's MPI handling using other mpi variants        */
/********************************** REMARK **********************/

struct mpifake_dtype {
  /* blocks */
  int count, length;
  /* bytes */
  int stride, tsize;
  struct mpifake_dtype *parent;
};

typedef struct mpifake_dtype *MPI_Datatype;
typedef void *MPI_Status;
typedef void *MPI_Comm;
typedef void *MPI_Errhandler;
typedef void *MPI_Request;

typedef void (MPI_User_function)(void *, void *, int *, MPI_Datatype *);
typedef void (MPI_Handler_function)(MPI_Comm *, int *, ...);

typedef MPI_User_function *MPI_Op;

void mpifake_copy(void *from, void *to, int *count, MPI_Datatype *dtype);
int mpifake_checked_copy(void *s, int scount, MPI_Datatype sdtype,
			 void *r, int rcount, MPI_Datatype rdtype);

#define MPI_LOR mpifake_copy
#define MPI_SUM mpifake_copy
#define MPI_MAX mpifake_copy
#define MPI_COPY mpifake_copy

#define MPI_SUCCESS 1

#define MPI_COMM_WORLD NULL

#define MPI_REQUEST_NULL NULL

extern struct mpifake_dtype mpifake_dtype_int;
extern struct mpifake_dtype mpifake_dtype_double;
extern struct mpifake_dtype mpifake_dtype_byte;
extern struct mpifake_dtype mpifake_dtype_long;
extern struct mpifake_dtype mpifake_dtype_char;

#define MPI_INT    (&mpifake_dtype_int)
#define MPI_DOUBLE (&mpifake_dtype_double)
#define MPI_BYTE   (&mpifake_dtype_byte)
#define MPI_LONG   (&mpifake_dtype_long)
#define MPI_CHAR   (&mpifake_dtype_char)

int MPI_Type_vector(int count, int length, int stride, MPI_Datatype oldtype, MPI_Datatype *newtype);
int MPI_Type_hvector(int count, int length, int stride, MPI_Datatype oldtype, MPI_Datatype *newtype);
int MPI_Type_contiguous(int count, MPI_Datatype oldtype, MPI_Datatype *newtype);

MDINLINE int MPI_Init(int *a, char ***b) { return MPI_SUCCESS; }
MDINLINE int MPI_Finalize(void) { return MPI_SUCCESS; }
MDINLINE int MPI_Comm_size(MPI_Comm comm, int *psize) { *psize = 1; return MPI_SUCCESS; }
MDINLINE int MPI_Comm_rank(MPI_Comm comm, int *rank) { *rank = 0; return MPI_SUCCESS; }
MDINLINE int MPI_Comm_split(MPI_Comm comm, int colour, int key, MPI_Comm *newcomm) { return MPI_SUCCESS; }
MDINLINE int MPI_Comm_free(MPI_Comm *comm) { return MPI_SUCCESS; }
MDINLINE int MPI_Type_commit(MPI_Datatype *dtype) { return MPI_SUCCESS; }
MDINLINE int MPI_Type_free(MPI_Datatype *dtype) { free(*dtype); *dtype = NULL; return MPI_SUCCESS; }
MDINLINE int MPI_Barrier(MPI_Comm comm) { return MPI_SUCCESS; }
MDINLINE int MPI_Waitall(int count, MPI_Request *reqs, MPI_Status *stats) { return MPI_SUCCESS; }
MDINLINE int MPI_Errhandler_create(MPI_Handler_function *errfunc, MPI_Errhandler *errhdl) { return MPI_SUCCESS; }
MDINLINE int MPI_Errhandler_set(MPI_Comm comm, MPI_Errhandler errhdl) { return MPI_SUCCESS; }
MDINLINE int MPI_Bcast(void *buff, int count, MPI_Datatype datatype, int root, MPI_Comm comm) { return MPI_SUCCESS; }
MDINLINE int MPI_Recv(void *buf, int count, MPI_Datatype dtype, int src, int tag, MPI_Comm comm, MPI_Status *stat) {
  fprintf(stderr, "MPI_Recv on a single node\n"); errexit(); return MPI_SUCCESS; }
MDINLINE int MPI_Irecv(void *buf, int count, MPI_Datatype dtype, int src, int tag, MPI_Comm comm, MPI_Request *req) {
  fprintf(stderr, "MPI_Recv on a single node\n"); errexit(); return MPI_SUCCESS; }
MDINLINE int MPI_Send(void *buf, int count, MPI_Datatype dtype, int dst, int tag, MPI_Comm comm) {
  fprintf(stderr, "MPI_Recv on a single node\n"); errexit(); return MPI_SUCCESS; }
MDINLINE int MPI_Isend(void *buf, int count, MPI_Datatype dtype, int dst, int tag, MPI_Comm comm, MPI_Request *req) {
  fprintf(stderr, "MPI_Recv on a single node\n"); errexit(); return MPI_SUCCESS; }
MDINLINE int MPI_Gather(void *sbuf, int scount, MPI_Datatype sdtype,
			void *rbuf, int rcount, MPI_Datatype rdtype,
			int root, MPI_Comm comm)
{ return mpifake_checked_copy(sbuf, scount, sdtype, rbuf, rcount, rdtype); }
MDINLINE int MPI_Allgather(void *sbuf, int scount, MPI_Datatype sdtype,
			   void *rbuf, int rcount, MPI_Datatype rdtype,
			   MPI_Comm comm)
{ return mpifake_checked_copy(sbuf, scount, sdtype, rbuf, rcount, rdtype); }
MDINLINE int MPI_Scatter(void *sbuf, int scount, MPI_Datatype sdtype,
			 void *rbuf, int rcount, MPI_Datatype rdtype,
			 int root, MPI_Comm comm)
{ return mpifake_checked_copy(sbuf, scount, sdtype, rbuf, rcount, rdtype); }
MDINLINE int MPI_Op_create(MPI_User_function func, int commute, MPI_Op *pop) { *pop = func; return MPI_SUCCESS; }
MDINLINE int MPI_Reduce(void *sbuf, void* rbuf, int count, MPI_Datatype dtype, MPI_Op op, int root, MPI_Comm comm)
{ op(sbuf, rbuf, &count, &dtype); return MPI_SUCCESS; }
MDINLINE int MPI_Allreduce(void *sbuf, void *rbuf, int count, MPI_Datatype dtype, MPI_Op op, MPI_Comm comm)
{ op(sbuf, rbuf, &count, &dtype); return MPI_SUCCESS; }

#endif
