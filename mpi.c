#include "mpi.h"

struct mpifake_dtype
  mpifake_dtype_int    = { 1, 1, sizeof(int), sizeof(int), NULL },
  mpifake_dtype_double = { 1, 1, sizeof(double), sizeof(double), NULL},
  mpifake_dtype_char   = { 1, 1, sizeof(char), sizeof(char), NULL },
  mpifake_dtype_byte   = { 1, 1, sizeof(char), sizeof(char), NULL },
  mpifake_dtype_long   = { 1, 1, sizeof(long), sizeof(long), NULL };


int MPI_Type_vector(int count, int length, int stride,
		    MPI_Datatype oldtype, MPI_Datatype *newtype)
{
  struct mpifake_dtype *dtype;
  dtype = *newtype = malloc(sizeof(struct mpifake_dtype));
  dtype->count  = count;
  dtype->length = length;
  dtype->stride = stride*oldtype->tsize;
  dtype->tsize  = count*dtype->stride;
  dtype->parent = oldtype;
  return MPI_SUCCESS;
}

int MPI_Type_hvector(int count, int length, int stride,
		     MPI_Datatype oldtype, MPI_Datatype *newtype)
{
  struct mpifake_dtype *dtype;
  dtype = *newtype = malloc(sizeof(struct mpifake_dtype));
  dtype->count  = count;
  dtype->length = length;
  dtype->stride = stride;
  dtype->tsize  = count*dtype->stride;
  dtype->parent = oldtype;
  return MPI_SUCCESS;
}

int MPI_Type_contiguous(int count, MPI_Datatype oldtype, MPI_Datatype *newtype)
{
  struct mpifake_dtype *dtype;
  dtype = *newtype = malloc(sizeof(struct mpifake_dtype));
  dtype->count  = count;
  dtype->length = 1;
  dtype->stride = oldtype->tsize;
  dtype->tsize = count*dtype->stride;
  dtype->parent = oldtype;
  return MPI_SUCCESS;
}

void mpifake_copy(void *from, void *to, int *count, MPI_Datatype *dtype)
{
  MPI_Datatype dt = *dtype, parent = dt->parent;
  int
    length = dt->length,
    stride = dt->stride,
    tsize  = dt->tsize,
    cnt    = dt->tsize**count;
  int c;
  for (c = 0; c < cnt; c += stride) {
    if (parent)
      mpifake_copy((char *)from + c, (char *)to + c, &length, &parent);
    else
      memcpy((char *)to + c, (char *)from + c, length*tsize);
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
