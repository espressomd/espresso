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
#ifndef _UTILS_H
#define _UTILS_H
/** \file utils.hpp
 *    Small functions that are useful not only for one modul.

 *  just some nice utilities... 
 *  and some constants...
 *
*/
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "config.hpp"
#include "debug.hpp"
#include "errorhandling.hpp"

/*************************************************************/
/** \name Mathematical, physical and chemical constants.     */
/*************************************************************/
/*@{*/
/** Pi. */
#define PI     3.14159265358979323846264338328 
/** Square root of Pi */
#define wupi   1.77245385090551602729816748334 
/** One over square root of Pi. */
#define wupii  0.56418958354775627928034964498 
/** Pi to the power 1/3. */
#define driwu2 1.25992104989487316476721060728 

/// error code if no error occured
#define ES_OK    0
/// error code if an error occured
#define ES_ERROR 1

/** space necessary for an (64-bit) integer with sprintf.
    Analog to Tcl
 */
#define ES_INTEGER_SPACE 24
/** space necessary for an double with sprintf. Precision
    is 17 digits, plus sign, dot, e, sign of exponent and
    3 digits exponent etc. Analog to Tcl
*/
#define ES_DOUBLE_SPACE 27

/*@}*/

/************************************************
 * data types
 ************************************************/

/** Integer list. 
    Use the functions specified in list operations. */
typedef struct {
  /** Dynamically allocated integer field. */
  int *e;
  /** number of used elements in the integer field. */
  int n;
  /** allocated size of the integer field. This value is ONLY changed
      in the routines specified in list operations ! */
  int max;
} IntList;

/** Double list.
    Use the functions specified in list operations. */
typedef struct {
  /** Dynamically allocated double field. */
  double *e;
  /** number of used elements in the double field. */
  int n;
  /** allocated size of the double field. This value is ONLY changed
      in the routines specified in list operations ! */
  int max;
} DoubleList;

/*************************************************************/
/** \name Dynamic memory allocation.                         */
/*************************************************************/
/*@{*/

/* to enable us to make sure that freed pointers are invalidated, we normally try to use realloc.
   Unfortunately allocating zero bytes (which should be avoided) actually allocates 16 bytes, and
   reallocating to 0 also. To avoid this, we use our own malloc and realloc procedures. */
#ifndef MEM_DEBUG

#ifdef realloc
#undef realloc
#endif

#ifdef malloc
#undef malloc
#endif

/** used instead of realloc.
    Makes sure that resizing to zero FREEs pointer */
inline void *prealloc(void *old, int size) {
  void *p;
  if (size <= 0) {
    free(old);
    return NULL;
  }
  p = (void *)realloc(old, size);
  if(p == NULL) {
    fprintf(stderr, "Could not allocate memory.\n");
    errexit();
  }
  return p;
}

/** used instead of malloc.
    Makes sure that a zero size allocation returns a NULL pointer */
inline void *pmalloc(int size)
{
  void *p;
  if (size <= 0) {
    return NULL;
  }
  p = (void *)malloc(size);
  if(p == NULL) {
    fprintf(stderr, "Could not allocate memory.\n");
    errexit();
  }
  return p;
}

/** use our own realloc which makes sure that realloc(0) is actually a free. */
#define realloc prealloc

/** use our own malloc which makes sure that malloc(0) returns NULL. */
#define malloc pmalloc

#endif
/*@}*/

/*************************************************************/
/* mass helper macro                                         */
/*************************************************************/

#ifdef MASS
/** macro for easy use of mass. If masses are not switched on, the particle mass is defined to 1,
    so it should be compiled out in most cases. */
#define PMASS(pt) (pt).p.mass
#else
#define PMASS(pt) 1
#endif

/*************************************************************/
/** \name List operations .                                  */
/*************************************************************/
/*@{*/

/** Initialize an \ref IntList.  */
inline void init_intlist(IntList *il)
{
  il->n   = 0;
  il->max = 0;
  il->e   = NULL;
}
extern int this_node;

/** Allocate an \ref IntList of size size. If you need an \ref IntList
    with variable size better use \ref realloc_intlist */
inline void alloc_intlist(IntList *il, int size)
{
  il->max = size;
  il->e = (int *) malloc(sizeof(int)*il->max);
}

/** Reallocate an \ref IntList */
inline void realloc_intlist(IntList *il, int size)
{
  if(size != il->max) {
    il->max = size;
    il->e = (int *) realloc(il->e, sizeof(int)*il->max);
  }
}

/** Allocate an \ref IntList, but only to multiples of grain. */
inline void alloc_grained_intlist(IntList *il, int size, int grain)
{
  il->max = grain*((size + grain - 1)/grain);
  il->e = (int *) malloc(sizeof(int)*il->max);
}

/** Reallocate an \ref IntList, but only to multiples of grain. */
inline void realloc_grained_intlist(IntList *il, int size, int grain)
{
  if(size >= il->max)
    il->max = grain*((size + grain - 1)/grain);
  else
    /* shrink not as fast, just lose half, rounded up */
    il->max = grain*(((il->max + size + 1)/2 +
		      grain - 1)/grain);

  il->e = (int *) realloc(il->e, sizeof(int)*il->max);
}

/** Check wether an \ref IntList contains the value c */
inline int intlist_contains(IntList *il, int c)
{
  int i;
  for (i = 0; i < il->n; i++)
    if (c == il->e[i]) return 1;
  return 0;
}

/** Initialize an \ref DoubleList.  */
inline void init_doublelist(DoubleList *il)
{
  il->n   = 0;
  il->max = 0;
  il->e   = NULL;
}

/** Allocate an \ref DoubleList of size size. If you need an \ref DoubleList
    with variable size better use \ref realloc_doublelist */
inline void alloc_doublelist(DoubleList *dl, int size)
{
  dl->max = size;
  dl->e = (double *) malloc(sizeof(double)*dl->max);
}

/** Reallocate an \ref DoubleList */
inline void realloc_doublelist(DoubleList *dl, int size)
{
  if(size != dl->max) {
    dl->max = size;
    dl->e = (double *) realloc(dl->e, sizeof(double)*dl->max);
  }
}

/** Allocate an \ref DoubleList, but only to multiples of grain. */
inline void alloc_grained_doublelist(DoubleList *dl, int size, int grain)
{
  dl->max = grain*((size + grain - 1)/grain);
  dl->e = (double *) malloc(sizeof(double)*dl->max);
}

/** Reallocate an \ref DoubleList, but only to multiples of grain. */
inline void realloc_grained_doublelist(DoubleList *dl, int size, int grain)
{
  if(size >= dl->max)
    dl->max = grain*((size + grain - 1)/grain);
  else
    /* shrink not as fast, just lose half, rounded up */
    dl->max = grain*(((dl->max + size + 1)/2 +
		      grain - 1)/grain);

  dl->e = (double *) realloc(dl->e, sizeof(double)*dl->max);
}
/*@}*/


/*************************************************************/
/** \name Mathematical functions.                            */
/*************************************************************/
/*@{*/

/** Calculates the maximum of 'double'-typed a and b, returning 'double'. */
inline double dmax(double a, double b) { return (a>b) ? a : b; }

/** Calculates the minimum of 'double'-typed a and b, returning 'double'. */
inline double dmin(double a, double b) { return (a<b) ? a : b; }

/** Calculates the maximum of 'int'-typed a and b, returning 'int'. */
inline int imax(int a, int b) { return (a>b) ? a : b; }

/** Calculates the minimum of 'int'-typed a and b, returning 'int'. */
inline int imin(int a, int b) { return (a<b) ? a : b; }

/** Check if a value is NaN. isnan() is only available in C++11 and C99, but not in C++98. **/
#ifndef isnan
  #define isnan(a) (a != a)
#endif

/** Calculates the remainder of a division */
inline double drem_down(double a, double b) { return a - floor(a/b)*b; }

/** vector difference */
inline void vector_subt(double res[3], double a[3], double b[3])
{
  int i;
  for (i=0;i<3;i++)
    res[i]=a[i]-b[i];
}

/** Very slow sort routine for small integer arrays. Sorts the values
    in decending order.  
    \param data   the integer array 
    \param size   size of the array
 */
inline void sort_int_array(int *data, int size)
{
  int i,j,tmp;
  for(i=0;i<size-1;i++)
    for(j=i+1;j<size;j++) {
      if(data[i]<data[j]) {
	tmp=data[i]; data[i]=data[j]; data[j]=tmp;
      }
    }
}

/** permute an interger array field of size size about permute positions. */
inline void permute_ifield(int *field, int size, int permute)
{
  int i,tmp;

  if(permute==0) return;
  if(permute<0) permute = (size + permute);
  while(permute>0) {
    tmp=field[0];
    for(i=1;i<size;i++) field[i-1] = field[i];
    field[size-1]=tmp;
    permute--;
  }
}

/** Mathematically rounds 'double'-typed x, returning 'double'. */
inline double dround(double x) { return floor(x+0.5); }

/** Calculates the SQuaRe of 'double' x, returning 'double'. */
inline double SQR(double x) { return x*x; }

/** approximates \f$ \exp(d^2) \mathrm{erfc}(d)\f$ by applying a formula from:
    Abramowitz/Stegun: Handbook of Mathematical Functions, Dover
    (9. ed.), chapter 7 */
inline double AS_erfc_part(double d)
{
#define AS_a1  0.254829592
#define AS_a2 -0.284496736
#define AS_a3  1.421413741
#define AS_a4 -1.453152027
#define AS_a5  1.061405429
#define AS_p   0.3275911
  double t;
  
  t = 1.0 / (1.0 + AS_p * d);
  
  return t * (AS_a1 + t * (AS_a2 + t * (AS_a3 + t * (AS_a4 + t * AS_a5) ) ) );
}

/** Calculates the sinc-function as sin(PI*x)/(PI*x).
 *
 * (same convention as in Hockney/Eastwood). In order to avoid
 * divisions by 0, arguments, whose modulus is smaller than epsi, will
 * be evaluated by an 8th order Taylor expansion of the sinc
 * function. Note that the difference between sinc(x) and this
 * expansion is smaller than 0.235e-12, if x is smaller than 0.1. (The
 * next term in the expansion is the 10th order contribution
 * PI^10/39916800 * x^10 = 0.2346...*x^12).  This expansion should
 * also save time, since it reduces the number of function calls to
 * sin().  
*/
inline double sinc(double d)
{
#define epsi 0.1

#define c2 -0.1666666666667e-0
#define c4  0.8333333333333e-2
#define c6 -0.1984126984127e-3
#define c8  0.2755731922399e-5

  double PId = PI*d, PId2;

  if (fabs(d)>epsi)
    return sin(PId)/PId;
  else {
    PId2 = SQR(PId);
    return 1.0 + PId2*(c2+PId2*(c4+PId2*(c6+PId2*c8)));
  }
}

/** factorizes small numbers up to a maximum of max factors. */
inline int calc_factors(int n, int *factors, int max)
{
  int f=2,i=0;
  while(n>1) {
    while(f<=n) {
      if(n%f==0) {
	if(i>=max) return 0;
	n /= f;
	factors[i]=f;
	i++;
	f=n;
      } 
      f++;
    }
    f=2;
  }
  return i;
}
/*@}*/


/*************************************************************/
/** \name Vector and matrix operations for three dimensons.  */
/*************************************************************/
/*@{*/

/** Subtracts vector v2 from vector v1 and stores resuld in vector dv */
inline void vecsub(double v1[3], double v2[3], double dv[3])
{
  dv[0] = v1[0] - v2[0];
  dv[1] = v1[1] - v2[1];
  dv[2] = v1[2] - v2[2];
}


/** calculates the length of a vector */
inline double normr(double v[3]) {
  double d2 = 0.0;
  int i;
  for(i=0;i<3;i++)
    d2 += SQR(v[i]);
  d2 = sqrtf(d2);
  return d2;
}

/** calculates the squared length of a vector */
inline double sqrlen(double v[3]) {
  double d2 = 0.0;
  int i;
  for(i=0;i<3;i++)
    d2 += SQR(v[i]);
  return d2;
}

/** calculates unit vector */
inline void unit_vector(double v[3],double y[3]) {
  double d = 0.0;
  int i;
  d=sqrt( sqrlen(v) );
  
  for(i=0;i<3;i++)
    y[i] = v[i]/d;
    
  return;
}

/** calculates the scalar product of two vectors a nd b */
inline double scalar(double a[3], double b[3]) {
  double d2 = 0.0;
  int i;
  for(i=0;i<3;i++)
    d2 += a[i]*b[i];
  return d2;
}

/** calculates the vector product c of two vectors a and b */
inline void vector_product(double a[3], double b[3], double c[3]) {
  c[0]=a[1]*b[2]-a[2]*b[1];
  c[1]=a[2]*b[0]-a[0]*b[2];
  c[2]=a[0]*b[1]-a[1]*b[0];
  return ;
}
 
/** rotates vector around axis by alpha */
inline void vec_rotate(double *axis, double alpha, double *vector, double *result){
  double sina,cosa,absa,a[3];
  sina=sin(alpha);
  cosa=cos(alpha);
  absa=sqrt(scalar(axis,axis));

  a[0]=axis[0]/absa;
  a[1]=axis[1]/absa;
  a[2]=axis[2]/absa;

  result[0]=(cosa+SQR(a[0])*(1-cosa))*vector[0]+(a[0]*a[1]*(1-cosa)-a[2]*sina)*vector[1]+(a[0]*a[2]*(1-cosa)+a[1]*sina)*vector[2];
  result[1]=(a[0]*a[1]*(1-cosa)+a[2]*sina)*vector[0]+(cosa+SQR(a[1])*(1-cosa))*vector[1]+(a[1]*a[2]*(1-cosa)-a[0]*sina)*vector[2];
  result[2]=(a[0]*a[2]*(1-cosa)-a[1]*sina)*vector[0]+(a[1]*a[2]*(1-cosa)+a[0]*sina)*vector[1]+(cosa+SQR(a[2])*(1-cosa))*vector[2];

  return;
}

/** Calc eigevalues of a 3x3 matrix stored in q as a 9x1 array*/
inline int calc_eigenvalues_3x3(double *q,  double *eva) {
  double q11,q22,q33,q12,q13,q23;
  double a,b,c;
  double QQ,R,R2,QQ3;
  double theta,root1,root2,x1,x2,x3,help;
  int anzdiff=3;

  /* copy tensor to local variables (This is really a fuc off code!) */
  q11 = *q++;    q12 = *q++;   q13 = *q;
  q22 = *(q+=2); q23 = *(++q); q33 = *(q+=3);
  
  /* solve the cubic equation of the eigenvalue problem */
  a = -(q11 + q22 + q33);
  b = q11*q22 + q11*q33 + q22*q33 - q12*q12 - q13*q13 - q23*q23;
  c = -(q11*q22*q33 + 2*q12*q23*q13 - q12*q12*q33 - q13*q13*q22 - q23*q23*q11);
  QQ = (a*a-3*b)/9.0;
  R = (2*a*a*a - 9*a*b + 27*c)/54;
  QQ3 = QQ*QQ*QQ;
  R2 = R*R;

  if (R2 < QQ3) {
    theta = acos(R/sqrt(QQ3));
    root1 = sqrt(QQ)*cos(theta/3.0);
    root2 = sqrt(QQ*3.0)*sin(theta/3.0);
    x1 = - 2*root1 - a/3.0;
    x2 =   root1 - root2 - a/3.0;
    x3 =   root1 + root2 - a/3.0;
  }
  else {
    double AA,BB,signum = 1.0;

    if (R < 0) { signum = -1.0; R = -R; }
    AA = - signum*exp( log(R + sqrt(R2-QQ3))/3.0);
    if (AA==0.0) BB = 0.0;
    else BB = QQ/AA;

    /* compute the second way of the diagonalization
     * remark : a diagonal matrix has to have real eigenvalues (=>x2=x3)
     * but (!), the general solution of this case can have in general 
     * imaginary solutions for x2,x3 (see numerical recipies p. 185) */
    x1 = (AA+BB) - a/3.0;
    x2 = 0.5*(AA+BB) - a/3.0;
    x3 = x2;
  }

  /* order the eigenvalues in decreasing order */
  if (x1<x2) { help = x1; x1 = x2; x2 = help; }
  if (x1<x3) { help = x1; x1 = x3; x3 = help; }
  if (x2<x3) { help = x2; x2 = x3; x3 = help; }
  eva[0] = x1; eva[1] = x2; eva[2] = x3;

  /* calculate number of different eigenvalues */
  if(x1==x2) anzdiff--;
  if(x2==x3) anzdiff--;
  
  return anzdiff;
}

/** Calc eigevectors of a 3x3 matrix stored in a as a 9x1 array*/
/** Given an eigenvalue (eva) returns the corresponding eigenvector (eve)*/
inline int calc_eigenvector_3x3(double *a,double eva,double *eve) {
  int i,j,ind1,ind2,ind3,row1,row2;
  double A_x1[3][3],coeff1,coeff2,norm;

  /* build (matrix - eva*unity) */ 
  for(i=0;i<3;i++) {
    for(j=0;j<3;j++) A_x1[i][j]=a[3*i+j];
    A_x1[i][i]-=eva;
  }

  /* solve the linear equation system */
  for(ind1=0;ind1<3;ind1++) {
    ind2=(ind1+1) % 3;    ind3=(ind1+2) % 3;
    row1=0;
    do {
      row1++;
    } while ((row1<3) && (A_x1[ind3][row1]==0));

    /* okay, if one whole column is empty then we can take this 
       direction as the eigenvector
       remember : {eigenvectors} = kernel(A_x1) */
    if (row1 == 3) {
      eve[ind3]=1.0;
      eve[(ind3+1) % 3]=0.0;
      eve[(ind3+2) % 3]=0.0;
      return 1;
    }

    for(i=1;i<3;i++) {
      row2=(row1+i) % 3;
      coeff1=A_x1[ind1][row1]*A_x1[ind3][row2]-
                    A_x1[ind1][row2]*A_x1[ind3][row1];
      coeff2=A_x1[ind2][row1]*A_x1[ind3][row2]-
                    A_x1[ind2][row2]*A_x1[ind3][row1];
      if (coeff1!=0.0) {
        eve[ind2]=1.0;
        eve[ind1]=-coeff2/coeff1;
        eve[ind3]=-( A_x1[ind2][row1]+eve[ind1]* A_x1[ind1][row1])
                              /A_x1[ind3][row1];
	norm = sqrt(eve[0]*eve[0] + eve[1]*eve[1] + eve[2]*eve[2]);
        eve[0]/=norm; eve[1]/=norm; eve[2]/=norm;
        return 1;
      }
      else {
        if (coeff2!=0.0) {
          eve[ind1]=1.0;
          eve[ind2]=-coeff1/coeff2;
          eve[ind3]=-( A_x1[ind1][row1]+eve[ind2]* A_x1[ind2][row1])
                            /A_x1[ind3][row1];
	  norm = sqrt(eve[0]*eve[0] + eve[1]*eve[1] + eve[2]*eve[2]);
	  eve[0]/=norm; eve[1]/=norm; eve[2]/=norm;
          return(1);
        }
      }
    } /* loop over the different rows */
  }   /* loop over the different columns */

  /* the try failed => not a singular matrix: only solution is (0,0,0) */                    
  return 0;
}
/*@}*/

/*************************************************************/
/** \name Linear algebra functions                           */
/*************************************************************/
/*@{*/

/** Calculate the LU decomposition of a matrix A. Uses Crout's method
 *  with partial implicit pivoting.  The original matrix A is replaced
 *  by its LU decomposition.  Due to the partial pivoting the result
 *  may contain row permutations which are recorded in perms.
 *  @return 0 for success, -1 otherwise (i.e. matrix is singular)
 *  @param A     Matrix to be decomposed (Input/Output) 
 *  @param n     Dimension of the matrix (Input) 
 *  @param perms Records row permutations effected by pivoting (Output)
 */
inline int lu_decompose_matrix(double **A, int n, int *perms) {
  int i, j, k, ip;
  double max, sum, tmp;

  double *scal = (double *)malloc(n*sizeof(double));

  /* loop over rows and store implicit scaling factors */
  for (i=0; i<n; i++) {
    max = 0.0;
    for (j=0; j<n; j++) {
      if ( (tmp=fabs(A[i][j])) > max) max=tmp;
    }
    if (max == 0.0) {
      /* matrix has a zero row and is singular */
      return -1;
    }
    scal[i] = 1.0/max;
  }

  /** Crout's algorithm: Calculate L and U columnwise from top to
   *  bottom. The diagonal elements of L are chosen to be 1. Only
   *  previously determined entries are used in the calculation. The
   *  original matrix elements are used only once and can be
   *  overwritten with the elements of L and U, the diagonal of L not
   *  being stored. Rows may be permuted according to the largest
   *  element (pivot) in the lower part, where rows are normalized to
   *  have the largest element scaled to unity.
   */

  /* loop over columns */
  for (j=0; j<n; j++) {

    /* calculate upper triangle part (without diagonal) */
    for (i=0; i<j; i++) {
      sum = A[i][j];
      for (k=0;k<=i-1;k++) sum -= A[i][k]*A[k][j];
      A[i][j] = sum;
    }
    
    /* calculate diagonal and lower triangle part */
    /* pivot is determined on the fly, but not yet divided by */
    ip = j;
    max = 0.0;
    for (i=j; i<n; i++) {
      sum = A[i][j];
      for (k=0; k<=j-1; k++) sum -= A[i][k]*A[k][j];
      A[i][j] = sum;
      if ((tmp=scal[i]*fabs(sum)) > max) {
	max = tmp;
	ip = i;
      }
    }

    /* swap rows according to pivot index */
    if (j != ip) {
      for (k=0; k<n; k++) {
	tmp = A[j][k];
	A[j][k] = A[ip][k];
	A[ip][k] = tmp;
      }
      scal[ip] = scal[j];
    }
    perms[j] = ip;

    if (A[j][j] == 0.0) {
      /* zero pivot indicates singular matrix */
      return -1;
    }

    /* now divide by pivot element */
    if (j != n) {
      tmp = 1.0/A[j][j];
      for (i=j+1; i<n; i++)
	A[i][j] *= tmp;
    }

  }

  free(scal);
  
  return 0;
}

/** Solve the linear equation system for a LU decomposed matrix A.
 *  Uses forward substitution for the lower triangle part and
 *  backsubstitution for the upper triangle part. Row permutations as
 *  indicated in perms are applied to b accordingly. The solution is
 *  written to b in place.
 *  @param A     Matrix in LU decomposed form (Input) 
 *  @param n     Dimension of the matrix (Input) 
 *  @param perms Indicates row permutations due to pivoting (Input)
 *  @param b     Right-hand side of equation system (Input).
 *               Is destroyed and contains the solution x afterwards (Output).
 */
inline void lu_solve_system(double **A, int n, int *perms, double *b) {
  int i, j;
  double sum;

  /* Step 1: Solve Ly=b */

  /* forward substitution for lower triangle part */
  for (i=0; i<n; i++) {
    /* take care of the correct permutations */
    sum=b[perms[i]];
    b[perms[i]]=b[i];
    for (j=0; j<=i-1; j++) sum -= A[i][j]*b[j];
    b[i] = sum;
  }

  /* Step 2: Solve Ux=y */

  /* backsubstitution for upper triangle part */
  for (i=n-1; i>=0; i--) {
    sum = b[i];
    for (j=i+1; j<n; j++) sum -= A[i][j]*b[j];
    b[i]=sum/A[i][i];
  }

}
/*@}*/

/*************************************************************/
/** \name Three dimensional grid operations                  */
/*************************************************************/
/*@{*/

/** get the linear index from the position (a,b,c) in a 3D grid
 *  of dimensions adim[]. returns linear index.
 *
 * @return        the linear index
 * @param a       x position 
 * @param b       y position 
 * @param c       z position 
 * @param adim    dimensions of the underlying grid  
 */
inline int get_linear_index(int a, int b, int c, int adim[3])
{
  return (a + adim[0]*(b + adim[1]*c));   
}

/** get the position a[] from the linear index in a 3D grid
 *  of dimensions adim[].
 *
 * @param i       linear index
 * @param a       x position (return value) 
 * @param b       y position (return value) 
 * @param c       z position (return value) 
 * @param adim    dimensions of the underlying grid  
 */
inline void get_grid_pos(int i, int *a, int *b, int *c, int adim[3])
{
  *a = i % adim[0];
  i /= adim[0];
  *b = i % adim[1];
  i /= adim[1];
  *c = i;
}

/** Malloc a 3d grid for doubles with dimension dim[3] . 
 * @param grid    pointer to grid.
 * @param dim  dimension of the grid.
*/
inline int malloc_3d_grid(double ****grid, int dim[3])
{
  int i,j;
  *grid = (double***)malloc(sizeof(double **)*dim[0]);
  if(*grid==NULL) return 0;
  for(i=0;i<dim[0];i++) {
    (*grid)[i] = (double**)malloc(sizeof(double *)*dim[1]);
    if((*grid)[i]==NULL) return 0;
    for(j=0;j<dim[1];j++) {
      (*grid)[i][j] = (double*)malloc(sizeof(double)*dim[2]);
      if((*grid)[i][j]==NULL) return 0;
    }
  }
  return 1;
}

/** print a block of a 3D array.
 *  @param data    3D array.
 *  @param start   start coordinate for the block.
 *  @param size    size of the block.
 *  @param dim     dimension of the array.
 *  @param element size of the elements in the array.
 *  @param num     number of element to print.

*/
inline void print_block(double *data, int start[3], int size[3], int dim[3], int element, int num)
{
  int i0,i1,i2,b=1;
  int divide=0,block1=0,start1;
  double tmp;

  while(divide==0) {
    if(b*size[2] > 7) {
      block1=b;
      divide = (int)ceil(size[1]/(double)block1);
    }
    b++;
  }
  fprintf(stderr,"?: print_block (%d of %d): (%d,%d,%d)+(%d,%d,%d) from grid (%d,%d,%d)\n",
	  num+1,element,
	  start[0],start[1],start[2],size[0],size[1],size[2],dim[0],dim[1],dim[2]);
  for(b=0;b<divide;b++) {
    start1 = b*block1+start[1];
    for(i0=start[0]+size[0]-1; i0>=start[0]; i0--) {
      for(i1=start1; i1<imin(start1+block1,start[1]+size[1]);i1++) {
	for(i2=start[2]; i2<start[2]+size[2];i2++) {
	  tmp=data[num+(element*(i2+dim[2]*(i1+dim[1]*i0)))];
	  if(tmp<0) fprintf(stderr,"%1.2e",tmp);
	  else      fprintf(stderr," %1.2e",tmp);
	}
	fprintf(stderr," | ");
      }
      fprintf(stderr,"\n");
    }
    fprintf(stderr,"\n");
  }
}
/*@}*/


/*************************************************************/
/** \name Distance calculations.  */
/*************************************************************/
/*@{*/

/** returns the distance between two position. 
 *  \param pos1 Position one.
 *  \param pos2 Position two.
*/
inline double distance(double pos1[3], double pos2[3])
{
  return sqrt( SQR(pos1[0]-pos2[0]) + SQR(pos1[1]-pos2[1]) + SQR(pos1[2]-pos2[2]) );
}

/** returns the distance between two positions squared.
 *  \param pos1 Position one.
 *  \param pos2 Position two.
*/
inline double distance2(double pos1[3], double pos2[3])
{
  return SQR(pos1[0]-pos2[0]) + SQR(pos1[1]-pos2[1]) + SQR(pos1[2]-pos2[2]);
}

/** Returns the distance between two positions squared and stores the
    distance vector pos1-pos2 in vec.
 *  \param pos1 Position one.
 *  \param pos2 Position two.
 *  \param vec  vecotr pos1-pos2.
 *  \return distance squared
*/
inline double distance2vec(double pos1[3], double pos2[3], double vec[3])
{
  vec[0] = pos1[0]-pos2[0];
  vec[1] = pos1[1]-pos2[1];
  vec[2] = pos1[2]-pos2[2];
  return SQR(vec[0]) + SQR(vec[1]) + SQR(vec[2]);
}

/** returns the distance between the unfolded coordintes of two particles. 
 *  \param pos1       Position of particle one.
 *  \param image_box1 simulation box index of particle one .
 *  \param pos2       Position of particle two.
 *  \param image_box2 simulation box index of particle two .
 *  \param box_l      size of simulation box.
*/
inline double unfolded_distance(double pos1[3], int image_box1[3], 
				  double pos2[3], int image_box2[3], double box_l[3])
{
  int i;
  double dist = 0;
  double lpos1[3],lpos2[3];
  for(i=0;i<3;i++){
    lpos1[i] = pos1[i];
    lpos2[i] = pos2[i];
    lpos1[i] += image_box1[i]*box_l[i];
    lpos2[i] += image_box2[i]*box_l[i];
    dist += SQR(lpos1[i]-lpos2[i]);
  }
  return sqrt(dist);
}
/*@}*/

/*************************************************************/
/** \name String helper functions                            */
/*************************************************************/
/*@{*/

/** extend a string with another one. Like strcat, just automatically
    increases the string space */
inline char *strcat_alloc(char *left, const char *right)
{
  if (!left) {
    char *res = (char *)malloc(strlen(right) + 1);
    strcpy(res, right);
    return res;
  }
  else {
    char *res = (char *)realloc(left, strlen(left) + strlen(right) + 1);
    strcat(res, right);
    return res;
  }
}

/*@}*/

/*************************************************************/
/** \name Object-in-fluid functions                          */
/*************************************************************/
/*@{*/

/** Computes the area of triangle between vectors P1 and P2, 
 *  by computing the crossproduct P1 x P2 and taking the half of its norm */
inline double area_triangle_new(double *P1, double *P2) {
 double area;
 double normal[3], n; //auxiliary variables
 vector_product(P1,P2,normal); 
 n=normr(normal);
 area = 0.5*n;
 return area;
}

/** Computes the area of triangle between vectors P1,P2,P3, 
 *  by computing the crossproduct P1P2 x P1P3 and taking the half of its norm */
inline double area_triangle(double *P1, double *P2, double *P3) {
	double area;
	double u[3],v[3],normal[3], n; //auxiliary variables
	u[0] = P2[0] - P1[0]; // u = P1P2
	u[1] = P2[1] - P1[1]; 
	u[2] = P2[2] - P1[2]; 
	v[0] = P3[0] - P1[0]; // v = P1P3
	v[1] = P3[1] - P1[1]; 
	v[2] = P3[2] - P1[2]; 
	vector_product(u,v,normal); 
	n=normr(normal);
	area = 0.5*n;
	return area;
}

/** Computes the normal vector to the plane given by points P1P2P3 */
inline void get_n_triangle(double* p1, double* p2, double* p3, double* n){
	n[0]=(p2[1]-p1[1])*(p3[2]-p1[2])-(p2[2]-p1[2])*(p3[1]-p1[1]);
	n[1]=(p2[2]-p1[2])*(p3[0]-p1[0])-(p2[0]-p1[0])*(p3[2]-p1[2]);
	n[2]=(p2[0]-p1[0])*(p3[1]-p1[1])-(p2[1]-p1[1])*(p3[0]-p1[0]);
}

/** This function returns the angle btw the triangle p1,p2,p3 and p2,p3,p4. 
 *  Be careful, the angle depends on the orientation of the trianlges! 
 *  You need to be sure that the orientation (direction of normal vector) 
 *  of p1p2p3 is given by the cross product p2p1 x p2p3. 
 *  The orientation of p2p3p4 must be given by p2p3 x p2p4. 
 * 
 *  Example: p1 = (0,0,1), p2 = (0,0,0), p3=(1,0,0), p4=(0,1,0). 
 *  The orientation of p1p2p3 should be in the direction (0,1,0) 
 *  and indeed: p2p1 x p2p3 = (0,0,1)x(1,0,0) = (0,1,0)
 *  This function is called in the beginning of the simulation when creating 
 *  bonds depending on the angle btw the triangles, the bending_force.
 *  Here, we determine the orientations by looping over the triangles 
 *  and checking the correct orientation. So when defining the bonds by tcl command
 *  "part p2 bond xxxx p1 p3 p4", we correctly input the particle id's.
 *  So if you have the access to the order of particles, you are safe to call this
 *  function with exactly this order. Otherwise you need to check the orientations. */
inline double angle_btw_triangles(double *P1, double *P2, double *P3, double *P4) {
	double phi;
	double u[3],v[3],normal1[3],normal2[3]; //auxiliary variables
	u[0] = P1[0] - P2[0]; // u = P2P1
	u[1] = P1[1] - P2[1]; 
	u[2] = P1[2] - P2[2]; 
	v[0] = P3[0] - P2[0]; // v = P2P3
	v[1] = P3[1] - P2[1]; 
	v[2] = P3[2] - P2[2]; 
	vector_product(u,v,normal1); 
	u[0] = P3[0] - P2[0]; // u = P2P3
	u[1] = P3[1] - P2[1]; 
	u[2] = P3[2] - P2[2]; 
	v[0] = P4[0] - P2[0]; // v = P2P4
	v[1] = P4[1] - P2[1]; 
	v[2] = P4[2] - P2[2]; 
	vector_product(u,v,normal2); 

	double tmp11,tmp22,tmp33;
	// Now we compute the scalar product of n1 and n2 divided by the norms of n1 and n2
	//tmp11 = dot(3,normal1,normal2);         // tmp11 = n1.n2
	tmp11 = scalar(normal1,normal2);         // tmp11 = n1.n2
	
	/*	
	tmp22 = normr(normal1);
	tmp33 = normr(normal2);
	tmp11 /= (tmp22*tmp33);  // tmp11 = n1.n2/(|n1||n2|)
*/
	tmp11 *= fabs(tmp11); // tmp11 = (n1.n2)^2
	tmp22 = sqrlen(normal1);  //tmp22 = |n1|^2
	tmp33 = sqrlen(normal2);  //tmp33 = |n1|^2
	tmp11 /= (tmp22*tmp33);  // tmp11 = (n1.n2/(|n1||n2|))^2
	if (tmp11 > 0 ) {
		tmp11 = sqrt(tmp11);
	} else {
		tmp11 = - sqrt(- tmp11);
	}		

	if(tmp11>=1.) { tmp11=0.0;}
	else if(tmp11<=-1.) { tmp11=M_PI;}
	phi = M_PI - acos(tmp11); 	// The angle between the faces (not considering the orientation, always less or equal to Pi) is
								// equal to Pi minus angle between the normals
	
	// Now we need to determine, if the angle btw two triangles is less than Pi or more than Pi. To do this we check, 
	// if the point P4 lies in the halfspace given by trianlge P1P2P3 and the normal to this triangle. If yes, we have 
	// angle less than Pi, if not, we have angle more than Pi.
	// General equation of the plane is n_x*x + n_y*y + n_z*z + d = 0 where (n_x,n_y,n_z) is the normal to the plane.
	// Point P1 lies in the plane, therefore d = -(n_x*P1_x + n_y*P1_y + n_z*P1_z)
	// Point P4 lies in the halfspace given by normal iff n_x*P4_x + n_y*P4_y + n_z*P4_z + d >= 0
	tmp11 = - (normal1[0]*P1[0] + normal1[1]*P1[1] + normal1[2]*P1[2]);
	if (normal1[0]*P4[0] + normal1[1]*P4[1] + normal1[2]*P4[2] + tmp11 < 0) phi = 2*M_PI - phi;
	return phi;
}

/*@}*/

#endif
