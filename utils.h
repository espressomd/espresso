// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
#ifndef UTILS_H
#define UTILS_H
/** \file utils.h
 *    Small functions that are useful not only for one modul.

 *  just some nice utilities... 
 *  and some constants...
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:limbach@mpip-mainz.mpg.de">Hanjo</a>
 *
 *  \todo General realloc function for dynamic arrays. So far these functions are spread all over the place and redone in nearly every modul.
 *  \todo General Send/Recv routine for two-step communications.
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "debug.h"

/************************************************
 * defines
 ************************************************/

#define PI     3.14159265358979323846264338328 /* Pi */
#define wupi   1.77245385090551602729816748334 /* root of PI */
#define wupii  0.56418958354775627928034964498 /* 1 over root of PI */
#define driwu2 1.25992104989487316476721060728 /* third root of 2 */

/************************************************
 * data types
 ************************************************/

/** Integer list. */
typedef struct {
  /* Dynamically allocated integer field. */
  int *e;
  /** number of used elements in the integer field. */
  int n;
  /** allocated size of the integer field. */
  int max;
} IntList;

/** Double list. */
typedef struct {
  /* Dynamically allocated double field. */
  double *e;
  /** number of used elements in the double field. */
  int n;
  /** allocated size of the double field. */
  int max;
} DoubleList;


/************************************************
 * functions: misc
 ************************************************/

/** exit ungracefully, core dump if switched on. Defined in main.c. */
void errexit();

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
MDINLINE void *prealloc(void *old, int size)
{
  if (size == 0) {
    free(old);
    return NULL;
  }
  return realloc(old, size);
}

/** used instead of malloc.
    Makes sure that a zero size allocation returns a NULL pointer */
MDINLINE void *pmalloc(int size)
{
  if (size <= 0) {
    return NULL;
  }
  return malloc(size);
}

/** use our own realloc which makes sure that realloc(0) is actually a free. */
#define realloc prealloc

/** use our own malloc which makes sure that malloc(0) returns NULL. */
#define malloc pmalloc

#endif

/************************************************
 * functions: lists
 ************************************************/


MDINLINE void init_intlist(IntList *il)
{
  il->n   = 0;
  il->max = 0;
  il->e   = NULL;
}
extern int this_node;

MDINLINE void alloc_intlist(IntList *il, int size)
{
  il->max = size;
  il->e = (int *) malloc(sizeof(int)*il->max);
}

MDINLINE void realloc_intlist(IntList *il, int size)
{
  if(size != il->max) {
    il->max = size;
    il->e = (int *) realloc(il->e, sizeof(int)*il->max);
  }
}

MDINLINE void alloc_grained_intlist(IntList *il, int size, int grain)
{
  il->max = grain*((size + grain - 1)/grain);
  il->e = (int *) malloc(sizeof(int)*il->max);
}

MDINLINE void realloc_grained_intlist(IntList *il, int size, int grain)
{
  if(size >= il->max)
    il->max = grain*((size + grain - 1)/grain);
  else
    /* shrink not as fast, just lose half, rounded up */
    il->max = grain*(((il->max + size + 1)/2 +
		      grain - 1)/grain);

  il->e = (int *) realloc(il->e, sizeof(int)*il->max);
}

MDINLINE int intlist_contains(IntList *il, int c)
{
  int i;
  for (i = 0; i < il->n; i++)
    if (c == il->e[i]) return 1;
  return 0;
}

MDINLINE void init_doublelist(DoubleList *il)
{
  il->n   = 0;
  il->max = 0;
  il->e   = NULL;
}

MDINLINE void alloc_doublelist(DoubleList *dl, int size)
{
  dl->max = size;
  dl->e = (double *) malloc(sizeof(double)*dl->max);
}

MDINLINE void realloc_doublelist(DoubleList *dl, int size)
{
  if(size != dl->max) {
    dl->max = size;
    dl->e = (double *) realloc(dl->e, sizeof(double)*dl->max);
  }
}

MDINLINE void alloc_grained_doublelist(DoubleList *dl, int size, int grain)
{
  dl->max = grain*((size + grain - 1)/grain);
  dl->e = (double *) malloc(sizeof(double)*dl->max);
}

MDINLINE void realloc_grained_doublelist(DoubleList *dl, int size, int grain)
{
  if(size >= dl->max)
    dl->max = grain*((size + grain - 1)/grain);
  else
    /* shrink not as fast, just lose half, rounded up */
    dl->max = grain*(((dl->max + size + 1)/2 +
		      grain - 1)/grain);

  dl->e = (double *) realloc(dl->e, sizeof(double)*dl->max);
}

/************************************************
 * functions: math
 ************************************************/

/** Calculates the maximum of 'double'-typed a and b, returning 'double'. */
MDINLINE double dmax(double a, double b) { return (a>b) ? a : b; }

/** Calculates the minimum of 'double'-typed a and b, returning 'double'. */
MDINLINE double dmin(double a, double b) { return (a<b) ? a : b; }

/** Calculates the maximum of 'int'-typed a and b, returning 'int'. */
MDINLINE int imax(int a, int b) { return (a>b) ? a : b; }

/** Calculates the minimum of 'int'-typed a and b, returning 'int'. */
MDINLINE int imin(int a, int b) { return (a<b) ? a : b; }

/** Mathematically rounds 'double'-typed x, returning 'double'. */
MDINLINE double dround(double x) { return floor(x+0.5); }

/** Calculates the SQuaRe of 'double' x, returning 'double'. */
MDINLINE double SQR(double x) { return x*x; }

/** calculates the squared length of a vector */
MDINLINE double sqrlen(double v[3]) {
  double d2 = 0.0;
  int i;
  for(i=0;i<3;i++)
    d2 += SQR(v[i]);
  return d2;
}

/** calculates unit vector */
MDINLINE void unit_vector(double v[3],double y[3]) {
  double d = 0.0;
  int i;
  for(i=0;i<3;i++) 
    d += SQR(v[i]);
	 
	d=sqrt(d);
  
  for(i=0;i<3;i++)
    y[i] = v[i]/d;
    
  return;
}

/** calculates the scalar product of two vectors */
MDINLINE double scalar(double a[3], double b[3]) {
  double d2 = 0.0;
  int i;
  for(i=0;i<3;i++)
    d2 += a[i]*b[i];
  return d2;
}

/** calculates the vector product of two vectors */
MDINLINE void vector_product(double a[3], double b[3], double c[3]) {
  c[0]=a[1]*b[2]-a[2]*b[1];
  c[1]=a[2]*b[0]-a[0]*b[2];
  c[2]=a[0]*b[1]-a[1]*b[0];
  return ;
}

/** approximates exp(d^2)*erfc(d) by applying a formula from:
    Abramowitz/Stegun: Handbook of Mathematical Functions, 
    Dover (9. ed.), chapter 7 */
MDINLINE double AS_erfc_part(double d)
{
  static double AS_a1 =  0.254829592;
  static double AS_a2 = -0.284496736;
  static double AS_a3 =  1.421413741;
  static double AS_a4 = -1.453152027;
  static double AS_a5 =  1.061405429;
  static double AS_p  =  0.3275911;
  static double t;
  
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
MDINLINE double sinc(double d)
{
  static double epsi =  0.1;

  static double   c2 = -0.1666666666667e-0;
  static double   c4 =  0.8333333333333e-2;
  static double   c6 = -0.1984126984127e-3;
  static double   c8 =  0.2755731922399e-5;

  double PId = PI*d, PId2;

  if (fabs(d)>epsi)
    return sin(PId)/PId;
  else {
    PId2 = SQR(PId);
    return 1.0 + PId2*(c2+PId2*(c4+PId2*(c6+PId2*c8)));
  }
  return 1.0;
}

/** factorizes small numbers up to a maximum of max factors. */
MDINLINE int calc_factors(int n, int *factors, int max)
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

/************************************************
 * functions: vectors and matrices
 ************************************************/
 
/** Calc eigevalues of a 3x3 matrix stored in q as a 9x1 array*/
MDINLINE int calc_eigenvalues_3x3(double *q,  double *eva) {
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
MDINLINE int calc_eigenvector_3x3(double *a,double eva,double *eve) {
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
    } while ((A_x1[ind3][row1]==0) && (row1<3));

    /* okay, if one whole column is empty then we can take this 
       direction as the eigenvector
       remember : {eigenvectors} = kernel(A_x1) */
    if (A_x1[ind3][row1]==0) {
      eve[row1]=1.0;
      eve[(row1+1) % 3]=0.0;
      eve[(row1+2) % 3]=0.0;
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


/** get the linear index from the position (a,b,c) in a 3D grid
 *  of dimensions adim[]. returns linear index.
 *
 * @return        the linear index
 * @param a       x position 
 * @param b       y position 
 * @param c       z position 
 * @param adim    dimensions of the underlying grid  
 */
MDINLINE int get_linear_index(int a, int b, int c, int adim[3])
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
MDINLINE void get_grid_pos(int i, int *a, int *b, int *c, int adim[3])
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
MDINLINE int malloc_3d_grid(double ****grid, int dim[3])
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

/* very slow sort routine for small integer arrays. */
MDINLINE void sort_int_array(int *data,int size)
{
  int i,j,tmp;
  for(i=0;i<size-1;i++)
    for(j=i+1;j<size;j++) {
      if(data[i]<data[j]) {
	tmp=data[i]; data[i]=data[j]; data[j]=tmp;
      }
    }
}

MDINLINE void permute_ifield(int *field, int size, int permute)
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

/************************************************
 * functions: Debugging
 ************************************************/

/** print a block of a 3D array.
 *  @param data    3D array.
 *  @param start   start coordinate for the block.
 *  @param size    size of the block.
 *  @param dim     dimension of the array.
 *  @param element size of the elements in the array.
 *  @param num     number of element to print.

*/
MDINLINE void print_block(double *data, int start[3], int size[3], int dim[3], int element, int num)
{
  int i0,i1,i2,b=1;
  int divide=0,block1,start1;
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



/************************************************
 * functions: Distance calculation
 ************************************************/

/** returns the distance between two position. 
 *  \param pos1 Position one.
 *  \param pos2 Position two.
*/
MDINLINE double distance(double pos1[3], double pos2[3])
{
  return sqrt( SQR(pos1[0]-pos2[0]) + SQR(pos1[1]-pos2[1]) + SQR(pos1[2]-pos2[2]) );
}

/** returns the distance between two positions squared.
 *  \param pos1 Position one.
 *  \param pos2 Position two.
*/
MDINLINE double distance2(double pos1[3], double pos2[3])
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
MDINLINE double distance2vec(double pos1[3], double pos2[3], double vec[3])
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
MDINLINE double unfolded_distance(double pos1[3], int image_box1[3], 
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
#endif
