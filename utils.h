#ifndef UTILS_H
#define UTILS_H
/** \file utils.h
 *    Small functions that are useful not only for one modul.
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:limbach@mpip-mainz.mpg.de">Hanjo</a>
 *
 *  \todo General realloc function for dynamic arrays. So far these functions are spread all over the place and redone in nearly every modul.
 *  \todo General Send/Recv routine for two-step communications.
*/
#include <math.h>

/* just some nice utilities... */
/* and some constants...*/

#define PI     3.14159265358979323846264338328 /* Pi */
#define wupi   1.77245385090551602729816748334 /* root of PI */
#define driwu2 1.25992104989487316476721060728 /* third root of 2 */

/** exit ungracefully, core dump if switched on. Defined in main.c. */
void errexit();

/** Maximum von  double a und  double b berechnen, gibt double zurueck. */
MDINLINE double dmax(double a, double b) { return (a>b) ? a : b; }

/** Maximum von  double a und  double b berechnen, gibt double zurueck. */
MDINLINE double dmin(double a, double b) { return (a<b) ? a : b; }

/** Maximum von  int a und  int b berechnen, gibt int zurueck. */
MDINLINE int imax(int a, int b) { return (a>b) ? a : b; }

/** Maximum von  int a und  int b berechnen, gibt int zurueck. */
MDINLINE int imin(int a, int b) { return (a<b) ? a : b; }

/** Wert von double x runden, gibt double zurueck. */
MDINLINE double dround(double x) { return floor(x+0.5); }

/** Berechnet das Quadrat von double x, gibt double zurueck. */
MDINLINE double SQR(double x) { return x*x; }

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
 * @param dim[3]  dimension of the grid.
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

/** print a block of a 3D array.
 *  @param data   3D array.
 *  @param start  start coordinate for the block.
 *  @param size   size of the block.
 *  @param dim    dimension of the array.
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
#endif
