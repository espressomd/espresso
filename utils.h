#ifndef UTILS_H
#define UTILS_H
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

#endif
