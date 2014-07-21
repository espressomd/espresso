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
/** \file uwerr.cpp
*/

#include <cmath>
#include <cstdlib>
#include <cfloat>
#include <cstdio>
#include "uwerr.hpp"

/*
enum UWerr_err_t {
  UW_NO_ERROR              = 0x000,
  UW_INVALID_DATA          = 0x001,
  UW_INVALID_ROWS          = 0x002,
  UW_INVALID_COLS          = 0x004,
  UW_INVALID_COL           = 0x008,
  UW_INVALID_FUNC          = 0x010,
  UW_INVALID_NREP          = 0x020,
  UW_INVALID_LEN           = 0x040,
  UW_NO_MEMORY             = 0x080,
  UW_WINDOWING_COND_FAILED = 0x100,
  UW_SIGMA_BIAS_CANCELLED  = 0x200
};
*/

/*
char * UWerr_err_to_str(struct UWerr_t s) {
  char * err_str[] = {"No error.",
		      "Invalid data.",
		      "Wrong number of rows for data.",
		      "Wrong number of columns for data.",
		      "Invalid function given.",
		      "Wrong vector for representations.",
		      "Wrong number for length of representations vector.",
		      "Out of memory.",
		      "Windowing condition failed at %d.",
		      "Sigma bias cancelled at %.3f."};

  return 0;
}
*/


static double gammaln(double a)
{
    int j;
    double x, y, tmp, ser;
    double cof[6] = { 76.18009172947146,
                     -86.50532032941677,
                      24.01409824083091,
                      -1.231739572450155,
                       0.1208650973866179e-2,
                      -0.5395239384953e-5};

    y = x = a;
    tmp = x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser = 1.000000000190015;
    for (j = 0; j < 6; ++j)
	ser += cof[j]/++y;

    return -tmp+log(2.5066282746310005*ser/x);
}

static void gammaser(double * gser, double a, double x)
{
    int n, ITMAX = 100;
    double sum, del, ap, gln, eps = DBL_EPSILON;

    gln = gammaln(a);

    if (x <= 0.0) {
	if (x < 0.0)
	    puts("uwerr: x less than 0 in gammaser.");
	*gser = 0.0;
	return;
    } else {
	ap = a;
	del = sum = 1.0/a;
	for (n = 0; n < ITMAX; ++n) {
	    ++ap;
	    del *= x/ap;
	    sum += del;
	    if (fabs(del) < fabs(sum)*eps) {
		*gser = sum*exp(-x+a*log(x)-gln);
		return;
	    }
	}
	puts("uwerr: a too large, ITMAX too small in gammaser.");
    }
}


static void gammacf(double * gcf, double a, double x)
{
    int i, ITMAX = 100;
    double an, b, c, d, del, h, gln, eps = DBL_EPSILON, FPMIN = DBL_MIN;

    gln=gammaln(a);
    b = x+1.0-a;
    c = 1.0/FPMIN;
    d = 1.0/b;
    h = d;
    for (i = 1; i <= ITMAX; ++i) {
	an = -i*(i-a);
	b += 2.0;
	d = an*d+b;
	if (fabs(d) < FPMIN)
	    d = FPMIN;
	c = b+an/c;
	if (fabs(c) < FPMIN)
	    c = FPMIN;
	d = 1.0/d;
	del = d*c;
	h *= del;
	if (fabs(del-1.0) <= eps)
	    break;
    }
    if (i > ITMAX)
	puts("uwerr: a too large, ITMAX too small in gammacf.");
    *gcf = exp(-x+a*log(x)-gln)*h;
}

/** The incomplete Gammafunction.

    This is the implementation of the incomplete Gammafunction as described in
    Nummerical Recepies in C++ The Art of Scientific Computing Second Edition.

    The incomplete Gammafunction is defined as
    \f[
       Q(a,x):=\frac{1}{\Gamma(a)}\int_x^\infty e^{-t}t^{a-1}\textrm{d}t\quad(a>0)
    \f]
 */

double gammaq(double a, double x)
{
    double retval=0;

    if (x < 0.0 || a <= 0.0) {
	puts("uwerr: Invalid arguments for gammaq.");
	return 0.0;
    }

    if (x < a+1.0) {
	gammaser(&retval, a, x);
	return 1.0-retval;
    }
    /* else */
    gammacf(&retval, a, x);
    return retval;
}

/** Sum up a vector.

    \param v A pointer to the vector.
    \param len The length of the vector.
    \return The sum of all vector elements.
 */
double UWerr_sum(double * v, int len)
{
  int i;
  double s = 0;
  
  if (v)
    for (i = 0; i < len; ++i)
      s += v[i];

  return s;
}

/** Sum up the product of two vectors.

    \f$$\sum_{k=1}^\textrm{len}v[k]w[k]$\f$

    \param v A pointer to the first vector.
    \param w A pointer to the second vector.
    \param len The length of the shortest vector.
    \return The sum of the products of vector elements.
 */
double UWerr_dsum_double(double * v, double * w, int len)
{
  int i;
  double s = 0;
  
  if (v && w)
    for (i = 0; i < len; ++i)
      s += v[i]*w[i];

  return s;
}

/** Sum up the product of two vectors.

    \f$$\sum_{k=1}^\textrm{len}v[k]w[k]$\f$

    \param v A pointer to the first vector.
    \param w A pointer to the second vector.
    \param len The length of the shortest vector.
    \return The sum of the products of vector elements.
 */
double UWerr_dsum_int(int * v, double * w, int len)
{
  int i;
  double s = 0;
  
  if (v && w)
    for (i = 0; i < len; ++i)
      s += v[i]*w[i];

  return s;
}
