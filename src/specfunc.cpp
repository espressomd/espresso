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
/** \file specfunc.cpp
    Special functions, see \ref specfunc.hpp "specfunc.h"
*/
#include <cmath>
#include "utils.hpp"
#include "specfunc.hpp"
#include "polynom.hpp"

/* Original gsl header
 * specfunc/bessel_K0.cpp
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/* Original Author:  G. Jungman */

/************************************************
 * chebychev expansions
 ************************************************/
/* Note that the first coefficient already includes the constant offsets */

/* based on SLATEC bk0(), bk0e() */

/* chebyshev expansions 

 series for bk0        on the interval  0.	    to  4.00000d+00
					with weighted error   3.57e-19
					 log weighted error  18.45
			       significant figures required  17.99
				    decimal places required  18.97

 series for ak0        on the interval  1.25000d-01 to  5.00000d-01
					with weighted error   5.34e-17
					 log weighted error  16.27
			       significant figures required  14.92
				    decimal places required  16.89

 series for ak02       on the interval  0.	    to  1.25000d-01
					with weighted error   2.34e-17
					 log weighted error  16.63
			       significant figures required  14.67
				    decimal places required  17.20
*/

static double bk0_data[11] = {
  -.5 -0.03532739323390276872,
   0.3442898999246284869, 
   0.03597993651536150163,
   0.00126461541144692592,
   0.00002286212103119451,
   0.00000025347910790261,
   0.00000000190451637722,
   0.00000000001034969525,
   0.00000000000004259816,
   0.00000000000000013744,
   0.00000000000000000035
};
static Polynom bk0_cs = { bk0_data, 11, 11 };

static double ak0_data[17] = {
  2.5 -0.07643947903327941,
  -0.02235652605699819,
   0.00077341811546938,
  -0.00004281006688886,
   0.00000308170017386,
  -0.00000026393672220,
   0.00000002563713036,
  -0.00000000274270554,
   0.00000000031694296,
  -0.00000000003902353,
   0.00000000000506804,
  -0.00000000000068895,
   0.00000000000009744,
  -0.00000000000001427,
   0.00000000000000215,
  -0.00000000000000033,
   0.00000000000000005
};
static Polynom ak0_cs = { ak0_data, 16, 16 };

static double ak02_data[14] = {
  2.5 -0.01201869826307592,
  -0.00917485269102569,
   0.00014445509317750,
  -0.00000401361417543,
   0.00000015678318108,
  -0.00000000777011043,
   0.00000000046111825,
  -0.00000000003158592,
   0.00000000000243501,
  -0.00000000000020743,
   0.00000000000001925,
  -0.00000000000000192,
   0.00000000000000020,
  -0.00000000000000002
};
static Polynom ak02_cs = { ak02_data, 13, 13 };

/* based on SLATEC besi0 */

/* chebyshev expansions

 series for bi0        on the interval  0.	    to  9.00000d+00
					with weighted error   2.46e-18
					 log weighted error  17.61
			       significant figures required  17.90
				    decimal places required  18.15

 series for ai0        on the interval  1.25000d-01 to  3.33333d-01
					with weighted error   7.87e-17
					 log weighted error  16.10
			       significant figures required  14.69
				    decimal places required  16.76


 series for ai02       on the interval  0.	    to  1.25000d-01
					with weighted error   3.79e-17
					 log weighted error  16.42
			       significant figures required  14.86
				    decimal places required  17.09
*/

static double bi0_data[12] = {
  5.5 -.07660547252839144951,
  1.92733795399380827000,
   .22826445869203013390, 
   .01304891466707290428,
   .00043442709008164874,
   .00000942265768600193,
   .00000014340062895106,
   .00000000161384906966,
   .00000000001396650044,
   .00000000000009579451,
   .00000000000000053339,
   .00000000000000000245
};
static Polynom bi0_cs = { bi0_data, 12, 12 };

static double ai0_data[21] = {
  .75 +.07575994494023796, 
   .00759138081082334,
   .00041531313389237,
   .00001070076463439,
  -.00000790117997921,
  -.00000078261435014,
   .00000027838499429,
   .00000000825247260,
  -.00000001204463945,
   .00000000155964859,
   .00000000022925563,
  -.00000000011916228,
   .00000000001757854,
   .00000000000112822,
  -.00000000000114684,
   .00000000000027155,
  -.00000000000002415,
  -.00000000000000608,
   .00000000000000314,
  -.00000000000000071,
   .00000000000000007
};
static Polynom ai0_cs = { ai0_data, 21, 21 };

static double ai02_data[22] = {
  .75 +.05449041101410882,
   .00336911647825569,
   .00006889758346918,
   .00000289137052082,
   .00000020489185893,
   .00000002266668991,
   .00000000339623203,
   .00000000049406022,
   .00000000001188914,
  -.00000000003149915,
  -.00000000001321580,
  -.00000000000179419,
   .00000000000071801,
   .00000000000038529,
   .00000000000001539,
  -.00000000000004151,
  -.00000000000000954,
   .00000000000000382,
   .00000000000000176,
  -.00000000000000034,
  -.00000000000000027,
   .00000000000000003
};
static Polynom ai02_cs = { ai02_data, 22, 22 };

/* based on SLATEC besk1(), besk1e() */

/* chebyshev expansions 

 series for bk1        on the interval  0.	    to  4.00000d+00
					with weighted error   7.02e-18
					 log weighted error  17.15
			       significant figures required  16.73
				    decimal places required  17.67

 series for ak1        on the interval  1.25000d-01 to  5.00000d-01
					with weighted error   6.06e-17
					 log weighted error  16.22
			       significant figures required  15.41
				    decimal places required  16.83

 series for ak12       on the interval  0.	    to  1.25000d-01
					with weighted error   2.58e-17
					 log weighted error  16.59
			       significant figures required  15.22
				    decimal places required  17.16
*/

static double bk1_data[11] = {
  1.5 +0.0253002273389477705,
  -0.3531559607765448760, 
  -0.1226111808226571480, 
  -0.0069757238596398643,
  -0.0001730288957513052,
  -0.0000024334061415659,
  -0.0000000221338763073,
  -0.0000000001411488392,
  -0.0000000000006666901,
  -0.0000000000000024274,
  -0.0000000000000000070
};

static Polynom bk1_cs = { bk1_data, 11, 11 };

static double ak1_data[17] = {
  2.5 +0.27443134069738830, 
   0.07571989953199368,
  -0.00144105155647540,
   0.00006650116955125,
  -0.00000436998470952,
   0.00000035402774997,
  -0.00000003311163779,
   0.00000000344597758,
  -0.00000000038989323,
   0.00000000004720819,
  -0.00000000000604783,
   0.00000000000081284,
  -0.00000000000011386,
   0.00000000000001654,
  -0.00000000000000248,
   0.00000000000000038,
  -0.00000000000000006
};
static Polynom ak1_cs = { ak1_data, 17, 17 };

static double ak12_data[14] = {
  2.5 +0.06379308343739001,
   0.02832887813049721,
  -0.00024753706739052,
   0.00000577197245160,
  -0.00000020689392195,
   0.00000000973998344,
  -0.00000000055853361,
   0.00000000003732996,
  -0.00000000000282505,
   0.00000000000023720,
  -0.00000000000002176,
   0.00000000000000215,
  -0.00000000000000022,
   0.00000000000000002
};
static Polynom ak12_cs = { ak12_data, 14, 14 };

/* based on SLATEC besi1(), besi1e() */

/* chebyshev expansions

 series for bi1        on the interval  0.	    to  9.00000d+00
					with weighted error   2.40e-17
					 log weighted error  16.62
			       significant figures required  16.23
				    decimal places required  17.14

 series for ai1        on the interval  1.25000d-01 to  3.33333d-01
					with weighted error   6.98e-17
					 log weighted error  16.16
			       significant figures required  14.53
				    decimal places required  16.82

 series for ai12       on the interval  0.	    to  1.25000d-01
				       with weighted error   3.55e-17
					log weighted error  16.45
			      significant figures required  14.69
				   decimal places required  17.12
*/

static double bi1_data[11] = {
  1.75 -0.001971713261099859,
   0.407348876675464810,
   0.034838994299959456,
   0.001545394556300123,
   0.000041888521098377,
   0.000000764902676483,
   0.000000010042493924,
   0.000000000099322077,
   0.000000000000766380,
   0.000000000000004741,
   0.000000000000000024
};
static Polynom bi1_cs = { bi1_data, 11, 11 };

static double ai1_data[21] = {
  .75 -0.02846744181881479,
  -0.01922953231443221,
  -0.00061151858579437,
  -0.00002069971253350,
   0.00000858561914581,
   0.00000104949824671,
  -0.00000029183389184,
  -0.00000001559378146,
   0.00000001318012367,
  -0.00000000144842341,
  -0.00000000029085122,
   0.00000000012663889,
  -0.00000000001664947,
  -0.00000000000166665,
   0.00000000000124260,
  -0.00000000000027315,
   0.00000000000002023,
   0.00000000000000730,
  -0.00000000000000333,
   0.00000000000000071,
  -0.00000000000000006
};
static Polynom ai1_cs = { ai1_data, 21, 21 };

static double ai12_data[22] = {
  .75 +0.02857623501828014,
  -0.00976109749136147,
  -0.00011058893876263,
  -0.00000388256480887,
  -0.00000025122362377,
  -0.00000002631468847,
  -0.00000000383538039,
  -0.00000000055897433,
  -0.00000000001897495,
   0.00000000003252602,
   0.00000000001412580,
   0.00000000000203564,
  -0.00000000000071985,
  -0.00000000000040836,
  -0.00000000000002101,
   0.00000000000004273,
   0.00000000000001041,
  -0.00000000000000382,
  -0.00000000000000186,
   0.00000000000000033,
   0.00000000000000028,
  -0.00000000000000003
};
static Polynom ai12_cs = { ai12_data, 22, 22 };

/************************************************
 * chebychev expansions
 ************************************************/

/* coefficients for Maclaurin summation in hzeta()
 * B_{2j}/(2j)!
 */
static double hzeta_c[15] = {
  1.00000000000000000000000000000,
  0.083333333333333333333333333333,
 -0.00138888888888888888888888888889,
  0.000033068783068783068783068783069,
 -8.2671957671957671957671957672e-07,
  2.0876756987868098979210090321e-08,
 -5.2841901386874931848476822022e-10,
  1.3382536530684678832826980975e-11,
 -3.3896802963225828668301953912e-13,
  8.5860620562778445641359054504e-15,
 -2.1748686985580618730415164239e-16,
  5.5090028283602295152026526089e-18,
 -1.3954464685812523340707686264e-19,
  3.5347070396294674716932299778e-21,
 -8.9535174270375468504026113181e-23
};

/************************************************
 * functions
 ************************************************/

double hzeta(double s, double q)
{
  double max_bits = 54.0;
  int jmax = 12, kmax = 10;
  int j, k;
  double pmax, scp, pcp, ans;

  if((s > max_bits && q < 1.0) || (s > 0.5*max_bits && q < 0.25))
    return pow(q, -s);
  if(s > 0.5*max_bits && q < 1.0) {
    double p1 = pow(q, -s);
    double p2 = pow(q/(1.0+q), s);
    double p3 = pow(q/(2.0+q), s);
    return p1 * (1.0 + p2 + p3);
  }
  /* Euler-Maclaurin summation formula 
   * [Moshier, p. 400, with several typo corrections]
   */
  pmax  = pow(kmax + q, -s);
  scp = s;
  pcp = pmax / (kmax + q);
  ans = pmax*((kmax+q)/(s-1.0) + 0.5);

  for(k=0; k<kmax; k++)
    ans += pow(k + q, -s);


  for(j=0; j<=jmax; j++) {
    double delta = hzeta_c[j+1] * scp * pcp;
    ans += delta;
    scp *= (s+2*j+1)*(s+2*j+2);
    pcp /= (kmax + q)*(kmax + q);
  }

  return ans;
}

double I0(double x)
{
  double c, y = fabs(x);
  if(y <= 3.0)
    return evaluateAsChebychevSeriesAt(&bi0_cs, y*y/4.5-1.0);

  c = (y <= 8.0) ?
    evaluateAsChebychevSeriesAt(&ai0_cs, (48.0/y-11.0)/5.0) :
    evaluateAsChebychevSeriesAt(&ai02_cs, 16.0/y-1.0);
  return exp(y) * c / sqrt(y);
}

double K0(double x)
{
  double c, I0;
  if(x <= 2.0) {
    c  = evaluateAsChebychevSeriesAt(&bk0_cs, 0.5*x*x-1.0);
    I0 = evaluateAsChebychevSeriesAt(&bi0_cs, x*x/4.5-1.0);
    return (-log(x) + M_LN2)*I0 + c;
  }
  c = (x <= 8.0) ?
    evaluateAsChebychevSeriesAt(&ak0_cs, (16.0/x-5.0)/3.0) :
    evaluateAsChebychevSeriesAt(&ak02_cs, 16.0/x-1.0);
  return exp(-x)*c/sqrt(x); 
}

double I1(double x)
{
  double c, y = fabs(x);
  if(y <= 3.0) {
    c = evaluateAsChebychevSeriesAt(&bi1_cs, y*y/4.5-1.0);
    return x*c;
  }
  c = (y <= 8.0) ?
    evaluateAsChebychevSeriesAt(&ai1_cs, (48.0/y-11.0)/5.0) :
    evaluateAsChebychevSeriesAt(&ai12_cs, 16.0/y-1.0);
  c = c / sqrt(y);
  if (x < 0)
    c = -c;
  return exp(y)*c;
}

double K1(double x)
{
  double c, I1;
  if(x <= 2.0) {
    c = evaluateAsChebychevSeriesAt(&bk1_cs, 0.5*x*x-1.0);
    I1 = x * evaluateAsChebychevSeriesAt(&bi1_cs, x*x/4.5-1.0);
    return (log(x) - M_LN2)*I1 + c/x;
  }
  c = (x <= 8.0) ?
    evaluateAsChebychevSeriesAt(&ak1_cs, (16.0/x-5.0)/3.0) :
    evaluateAsChebychevSeriesAt(&ak12_cs, 16.0/x-1.0);
  return exp(-x)*c/sqrt(x);
}

/***********************************************************
 * optimized K0/1 implementations for 10^(-14) precision
 ***********************************************************/

/** necessary orders for K0/1 from 2 up to 22 for 10^-14 precision. Note that at 8
    the expansion changes. From 23 to 26 order 2 is used, above order 1. For the
    latter cases separate implementations are necessary. */
static int ak01_orders[] = {
  /* 2 - 8 */
     11, 11, 10,
  10, 9, 9,
  /* 8 - 26 */
           6, 6,
  5, 5, 5, 4, 4,
  4, 3, 3, 2, 2,
  2, 2, 2
};

double LPK0(double x)
{
  if (x >= 27.) {
    double tmp = .5*exp(-x)/sqrt(x);
    return tmp*ak0_data[0];
  }
  if (x >= 23.) {
    double tmp = exp(-x)/sqrt(x), xx = (16./3.)/x - 5./3.;
    return tmp*(xx*ak0_data[1] + 0.5*ak0_data[0]);
  }
  if (x > 2) {
    int j = ak01_orders[((int)x) - 2];
    double tmp, x2;
    double dd0, d0;
    double *s0;
    if (x <= 8) {
      s0 = ak0_data;
      x2 = (2.*16./3.)/x - 2.*5./3.;
    } else {
      s0 = ak02_data;
      x2 = (2.*16.)/x - 2.;
    }
    dd0 = s0[j];
    d0  = x2*dd0 + s0[j - 1];
    for(j -= 2; j >= 1; j--) {
      double tmp0 = d0;
      d0 = x2*d0 - dd0 + s0[j];
      dd0 = tmp0;
    }
    tmp = exp(-x)/sqrt(x);
    return tmp*(0.5*(s0[0] + x2*d0) - dd0);
  }
  /* x <= 2 */
  {
    /* I0/1 series */
    int j = 10;
    double ret, tmp, x2 = (2./4.5)*x*x - 2.;
    double dd0, d0;
    dd0 = bi0_data[j];
    d0  = x2*dd0 + bi0_data[j - 1];
    for(j -= 2; j >= 1; j--) {
      double tmp0 = d0;
      d0 = x2*d0 - dd0 + bi0_data[j];
      dd0 = tmp0;
    }
    tmp = log(x) - M_LN2;
    ret =  -tmp*(0.5*(bi0_data[0] + x2*d0)- dd0);

    /* K0/K1 correction */
    j = 9;
    x2 = x*x - 2.;
    dd0 = bk0_data[j];
    d0  = x2*dd0 + bk0_data[j - 1];
    for(j -= 2; j >= 1; j--) {
      double tmp0 = d0;
      d0 = x2*d0 - dd0 + bk0_data[j];
      dd0 = tmp0;
    }
    return ret + (0.5*(x2*d0 + bk0_data[0]) - dd0);
  }
}

double LPK1(double x)
{
  if (x >= 27.) {
    double tmp = .5*exp(-x)/sqrt(x);
    return tmp*ak1_data[0];
  }
  if (x >= 23.) {
    double tmp = exp(-x)/sqrt(x), xx = (16./3.)/x - 5./3.;
    return tmp*(xx*ak1_data[1] + 0.5*ak1_data[0]);
  }
  if (x > 2) {
    int j = ak01_orders[((int)x) - 2];
    double tmp, x2;
    double dd1, d1;
    double *s1;
    if (x <= 8) {
      s1 = ak1_data;
      x2 = (2.*16./3.)/x - 2.*5./3.;
    } else {
      s1 = ak12_data;
      x2 = (2.*16.)/x - 2.;
    }
    dd1 = s1[j];
    d1  = x2*dd1 + s1[j - 1];
    for(j -= 2; j >= 1; j--) {
      double tmp1 = d1;
      d1 = x2*d1 - dd1 + s1[j];
      dd1 = tmp1;      
    }
    tmp = exp(-x)/sqrt(x);
    return tmp*(0.5*(s1[0] + x2*d1) - dd1);
  }
  /* x <= 2 */
  {
    /* I0/1 series */
    int j = 10;
    double ret, tmp, x2 = (2./4.5)*x*x - 2.;
    double dd1, d1;
    dd1 = bi1_data[j];
    d1  = x2*dd1 + bi1_data[j - 1];
    for(j -= 2; j >= 1; j--) {
      double tmp1 = d1;
      d1 = x2*d1 - dd1 + bi1_data[j];
      dd1 = tmp1;      
    }
    tmp = log(x) - M_LN2;
    ret = x*tmp*(0.5*(bi1_data[0] + x2*d1) - dd1);

    /* K0/K1 correction */
    j = 9;
    x2 = x*x - 2.;
    dd1 = bk1_data[j];
    d1  = x2*dd1 + bk1_data[j - 1];
    for(j -= 2; j >= 1; j--) {
      double tmp1 = d1;
      d1 = x2*d1 - dd1 + bk1_data[j];
      dd1 = tmp1;      
    }
    return ret + (0.5*(x2*d1 + bk1_data[0]) - dd1)/x;
  }
}

void LPK01(double x, double *K0, double *K1)
{
  if (x >= 27.) {
    double tmp = .5*exp(-x)/sqrt(x);
    *K0 = tmp*ak0_data[0];
    *K1 = tmp*ak1_data[0];
    return;
  }
  if (x >= 23.) {
    double tmp = exp(-x)/sqrt(x), xx = (16./3.)/x - 5./3.;
    *K0 = tmp*(xx*ak0_data[1] + 0.5*ak0_data[0]);
    *K1 = tmp*(xx*ak1_data[1] + 0.5*ak1_data[0]);
    return;    
  }
  if (x > 2) {
    int j = ak01_orders[((int)x) - 2];
    double tmp, x2;
    double dd0, dd1, d0, d1;
    double *s0, *s1;
    if (x <= 8) {
      s0 = ak0_data; s1 = ak1_data;
      x2 = (2.*16./3.)/x - 2.*5./3.;
    } else {
      s0 = ak02_data; s1 = ak12_data;
      x2 = (2.*16.)/x - 2.;
    }
    dd0 = s0[j];
    dd1 = s1[j];
    d0  = x2*dd0 + s0[j - 1];
    d1  = x2*dd1 + s1[j - 1];
    for(j -= 2; j >= 1; j--) {
      double tmp0 = d0, tmp1 = d1;
      d0 = x2*d0 - dd0 + s0[j];
      d1 = x2*d1 - dd1 + s1[j];
      dd0 = tmp0;
      dd1 = tmp1;      
    }
    tmp = exp(-x)/sqrt(x);
    *K0 = tmp*(0.5*(s0[0] + x2*d0) - dd0);
    *K1 = tmp*(0.5*(s1[0] + x2*d1) - dd1);
    return;
  }
  /* x <= 2 */
  {
    /* I0/1 series */
    int j = 10;
    double tmp, x2 = (2./4.5)*x*x - 2.;
    double dd0, dd1, d0, d1;
    dd0 = bi0_data[j];
    dd1 = bi1_data[j];
    d0  = x2*dd0 + bi0_data[j - 1];
    d1  = x2*dd1 + bi1_data[j - 1];
    for(j -= 2; j >= 1; j--) {
      double tmp0 = d0, tmp1 = d1;
      d0 = x2*d0 - dd0 + bi0_data[j];
      d1 = x2*d1 - dd1 + bi1_data[j];
      dd0 = tmp0;
      dd1 = tmp1;      
    }
    tmp = log(x) - M_LN2;
    *K0 =  -tmp*(0.5*(bi0_data[0] + x2*d0)- dd0);
    *K1 = x*tmp*(0.5*(bi1_data[0] + x2*d1) - dd1);

    /* K0/K1 correction */
    j = 9;
    x2 = x*x - 2.;
    dd0 = bk0_data[j];
    dd1 = bk1_data[j];
    d0  = x2*dd0 + bk0_data[j - 1];
    d1  = x2*dd1 + bk1_data[j - 1];
    for(j -= 2; j >= 1; j--) {
      double tmp0 = d0, tmp1 = d1;
      d0 = x2*d0 - dd0 + bk0_data[j];
      d1 = x2*d1 - dd1 + bk1_data[j];
      dd0 = tmp0;
      dd1 = tmp1;      
    }
    *K0 += (0.5*(x2*d0 + bk0_data[0]) - dd0);
    *K1 += (0.5*(x2*d1 + bk1_data[0]) - dd1)/x;
    return;
  }
}
