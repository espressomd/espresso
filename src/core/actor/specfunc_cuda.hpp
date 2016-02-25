/*
  Copyright (C) 2014,2015,2016 The ESPResSo project
  
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
#include "config.hpp"

const mmm1dgpu_real M_LN2f = M_LN2;

// Adapted from specfunc.c and polynom.h

__constant__ static mmm1dgpu_real bk0_data[11] = {
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
__constant__ static int bk0_size = 11;

__constant__ static mmm1dgpu_real ak0_data[17] = {
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
__constant__ static int ak0_size = 16;

__constant__ static mmm1dgpu_real ak02_data[14] = {
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
__constant__ static int ak02_size = 13;

__constant__ static mmm1dgpu_real bi0_data[12] = {
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
__constant__ static int bi0_size = 12;

__constant__ static mmm1dgpu_real bk1_data[11] = {
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
__constant__ static int bk1_size = 11;

__constant__ static mmm1dgpu_real ak1_data[17] = {
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
__constant__ static int ak1_size = 17;

__constant__ static mmm1dgpu_real ak12_data[14] = {
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
__constant__ static int ak12_size = 14;

__constant__ static mmm1dgpu_real bi1_data[11] = {
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
__constant__ static int bi1_size = 11 ;

__device__ mmm1dgpu_real evaluateAsChebychevSeriesAt(mmm1dgpu_real *c, int n, mmm1dgpu_real x)
{
  int j;
  mmm1dgpu_real x2 = 2 * x;
  mmm1dgpu_real dd = c[n - 1];
  mmm1dgpu_real d  = x2*dd + c[n - 2];
  for(j = n - 3; j >= 1; j--) {
    mmm1dgpu_real tmp = d;
    d = x2*d - dd + c[j];
    dd = tmp;
  }
  return x*d - dd + c[0]/2;
}

__device__ mmm1dgpu_real evaluateAsTaylorSeriesAt(mmm1dgpu_real *c, int n, mmm1dgpu_real x)
{
  int cnt = n - 1;
  mmm1dgpu_real r = c[cnt];
  while (--cnt >= 0)
    r = r*x + c[cnt];
  return r;
}

__device__ mmm1dgpu_real dev_K0(mmm1dgpu_real x)
{
	mmm1dgpu_real c = evaluateAsChebychevSeriesAt(
    x<=2?	bk0_data	:x<=8?	ak0_data	:	ak02_data,
    x<=2?	bk0_size	:x<=8?	ak0_size	:	ak02_size,
    x<=2?	x*x/2-1.0f	:x<=8?	(16/x-5.0f)/3.0f	:	(16/x-1.0f)
    );
  if (x <= 2) {
    mmm1dgpu_real I0 = evaluateAsChebychevSeriesAt(bi0_data, bi0_size, x*x/4.5f-1.0f);
    return (-log(x) + M_LN2f)*I0 + c;
  }
	return exp(-x)*c*rsqrt(x);
}

__device__ mmm1dgpu_real dev_K1(mmm1dgpu_real x)
{
  mmm1dgpu_real c = evaluateAsChebychevSeriesAt(
    x<=2? bk1_data  :x<=8?  ak1_data  : ak12_data,
    x<=2? bk1_size  :x<=8?  ak1_size  : ak12_size,
    x<=2? x*x/2-1.0f  :x<=8?  (16/x-5.0f)/3.0f  : (16/x-1.0f)
    );
	if(x <= 2) {
    mmm1dgpu_real I1 = x * evaluateAsChebychevSeriesAt(bi1_data, bi1_size, x*x/4.5f-1.0f);
		return (log(x) - M_LN2f)*I1 + c/x;
	}
	return exp(-x)*c*rsqrt(x);
}
