/*
  Copyright (C) 2016 The ESPResSo project
  Copyright (C) 2012 Alexander (Polyakov) Peletskyi

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

#ifndef DIPOLARBARNESHUT_CUH_
#define DIPOLARBARNESHUT_CUH_

#include "config.hpp"

#ifdef BARNES_HUT

/**
 *
 */

//#include <cuda.h>
//#include <curand.h>
//#include <curand_kernel.h>
#include <cstdio>
#include <iostream>
#include <time.h>

typedef float dds_float ;

#define SHARED_ARRAY_BH 512

struct BHBox{
	float *minx;
	float *maxx;
	float *miny;
	float *maxy;
	float *minz;
	float *maxz;
};

struct BHArrays{
	int *err;
	int *sort;
	int *child;
	int *count;
	int *start;
};

// thread count
#define THREADS0 512	/* must be a power of 2 */
#define THREADS1 512	/* must be a power of 2 */
#define THREADS2 1024
#define THREADS3 1024
#define THREADS4 256
#define THREADS5 128	//256
#define THREADS6 512

// block count = factor * #SMs
#define FACTOR0 1
#define FACTOR1 3
#define FACTOR2 1
#define FACTOR3 1	/* must all be resident at the same time */
#define FACTOR4 1	/* must all be resident at the same time */
#define FACTOR5 5
#define FACTOR6 3

#define WARPSIZE 32
#define MAXDEPTH 32


void fillConstantPointers(float* rx, float* ry, float* rz,
		float* dipx, float* dipy, float* dipz,
		int nbodies, int nnodes, BHArrays arrl, BHBox boxl, float* mass);

void initBH(int blocks);

//void calcU(int blocks);

void buildBoxBH(int blocks);

void buildTreeBH(int blocks);

void summarizeBH(int blocks);

void sortBH(int blocks);

void forceBH(int blocks, dds_float k, float* f, float* torque, dds_float box_l[3],int periodic[3]);

void energyBH(int blocks, dds_float k, dds_float box_l[3],int periodic[3],float* E);

//void integrate(int blocks);

//void fillGloabalChangableBH(ChangableParams* chPar);

//void fillGloabalConstantBH(ConstParams* cPar);

//float3 averageU(int blocks);

//void getCurrentPhysBH(int blocks, float* devOldPhy);

//void getRingStatBH(int blocks, int nRings, int nPart, float** devOldPhy, float* dphy, int* pd);

#endif /* FFLCOMPUTEBH_CUH_ */
#endif // DIPOLAR_DIRECT_SUM
