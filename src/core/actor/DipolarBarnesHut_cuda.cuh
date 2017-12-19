/*
  Copyright (C) 2016,2017 The ESPResSo project
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

#ifdef DIPOLAR_BARNES_HUT

#include <cstdio>
#include <iostream>
#include <time.h>

typedef float dds_float ;

#define SHARED_ARRAY_BH 512

struct BHBox{
    // min positions' coordinates of the BH box.
	float *minp;
	// max positions' coordinates of the BH box.
	float *maxp;
};

struct BHArrays{
    // Error report.
	int *err;
	// Indices of particles sorted according to the tree linear representation.
	int *sort;
	// The tree linear representation.
	int *child;
	// Supplementary array: a tree nodes (division octant cells/particles inside) counting.
	int *count;
	// Start indices for the per-cell sorting.
	int *start;
};

// thread count for different kernels (see kernel calls from below functions).
#define THREADS1 512	/* must be a power of 2 */
#define THREADS2 1024
#define THREADS3 1024
#define THREADS4 256
#define THREADS5 128	//256

// block count = factor * #SMs
// for different kernels (see kernel calls from below functions).
#define FACTOR1 3
#define FACTOR2 1
#define FACTOR3 1	/* must all be resident at the same time */
#define FACTOR4 1	/* must all be resident at the same time */
#define FACTOR5 5

// Warp size.
#define WARPSIZE 32
// Max possible depth of the Barnes-Hut tree branching.
#define MAXDEPTH 32

// Function to set the BH method parameters.
void setBHPrecision(float epssq, float itolsq);

// Populating of array pointers allocated in GPU device from .cu part of the Espresso interface
void fillConstantPointers(float* r, float* dip,
		int nbodies, int nnodes, BHArrays arrl, BHBox boxl, float* mass);

// Required BH CUDA init.
void initBHgpu(int blocks);

// Building Barnes-Hut spatial min/max position box
void buildBoxBH(int blocks);

// Building Barnes-Hut tree in a linear childd array representation
// of octant cells and particles inside.
void buildTreeBH(int blocks);

// Calculate octant cells masses and cell index counts.
// Determine cells centers of mass and total dipole moments
// on all possible levels of the BH tree.
void summarizeBH(int blocks);

// Sort particle indexes according to the BH tree representation.
// Crucial for the per-warp perfomance tuning of forceCalculationKernel and energyCalculationKernel.
void sortBH(int blocks);

// Force calculation.
void forceBH(int blocks, dds_float k, float* f, float* torque, dds_float box_l[3],int periodic[3]);

// Energy calculation.
void energyBH(int blocks, dds_float k, dds_float box_l[3],int periodic[3],float* E);

#endif /* DIPOLARBARNESHUT_CUH_ */
#endif // BARNES_HUT
