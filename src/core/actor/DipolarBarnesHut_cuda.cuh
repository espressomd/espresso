/*
 * Copyright (C) 2016-2019 The ESPResSo project
 * Copyright (C) 2012 Alexander (Polyakov) Peletskyi
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef DIPOLARBARNESHUT_CUH_
#define DIPOLARBARNESHUT_CUH_

#include "config.hpp"

#ifdef DIPOLAR_BARNES_HUT

typedef float dds_float;

typedef struct {
  // CUDA blocks
  int blocks;
  // each node corresponds to a split of the cubic box in 3D space to equal
  // cubic boxes hence, 8 octant nodes per particle is a theoretical octree
  // limit: a maximal number of octree nodes is "nnodesd" and a number of
  // particles "nbodiesd" respectively.
  int nbodies;
  int nnodes;
  // particle positions on the device:
  float *r;
  // particle dipole moments on the device:
  float *u;
  // Not a real mass. Just a node weight coefficient.
  float *mass;
  // min positions' coordinates of the BH box.
  float *minp;
  // max positions' coordinates of the BH box.
  float *maxp;
  // Error report.
  int *err;
  // Indices of particles sorted according to the tree linear representation.
  int *sort;
  // The tree linear representation.
  int *child;
  // Supplementary array: a tree nodes (division octant cells/particles inside)
  // counting.
  int *count;
  // Start indices for the per-cell sorting.
  int *start;
  // trace the max loops for a threads' sync
  int *max_lps;
} BHData;

// thread count for different kernels (see kernel calls from below functions).
#define THREADS1 512
#define THREADS2 1024
#define THREADS3 1024
#define THREADS4 1024
#define THREADS5 256

// block count = factor * #SMs
// for different kernels (see kernel calls from below functions).
#define FACTOR1 2
#define FACTOR2 1
#define FACTOR3 1 /* must all be resident at the same time */
#define FACTOR4 1 /* must all be resident at the same time */
#define FACTOR5 4

// Warp size.
#define WARPSIZE 32
// Max possible depth of the Barnes-Hut tree branching.
#define MAXDEPTH 32

// Function to set the BH method parameters.
void setBHPrecision(float *epssq, float *itolsq);

// An allocation of the GPU device memory and an initialization where it is
// needed.
void allocBHmemCopy(int nbodies, BHData *bh_data);

// Populating of array pointers allocated in GPU device before.
// Copy the particle data to the Barnes-Hut related arrays.
void fillConstantPointers(float *r, float *dip, BHData bh_data);

// Required BH CUDA init.
void initBHgpu(int blocks);

// Building Barnes-Hut spatial min/max position box
void buildBoxBH(int blocks);

// Building Barnes-Hut tree in a linear child array representation
// of octant cells and particles inside.
void buildTreeBH(int blocks);

// Calculate octant cells masses and cell index counts.
// Determine cells centers of mass and total dipole moments
// on all possible levels of the BH tree.
void summarizeBH(int blocks);

// Sort particle indexes according to the BH tree representation.
// Crucial for the per-warp performance tuning of forceCalculationKernel and
// energyCalculationKernel.
void sortBH(int blocks);

// Force calculation.
int forceBH(BHData *bh_data, dds_float k, float *f, float *torque);

// Energy calculation.
int energyBH(BHData *bh_data, dds_float k, float *E);

#endif // DIPOLAR_BARNES_HUT
#endif /* DIPOLARBARNESHUT_CUH_ */
