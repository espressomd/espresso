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

#include "config.hpp"
#include "thrust/reduce.h"
#include "thrust/device_ptr.h"
#include "cuda.h"
#include <curand.h>
#include <curand_kernel.h>
#include "../cuda_utils.hpp"

typedef float dds_float ;

#ifdef BARNES_HUT

#include "DipolarBarnesHut_cuda.cuh"

#define SQ(x) ((x)*(x))
#define IND (blockDim.x * blockIdx.x + threadIdx.x)

__device__ float sgn(float x){
	return (x > 0.0f) ? 1.0f : ((x < 0.0f) ? -1.0f : 0.0f);
}

using namespace std;

// each node corresponds to a split of the cubic box in 3D space to equal cubic boxes
// hence, 8 nodes per particle is a theoretical octree limit:
// a number of nodes is "nnodesd" and a number of particles "nbodiesd" respectively
__constant__ int nnodesd, nbodiesd;
__constant__ volatile float epssqd, itolsqd;
// blkcntd is a factual blocks' count
// bottomd is a bottom node
// maxdepthd is a largest length of the octree "branch" till the "leaf"
__device__ volatile int bottomd, maxdepthd, blkcntd;
// hald edge of the BH box
__device__ volatile float radiusd;
__device__ __constant__ volatile float* xd;
__constant__ volatile float* uxd;
__constant__ volatile float* massd; // not a real mass. Just a node weight coefficient.
__constant__ volatile float *mind;
__constant__ volatile float *maxd;
__constant__ volatile int *errd;
__constant__ volatile int *sortd;
__constant__ volatile int *childd;
__constant__ volatile int *countd;
__constant__ volatile int *startd;
__constant__ volatile float *box_ld;

__device__ void dds_sumReduction_BH(dds_float *input, dds_float *sum)
{
	int tid = threadIdx.x;
	for (int i = blockDim.x; i > 1; i /= 2)
	{
		__syncthreads();
		if (tid < i/2)
			input[tid] += input[i/2+tid];
		if ((i%2==1) && (tid ==0))
		  input[tid]+=input[i-1];
	}
	__syncthreads();
	if (threadIdx.x==0)
	{
 	  sum[0] = input[0];
    }

}

/******************************************************************************/
/*** initialize memory ********************************************************/
/******************************************************************************/

__global__ void initializationKernel()
{
	register int ind;
	ind = IND;
	if (ind == 0) {
		*errd = 0;
		maxdepthd = 1;
		blkcntd = 0;
	}
}

/******************************************************************************/
/*** compute center and radius ************************************************/
/******************************************************************************/

__global__
__launch_bounds__(THREADS1, FACTOR1)
void boundingBoxKernel()
{
	register int i, j, k, l, inc;
	register float val;
	float minp[3], maxp[3];
	__shared__ volatile float smin[3*THREADS1], smax[3*THREADS1];
	for (l = 0; l < 3; l++) {
	  minp[l] = maxp[l] = xd[l];
	}

	// scan all bodies
	i = threadIdx.x;
	inc = THREADS1 * gridDim.x;
	for (j = i + blockIdx.x * THREADS1; j < nbodiesd; j += inc)
	  for (l = 0; l < 3; l++) {
	    val = xd[3 * j + l];
	    minp[l] = min(minp[l], val);
	    maxp[l] = max(maxp[l], val);
	  }

	// reduction in shared memory:
	for (l = 0; l < 3; l++) {
	  smin[3 * i + l] = minp[l];
	  smax[3 * i + l] = maxp[l];
	}

	for (j = THREADS1 / 2; j > 0; j /= 2) {
		__syncthreads();
		if (i < j) {
			k = i + j;
			for (l = 0; l < 3; l++) {
			  smin[3 * i + l] = minp[l] = min(minp[l], smin[3 * k + l]);
			  smax[3 * i + l] = maxp[l] = max(maxp[l], smax[3 * k + l]);
			}
		}
	}

	// write a given block result to the global memory
	if (i == 0) {
		k = blockIdx.x;
		for (l = 0; l < 3; l++) {
		  mind[3 * k + l] = minp[l];
		  maxd[3 * k + l] = maxp[l];
		}

		inc = gridDim.x - 1;
		if (inc == atomicInc((unsigned int *)&blkcntd, inc)) {
		// I'm the last block, so combine all block results over the index j
			for (j = 0; j <= inc; j++)
			  for (l = 0; l < 3; l++) {
			    minp[l] = min(minp[l], mind[3 * j + l]);
			    maxp[l] = max(maxp[l], maxd[3 * j + l]);
			  }

			// compute 'radius'
			val = max(maxp[0] - minp[0], maxp[1] - minp[1]);
			radiusd = max(val, maxp[2] - minp[2]) * 0.5f;

			k = nnodesd;
			// create the root node of the Barnes-Hut octree
			// bottom node is defined with max possible index just to start
			// it will be updated within further tree building in
			// corresponding kernel
			bottomd = k;
			// mass of the root node has no physical sense cause
			// it does not correspond to any real particle
			massd[k] = -1.0f;
			// sorting init
			startd[k] = 0;
			// position of the root node should be in the center of the just defined BH box
			for (l = 0; l < 3; l++) xd[3 * k + l] = (minp[l] + maxp[l]) * 0.5f;
			k *= 8;
			// init further tree building:
			for (i = 0; i < 8; i++) childd[k + i] = -1;
		}
	}
}


/******************************************************************************/
/*** build tree ***************************************************************/
/******************************************************************************/

__global__
__launch_bounds__(THREADS2, FACTOR2)
void treeBuildingKernel()
{
	register int i, j, k, depth, localmaxdepth, skip, inc;
	register float x, y, z, r;
	register float px, py, pz;
	register int ch, n, cell, locked, patch;
	register float radius, rootx, rooty, rootz;

	// cache root data
	radius = radiusd;
	rootx = xd[3 * nnodesd];
	rooty = xd[3 * nnodesd + 1];
	rootz = xd[3 * nnodesd + 2];

	localmaxdepth = 1;
	skip = 1;
	inc = blockDim.x * gridDim.x;
	i = threadIdx.x + blockIdx.x * blockDim.x;

	// iterate over all bodies assigned to thread
	while (i < nbodiesd) {
		if (skip != 0) {
			// new body, so start traversing at root
			skip = 0;
			px = xd[3 * i];
			py = xd[3 * i + 1];
			pz = xd[3 * i + 2];
			n = nnodesd;
			depth = 1;
			r = radius;
			j = 0;
			// determine which child to follow
			if (rootx < px) j = 1;
			if (rooty < py) j += 2;
			if (rootz < pz) j += 4;
		}

		// follow path to leaf cell
		ch = childd[n * 8 + j];
		while (ch >= nbodiesd) {
			n = ch;
			depth++;
			r *= 0.5f;
			j = 0;
			// determine which child to follow
			if (xd[3 * n] < px) j = 1;
			if (xd[3 * n + 1] < py) j += 2;
			if (xd[3 * n + 2] < pz) j += 4;
			ch = childd[n * 8 + j];
		}

		if (ch != -2) {	// skip if child pointer is locked and try again later
			locked = n  * 8 + j;
			if (ch == atomicCAS((int *)&childd[locked], ch, -2)) {	// try to lock
				if (ch == -2) {
					printf("Error: ch = -2\n");
					break;
				}
				if (ch == -1) {
					// if null, just insert the new body
					childd[locked] = i;
				} else {	// there already is a body in this position
					patch = -1;
					// create new cell(s) and insert the old and new body
					do {
						depth++;

						cell = atomicSub((int *)&bottomd, 1) - 1;
						if (cell <= nbodiesd) {
							*errd = 1;
							bottomd = nnodesd;
						}
						patch = max(patch, cell);

						x = (j & 1) * r;
						y = ((j >> 1) & 1) * r;
						z = ((j >> 2) & 1) * r;
						r *= 0.5f;

						massd[cell] = -1.0f;
						startd[cell] = -1;
						x = xd[3 * cell] = xd[3 * n] - r + x;
						y = xd[3 * cell + 1] = xd[3 * n + 1] - r + y;
						z = xd[3 * cell + 2] = xd[3 * n + 2] - r + z;
						for (k = 0; k < 8; k++)
							childd[cell* 8 + k] = -1;

						if (patch != cell) {
							childd[n * 8 + j] = cell;
						}

						j = 0;
						if (x < xd[3 * ch]) j = 1;
						if (y < xd[3 * ch + 1]) j += 2;
						if (z < xd[3 * ch + 2]) j += 4;
						childd[cell * 8 + j] = ch;

						n = cell;
						j = 0;
						if (x < px) j = 1;
						if (y < py) j += 2;
						if (z < pz) j += 4;

						ch = childd[n * 8 + j];
						// repeat until the two bodies are different children
					} while (ch >= 0);

					childd[n * 8 + j] = i;
					__threadfence();	// push out subtree
					childd[locked] = patch;
				}

				localmaxdepth = max(depth, localmaxdepth);
				i += inc;	// move on to next body
				skip = 1;
			}
		}
		__syncthreads();	// throttle
	}
	// record maximum tree depth
	atomicMax((int *)&maxdepthd, localmaxdepth);
#ifdef DEBUG
//	if ( maxdepthd == localmaxdepth ) printf ("Max tree depth: %d\n", maxdepthd);
#endif
}

/******************************************************************************/
/*** compute center of mass ***************************************************/
/******************************************************************************/

__global__
__launch_bounds__(THREADS3, FACTOR3)
void summarizationKernel()
{
	register int i, j, k, ch, inc, missing, cnt, bottom;
	register float m, cm, px, py, pz, ux, uy, uz;
	__shared__ volatile int child[THREADS3 * 8];

	bottom = bottomd;
	inc = blockDim.x * gridDim.x;
	k = (bottom & (-WARPSIZE)) + threadIdx.x + blockIdx.x * blockDim.x;	// align to warp size
	if (k < bottom) k += inc;

	missing = 0;
	// iterate over all cells assigned to thread
	while (k <= nnodesd) {
		//iteration++;
		if (missing == 0) {
			// new cell, so initialize
			cm = 0.0f;
			px = 0.0f;
			py = 0.0f;
			pz = 0.0f;
			ux = 0.0f;
			uy = 0.0f;
			uz = 0.0f;
			cnt = 0;
			j = 0;
			for (i = 0; i < 8; i++) {
				ch = childd[k * 8 + i];
				if (ch >= 0) {
			if (i != j) {
				// move children to front (needed later for speed)
				childd[k* 8 + i] = -1;
				childd[k* 8 + j] = ch;
			}
			child[missing * THREADS3 + threadIdx.x] = ch;	// cache missing children
			m = massd[ch];
			missing++;
			if (m >= 0.0f) {
				// child is ready
				missing--;
				if (ch >= nbodiesd) {	// count bodies (needed later)
					cnt += countd[ch] - 1;
				}
				// add child's contribution
				cm += m;
				px += xd[3 * ch] * m;
				py += xd[3 * ch + 1] * m;
				pz += xd[3 * ch + 2] * m;
				ux += uxd[3 * ch];
				uy += uxd[3 * ch + 1];
				uz += uxd[3 * ch + 2];
			}
			j++;
				}
			}
			cnt += j;
		}

		if (missing != 0) {
			do {
				// poll missing child
				ch = child[(missing - 1) * THREADS3 + threadIdx.x];
				m = massd[ch];
				if (m >= 0.0f) {
					// child is now ready
					missing--;
					if (ch >= nbodiesd) {
						// count bodies (needed later)
						cnt += countd[ch] - 1;
					}
					// add child's contribution
					cm += m;
					px += xd[3 * ch] * m;
					py += xd[3 * ch + 1] * m;
					pz += xd[3 * ch + 2] * m;
					ux += uxd[3 * ch];
					uy += uxd[3 * ch + 1];
					uz += uxd[3 * ch + 2];
				}
				// repeat until we are done or child is not ready
			} while ((m >= 0.0f) && (missing != 0));
		}

		if (missing == 0) {
			// all children are ready, so store computed information
				countd[k] = cnt;
				m = 1.0f / cm;
				xd[3 * k] = px * m;
				xd[3 * k + 1] = py * m;
				xd[3 * k + 2] = pz * m;

				uxd[3 * k] = ux;
				uxd[3 * k + 1] = uy;
				uxd[3 * k + 2] = uz;
			__threadfence();	// make sure data are visible before setting mass
			massd[k] = cm;
			k += inc;	// move on to next cell
		}
	}	//while
}

/******************************************************************************/
/*** sort bodies **************************************************************/
/******************************************************************************/

__global__
__launch_bounds__(THREADS4, FACTOR4)
void sortKernel()
{
	register int i, k, ch, dec, start, bottom;

	bottom = bottomd;
	dec = blockDim.x * gridDim.x;
	k = nnodesd + 1 - dec + threadIdx.x + blockIdx.x * blockDim.x;

	// iterate over all cells assigned to thread
	while (k >= bottom) {
		start = startd[k];
		if (start >= 0) {
			for (i = 0; i < 8; i++) {
				ch = childd[k* 8 + i];
				if (ch >= nbodiesd) {
			// child is a cell
			startd[ch] = start;	// set start ID of child
			start += countd[ch];	// add #bodies in subtree
				} else if (ch >= 0) {
			// child is a body
			sortd[start] = ch;	// record body in 'sorted' array
			start++;
				}
			}
			k -= dec;	// move on to next cell
		}
		__syncthreads();	// throttle
	}
}


/******************************************************************************/
/*** compute force ************************************************************/
/******************************************************************************/

__global__
__launch_bounds__(THREADS5, FACTOR5)
void forceCalculationKernel(dds_float pf,
	     float *f, float* torque, dds_float box_l[3], int periodic[3])
{
	register int i, j, k, n, depth, base, sbase, diff, t;
	register float px, py, pz, dx, dy, dz, tmp, fx, fy, fz, hx, hy, hz, ux, uy, uz;
	register float ucx, ucy, ucz;
	__shared__ volatile int pos[MAXDEPTH * THREADS5/WARPSIZE], node[MAXDEPTH * THREADS5/WARPSIZE];
	__shared__ float dq[MAXDEPTH * THREADS5/WARPSIZE];
	float b, b2, d1, dd5;
	float bb2d7, umd5;
	float pow3s2d2, flj;

	if (0 == threadIdx.x) {
		tmp = radiusd;
		// precompute values that depend only on tree level
		dq[0] = tmp * tmp * itolsqd;
		for (i = 1; i < maxdepthd; i++) {
			dq[i] = dq[i - 1] * 0.25f;
			dq[i - 1] += epssqd;
		}
		dq[i - 1] += epssqd;

		if (maxdepthd > MAXDEPTH) {
			*errd = maxdepthd;
		}
	}
	__syncthreads();

	if (maxdepthd <= MAXDEPTH) {
		// figure out first thread in each warp (lane 0)
		base = threadIdx.x / WARPSIZE;
		sbase = base * WARPSIZE;
		j = base * MAXDEPTH;

		diff = threadIdx.x - sbase;
		// make multiple copies to avoid index calculations later
		if (diff < MAXDEPTH) {
			dq[diff+j] = dq[diff];
		}
		__syncthreads();

		// iterate over all bodies assigned to thread
		for (k = threadIdx.x + blockIdx.x * blockDim.x; k < nbodiesd; k += blockDim.x * gridDim.x) {
			i = sortd[k];	// get permuted/sorted index
			// cache position info
			px = xd[3 * i];
			py = xd[3 * i + 1];
			pz = xd[3 * i + 2];

			ux = uxd[3 * i];
			uy = uxd[3 * i + 1];
			uz = uxd[3 * i + 2];

			hx = 0.0f;
			hy = 0.0f;
			hz = 0.0f;

			fx = 0.0f;
			fy = 0.0f;
			fz = 0.0f;

			// initialize iteration stack, i.e., push root node onto stack
			depth = j;
			if (sbase == threadIdx.x) {
				node[j] = nnodesd;
				pos[j] = 0;
			}

			while (depth >= j) {
				// stack is not empty
				while ((t = pos[depth]) < 8) {
					// node on top of stack has more children to process
					n = childd[node[depth] * 8 + t];	// load child pointer
					if (sbase == threadIdx.x) {
						// I'm the first thread in the warp
						pos[depth] = t + 1;
					}
					if (n >= 0) {
						dx = -(xd[3 * n] - px);
						dy = -(xd[3 * n + 1] - py);
						dz = -(xd[3 * n + 2] - pz);
						tmp = dx * dx + (dy * dy + dz * dz);	// compute distance squared
						if ((n < nbodiesd) || __all(tmp >= dq[depth])) {	// check if all threads agree that cell is far enough away (or is a body)
							if (n != i) {

								ucx = uxd[3 * n];
								ucy = uxd[3 * n + 1];
								ucz = uxd[3 * n + 2];

								b = ucx * dx + ucy * dy + ucz * dz;
								b2 = ux * dx + uy * dy + uz * dz;

								d1 = sqrtf(tmp/*, 0.5f*/);
								dd5 = __fdividef(1.0f, tmp * tmp * d1);
								bb2d7 = 15.0f * b * b2 * __fdividef(dd5, tmp);
								umd5 = - 3.0f * (ux*ucx + uy*ucy + uz*ucz) * dd5;

								hx += (b * 3.0f * dx - tmp * ucx) * dd5;
								hy += (b * 3.0f * dy - tmp * ucy) * dd5;
								hz += (b * 3.0f * dz - tmp * ucz) * dd5;


								fx += -dx * (umd5 + bb2d7)
								+ 3.0f * (b * ux + b2 * ucx) * dd5;

								fy += -dy * (umd5  +  bb2d7)
								+ 3.0f * (b * uy + b2 * ucy) * dd5;

								fz += -dz * (umd5  +  bb2d7)
								+ 3.0f * (b * uz + b2 * ucz) * dd5;
							}
						} else {
							// push cell onto stack
							depth++;
							if (sbase == threadIdx.x) {
								node[depth] = n;
								pos[depth] = 0;
							}
						}
					} else {
						depth = max(j, depth - 1);	// early out because all remaining children are also zero
					}
				}
				depth--;	// done with this level
			}

#define	NX (uy * hz - uz * hy)
#define	NY (uz * hx - ux * hz)
#define	NZ (ux * hy - uy * hx)

#define MAX_DELTA 0.1f

			if (fx != fx || fy != fy || fz != fz ||
				hx != hx || hy != hy || hz != hz) {	//nan
				printf("Force Kernel: NAN in particle[%d]\n", i);
				printf("x = %f, y = %f, z = %f,\n", px, py, pz);
				printf("fx = %f, fy = %f, fz = %f,\n", fx, fy, fz);
				printf("hx = %f, hy = %f, hz = %f,\n", hx, hy, hz);

				fx = 1E10;
				fy = 1E10;
				fz = 1E10;
				hx = 0.0f;
				hy = 0.0f;
				hz = 0.0f;
			}

			atomicAdd(f+3*i+0, fx * pf);
			atomicAdd(torque+3*i+0, NX * pf);
			atomicAdd(f+3*i+1, fy * pf);
			atomicAdd(torque+3*i+1, NY * pf);
			atomicAdd(f+3*i+2, fz * pf);
			atomicAdd(torque+3*i+2, NZ * pf);
		}
	}
}

/******************************************************************************/
/*** compute energy ************************************************************/
/******************************************************************************/

__global__
__launch_bounds__(THREADS5, FACTOR5)
void energyCalculationKernel(dds_float pf,
	     dds_float box_l[3], int periodic[3],dds_float* energySum)
{
	register int i, j, k, n, depth, base, sbase, diff, t;
	register float px, py, pz, dx, dy, dz, tmp, fx, fy, fz, hx, hy, hz, ux, uy, uz;
	register float ucx, ucy, ucz;
	__shared__ volatile int pos[MAXDEPTH * THREADS5/WARPSIZE], node[MAXDEPTH * THREADS5/WARPSIZE];
	__shared__ float dq[MAXDEPTH * THREADS5/WARPSIZE];
	dds_float sum=0.0;
	extern __shared__ dds_float res[];

	float b, b2, d1, dd5;
	float bb2d7, umd5;
	float pow3s2d2, flj;

	if (0 == threadIdx.x) {
		tmp = radiusd;
		// precompute values that depend only on tree level
		dq[0] = tmp * tmp * itolsqd;
		for (i = 1; i < maxdepthd; i++) {
			dq[i] = dq[i - 1] * 0.25f;
			dq[i - 1] += epssqd;
		}
		dq[i - 1] += epssqd;

		if (maxdepthd > MAXDEPTH) {
			*errd = maxdepthd;
		}
	}
	__syncthreads();

	if (maxdepthd <= MAXDEPTH) {
		// figure out first thread in each warp (lane 0)
		base = threadIdx.x / WARPSIZE;
		sbase = base * WARPSIZE;
		j = base * MAXDEPTH;

		diff = threadIdx.x - sbase;
		// make multiple copies to avoid index calculations later
		if (diff < MAXDEPTH) {
			dq[diff+j] = dq[diff];
		}
		__syncthreads();

		// iterate over all bodies assigned to thread
		for (k = threadIdx.x + blockIdx.x * blockDim.x; k < nbodiesd; k += blockDim.x * gridDim.x) {
			i = sortd[k];	// get permuted/sorted index
			// cache position info
			px = xd[3 * i];
			py = xd[3 * i + 1];
			pz = xd[3 * i + 2];

			ux = uxd[3 * i];
			uy = uxd[3 * i + 1];
			uz = uxd[3 * i + 2];

			hx = 0.0f;
			hy = 0.0f;
			hz = 0.0f;

			fx = 0.0f;
			fy = 0.0f;
			fz = 0.0f;

			// initialize iteration stack, i.e., push root node onto stack
			depth = j;
			if (sbase == threadIdx.x) {
				node[j] = nnodesd;
				pos[j] = 0;
			}

			while (depth >= j) {
				// stack is not empty
				while ((t = pos[depth]) < 8) {
					// node on top of stack has more children to process
					n = childd[node[depth] * 8 + t];	// load child pointer
					if (sbase == threadIdx.x) {
						// I'm the first thread in the warp
						pos[depth] = t + 1;
					}
					if (n >= 0) {
						dx = -(xd[3 * n] - px);
						dy = -(xd[3 * n + 1] - py);
						dz = -(xd[3 * n + 2] - pz);
						tmp = dx * dx + (dy * dy + dz * dz);	// compute distance squared
						if ((n < nbodiesd) || __all(tmp >= dq[depth])) {	// check if all threads agree that cell is far enough away (or is a body)
							if (n != i) {

								ucx = uxd[3 * n];
								ucy = uxd[3 * n + 1];
								ucz = uxd[3 * n + 2];

								b = ucx * dx + ucy * dy + ucz * dz;
								b2 = ux * dx + uy * dy + uz * dz;

								d1 = sqrtf(tmp/*, 0.5f*/);
								dd5 = __fdividef(1.0f, tmp * tmp * d1);
								bb2d7 = 15.0f * b * b2 * __fdividef(dd5, tmp);
								umd5 = - 3.0f * (ux*ucx + uy*ucy + uz*ucz) * dd5;

								hx += (b * 3.0f * dx - tmp * ucx) * dd5;
								hy += (b * 3.0f * dy - tmp * ucy) * dd5;
								hz += (b * 3.0f * dz - tmp * ucz) * dd5;
							}
						} else {
							// push cell onto stack
							depth++;
							if (sbase == threadIdx.x) {
								node[depth] = n;
								pos[depth] = 0;
							}
						}
					} else {
						depth = max(j, depth - 1);	// early out because all remaining children are also zero
					}
				}
				depth--;	// done with this level
			}

			sum += - (ux * hx + uy * hy + uz * hz);

#define	NX (uy * hz - uz * hy)
#define	NY (uz * hx - ux * hz)
#define	NZ (ux * hy - uy * hx)

#define MAX_DELTA 0.1f


			if (fx != fx || fy != fy || fz != fz ||
				hx != hx || hy != hy || hz != hz) {	//nan
				printf("Force Kernel: NAN in particle[%d]\n", i);
				printf("x = %f, y = %f, z = %f,\n", px, py, pz);
				printf("fx = %f, fy = %f, fz = %f,\n", fx, fy, fz);
				printf("hx = %f, hy = %f, hz = %f,\n", hx, hy, hz);

				fx = 1E10;
				fy = 1E10;
				fz = 1E10;
				hx = 0.0f;
				hy = 0.0f;
				hz = 0.0f;
			}
		}
	}

	sum *= pf;
    sum /= 2.0f; // the self-consistent field energy
    // Save per thread result into block shared mem
    res[threadIdx.x] =sum;
    // Sum results within a block
    __syncthreads(); // Wait til all threads in block are done
    dds_sumReduction_BH(res,&(energySum[blockIdx.x]));
}

static void CudaTest(char *msg) {
	cudaError_t e;

	cudaThreadSynchronize();
	if (cudaSuccess != (e = cudaGetLastError())) {
		fprintf(stderr, "%s: %d\n", msg, e);
		fprintf(stderr, "%s\n", cudaGetErrorString(e));
		exit(-1);
	}
}

void initBHgpu(int blocks) {
	initializationKernel<<<blocks * FACTOR5, THREADS5>>>();
	CudaTest("init kernel launch failed");

	cudaFuncSetCacheConfig(boundingBoxKernel, cudaFuncCachePreferShared);
	cudaFuncSetCacheConfig(boundingBoxKernel, cudaFuncCachePreferShared);
	cudaFuncSetCacheConfig(treeBuildingKernel, cudaFuncCachePreferL1);
	cudaFuncSetCacheConfig(summarizationKernel, cudaFuncCachePreferShared);
	cudaFuncSetCacheConfig(sortKernel, cudaFuncCachePreferL1);
	cudaFuncSetCacheConfig(forceCalculationKernel, cudaFuncCachePreferL1);

	cudaGetLastError();	// reset error value
}

void buildBoxBH(int blocks) {
	cudaThreadSynchronize();
	boundingBoxKernel<<<blocks * FACTOR1, THREADS1>>>();
	CudaTest("kernel 1 launch failed");
	cudaThreadSynchronize();
}

void buildTreeBH(int blocks) {
	treeBuildingKernel<<<blocks * FACTOR2, THREADS2>>>();
	CudaTest("kernel 2 launch failed");
	cudaThreadSynchronize();
}

void summarizeBH(int blocks) {
	summarizationKernel<<<blocks * FACTOR3, THREADS3>>>();
	CudaTest("kernel 3 launch failed");
	cudaThreadSynchronize();
}

void sortBH(int blocks) {
	sortKernel<<<blocks * FACTOR4, THREADS4>>>();
	CudaTest("kernel 4 launch failed");
	cudaThreadSynchronize();
}

void forceBH(int blocks, dds_float k, float* f, float* torque, dds_float box_l[3],int periodic[3]) {
	dds_float* box_l_gpu;
	int* periodic_gpu;
	cuda_safe_mem(cudaMalloc((void**)&box_l_gpu,3*sizeof(dds_float)));
	cuda_safe_mem(cudaMalloc((void**)&periodic_gpu,3*sizeof(int)));
	cuda_safe_mem(cudaMemcpy(box_l_gpu,box_l,3*sizeof(dds_float),cudaMemcpyHostToDevice));
	cuda_safe_mem(cudaMemcpy(periodic_gpu,periodic,3*sizeof(int),cudaMemcpyHostToDevice));

	forceCalculationKernel<<<blocks * FACTOR5, THREADS5>>>(k, f, torque, box_l_gpu, periodic_gpu);
	CudaTest("kernel 5 launch failed");
	cudaThreadSynchronize();

	cudaFree(box_l_gpu);
	cudaFree(periodic_gpu);
}

void energyBH(int blocks, dds_float k, dds_float box_l[3],int periodic[3],float* E) {
	dds_float* box_l_gpu;
	int* periodic_gpu;
	dim3 grid(1,1,1);
	dim3 block(1,1,1);

	grid.x = blocks * FACTOR5;
	block.x = THREADS5;

	cuda_safe_mem(cudaMalloc((void**)&box_l_gpu,3*sizeof(dds_float)));
	cuda_safe_mem(cudaMalloc((void**)&periodic_gpu,3*sizeof(int)));
	cuda_safe_mem(cudaMemcpy(box_l_gpu,box_l,3*sizeof(float),cudaMemcpyHostToDevice));
	cuda_safe_mem(cudaMemcpy(periodic_gpu,periodic,3*sizeof(int),cudaMemcpyHostToDevice));

	dds_float *energySum;
	cuda_safe_mem(cudaMalloc(&energySum,(int)(sizeof(dds_float)*grid.x)));

	energyCalculationKernel<<<grid, block, THREADS5*sizeof(dds_float)>>>(k, box_l_gpu, periodic_gpu,energySum);
	CudaTest("kernel 6 launch failed");
	cudaThreadSynchronize();

	// Sum the results of all blocks
	// One energy part per block in the prev kernel
	thrust::device_ptr<dds_float> t(energySum);
	float x=thrust::reduce(t,t+grid.x);
	cuda_safe_mem(cudaMemcpy(E,&x,sizeof(float),cudaMemcpyHostToDevice));

	cuda_safe_mem(cudaFree(energySum));
	cuda_safe_mem(cudaFree(box_l_gpu));
	cuda_safe_mem(cudaFree(periodic_gpu));
}

void setBHPrecision(float epssq, float itolsq) {
    float epssq_loc, itolsq_loc;
    epssq_loc = epssq;
    itolsq_loc = itolsq;
    cuda_safe_mem(cudaMemcpyToSymbol(epssqd, &epssq_loc, sizeof(float), 0, cudaMemcpyHostToDevice));
    cuda_safe_mem(cudaMemcpyToSymbol(itolsqd, &itolsq_loc, sizeof(float), 0, cudaMemcpyHostToDevice));
}

// Populating of array pointers allocated in GPU device from .cu part of the Espresso interface
void fillConstantPointers(float* r, float* dip, int nbodies, int nnodes, BHArrays arrl, BHBox boxl, float* mass) {
	cuda_safe_mem(cudaMemcpyToSymbol(nnodesd, &nnodes, sizeof(int), 0, cudaMemcpyHostToDevice));
	cuda_safe_mem(cudaMemcpyToSymbol(nbodiesd, &nbodies, sizeof(int), 0, cudaMemcpyHostToDevice));
	cuda_safe_mem(cudaMemcpyToSymbol(errd, &(arrl.err), sizeof(void*), 0, cudaMemcpyHostToDevice));
	cuda_safe_mem(cudaMemcpyToSymbol(sortd, &(arrl.sort), sizeof(void*), 0, cudaMemcpyHostToDevice));
	cuda_safe_mem(cudaMemcpyToSymbol(childd, &(arrl.child), sizeof(void*), 0, cudaMemcpyHostToDevice));
	cuda_safe_mem(cudaMemcpyToSymbol(countd, &(arrl.count), sizeof(void*), 0, cudaMemcpyHostToDevice));
	cuda_safe_mem(cudaMemcpyToSymbol(startd, &(arrl.start), sizeof(void*), 0, cudaMemcpyHostToDevice));
	cuda_safe_mem(cudaMemcpyToSymbol(xd, &r, sizeof(float*), 0, cudaMemcpyHostToDevice));
	cuda_safe_mem(cudaMemcpyToSymbol(uxd, &dip, sizeof(float*), 0, cudaMemcpyHostToDevice));
	cuda_safe_mem(cudaMemcpyToSymbol(massd, &(mass), sizeof(void*), 0, cudaMemcpyHostToDevice));
	cuda_safe_mem(cudaMemcpyToSymbol(maxd, &(boxl.maxp), sizeof(void*), 0, cudaMemcpyHostToDevice));
	cuda_safe_mem(cudaMemcpyToSymbol(mind, &(boxl.minp), sizeof(void*), 0, cudaMemcpyHostToDevice));
}

#endif // BARNES_HUT
