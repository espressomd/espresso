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
	register float minp[3], maxp[3];
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
	register int i, j, k, l, depth, localmaxdepth, skip, inc;
	register float r;
	register float pos[3];
	register float p[3];
	register int ch, n, cell, locked, patch;
	register float radius;
	register float root[3];

	// cache root data
	radius = radiusd;
	for (l = 0; l < 3; l++) root[l] = xd[3 * nnodesd + l];

	localmaxdepth = 1;
	skip = 1;
	// increment to move among the bodies assigned to the given thread
	// hence, one should step over all other threads in GPU with
	// a quantity of blockDim.x * gridDim.x
	inc = blockDim.x * gridDim.x;
	// just a regular 1D GPU index
	i = threadIdx.x + blockIdx.x * blockDim.x;

	// iterate over all bodies assigned to thread
	while (i < nbodiesd) {
		if (skip != 0) {
			// new body, so start traversing at root
			skip = 0;
			for (l = 0; l < 3; l++) p[l] = xd[3 * i + l];
			n = nnodesd;
			depth = 1;
			r = radius;
			j = 0;
			// determine which child to follow
			for (l = 0; l < 3; l++) if (root[l] < p[l]) j += pow(2, l);
		}

		// follow path to leaf cell
		ch = childd[n * 8 + j];
		while (ch >= nbodiesd) {
			n = ch;
			depth++;
			r *= 0.5f;
			j = 0;
			// determine which child to follow
			for (l = 0; l < 3; l++) if (xd[3 * n + l] < p[l]) j += pow(2, l);
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

						for (l = 0; l < 3; l++) pos[l] = ((j >> l) & 1) * r;
						r *= 0.5f;

						massd[cell] = -1.0f;
						startd[cell] = -1;
						for (l = 0; l < 3; l++) pos[l] = xd[3 * cell + l] = xd[3 * n + l] - r + pos[l];
						for (k = 0; k < 8; k++)
							childd[cell* 8 + k] = -1;

						if (patch != cell) {
							childd[n * 8 + j] = cell;
						}

						j = 0;
						for (l = 0; l < 3; l++) if (pos[l] < xd[3 * ch + l]) j += pow(2, l);
						childd[cell * 8 + j] = ch;

						n = cell;
						j = 0;
						for (l = 0; l < 3; l++) if (pos[l] < p[l]) j += pow(2, l);

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
}

/******************************************************************************/
/*** compute center of mass ***************************************************/
/******************************************************************************/

__global__
__launch_bounds__(THREADS3, FACTOR3)
void summarizationKernel()
{
	register int i, j, k, l, ch, inc, missing, cnt, bottom;
	// the node "mass" and its count respectively:
	register float m, cm;
	// position of equivalent total dipole and its magnitude:
	// (like a mass and the center of mass)
	register float p[3], u[3];
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
			for (l = 0; l < 3; l++) {
			  p[l] = 0.0f;
			  u[l] = 0.0f;
			}
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
				for (l = 0; l < 3; l++) {
				  p[l] += xd[3 * ch + l] * m;
				  u[l] += uxd[3 * ch + l];
				}
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
					for (l = 0; l < 3; l++) {
					    p[l] += xd[3 * ch + l] * m;
					    u[l] += uxd[3 * ch + l];
					  }
				}
				// repeat until we are done or child is not ready
			} while ((m >= 0.0f) && (missing != 0));
		}

		if (missing == 0) {
			// all children are ready, so store computed information
				countd[k] = cnt;
				m = 1.0f / cm;
				for (l = 0; l < 3; l++) {
				  xd[3 * k + l] = p[l] * m;
				  uxd[3 * k + l] = u[l];
				}
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
/*** compute dipole-dipole force and torque ***********************************/
/******************************************************************************/

__global__
__launch_bounds__(THREADS5, FACTOR5)
void forceCalculationKernel(dds_float pf,
	     float *force, float* torque, dds_float box_l[3], int periodic[3])
{
	register int i, j, k, l, n, depth, base, sbase, diff, t;
	register float tmp;
	register float dr[3], f[3], h[3], u[3], uc[3], N[3];
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
			for (l = 0; l < 3; l++) {
              u[l] = uxd[3 * i + l];
              h[l] = 0.0f;
              f[l] = 0.0f;
			}

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
					    tmp = 0.0f;    // compute distance squared
					    for (l = 0; l < 3; l++) {
					      dr[l] = -xd[3 * n + l] + xd[3 * i + l];
					      tmp += dr[l] * dr[l];
					    }

						if ((n < nbodiesd) || __all(tmp >= dq[depth])) {	// check if all threads agree that cell is far enough away (or is a body)
							if (n != i) {

							    d1 = sqrtf(tmp/*, 0.5f*/);
							    dd5 = __fdividef(1.0f, tmp * tmp * d1);
							    b = 0.0f;
							    b2 = 0.0f;
							    umd5 = 0.0f;
							    for (l = 0; l < 3; l++) {
							      uc[l] = uxd[3 * n + l];
							      b += uc[l] * dr[l];
							      b2 += u[l] * dr[l];
							      umd5 += u[l] * uc[l];
							    }
							    umd5 *= - 3.0f * dd5;

								bb2d7 = 15.0f * b * b2 * __fdividef(dd5, tmp);

								for (l = 0; l < 3; l++) {
								  h[l] += (b * 3.0f * dr[l] - tmp * uc[l]) * dd5;
								  f[l] += -dr[l] * (umd5 + bb2d7)
								     + 3.0f * (b * u[l] + b2 * uc[l]) * dd5;
								}
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

			N[0] = u[1] * h[2] - u[2] * h[1];
			N[1] = u[2] * h[0] - u[0] * h[2];
			N[2] = u[0] * h[1] - u[1] * h[0];

			for (l = 0; l < 3; l++)
			{
			  if (f[l] != f[l] || h[l] != h[l]) {	//nan
			  	printf("Force Kernel: NAN in particle[%d]\n", i);
			  	printf("x = %f, y = %f, z = %f,\n", uxd[3 * i + 0], uxd[3 * i + 1], uxd[3 * i + 2]);
			  	printf("fx = %f, fy = %f, fz = %f,\n", f[0], f[1], f[2]);
			  	printf("hx = %f, hy = %f, hz = %f,\n", h[0], h[1], h[2]);
			  }
			  atomicAdd(force + 3 * i + l, f[l] * pf);
			  atomicAdd(torque + 3 * i + l, N[l] * pf);
			}
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
	register int i, j, k, l, n, depth, base, sbase, diff, t;
	register float tmp;
	register float dr[3], h[3], u[3], uc[3];
	__shared__ volatile int pos[MAXDEPTH * THREADS5/WARPSIZE], node[MAXDEPTH * THREADS5/WARPSIZE];
	__shared__ float dq[MAXDEPTH * THREADS5/WARPSIZE];
	dds_float sum=0.0;
	extern __shared__ dds_float res[];

	float b, d1, dd5;
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
			for (l = 0; l < 3; l++) {
			  u[l] = uxd[3 * i + l];
			  h[l] = 0.0f;
			}

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
					    tmp = 0.0f;
					    for (l = 0; l < 3; l++) {
					      dr[l] = -xd[3 * n + l] + xd[3 * i + l];
					      tmp += dr[l] * dr[l];
					    }
						if ((n < nbodiesd) || __all(tmp >= dq[depth])) {	// check if all threads agree that cell is far enough away (or is a body)
							if (n != i) {
								d1 = sqrtf(tmp/*, 0.5f*/);
								dd5 = __fdividef(1.0f, tmp * tmp * d1);
								b = 0.0f;
								for (l = 0; l < 3; l++) {
								  uc[l] = uxd[3 * n + l];
								  b += uc[l] * dr[l];
								}

								for (l = 0; l < 3; l++) h[l] += (b * 3.0f * dr[l] - tmp * uc[l]) * dd5;
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

			for (l = 0; l < 3; l++)
			{
			  sum += -u[l] * h[l];
			  if (h[l] != h[l]) {   //nan
			    printf("Energy Kernel: NAN in particle[%d]\n", i);
			    printf("x = %f, y = %f, z = %f,\n", uxd[3 * i + 0], uxd[3 * i + 1], uxd[3 * i + 2]);
			    printf("hx = %f, hy = %f, hz = %f,\n", h[0], h[1], h[2]);
			  }
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
