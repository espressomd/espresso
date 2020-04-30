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
/** @file
 *  The method concept is revealed within @cite burtscher11a
 */

#include "cuda_wrapper.hpp"

#include "../cuda_init.hpp"
#include "../cuda_utils.hpp"
#include "config.hpp"
#include <thrust/device_ptr.h>
#include <thrust/reduce.h>

typedef float dds_float;

#ifdef DIPOLAR_BARNES_HUT

#include "DipolarBarnesHut_cuda.cuh"

#define IND (blockDim.x * blockIdx.x + threadIdx.x)

using namespace std;

// Method performance/accuracy parameters
__constant__ float epssqd[1], itolsqd[1];
// blkcntd is a factual blocks' count.
// bottomd is a bottom Barnes-Hut node (the division octant cell) in a linear
// array representation. maxdepthd is a largest length of the octree "branch"
// till the "leaf".
__device__ volatile int bottomd, maxdepthd, blkcntd;
// half edge of the BH box
__device__ volatile float radiusd;
// the struct containing all the device pointers
__device__ __constant__ volatile BHData bhpara[1];

// The "half-convolution" multi-thread reduction.
// The thread with a lower index will operate longer and
// final result (here: the sum) is flowing down towards zero thread.
__device__ void dds_sumReduction_BH(dds_float *input, dds_float *sum) {
  auto tid = static_cast<int>(threadIdx.x);
  for (auto i = static_cast<int>(blockDim.x); i > 1; i /= 2) {
    __syncthreads();
    if (tid < i / 2)
      input[tid] += input[i / 2 + tid];
    if ((i % 2 == 1) && (tid == 0))
      input[tid] += input[i - 1];
  }
  __syncthreads();
  if (threadIdx.x == 0) {
    sum[0] = input[0];
  }
}

/******************************************************************************/
/*** initialize memory ********************************************************/
/******************************************************************************/

__global__ void initializationKernel() {
  int ind;
  ind = IND;
  if (ind == 0) {
    *bhpara->err = 0;
    *bhpara->max_lps = 0;
    maxdepthd = 1;
    blkcntd = 0;
  }
}

/******************************************************************************/
/*** compute center and radius ************************************************/
/******************************************************************************/

__global__ __launch_bounds__(THREADS1, FACTOR1) void boundingBoxKernel() {
  int i, j, k, l, t, inc, n_blocks;
  float val;
  // min/max positions per the thread:
  float minp[3], maxp[3];
  // min/max positions per block:
  __shared__ float smin[3 * THREADS1], smax[3 * THREADS1];
  for (l = 0; l < 3; l++) {
    minp[l] = maxp[l] = bhpara->r[l];
  }

  // Scan all bodies.
  // In order to iterate over all bodies assigned to thread,
  // it is necessary to step over all threads in the GPU:
  // inc = [number of blocks: gridDim.x] * [THREADS1 per block within given
  // kernel]. Hence, this approach could handle an infinite number of bodies
  // (particles)
  i = static_cast<int>(threadIdx.x);
  inc = THREADS1 * gridDim.x;
  // j is an absolute index of the particle.
  // It is shifted over a count of the passed block threads behind: blockIdx.x *
  // THREADS1. NOTE: this loop is extrema search among all particles of the
  // given thread in the present block. However, one is not among all threads of
  // this block.
  for (j = i + static_cast<int>(blockIdx.x) * THREADS1; j < bhpara->nbodies;
       j += inc)
    for (l = 0; l < 3; l++) {
      val = bhpara->r[3 * j + l];
      minp[l] = min(minp[l], val);
      maxp[l] = max(maxp[l], val);
    }

  // For a start point of a reduction in the given block shared memory
  // of the i-th thread extrema:
  for (l = 0; l < 3; l++) {
    smin[3 * i + l] = minp[l];
    smax[3 * i + l] = maxp[l];
  }

  // Now it's time to (min)maximize among all threads of the given block.
  // Each mim/max operation will be applied between
  // the shared memory smin/smax and the given thread.
  // The "half-convolution" multi-thread reduction
  // the thread with a lower index will operate longer and
  // final result (here: the shared memory extrema)
  // is flowing down towards zero thread.
  for (t = THREADS1 / 2; t > 0; t /= 2) {
    __syncthreads();
    if (i < t) {
      k = i + t;
      for (l = 0; l < 3; l++) {
        smin[3 * i + l] = minp[l] = min(minp[l], smin[3 * k + l]);
        smax[3 * i + l] = maxp[l] = max(maxp[l], smax[3 * k + l]);
        // very last minp/maxp assignment will be made by zero thread (i == 0)
      }
    }
  }

  // Thread i == 0 is responsible for a writing
  // of a given block result into the global memory
  // and other per-block operations.
  if (i == 0) {
    // per k-th block
    k = static_cast<int>(blockIdx.x);
    for (l = 0; l < 3; l++) {
      // global memory storage of the per-block extrema
      bhpara->minp[3 * k + l] = minp[l];
      bhpara->maxp[3 * k + l] = maxp[l];
      // note, that we are in zero thread and its variables minp/maxp
      // contain de facto already reduced (see above) shared extrema smin/smax
    }

    n_blocks = static_cast<int>(gridDim.x) - 1;
    // The block increment is performing by its zero thread.
    if (n_blocks == atomicInc((unsigned int *)&blkcntd, n_blocks)) {
      // I'm the (randomly) last block, so combine all other blocks' results
      // over the index j:
      for (j = 0; j <= n_blocks; j++)
        for (l = 0; l < 3; l++) {
          minp[l] = min(minp[l], bhpara->minp[3 * j + l]);
          maxp[l] = max(maxp[l], bhpara->maxp[3 * j + l]);
        }

      // Compute 'radius':
      val = max(maxp[0] - minp[0], maxp[1] - minp[1]);
      radiusd = max(val, maxp[2] - minp[2]) * 0.5f;

      // NOTE: now the extrema are global.
      // Present code fragment will be executed once: in zero thread of the last
      // block.

      k = bhpara->nnodes;
      // Create the root node of the Barnes-Hut octree.
      // Bottom node is defined with max possible index just to start
      // It will be updated within further tree building in
      // corresponding kernel.
      bottomd = k;
      // Weight of the root node init.
      bhpara->mass[k] = -1.0f;
      // Sorting init for the tree root.
      bhpara->start[k] = 0;
      // Position of the root node should be in the center of just defined BH
      // box:
      for (l = 0; l < 3; l++)
        bhpara->r[3 * k + l] = (minp[l] + maxp[l]) * 0.5f;
      // Init further tree building octo- meaning their absence at the
      // beginning:
      for (i = 0; i < 8; i++)
        bhpara->child[8 * k + i] = -1;
    }
  }
}

/******************************************************************************/
/*** build tree ***************************************************************/
/******************************************************************************/

__global__ __launch_bounds__(THREADS2, FACTOR2) void treeBuildingKernel() {
  int i, j, k, l, depth, localmaxdepth, skip, inc, lps;
  float r;
  float pos[3];
  float p[3];
  int ch, n, cell, locked, patch;
  float radius;
  float root[3];

  // Radius is determined in boundingBoxKernel
  radius = radiusd;
  // The root node has been created at the end of the boundingBoxKernel.
  // Cache the root data:
  for (l = 0; l < 3; l++)
    root[l] = bhpara->r[3 * bhpara->nnodes + l];
  // Maximum tree depth within the given thread.
  localmaxdepth = 1;
  // Skip the branch following and start from the root.
  skip = 1;
  // Number of loops for the threads sync algorithm
  lps = 0;
  // Increment to move among the bodies assigned to the given thread.
  // Hence, one should step over all other threads in GPU with
  // a quantity of blockDim.x * gridDim.x.
  inc = static_cast<int>(blockDim.x * gridDim.x);
  // Just a regular 1D GPU index
  i = static_cast<int>(threadIdx.x + blockIdx.x * blockDim.x);
  // AMD-specific threads sync
#if defined(__HIPCC__) and not defined(__CUDACC__)
  __syncthreads();
#endif

  // Iterate over all bodies assigned to thread.
  while (i < bhpara->nbodies) {
    if (skip != 0) {
      // New body, so start traversing at root. Skip it further.
      skip = 0;
      // Particle position corresponding to the given thread and block:
      for (l = 0; l < 3; l++)
        p[l] = bhpara->r[3 * i + l];
      // Let's start a moving via the tree from the root node 8 * nbodiesd:
      n = bhpara->nnodes;
      depth = 1;
      r = radius;
      // Determine which child to follow.
      // j=0..7 determines the octant in a binary representations
      j = 0;
      for (l = 0; l < 3; l++)
        if (root[l] < p[l])
          j += pow(2, l);
    }

    // Follow path to leaf cell. Should not happen at the first iteration of
    // this loop.
    ch = bhpara->child[n * 8 + j];

    // Global memory writing related threads sync
    if (lps++ > THREADS2) {
      // AMD-specific threads sync
#if defined(__HIPCC__) and not defined(__CUDACC__)
      atomicInc((unsigned int *)bhpara->max_lps, 0);
#else
      *bhpara->max_lps = lps;
#endif
    }
    //.. now wait for global memory updates. This impacts on race conditions and
    // frameworks level optimizations. Further kernels contain a similar code
    // fragment.
    __threadfence();

    // The child with the index higher than nbodiesd (number of particles) means
    // that it is an octant cell, not a body.
    // Actually, we need nnodesd == 8 * nbodiesd nodes for the cells storage.
    // Let's iterate this "while" loop before we will reach the particle
    // or an absence of a child:
    while (ch >= bhpara->nbodies) {
      n = ch;
      depth++;
      // Going down through the octree depth: radius split corresponds to the
      // cube split. Corresponding octants are cells, not bodies. Smaller radius
      // will be used until the next skip == 1
      r *= 0.5f;
      j = 0;
      // Determine which child octant to follow based on body coordinates.
      // j=0..7 determines the octant in a binary representations.
      for (l = 0; l < 3; l++)
        if (bhpara->r[3 * n + l] < p[l])
          j += pow(2, l);
      ch = bhpara->child[n * 8 + j];
    }
    // Now we are deep enough in the tree, passed all levels of cells and
    // reached the body (particle).
    if (ch != -2) { // Skip if child pointer is locked (-2) and try again later.
      locked = n * 8 + j;
      // Try to lock and iterate towards next body:
      if (ch == atomicCAS((int *)&bhpara->child[locked], ch, -2)) {
        // If we are here then child[locked] was equal to "ch" and now is
        // assigned to -2, it is locked. We will not came here in a case if this
        // child will be already locked by other threads because other particle
        // could already do it in other thread. Several particles simultaneously
        // could force a new octants' split. In this case we will not came here
        // and "i += inc" will be not executed below. Hence, the present loop
        // "while (i < nbodiesd)" will mandatory repeat in the given thread for
        // the given i-th particle until other threads will unlock this cell for
        // either a body insertion and/or new octant level cells creation.
        if (ch == -2) {
          // Cannot be here..
          *bhpara->err = 1;
          break;
        }
        if (ch == -1) {
          // If -1 (i.e. no child index) then just insert a new body
          bhpara->child[locked] = i;
        } else {
          patch = -1;
          // There already is a body and/or cell(s) in this position.
          // We should start from a loop of the deepest child cell
          // determination. Then, we need to create new cells and to distribute
          // existing and new bodies between these new cells.
          do {
            // If we are here then the tree bottom level is moving further:
            // Note, that the bottomd is not a global tree depth.
            // It is rather a tree size in 1D representation of nnodesd == 8 *
            // nbodiesd nodes. These lists will be correctly handled by the
            // sortKernel later.
            depth++;
            cell = atomicSub((int *)&bottomd, 1) - 1;
            if (cell <= bhpara->nbodies) {
              // This should not happen. A cell cannot have such index. Error.
              *bhpara->err = 1;
              bottomd = bhpara->nnodes;
            }
            // The "patch" is saving the information about a first cell created
            // in the current thread before it continues to dive into the tree
            // depth. Hence, the "patch" defines the new branch inception index.
            patch = max(patch, cell);

            // Note: "j" is defined by the body position against the parent cell
            // within above "while (ch >= nbodiesd)" loop.
            // The new cell will be placed below relatively the center of
            // corresponding j-th octant:
            for (l = 0; l < 3; l++)
              pos[l] = static_cast<float>((j >> l) & 1) * r;
            // Note, that negative octants correspond to pos[l] == 0 and
            // positive octants correspond to pos[l] == r.

            // Going down through the octree depth: radius split corresponds to
            // the cube split. Corresponding octants are cells, not bodies.
            // Smaller radius will be used until the next skip == 1
            r *= 0.5f;

            // Init the node weight coefficients.
            // Note: particles has mass=1.0 is defined in allocBHmemCopy().
            bhpara->mass[cell] = -1.0f;

            // The startd array is crucial for the sortKernel.
            // The original root node already has startd = 0.
            // Here, we assign -1 to the cells in contrast to bodies.
            // Bodies do not need this array. They need only array sortd,
            // which will be defined in the sortKernel for a usage by the force
            // and energy calculation kernels.
            bhpara->start[cell] = -1;

            // Now, let's save the cell coordinates locally (pos[l]) and
            // globally (xd[3 * cell + l]). This location should be shifted from
            // the octant center defined above (pos[l] before this assignment).
            // Parent cell coordinates bhpara->r[3 * n + l] will be added.
            // Parent radius now is equal to 2 * r, where "r" is already updated
            // above: r *= 0.5f. Hence, the negative octant is defined above by
            // pos[l] == 0 and positive - by pos[l] == 2 * r. In order to
            // transform these coordinates into relative octant positions, we
            // need to add -r to obtain -r and r for negative and positive
            // octants. Now, the child (cell) octants centers are deriving from
            // the parent (n) octant center:
            for (l = 0; l < 3; l++)
              pos[l] = bhpara->r[3 * cell + l] =
                  bhpara->r[3 * n + l] - r + pos[l];

            // By default, the new cell has no children in all k-th octants:
            for (k = 0; k < 8; k++)
              bhpara->child[cell * 8 + k] = -1;

            // This condition should always be true cause "patch" is -1 at the
            // beginning and the bottomd/cell reduces further.
            if (patch != cell) {
              // New cell is assigned as a child of previous "n" parent:
              bhpara->child[n * 8 + j] = cell;
            }

            // pos[l] already contains the child cell coordinates.
            // Let's assign "child" then. First the octant should be selected:
            j = 0;
            for (l = 0; l < 3; l++)
              if (pos[l] < bhpara->r[3 * ch + l])
                j += pow(2, l);
            // New element just appeared in the chain of cells. Hence, that what
            // supposed to be a child ("ch") before entering the present
            // iteration, now will be a child of the new cell (after this
            // smallest octant split into new octants):
            bhpara->child[cell * 8 + j] = ch;

            // Now cell is claimed to be a parent of further iteration of the
            // present loop.
            n = cell;
            j = 0;
            __threadfence();
            // Let's handle the particle position (p[l]) corresponding to the
            // given thread and block against new octant cell (pos[l]):
            for (l = 0; l < 3; l++)
              if (pos[l] < p[l])
                j += pow(2, l);

            // Now the current cell's child should be considering in the new
            // particle new octant:
            ch = bhpara->child[n * 8 + j];
            // Repeat until the two bodies will be different children.
            // Hence, the current "child" should have no children.
            // It is equivalent to an absence of other particles
            // in the i-th particle new smallest octant, otherwise we should
            // split octants further until these two particles will come to
            // different octants:
          } while (ch >= 0);

          // i-th particle assignment as a child to the last created cell:
          bhpara->child[n * 8 + j] = i;
          // Push out the subtree among the whole grid.
          // Data setting must be completed after this point.
          __threadfence();
          // The final locked child index is defined by a maximal cell index,
          // i.e. by a beginning of the new tree of cells created within
          // the loop "while (ch >= 0)".
          // The "patch" defines the new just created branch inception index:
          bhpara->child[locked] = patch;
        }

        localmaxdepth = max(depth, localmaxdepth);
        // Each thread started from the skip=1 and made the above
        // tree building loop procedure (while (ch >= 0)).
        // They should do the same for remaining (each) particles.
        // Note, that bottomd is a global variable and it is already updated.
        // This is already taken into account in further sortKernel logic.
        // Hence, move on to the next body assigned to the given thread:
        i += inc;
        skip = 1;
        lps = 0;
      }
    }
    // AMD-specific threads sync
#if defined(__HIPCC__) and not defined(__CUDACC__)
    __syncthreads();
#endif
  }
  // Record maximum tree depth:
  atomicMax((int *)&maxdepthd, localmaxdepth);
#if defined(__HIPCC__) and not defined(__CUDACC__)
  __syncthreads();
#endif
}

/******************************************************************************/
/*** compute center of mass ***************************************************/
/******************************************************************************/

__global__ __launch_bounds__(THREADS3, FACTOR3) void summarizationKernel() {
  int i, j, k, l, im, ch, inc, missing, missing_max, cnt, bottom, lps;
  int iteration, repeat_flag;
  // the node "mass" and its count respectively:
  float m, cm;
  // position of equivalent total dipole and its magnitude:
  // (like a mass and the center of mass)
  float p[3], u[3];
  // Per-block BH tree caching:
  __shared__ int child[THREADS3 * 8];

  // no children by default:
  for (i = 0; i < 8; i++)
    child[i * THREADS3 + threadIdx.x] = -1;
  bottom = bottomd;
  // Increment towards other particles assigned to the given thread:
  inc = static_cast<int>(blockDim.x * gridDim.x);
  // Nodes iteration "k" should start from the "bottomd" level of the cells,
  // which is a minimal index of the last created cell.
  // Starting "k" value should be aligned using the warp size
  // according to the designed threads performance.
  // k = (bottom & (-WARPSIZE)) + threadIdx.x + blockIdx.x * blockDim.x;
  k = bottom + static_cast<int>(threadIdx.x + blockIdx.x * blockDim.x);
  // Threads below the bottom line could proceed to their next cells.
  // if (k < bottom) k += inc;

  // Assume no missing children:
  missing = 0;
  iteration = 0;
  repeat_flag = 0;
  __syncthreads(); // throttle
  // threads sync related
  lps = 0;
  // Iterate over all cells (not particles) assigned to the thread:
  while (k <= bhpara->nnodes) {
    if (lps++ > THREADS3) {
      *bhpara->max_lps = lps;
      __threadfence();
    }
    if (bhpara->mass[k] < 0.) {
      iteration++;
      if (missing == 0) {
        // New cell, so initialize:
        cm = 0.0f;
        for (l = 0; l < 3; l++) {
          p[l] = 0.0f;
          u[l] = 0.0f;
        }
        cnt = 0;
        j = 0;
        for (i = 0; i < 8; i++) {
          ch = bhpara->child[k * 8 + i];
          if (ch >= 0) {
            if (i != j) {
              // Move child to front (needed later for a speed only).
              // The child's octant change is incorrect from
              // a tree organization perspective. However, the sum
              // will be the same.
              bhpara->child[k * 8 + i] = -1;
              bhpara->child[k * 8 + j] = ch;
            }
            // Cache a missing child in the block shared memory:
            child[missing * THREADS3 + threadIdx.x] = ch;
            m = bhpara->mass[ch];
            // Is a child the particle? Only particles have non-negative mass
            // initialized originally. Another option: a cell which already
            // aggregated masses of other cells and particles. "missing" means
            // that a non-zero contribution of such kind is missing:
            missing++;
            if (m >= 0.0f) {
              // child is ready
              missing--;
              // The child is a cell, not a body (ch >= nbodiesd).
              // Also, the previous condition (m >= 0.0f) reveals
              // that its' children total mass is already calculated.
              // Hence, below command "countd[k] = cnt" is already executed by
              // other threads/blocks and we can add this count
              if (ch >= bhpara->nbodies) { // count bodies (needed later)
                // As far as a child is a cell, its "countd" was already
                // calculated.
                cnt += bhpara->count[ch] - 1;
              }
              // add child's contribution
              cm += m;
              for (l = 0; l < 3; l++) {
                p[l] += bhpara->r[3 * ch + l] * m;
                u[l] += bhpara->u[3 * ch + l];
              }
            }
            j++;
          } // if (ch >= 0)
        }
        missing_max = missing;
        // Count of childs:
        cnt += j;
      }

      //__syncthreads();    // throttle

      if (missing != 0) {
        for (im = 0; im < missing_max; im++) {
          // poll missing child
          ch = child[im * THREADS3 + threadIdx.x];
          if (ch >= 0) {
            m = bhpara->mass[ch];
            // Is a child the particle? Only particles have non-negative mass
            // initialized originally. Another option: a cell which already
            // aggregated masses of other cells and particles.
            if (m >= 0.0f) {
              // child is now ready
              missing--;
              child[im * THREADS3 + threadIdx.x] = -1;
              // The child is a cell, not a body (ch >= nbodiesd).
              if (ch >= bhpara->nbodies) {
                // count bodies (needed later)
                cnt += bhpara->count[ch] - 1;
              }
              // add child's contribution
              cm += m;
              for (l = 0; l < 3; l++) {
                p[l] += bhpara->r[3 * ch + l] * m;
                u[l] += bhpara->u[3 * ch + l];
              }
            } // m >= 0.0f
          }   // ch >= 0
        }     // missing_max
        // repeat until we are done or child is not ready
      }

      //__syncthreads(); // throttle

      // (missing == 0) could be true and threads will move to next particles (k
      // += inc) only if previous conditions (m >= 0.0f) will be true. It can
      // happen only if cell will obtain the mass (only here below: "massd[k] =
      // cm") or they will find the very last childs: particles. Before that:
      // do/while loop will continue.
      if (missing == 0) {
        // all children are ready, so store computed information
        bhpara->count[k] = cnt;
        m = 1.0f / cm;
        for (l = 0; l < 3; l++) {
          bhpara->r[3 * k + l] = p[l] * m;
          bhpara->u[3 * k + l] = u[l];
        }
        __threadfence(); // make sure data are visible before setting mass
        bhpara->mass[k] = cm;
        __threadfence();
        k += inc;
        iteration = 0;
        lps = 0;
      }
      //__syncthreads(); // throttle
      if (iteration > THREADS3 + 1) {
        k += inc;
        repeat_flag = 1;
        iteration = 0;
        missing = 0;
      }
    } else {
      k += inc;
    }
    if ((k > bhpara->nnodes) && (repeat_flag)) {
      repeat_flag = 0;
      missing = 0;
      k = bottom + static_cast<int>(threadIdx.x + blockIdx.x * blockDim.x);
    }
  } // while
}

/******************************************************************************/
/*** sort bodies **************************************************************/
/******************************************************************************/

// This kernel concurrently places the bodies into an array such that the bodies
// appear in the same order in the array as they would during an in-order
// traversal of the octree. This sorting groups spatially close bodies (in the
// same octant cells) together, and these grouped bodies are crucial to speed up
// forceCalculationKernel and energyCalculationKernel
__global__ __launch_bounds__(THREADS4, FACTOR4) void sortKernel() {
  int i, k, ch, dec, start, bottom, lps;

  bottom = bottomd;
  dec = static_cast<int>(blockDim.x * gridDim.x);
  // Start from the end of the nnodesd == 8 * nbodiesd.
  // Reverse order is required now cause octant cells which are more close
  // to the root have a larger count of entities inside (countd[k]).
  // Particles should be sorted over all entities count in the tree array
  // representation made by treeBuildingKernel.
  k = bhpara->nnodes + 1 - dec +
      static_cast<int>(threadIdx.x + blockIdx.x * blockDim.x);
  // threads sync related
  lps = 0;

  // iterate over all cells assigned to thread
  while (k >= bottom) {
    start = bhpara->start[k];
    // Threads sync related
    if (lps++ > THREADS4) {
      *bhpara->max_lps = lps;
      __threadfence();
    }
    // Let's start from the root which has only startd=0 defined
    // in boundingBoxKernel. All other bodies and cells have -1.
    if (start >= 0) {
      for (i = 0; i < 8; i++) {
        ch = bhpara->child[k * 8 + i];
        if (ch >= bhpara->nbodies) {
          // child is a cell
          bhpara->start[ch] = start;  // set start ID of child
          start += bhpara->count[ch]; // add # of bodies in subtree
        } else if (ch >= 0) {
          // Child is a body.
          // This particle should be saved with a stepping over
          // a count of particles in the cells.
          // treeBuildingKernel already has ordered cells in a
          // linear array way. The sortKernel just order random particle
          // indices in the same order. Hence, they will be much faster accessed
          // by forceCalculationKernel and energyCalculationKernel.
          bhpara->sort[start] = ch; // record body in 'sorted' array
          start++;
        }
      }
      k -= dec; // move on to next cell
      // Threads sync related
      lps = 0;
    }
    //__syncthreads(); // throttle
  }
}

/******************************************************************************/
/*** compute dipole-dipole force and torque ***********************************/
/******************************************************************************/

__global__ __launch_bounds__(THREADS5, FACTOR5) void forceCalculationKernel(
    dds_float pf, float *force, float *torque) {
  int i, j, k, l, n, depth, base, sbase, diff, t;
  float tmp;
  // dr is a distance between particles.
  // f,h, and N are a target force, field, and torque respectively.
  // u and uc are dipole moments of two particles.
  float dr[3], f[3], h[3], u[3], uc[3], N[3];
  // Shared memory aggregates of each warp-specific stacks
  // with the MAXDEPTH size per each warp:
  // "node" is the BH octant sub-cell in the stack.
  // "pos"=0..7 - which octant we are examining now in the stack.
  // dq is an array used to determine that the given BH cell is far enough.
  __shared__ int pos[MAXDEPTH * THREADS5 / WARPSIZE],
      node[MAXDEPTH * THREADS5 / WARPSIZE];
  __shared__ float dq[MAXDEPTH * THREADS5 / WARPSIZE];
  float b, b2, d1, dd5;
  float bb2d7, umd5;

  // Zero thread of the block initialize shared data for all warps.
  if (0 == threadIdx.x) {
    // Precompute values that depend only on tree level.
    // The method's parameters (a trade-off accuracy/performance)
    // which determine that the
    // cell is far enough are "itolsqd" and
    // "epssqd" which define a fraction of the octant cell and
    // an additive distance respectively.
    // Their joint contribution for the given tree depth are
    // calculated within the array dq[i], which will
    // be compared later with squared distance between the particle
    // and the cell depending on a cell level.
    // Original tree box edge (2*radiusd) should be divided *0.5
    // as much as the tree depth takes place.
    tmp = radiusd;
    dq[0] = tmp * tmp * *itolsqd;
    for (i = 1; i < maxdepthd; i++) {
      dq[i] = dq[i - 1] * 0.25f;
      dq[i - 1] += *epssqd;
    }
    dq[i - 1] += *epssqd;

    // Only maximal Barnes-Hut tree depth is allowed.
    // This error is technically possible, however, most applications
    // are far from the 1/2^32 particles' convergence.
    if (maxdepthd > MAXDEPTH) {
      *bhpara->err = maxdepthd;
    }
  }
  __syncthreads();

  // Only maximal Barnes-Hut tree depth is allowed.
  if (maxdepthd <= MAXDEPTH) {
    // How many warps are behind the current thread (?):
    base = static_cast<int>(threadIdx.x) / WARPSIZE;
    // Figure out first thread in each warp (lane 0):
    sbase = base * WARPSIZE;
    // Initial stack index is its MAXDEPTH portion start for the given warp
    // count base:
    j = base * MAXDEPTH;

    // How far the thread is from the warp beginning (?):
    diff = static_cast<int>(threadIdx.x) - sbase;
    // Make multiple copies to avoid index calculations later:
    if (diff < MAXDEPTH) {
      // Each thread copies its own dq[] element to a part of
      // dq array dedicated to the given warp:
      dq[diff + j] = dq[diff];
    }
    __syncthreads();

    // Iterate over all bodies assigned to thread:
    for (k = static_cast<int>(threadIdx.x + blockIdx.x * blockDim.x);
         k < bhpara->nbodies; k += static_cast<int>(blockDim.x * gridDim.x)) {
      // Sorted body indexes assigned to me:
      i = bhpara->sort[k]; // get permuted/sorted index
      // Cache the particle position info:
      for (l = 0; l < 3; l++) {
        u[l] = bhpara->u[3 * i + l];
        h[l] = 0.0f;
        f[l] = 0.0f;
      }

      // Force will be calculated for i-th particle.
      // All other space is of interest, hence all cells will be considered.
      // Hence, we should start from the root node (whole Barnes-Hut cube).
      // Initialize iteration stack, i.e., push root node onto stack.
      // Let's start from zero octant.
      depth = j;
      if (sbase == threadIdx.x) {
        node[j] = bhpara->nnodes;
        pos[j] = 0;
      }

      while (depth >= j) {
        // Stack is not empty (depth is still higher than j-level).
        // Hence, there are still some children to consider.
        while ((t = pos[depth]) < 8) {
          // Node on top of stack has more children to process:
          n = bhpara->child[node[depth] * 8 + t]; // load child pointer
          if (sbase == threadIdx.x) {
            // I'm the first thread in the warp.
            // Let me check for the current depth level of the tree
            // whether we have more children among 8 octant cells.
            // Hence, let's go to the next octant in the next iteration only:
            pos[depth] = t + 1;
          }
          __threadfence_block();
          // There is a child (octant cell) with a dipole moment uxd[3 * n + l]
          // and the center position bhpara->r[3 * n + l]:
          if (n >= 0) {
            tmp = 0.0f; // compute distance squared
            for (l = 0; l < 3; l++) {
              dr[l] = -bhpara->r[3 * n + l] + bhpara->r[3 * i + l];
              tmp += dr[l] * dr[l];
            }

            // NOTE: i-th particle is specific for the given thread.
            // However, the above-mentioned index "i" is already sorted by the
            // sortKernel. Hence, below stack-diving operations will be
            // performed by different threads of the same warp almost
            // synchronously because adjacent threads have particles located in
            // the space not far from each other (=adjacent octants of
            // corresponding local maxdepth). Hence, such particles will have a
            // similar stack of the below execution from a perspective of other
            // Barnes-Hut tree branches distancing against them. Same global
            // memory fragments (child[]) will be loaded to the same warp
            // access and no corresponding threads' divergence will take place.
            // Only in this way, zero thread could control others: see "if
            // (sbase == threadIdx.x)" above and below.. and this control will
            // be correct. It will be even with a single stack arrays per the
            // single warp: The pos, node and dq array fragments are shared
            // between all threads of the whole warp.

            // Check if all threads agree that cell is far enough away (or is a
            // body, i.e. n < nbodiesd).
#if defined(__CUDACC__) && CUDA_VERSION >= 9000
            if ((n < bhpara->nbodies) ||
                __all_sync(__activemask(), tmp >= dq[depth])) {
#else
            if ((n < bhpara->nbodies) || __all(tmp >= dq[depth])) {
#endif
              if (n != i) {

                d1 = sqrtf(tmp /*, 0.5f*/);
                dd5 = __fdividef(1.0f, tmp * tmp * d1);
                b = 0.0f;
                b2 = 0.0f;
                umd5 = 0.0f;
                for (l = 0; l < 3; l++) {
                  uc[l] = bhpara->u[3 * n + l];
                  b += uc[l] * dr[l];
                  b2 += u[l] * dr[l];
                  umd5 += u[l] * uc[l];
                }
                umd5 *= -3.0f * dd5;

                bb2d7 = 15.0f * b * b2 * __fdividef(dd5, tmp);

                for (l = 0; l < 3; l++) {
                  h[l] += (b * 3.0f * dr[l] - tmp * uc[l]) * dd5;
                  f[l] += -dr[l] * (umd5 + bb2d7) +
                          3.0f * (b * u[l] + b2 * uc[l]) * dd5;
                }
              }
            } else {
              // If it is not then let's split octants further (more depth).
              // Push the cell onto the stack:
              depth++;
              if (sbase == threadIdx.x) {
                // The given warp should start from this child as a root
                // further:
                node[depth] = n;
                // Let's start from it zero octant:
                pos[depth] = 0;
              }
              __threadfence_block();
            }
          } else {
            // Early out because all remaining children are also zero.
            // We should move to the next octant or to the next depth if other
            // threads already checked other octants:
            depth = max(j, depth - 1);
          }
        }
        depth--; // Done with this level
      }

      // Torque:
      N[0] = u[1] * h[2] - u[2] * h[1];
      N[1] = u[2] * h[0] - u[0] * h[2];
      N[2] = u[0] * h[1] - u[1] * h[0];

      for (l = 0; l < 3; l++) {
        if (f[l] != f[l] || h[l] != h[l]) { // nan
          printf("Force Kernel: NAN in particle[%d]\n", i);
          printf("x = %f, y = %f, z = %f,\n", bhpara->u[3 * i + 0],
                 bhpara->u[3 * i + 1], bhpara->u[3 * i + 2]);
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
/*** compute energy
 * ************************************************************/
/******************************************************************************/

__global__ __launch_bounds__(THREADS5, FACTOR5) void energyCalculationKernel(
    dds_float pf, dds_float *energySum) {
  // NOTE: the algorithm of this kernel is almost identical to
  // forceCalculationKernel. See comments there.

  int i, j, k, l, n, depth, base, sbase, diff, t;
  float tmp;
  float dr[3], h[3], u[3], uc[3];
  __shared__ int pos[MAXDEPTH * THREADS5 / WARPSIZE],
      node[MAXDEPTH * THREADS5 / WARPSIZE];
  __shared__ float dq[MAXDEPTH * THREADS5 / WARPSIZE];
  dds_float sum = 0.0;
  HIP_DYNAMIC_SHARED(dds_float, res)

  float b, d1, dd5;

  if (0 == threadIdx.x) {
    tmp = radiusd;
    // precompute values that depend only on tree level
    dq[0] = tmp * tmp * *itolsqd;
    for (i = 1; i < maxdepthd; i++) {
      dq[i] = dq[i - 1] * 0.25f;
      dq[i - 1] += *epssqd;
    }
    dq[i - 1] += *epssqd;

    if (maxdepthd > MAXDEPTH) {
      *bhpara->err = maxdepthd;
    }
  }
  __syncthreads();

  if (maxdepthd <= MAXDEPTH) {
    // figure out first thread in each warp (lane 0)
    base = static_cast<int>(threadIdx.x) / WARPSIZE;
    sbase = base * WARPSIZE;
    j = base * MAXDEPTH;

    diff = static_cast<int>(threadIdx.x) - sbase;
    // make multiple copies to avoid index calculations later
    if (diff < MAXDEPTH) {
      dq[diff + j] = dq[diff];
    }
    __syncthreads();

    // iterate over all bodies assigned to thread
    for (k = static_cast<int>(threadIdx.x + blockIdx.x * blockDim.x);
         k < bhpara->nbodies; k += static_cast<int>(blockDim.x * gridDim.x)) {
      i = bhpara->sort[k]; // get permuted/sorted index
      // cache position info
      for (l = 0; l < 3; l++) {
        u[l] = bhpara->u[3 * i + l];
        h[l] = 0.0f;
      }

      // initialize iteration stack, i.e., push root node onto stack
      depth = j;
      if (sbase == threadIdx.x) {
        node[j] = bhpara->nnodes;
        pos[j] = 0;
      }

      while (depth >= j) {
        // stack is not empty
        while ((t = pos[depth]) < 8) {
          // node on top of stack has more children to process
          n = bhpara->child[node[depth] * 8 + t]; // load child pointer
          if (sbase == threadIdx.x) {
            // I'm the first thread in the warp
            pos[depth] = t + 1;
          }
          __threadfence_block();
          if (n >= 0) {
            tmp = 0.0f;
            for (l = 0; l < 3; l++) {
              dr[l] = -bhpara->r[3 * n + l] + bhpara->r[3 * i + l];
              tmp += dr[l] * dr[l];
            }
#if defined(__CUDACC__) && CUDA_VERSION >= 9000
            if ((n < bhpara->nbodies) ||
                __all_sync(
                    __activemask(),
                    tmp >= dq[depth])) { // check if all threads agree that cell
                                         // is far enough away (or is a body)
#else
            if ((n < bhpara->nbodies) ||
                __all(tmp >=
                      dq[depth])) { // check if all threads agree that cell
                                    // is far enough away (or is a body)
#endif
              if (n != i) {
                d1 = sqrtf(tmp /*, 0.5f*/);
                dd5 = __fdividef(1.0f, tmp * tmp * d1);
                b = 0.0f;
                for (l = 0; l < 3; l++) {
                  uc[l] = bhpara->u[3 * n + l];
                  b += uc[l] * dr[l];
                }

                for (l = 0; l < 3; l++)
                  h[l] += (b * 3.0f * dr[l] - tmp * uc[l]) * dd5;
              }
            } else {
              // push cell onto stack
              depth++;
              if (sbase == threadIdx.x) {
                node[depth] = n;
                pos[depth] = 0;
              }
              __threadfence_block();
            }
          } else {
            depth = max(j, depth - 1); // early out because all remaining
                                       // children are also zero
          }
        }
        depth--; // done with this level
      }

      for (l = 0; l < 3; l++) {
        sum += -u[l] * h[l];
        if (h[l] != h[l]) { // nan
          printf("Energy Kernel: NAN in particle[%d]\n", i);
          printf("x = %f, y = %f, z = %f,\n", bhpara->u[3 * i + 0],
                 bhpara->u[3 * i + 1], bhpara->u[3 * i + 2]);
          printf("hx = %f, hy = %f, hz = %f,\n", h[0], h[1], h[2]);
        }
      }
    }
  }

  sum *= pf;
  // the self-consistent field energy;
  // the Barnes-Hut algorithm, probably, does not allow to avoid this /2 cause
  // it is not symmetric:
  sum /= 2.0f;
  // Save per thread result into block shared mem
  res[threadIdx.x] = sum;
  // Sum results within a block
  __syncthreads(); // Wait til all threads in block are done
  dds_sumReduction_BH(res, &(energySum[blockIdx.x]));
}

// Required BH CUDA init.
void initBHgpu(int blocks) {
  dim3 grid(1, 1, 1);
  dim3 block(1, 1, 1);

  grid.x = blocks * FACTOR5;
  block.x = THREADS5;

  KERNELCALL(initializationKernel, grid, block);

  // According to the experimental performance optimization:
  cudaFuncSetCacheConfig(boundingBoxKernel, cudaFuncCachePreferShared);
  cudaFuncSetCacheConfig(treeBuildingKernel, cudaFuncCachePreferL1);
  cudaFuncSetCacheConfig(summarizationKernel, cudaFuncCachePreferShared);
  cudaFuncSetCacheConfig(sortKernel, cudaFuncCachePreferL1);
  cudaFuncSetCacheConfig(forceCalculationKernel, cudaFuncCachePreferL1);
  cudaFuncSetCacheConfig(energyCalculationKernel, cudaFuncCachePreferL1);

  cudaGetLastError(); // reset error value
}

// Building Barnes-Hut spatial min/max position box
void buildBoxBH(int blocks) {
  dim3 grid(1, 1, 1);
  dim3 block(1, 1, 1);

  grid.x = blocks * FACTOR1;
  block.x = THREADS1;

  cudaDeviceSynchronize();
  KERNELCALL(boundingBoxKernel, grid, block);
  cuda_safe_mem(cudaDeviceSynchronize());
}

// Building Barnes-Hut tree in a linear child array representation
// of octant cells and particles inside.
void buildTreeBH(int blocks) {
  dim3 grid(1, 1, 1);
  dim3 block(1, 1, 1);

  grid.x = blocks * FACTOR2;
  block.x = THREADS2;

  KERNELCALL(treeBuildingKernel, grid, block);
  cuda_safe_mem(cudaDeviceSynchronize());
}

// Calculate octant cells masses and cell index counts.
// Determine cells centers of mass and total dipole moments
// on all possible levels of the BH tree.
void summarizeBH(int blocks) {
  dim3 grid(1, 1, 1);
  dim3 block(1, 1, 1);

  grid.x = blocks * FACTOR3;
  block.x = THREADS3;

  KERNELCALL(summarizationKernel, grid, block);
  cuda_safe_mem(cudaDeviceSynchronize());
}

// Sort particle indexes according to the BH tree representation.
// Crucial for the per-warp performance tuning of forceCalculationKernel and
// energyCalculationKernel.
void sortBH(int blocks) {
  dim3 grid(1, 1, 1);
  dim3 block(1, 1, 1);

  grid.x = blocks * FACTOR4;
  block.x = THREADS4;

  KERNELCALL(sortKernel, grid, block);
  cuda_safe_mem(cudaDeviceSynchronize());
}

// Force calculation.
int forceBH(BHData *bh_data, dds_float k, float *f, float *torque) {
  int error_code = 0;
  dim3 grid(1, 1, 1);
  dim3 block(1, 1, 1);

  grid.x = bh_data->blocks * FACTOR5;
  block.x = THREADS5;

  KERNELCALL(forceCalculationKernel, grid, block, k, f, torque);
  cuda_safe_mem(cudaDeviceSynchronize());

  cuda_safe_mem(cudaMemcpy(&error_code, bh_data->err, sizeof(int),
                           cudaMemcpyDeviceToHost));
  return error_code;
}

// Energy calculation.
int energyBH(BHData *bh_data, dds_float k, float *E) {
  int error_code = 0;
  dim3 grid(1, 1, 1);
  dim3 block(1, 1, 1);

  grid.x = bh_data->blocks * FACTOR5;
  block.x = THREADS5;

  dds_float *energySum;
  cuda_safe_mem(cudaMalloc(&energySum, (int)(sizeof(dds_float) * grid.x)));
  // cleanup the memory for the energy sum
  cuda_safe_mem(cudaMemset(energySum, 0, (int)(sizeof(dds_float) * grid.x)));

  KERNELCALL_shared(energyCalculationKernel, grid, block,
                    block.x * sizeof(dds_float), k, energySum);
  cuda_safe_mem(cudaDeviceSynchronize());

  // Sum the results of all blocks
  // One energy part per block in the prev kernel
  thrust::device_ptr<dds_float> t(energySum);
  float x = thrust::reduce(t, t + grid.x);
  cuda_safe_mem(cudaMemcpy(E, &x, sizeof(float), cudaMemcpyHostToDevice));

  cuda_safe_mem(cudaFree(energySum));
  cuda_safe_mem(cudaMemcpy(&error_code, bh_data->err, sizeof(int),
                           cudaMemcpyDeviceToHost));
  return error_code;
}

// Function to set the BH method parameters.
void setBHPrecision(float *epssq, float *itolsq) {
  cuda_safe_mem(cudaMemcpyToSymbol(HIP_SYMBOL(epssqd), epssq, sizeof(float), 0,
                                   cudaMemcpyHostToDevice));
  cuda_safe_mem(cudaMemcpyToSymbol(HIP_SYMBOL(itolsqd), itolsq, sizeof(float),
                                   0, cudaMemcpyHostToDevice));
}

// An allocation of the GPU device memory and an initialization where it is
// needed.
void allocBHmemCopy(int nbodies, BHData *bh_data) {
  if (bh_data->nbodies == nbodies)
    return;

  bh_data->nbodies = nbodies;

  int devID = -1;
  EspressoGpuDevice dev;

  devID = cuda_get_device();
  cuda_get_device_props(devID, dev);

  bh_data->blocks = dev.n_cores;
  // Each node corresponds to a split of the cubic box in 3D space to equal
  // cubic boxes hence, 8 nodes per particle is a theoretical octree limit:
  bh_data->nnodes = bh_data->nbodies * 8;

  int n_total_threads = 1024 * bh_data->blocks;
  if (bh_data->nnodes < n_total_threads)
    bh_data->nnodes = n_total_threads;
  else
    bh_data->nnodes = (bh_data->nnodes / n_total_threads) * n_total_threads;

  if (bh_data->err != nullptr)
    cuda_safe_mem(cudaFree(bh_data->err));
  cuda_safe_mem(cudaMalloc((void **)&(bh_data->err), sizeof(int)));

  if (bh_data->max_lps != nullptr)
    cuda_safe_mem(cudaFree(bh_data->max_lps));
  cuda_safe_mem(cudaMalloc((void **)&(bh_data->max_lps), sizeof(int)));

  if (bh_data->child != nullptr)
    cuda_safe_mem(cudaFree(bh_data->child));
  cuda_safe_mem(cudaMalloc((void **)&(bh_data->child),
                           sizeof(int) * (bh_data->nnodes + 1) * 8));

  if (bh_data->count != nullptr)
    cuda_safe_mem(cudaFree(bh_data->count));
  cuda_safe_mem(cudaMalloc((void **)&(bh_data->count),
                           sizeof(int) * (bh_data->nnodes + 1)));

  if (bh_data->start != nullptr)
    cuda_safe_mem(cudaFree(bh_data->start));
  cuda_safe_mem(cudaMalloc((void **)&(bh_data->start),
                           sizeof(int) * (bh_data->nnodes + 1)));

  if (bh_data->sort != nullptr)
    cuda_safe_mem(cudaFree(bh_data->sort));
  cuda_safe_mem(cudaMalloc((void **)&(bh_data->sort),
                           sizeof(int) * (bh_data->nnodes + 1)));

  // Weight coefficients of m_bhnnodes nodes: both particles and octant cells
  if (bh_data->mass != nullptr)
    cuda_safe_mem(cudaFree(bh_data->mass));
  cuda_safe_mem(cudaMalloc((void **)&(bh_data->mass),
                           sizeof(float) * (bh_data->nnodes + 1)));

  // n particles have unitary weight coefficients.
  // Cells will be defined with -1 later.
  auto *mass_tmp = new float[bh_data->nbodies];
  for (int i = 0; i < bh_data->nbodies; i++) {
    mass_tmp[i] = 1.0f;
  }
  cuda_safe_mem(cudaMemcpy(bh_data->mass, mass_tmp,
                           sizeof(float) * bh_data->nbodies,
                           cudaMemcpyHostToDevice));
  delete[] mass_tmp;
  // (max[3*i], max[3*i+1], max[3*i+2])
  // are the octree box dynamical spatial constraints
  // this array is updating per each block at each interaction calculation
  // within the boundingBoxKernel
  if (bh_data->maxp != nullptr)
    cuda_safe_mem(cudaFree(bh_data->maxp));
  cuda_safe_mem(cudaMalloc((void **)&(bh_data->maxp),
                           sizeof(float) * bh_data->blocks * FACTOR1 * 3));
  // (min[3*i], min[3*i+1], min[3*i+2])
  // are the octree box dynamical spatial constraints
  // this array is updating per each block at each interaction calculation
  // within the boundingBoxKernel
  if (bh_data->minp != nullptr)
    cuda_safe_mem(cudaFree(bh_data->minp));
  cuda_safe_mem(cudaMalloc((void **)&(bh_data->minp),
                           sizeof(float) * bh_data->blocks * FACTOR1 * 3));

  if (bh_data->r != nullptr)
    cuda_safe_mem(cudaFree(bh_data->r));
  cuda_safe_mem(
      cudaMalloc(&(bh_data->r), 3 * (bh_data->nnodes + 1) * sizeof(float)));

  if (bh_data->u != nullptr)
    cuda_safe_mem(cudaFree(bh_data->u));
  cuda_safe_mem(
      cudaMalloc(&(bh_data->u), 3 * (bh_data->nnodes + 1) * sizeof(float)));
}

// Populating of array pointers allocated in GPU device before.
// Copy the particle data to the Barnes-Hut related arrays.
void fillConstantPointers(float *r, float *dip, BHData bh_data) {
  cuda_safe_mem(cudaMemcpyToSymbol(HIP_SYMBOL(bhpara), &bh_data, sizeof(BHData),
                                   0, cudaMemcpyHostToDevice));
  cuda_safe_mem(cudaMemcpy(bh_data.r, r, 3 * bh_data.nbodies * sizeof(float),
                           cudaMemcpyDeviceToDevice));
  cuda_safe_mem(cudaMemcpy(bh_data.u, dip, 3 * bh_data.nbodies * sizeof(float),
                           cudaMemcpyDeviceToDevice));
}

#endif // BARNES_HUT
