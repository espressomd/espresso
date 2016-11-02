/*
  Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
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

#ifndef __INTEGRATE_SD_CUDA_KERNEL_HPP
#define __INTEGRATE_SD_CUDA_KERNEL_HPP

#include "integrate_sd.hpp"

#ifdef SD

// This computes the farfield contribution.
// r is the vector of [x_1, y_1, z_1, x_2, y_2, z_2, ...]
// N is the number of particles
// self_mobility is 1./(6.*PI*viscosity*radius)
// a is the particle radius
// mobility is the mobility matrix which will be retruned
// L is the boxlength
__global__ void sd_compute_mobility_matrix(const real * r, int N, real self_mobility, real a, real * mobility);



/// This computes the farfield contribution of the mobility with ewald summation
/// this version assumes that the real-space cutoff is smaller than the boxlength
// r is the vector of [x_1, y_1, z_1, x_2, y_2, z_2, ...]
// N is the number of particles
// self_mobility is 1./(6.*PI*eta*a)
// a is the particle radius
// mobility is the mobility matrix which will be retruned
// L_d is the boxlength
// cutoff the realspace cutoff distance
// xi the splitting parameter as defined by Beenakker 1986
// xa  = xi * a
// xa3 = xa * xa * xa
__global__ void sd_compute_mobility_matrix_real_short(const real * r, int N, real self_mobility, real a, const real * L_g, real * mobility,
						      real cutoff, real xi, real xa, real xa3);

/// This computes the farfield contribution of the mobility with ewald summation
/// this kernel computes the sines and cosinus.
// r is the vector of [x_1, y_1, z_1, x_2, y_2, z_2, ...]
// N is the number of particles
// vecs are the k vectors
// num is the number of k vectors
// cosines is the pointer where the cosines : cos( position \times k-vector) are saved
// sines the same for sines
// ldd the rounded number of particles
__global__ void sd_compute_mobility_sines(const real * r, int N, const real * vecs, int num, real * sines, real * cosines, int ldd);

/// adds to each of the diagonal elemnts of the size*size matrix \param matrix
/// with lda \param lda 1
__global__ void sd_add_identity_matrix(real * matrix, int size, int lda);

void _cudaCheckError(const char *msg, const char * file, const int line);
// this computes the near field
// it calculates the ResistanceMatrix
__global__ void sd_compute_resistance_matrix(const real * r, int N, real self_mobility, real a, const real * L, real * resistance, int * myInfo);
// this computes the near field
// it calculates the Sparse ResistanceMatrix

__global__ void sd_compute_resistance_matrix_sparse(const real * pos, const int N, const real self_mobility,const real a,const real * _L,
  						    real * resistance,const int * col_idx,const int * row_l, int * myInfo);
// TODO: make the order of arguments uniform (and logical?)
// TODO: description here
__global__ void sd_compute_brownian_force_nearfield(const real * r,const real * gaussian_nf,int N,const real * L, real a, real self_mobility,real * brownian_force_nf);


/// This computes the farfield contribution of the mobility with ewald summation
/// this kernel computes the sines and cosinus.
// forces is the vector of [fx_1, fy_1, fz_1, fx_2, fy_2, fz_2, ...]
// a is the particle radius
// mobility is the mobility matrix which will be retruned (in/out)
// L_d is the boxlength
// cutoff the wavespace cutoff distance
__global__ void sd_wavepart_sum(const real * const forces, const real * const vecs,    int num, const  real * const matrices_d,
                                const real * const sines,  const real * const cosines, int ldd, int N, real * sin_sum, real * cos_sum, real max);


__global__ void sd_wavepart_assemble(const int num, const real * const sines, const real * const cosines, const real * const sin_sum,
				     const real * const cos_sum, const int ldd, real * out, real max, int N, const real factor);


/// Add the wavepart contribution of the mobility to the matrix
__global__ void sd_wavepart_addto_matrix(const int num, const  real * const matrices_d, const real * const sines,  const real * const cosines,
					 const int ldd, const int N, real * const mobility_d, const int mat_ldd);

// make sure to have one thread per particle
__global__ void sd_real_integrate_prepare( const real * r_d , real * disp_d, const real * L, real a, int N);
__global__ void sd_real_integrate( real * r_d , const real * disp_d, const real * L, real a, int N);


// this sets a block to zero
// matrix: pointer to the given matrix
// size  : the size of the matrix (in the example below 3N)
// ldd   : the leading dimension of the matrix
__global__ void sd_set_zero_matrix(real * matrix, int size, int ldd);


// this sets a block to zero
// data  : pointer to the given data
// size  : the size of the data
__global__ void sd_set_zero(real * data, int size);

// this sets a block to an given integer
// data  : pointer to the given data
// size  : the size of the data
// value : the value written to the data block
__global__ void sd_set_value(int * data, int size, int value);

// this sets a block to an given real
// data  : pointer to the given data
// size  : the size of the data
// value : the value written to the data block
__global__ void sd_set_value(real * data, int size, real value);


// This function reads the diagonal element entries of a matrix
// and stores them in the vector diag
// and the inverse of them in diag_i
__global__ void sd_get_diag(int size, const real * mat_a,int lda,real * diag,real * diag_i);



// implementation of a bucket sort algorithm
// puts all the N particles with given position pos 
// and particle radius a within the periodic boundary 
// conditions of boxSize L_i = bucketSize_i * bucketNum_i
// puts them in the list particleList
// pos                device array of particle position xyz
// bucketSize         device array with the number of buckets in x y and z direction
// bucketNum          device array with the size of a bucket in x y and z direction
// N                  number of particles
// particleCount      device array of the numbers of particles per bucket. must be initalized to zero
// particleList       device array of the partilces in each bucket
// maxParticlePerCell maximum particles per cell
// totalBucketNUm     bucketNum[0]*bucketNum[1]*bucketNum[2] - the total number of buckets
__global__ void sd_bucket_sort( const real * pos , const real * bucketSize, const int * bucketNum, int N, int * particleCount,
				int * particleList, int maxParticlePerCell, int totalBucketNum, int * particleToBucketList);


// matrix multiplication for sparse matrizes
// never use direct, use always wrapper
__global__ void sd_multiply_sparse_matrix_vector(const int size, const real factor, const real * matrix, const int ldd, const int ldd_short,
						 const int * col_idx, const int * row_l, const real * vec_in, real * vec_out );

// finding all interacting particles, given a certain iteraction range.
// this requires bucket sort and is therefor O(N)
// if interaction_range is > min(L)/2 particles can appeare more often in the list.
// self_interaction is excluded.
__global__ void sd_find_interacting_particles(const real * pos,const real * _L, const int N, int * interacting, int * num_interacting,
					      const int * particle_count, const int * particle_sorted_list, const real * bucket_size_,
					      const int * bucket_num_, const int * particle_to_bucket_list, const real interaction_range,
					      const int total_bucket_num);

// finding all interacting particles, given a certain iteraction range.
// this works without bucket sort and is O(N^2)
__global__ void sd_find_interacting_particles(const real * pos,const real * L_g, const int N, int * interacting, int * num_interacting,
					      const real interaction_range);



/// computing the inf norm of a given vector (or in other words the largest value)
/// if more than one block is used, this kernel does not perform the final reduction
/// size     : size of the vector
/// vec      : data of the vector
/// res      : the 1 norm (OUT)
__global__ void sd_nrm_inf(const int size, const int * const vec, int * res);



__global__ void sd_nrm1(const int size,const real * const vec, real * erg);

#endif /* SD */

#endif
