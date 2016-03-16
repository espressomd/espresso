
/*
   Copyright (C) 2010,2011,2012,2016 The ESPResSo project

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
#ifdef SD /* Terminates at end of file */
//#define SD_PC

/** \file integrate_sd_cuda.cu    Stokes dynamics integrator.
 *
 *  Here the calculation of the displacements is implemented.
*/


#include <stdio.h>
#include <iostream>
#include "cuda_runtime.h"
#include "curand.h"
#include <device_functions.h>
// for std::swap:
// C++98:
#include <algorithm>
// C++11:
//#include <utility>

#include "assert.h"
#include "integrate_sd_cuda.hpp"
#include "integrate_sd_cuda_debug.hpp"
#include "integrate_sd_cuda_kernel.hpp"
#include "integrate_sd.hpp"
#include <cublas_v2.h>
//#include "errorhandling.hpp"
#include "global.hpp"


// global variables for usage in this file
cublasHandle_t     cublas    = NULL;
curandGenerator_t  generator = NULL;
int reprint=-1;
/* *************************************************************************************************************** *
 * ********************************************     implementation    ******************************************** *
 * *************************************************************************************************************** */
/* *************************************************************************************************************** *
 * *******     III MM   MM PPP  L     EEEE MM   MM EEEE NN    N TTTTTTT  AAA  TTTTTTT III  OOO  NN    N    ******* *
 * *******      I  M M M M P  P L     E    M M M M E    N N   N    T    A   A    T     I  O   O N N   N    ******* *
 * *******      I  M  M  M PPP  L     EEE  M  M  M EEE  N  N  N    T    AAAAA    T     I  O   O N  N  N    ******* *
 * *******      I  M     M P    L     E    M     M E    N   N N    T    A   A    T     I  O   O N   N N    ******* *
 * *******     III M     M P    LLLL  EEEE M     M EEEE N    NN    T    A   A    T    III  OOO  N    NN    ******* *
 * *************************************************************************************************************** */
/* *************************************************************************************************************** */

/* *************************************************************************************************************** *
 * ********************************************     HOST-Functions    ******************************************** *
 * *************************************************************************************************************** */

// this calls all the functions to:
//  * generate the mobility matrix (farfield and nearfield)
//  * compute the displacements
//  * add the displacements to the positions
// PARAMTERS:
// box_l_h : the size of the box in x,y and z-direction, on the host (in)
// N       : Number of particles (in)
// pos_h   : position of the particles, simple* array on host (in and out)
// force_h : forces on the particles, simple* array on host (in)
// velo_h  : velocities of the particles, simple* array on host (in and out)
// * : a simple array is e.g. [x_1, y_1, z_1, x_2, y_2, z_2, ...]
void propagate_pos_sd_cuda(real * box_l_h, int N,real * pos_h, real * force_h, real * velo_h){
  //printVectorHost(pos_h,3*N,"pos after call");
  real viscosity=sd_viscosity;
  real radius   =sd_radius;
  if (viscosity  < 0){
    std::cerr << "The viscosity for SD was not set\n";
    errexit();
  }
  if (radius  < 0){
    std::cerr << "The particle radius for SD was not set\n";
    errexit();
  }
  if (time_step < 0){
    std::cerr << "The timestep was not set\n";
    errexit();
  }
  
  int lda=((3*N+31)/32)*32;
  
  //static cublasHandle_t cublas=NULL;
  if (cublas==NULL){
    cublasCall(cublasCreate(&cublas));
    //magma_init();
  }

  real * box_l_d=NULL;
  cuda_safe_mem(cudaMalloc((void**)&box_l_d, 3*sizeof(real)));
  cuda_safe_mem(cudaMemcpy(box_l_d,box_l_h,3*sizeof(real),cudaMemcpyHostToDevice));
  real * pos_d=NULL;
  cuda_safe_mem(cudaMalloc((void**)&pos_d, (DIM)*(((N+31)/32)*32)*sizeof(real)));
  cuda_safe_mem(cudaMemcpy(pos_d,pos_h,N*DIM*sizeof(real),cudaMemcpyHostToDevice));
  //printVectorDev(pos_d,3*N,"pos after copy");
  real * force_d=NULL;
  cuda_safe_mem(cudaMalloc((void**)&force_d, DIM*(((N+31)/32)*32)*sizeof(real)));
  cuda_safe_mem(cudaMemcpy(force_d,force_h,N*DIM*sizeof(real),cudaMemcpyHostToDevice));
  real * mobility_d=NULL;
  cuda_safe_mem(cudaMalloc((void**)&mobility_d, lda*N*3*sizeof(real)));
  real * disp_d=NULL;
  cuda_safe_mem(cudaMalloc((void**)&disp_d, DIM*(((N+31)/32)*32)*sizeof(real)));
  int myInfo_h[]={0,0,0};
  int * myInfo_d=NULL;
  cuda_safe_mem(cudaMalloc((void**)&myInfo_d, 3*sizeof(int)));
  cuda_safe_mem(cudaMemcpy(myInfo_d,myInfo_h,3*sizeof(int),cudaMemcpyHostToDevice));
  
  sd_compute_displacement( pos_d, N, viscosity, radius, box_l_d, box_l_h, mobility_d, force_d, disp_d, myInfo_d);
  
  int numBlocks = (N+numThreadsPerBlock-1)/numThreadsPerBlock;

  // copy before rescaling to get velocities
  cuda_safe_mem(cudaMemcpy(velo_h,disp_d,N*DIM*sizeof(real),cudaMemcpyDeviceToHost));
  // rescale displacements
  real alpha=time_step;
  cublasCall(cublasRscal( cublas, 3*N, &alpha, disp_d, 1));  
#ifdef SD_MAX_STEP
  sd_real_integrate_prepare<<< numBlocks , numThreadsPerBlock  >>>(pos_d , disp_d, box_l_d, sd_radius, N);
#endif
  sd_real_integrate<<< numBlocks , numThreadsPerBlock  >>>(pos_d , disp_d, box_l_d, sd_radius, N);
  
  // copy back the positions
  cuda_safe_mem(cudaMemcpy(pos_h,pos_d,N*DIM*sizeof(real),cudaMemcpyDeviceToHost));
  // save the displacements as velocities (maybe somebody is interested)
  //alpha=1/time_step;
  //cublasCall(cublasRscal(cublas, DIM*N, &alpha, disp_d, 1));
  
  
  cuda_safe_mem(cudaFree((void*)box_l_d));
  cuda_safe_mem(cudaFree((void*)pos_d));
  cuda_safe_mem(cudaFree((void*)force_d));
  cuda_safe_mem(cudaFree((void*)mobility_d));
  cuda_safe_mem(cudaFree((void*)disp_d));
  cuda_safe_mem(cudaFree((void*)myInfo_d));
}



// calculate the farfield and the nearfield and add them
// PARAMETERS:
// cublas    : a valid handle of cublas (in)
// r_d       : position of the particles on the device, size 3*N (in)
//             the form has to be [x_1, y_1, z_1, x_2, y_2, z_2, ...]
// N         : Number of particles (in)
// viscosity : viscositiy of the fluid (in)
// radius    : Particle radius (in)
// L_d       : boxsize in x y and z-directions (in)
// total_mobility_d: matrix of the computed total mobility, size 3*3*N*N (in/out, is overwritten)
void sd_compute_displacement(const real * r_d, int N, real viscosity, real radius,const real * L_d,const real * L_h,
			     const real * total_mobility_d, const real * force_d, real * disp_d, int * myInfo_d)
{
#ifndef SD_FF_ONLY
  real det_prec=1e-3;
#endif
  real ran_prec=sd_random_precision;
  static real ran_prec_last=1e-3;
  if ( ran_prec_last != ran_prec){
    ran_prec_last=ran_prec;
    printf("\nSetting the precision for the random part to %e\n",ran_prec);
  }
  cudaThreadSynchronize(); // just for debugging
  cudaCheckError("START");
  const int lda=((3*N+31)/32)*32;
  const int numBlocks = (N+numThreadsPerBlock-1)/numThreadsPerBlock;
  // sort particles in buckets
#ifndef SD_FF_ONLY
  bool _use_buckets=false;
  for (int d=0;d<3;d++){
    if (L_h[d] > 4*4*radius){ // check if we can make at least 4 boxes in any direction
      //_use_buckets=true;
      //#warning "buckets disabled - they are still buggy"
    }
  }
  const bool use_buckets=_use_buckets;
  // variables for the bucket version
  real * bucket_size=NULL;
  int * bucket_num=NULL;
  int * particle_count=NULL;
  int * particle_sorted_list=NULL;
  int * particle_to_bucket_list=NULL;
  int total_bucket_num=1;
  int max_particles_per_bucket=N;
  //const int bucket_size_factor=3; // is global
  if ( use_buckets ){
    int bucket_num_h[3];
    real bucket_size_h[3];
    real bucket_volume=1;
    real bucket_volume_with_border=1;
    for (int d=0;d<3;d++){
      bucket_num_h[d]=L_h[d]/4/radius*bucket_size_factor;
      //bucket_num_h[d]=bucket_num_h[d]>(2*bucket_size_factor+1)?bucket_num_h[d]:1; // buckets are propably only usefull if there are more than 3 ...
      bucket_num_h[d]=bucket_num_h[d]>N?N:bucket_num_h[d]; // never more than N**3 buckets ...
      total_bucket_num *= bucket_num_h[d];
      bucket_size_h[d]=L_h[d]/bucket_num_h[d];
      bucket_volume*=bucket_size_h[d];
      bucket_volume_with_border*=bucket_size_h[d]+2*radius;
    }
    max_particles_per_bucket=bucket_volume_with_border/(4*M_PI*radius*radius*radius/3)*0.74048;
    max_particles_per_bucket=max_particles_per_bucket > N || max_particles_per_bucket < 0?N:max_particles_per_bucket; // there cannot be more than all particles in one bucket
    cuda_safe_mem(cudaMalloc((void **) &bucket_size,             sizeof(real)*3));
    cuda_safe_mem(cudaMemcpy(bucket_size,bucket_size_h,          sizeof(real)*3, cudaMemcpyHostToDevice));
    cuda_safe_mem(cudaMalloc((void **) &bucket_num,              sizeof(int) *3));
    cuda_safe_mem(cudaMemcpy(bucket_num,bucket_num_h,            sizeof(int) *3, cudaMemcpyHostToDevice));
    cuda_safe_mem(cudaMalloc((void **) &particle_count,          sizeof(int) *total_bucket_num));
    sd_set_value<<<32,32>>>(particle_count,total_bucket_num,0);
    cuda_safe_mem(cudaMalloc((void **) &particle_to_bucket_list, sizeof(int) *N));
    cuda_safe_mem(cudaMalloc((void **) &particle_sorted_list,    sizeof(int) *total_bucket_num*max_particles_per_bucket));
    sd_bucket_sort<<<numBlocks, numThreadsPerBlock >>>(r_d, bucket_size, bucket_num, N, particle_count, particle_sorted_list,
						       max_particles_per_bucket, total_bucket_num, particle_to_bucket_list);
  }
#endif  
  
  // compute the mobility Matrix
  matrix mobility;
  mobility.size      = N*3;
  mobility.ldd       = ((N*3+31)/32)*32;
  mobility.ldd_short = ((N+31)/32)*32;
  real * &mobility_d=mobility.data;
  cuda_safe_mem(cudaMalloc( (void**)&mobility_d, lda*DIM*N*sizeof(real) ));
  assert(mobility_d);
  sd_set_zero_matrix<<<numBlocks, numThreadsPerBlock >>>(mobility_d,mobility.size, mobility.ldd);
  cudaThreadSynchronize(); // just for debugging
  cudaCheckError("");
#ifdef SD_NOT_PERIODIC
  sd_compute_mobility_matrix<<< numBlocks , numThreadsPerBlock  >>>(r_d,N,1./(6.*M_PI*viscosity*radius), radius, mobility_d);
  cudaThreadSynchronize(); // just for debugging
  cudaCheckError("compute mobility error");
#else
  real ew_prec=1e-8;
  real L_min=L_h[0];
  L_min=min(L_min,L_h[1]);
  L_min=min(L_min,L_h[2]);
  real L_max=L_h[0];
  L_max=max(L_max,L_h[1]);
  L_max=max(L_max,L_h[2]);
  assert(ew_prec < 1);
  real xi=2./L_min*sqrt(log(1./ew_prec));
  real kmax2=-4*xi*xi*log(ew_prec);
  kmax2*=SQR(L_max/2./M_PI);
  {
    static int printed=0;
    if ((printed%reprint == 1)){
      fprintf(stderr,"\nL_min is %f\n",L_min);
      fprintf(stderr,"xi is %f\n",xi);
      fprintf(stderr,"kmax is %e\n",sqrt(kmax2));
      fprintf(stderr,"selfmobility is %e\n",1./(6.*M_PI*viscosity*radius));
      //printed=1;
    }
    printed++;
  }
  //fprintf(stderr,"kmax is %e\n",sqrt(kmax2));
  sd_compute_mobility_matrix_real_short<<<numBlocks, numThreadsPerBlock>>>(r_d, N, 1./(6.*M_PI*viscosity*radius), radius, L_d, 
									   mobility.data, L_min/2., xi, xi*radius, xi*radius*xi*radius*xi*radius);
  //std::cerr << "kmax is " << ceil(sqrt(kmax2)) << std::endl;
  sd_compute_mobility_matrix_wave_cpu(ceil(sqrt(kmax2))+1, (ceil(sqrt(kmax2))+1)*(ceil(sqrt(kmax2))+1), mobility,
				      radius, L_h, xi, 1./(6.*M_PI*viscosity*radius), ew_prec);
  
  //cudaCheckError("");
  //mobility.hash_data();
  if (mobility.wavespace){
    cudaThreadSynchronize(); 
    cudaCheckError("sd_compute_mobility_matrix_real_short");
    sd_compute_mobility_sines<<<numBlocks, numThreadsPerBlock>>>(r_d, N, mobility.wavespace->vecs, mobility.wavespace->num,
								 mobility.wavespace->sines, mobility.wavespace->cosines, 
								 mobility.ldd_short);
    cudaThreadSynchronize(); 
    cudaCheckError("sd_compute_mobility_sines");

    //cuda_safe_mem(cudaMalloc( (void**)&mobility.dense,    mobility.ldd*mobility.size*sizeof(real) ));
    //cublasCall(cublasRcopy(cublas, mobility.size*mobility.ldd,mobility.data, 1, mobility.dense,1));
    if (false) {
      sd_wavepart_addto_matrix<<<dim3(numBlocks,N,1), dim3(numThreadsPerBlock,1,1),
                                 numThreadsPerBlock*3*sizeof(real)>>>(mobility.wavespace->num,   mobility.wavespace->matrices,
								      mobility.wavespace->sines, mobility.wavespace->cosines,
								      mobility.ldd_short,        N,
								      mobility.data,            mobility.ldd);
      cudaThreadSynchronize(); 
      cudaCheckError("add wavepart to matrix");
      mobility._free_wavespace();
    }
#ifdef SD_DEBUG    
    /* *********************************
     * **   Test to check wavespace   **
     * ********************************* *
    real * testvec=NULL;
    cuda_safe_mem(cudaMalloc( (void**)&testvec,    mobility.ldd*3*sizeof(real) ));
    if (true){
      if (generator == NULL){
	curandCall(curandCreateGenerator(&generator, CURAND_RNG_PSEUDO_DEFAULT));
      }
      curandCall(curandGenerateNormalReal( generator, testvec, ((mobility.size+1)/2)*2, 0, sqrt(2.)));
    } else {
      sd_set_zero<<<numBlocks, numThreadsPerBlock>>>(testvec, mobility.size);
      //real val=1e-20;
      //for (int i = 0 ; i < mobility.size;i++){
      //if (i > 30)
      //sd_set_value<<<1,1>>>(testvec+i,1,val);
      //val*=10;
      //}
      sd_set_value<<<1,1>>>(testvec+6,3,1);
    }
    std::cout << "called ";
    real * testout1 = testvec+mobility.ldd*1;
    real * testout2 = testvec+mobility.ldd*2;
    const real one=1;
    sd_Rgemv(&one, mobility, testvec, testout1);
    if (warnings > 2){
      mobility.print();
      }
#endif
    //sd_wavepart_addto_matrix<<<dim3(3,N,1),dim3(2,1,1),numThreadsPerBlock*3*sizeof(real)>>>(mobility.wavespace->num, 
    sd_wavepart_addto_matrix<<<dim3(numBlocks,N,1), dim3(numThreadsPerBlock,1,1),
                               numThreadsPerBlock*3*sizeof(real)>>>(mobility.wavespace->num,   mobility.wavespace->matrices,
								    mobility.wavespace->sines, mobility.wavespace->cosines,
								    mobility.ldd_short,        N,
								    mobility.data,             mobility.ldd);
    
    cudaCheckError("");
    mobility._free_wavespace();
    cudaCheckError("");
#ifdef SD_DEBUG
    sd_Rgemv(&one, mobility, testvec, testout2);
    const real minusone=-1;
    cublasCall(cublasRaxpy(cublas, mobility.size, &minusone, testout2, 1, testout1, 1));
    real erg;
    cublasCall(cublasRnrm2(cublas, mobility.size, testout1, 1, &erg));
    if (erg > 0.01){
      printVectorDev(testvec,  mobility.size, "in         ");
      printVectorDev(testout2, mobility.size, "correct out");
      printVectorDev(testout1, mobility.size, "difference ");
      mobility.print();
      fflush(stdout);
    }
    assert(erg < 0.01);
    */
#endif
  }
  //mobility.printStats();
  if (!isSymmetricDev(mobility)){
    if (warnings > 2){
      printMatrixDev(mobility);
    }
  }
  assert(isSymmetricDev(mobility));
  if (warnings > 20){
    //printMatrixDev(mobility);
    //mobility.printWavespace();
  }
#endif
    
#ifndef SD_FF_ONLY
  // compute the resistance matrix
  matrix resistance;
  resistance.size      = N*3;
  resistance.ldd       = ((N*3+31)/32)*32;
  resistance.ldd_short = ((N+31)/32)*32;
  int lda_short=resistance.ldd_short;
  resistance.is_sparse=true;
  const int max_part_in_lc=125;
  cuda_safe_mem(cudaMalloc( (void**)&resistance.data,    lda* 3*    max_part_in_lc* sizeof(real) ));
  cuda_safe_mem(cudaMalloc( (void**)&resistance.col_idx, lda_short* max_part_in_lc* sizeof(int) ));
  cuda_safe_mem(cudaMalloc( (void**)&resistance.row_l,   lda_short*sizeof(int) ));
  
  if (use_buckets){
    cudaThreadSynchronize(); 
    cudaCheckError("finding interaction error");
    sd_find_interacting_particles<<<numBlocks, numThreadsPerBlock>>>(r_d, L_d, N, resistance.col_idx, resistance.row_l, 
								     particle_count, particle_sorted_list, bucket_size,
								     bucket_num, particle_to_bucket_list, 4*radius,total_bucket_num);
    cudaThreadSynchronize(); 
    cudaCheckError("finding interaction error");
  } else {
    cudaThreadSynchronize(); 
    cudaCheckError("finding interaction error");
    sd_find_interacting_particles<<<numBlocks, numThreadsPerBlock>>>(r_d, L_d, N, resistance.col_idx, resistance.row_l, 
								     4*radius);
    cudaThreadSynchronize(); 
    cudaCheckError("finding interaction error");
    
  }
  //#ifdef SD_DEBUG
  {
    assert(sd_find_max(resistance.row_l,N) < max_part_in_lc );
  }
  //#endif // SD_DEBUG
  sd_compute_resistance_matrix_sparse<<<numBlocks, numThreadsPerBlock>>>(r_d, N, -1./(6.*M_PI*viscosity*radius), radius, L_d,
									   resistance.data, resistance.col_idx, resistance.row_l, myInfo_d);
    /*} else {
  resistance.is_sparse=false;
    cuda_safe_mem(cudaMalloc( (void**)&resistance_d, lda*N*3*sizeof(real) ));
    sd_set_zero_matrix<<<numBlocks, numThreadsPerBlock >>>(resistance.data,resistance.size,resistance.ldd);
    cudaThreadSynchronize(); // debug
    cudaCheckError("sd_set_zero");
    sd_compute_resistance_matrix<<< numBlocks , numThreadsPerBlock  >>>(r_d,N,-1./(6.*M_PI*viscosity*radius), radius, L_d, resistance_d, myInfo_d);
    //#warning "bug"
    }*/
  cudaThreadSynchronize(); // we need both matrices to continue;
  cudaCheckError("compute resistance error");
#endif
  //printVectorDev(resistance.data,resistance.size,"resis first line");
#ifdef SD_DEBUG
  assert(!hasAnyNanDev(mobility_d,N*3*lda));
#ifndef SD_FF_ONLY
  if (hasAnyNanDev(resistance_d,N*3*lda)){
    printPosDev(r_d,N,"positions");
    //printMatrixDev(resistance_d,lda, 3*N,"Resistance with nans ...");
  }
  assert(!hasAnyNanDev(resistance_d,N*3*lda));
  assert(isSymmetricDev(resistance_d,lda,N*3));
#endif
#endif
  // initialize displacement
  static int last_N=-1;
  static real * last_det_disp=NULL;
  if (N == last_N && last_det_disp!=NULL){
    cublasCall(cublasRcopy(cublas, N*3, last_det_disp,1, disp_d,1));
  } else {
    sd_set_zero<<<numBlocks, numThreadsPerBlock>>>(disp_d,N*3);
  }
  //printVectorDev(r_d, 3*N,"pos");
#ifdef SD_FF_ONLY
  const real one=1;
  //cublasCall(cublasRgemv(cublas, CUBLAS_OP_N, 3*N, 3*N, &one, mobility_d, lda, force_d, 1, &zero, disp_d, 1));
  sd_Rgemv(&one, mobility,force_d,disp_d);
#else
  real err_force = sd_iterative_solver(mobility, resistance, force_d, disp_d,det_prec,0, true);
#endif
#ifdef SD_DEBUG
  if (hasAnyNanDev(disp_d,N*3)){
    printVectorDev(disp_d,N*3,"disp");
    printVectorDev(force_d,N*3,"force");
#ifndef SD_FF_ONLY
    //printMatrixDev(resistance_d,lda,N*3,"resistance produces nans?");
#endif
  }
  assert(!hasAnyNanDev(disp_d,N*3));
#endif
  // save deterministic displacement
  if (N > last_N){
    if (last_det_disp!= NULL){
      cuda_safe_mem(cudaFree((void*)last_det_disp));
      last_det_disp=NULL;
    }
  }
  if (last_det_disp==NULL){
    cuda_safe_mem(cudaMalloc( (void **) &last_det_disp, N*3*sizeof(real)));        assert(last_det_disp !=NULL);
  }
  cublasCall(cublasRcopy(cublas, N*3, disp_d, 1, last_det_disp,1));
  last_N=N;
  {
    static int printed=0;
    if (!printed){
      ;//printMatrixDev(mobility);
    }
  }
  // brownian part
  if (temperature > 0 ){
#ifndef SD_FF_ONLY
    int myInfo_h[3];
    cuda_safe_mem(cudaMemcpy(myInfo_h,myInfo_d,3*sizeof(int),cudaMemcpyDeviceToHost));
#else
    const int myInfo_h[3]={0,0,0};
#endif
    int N_ldd = ((N+31)/32)*32;
    int num_of_rands = N_ldd*myInfo_h[2]*2*DIM+N_ldd*DIM; // has to be even!
    if (myInfo_h[0]){
      fprintf(stderr,"We had %d overlapping particles!\n",myInfo_h[0]);
    }
    real * brownian_force_nf = NULL;
    if (myInfo_h[2]){
      cuda_safe_mem(cudaMalloc( (void**)&brownian_force_nf, (3*N)*sizeof(real) ));     assert(brownian_force_nf != NULL);
    }
    real * brownian_force_ff = NULL;
    cuda_safe_mem(cudaMalloc( (void**)&brownian_force_ff, (3*N)*sizeof(real) ));     assert(brownian_force_ff != NULL);
    real * gaussian = NULL;
    cuda_safe_mem(cudaMalloc( (void**)&gaussian, (num_of_rands)*sizeof(real) ));     assert(gaussian != NULL);
    real * gaussian_ff = gaussian;
    real * gaussian_nf = gaussian+N_ldd*DIM;
    unsigned long long *      sd_random_generator_offset     = (unsigned long long *)sd_random_state;
    unsigned long long *      sd_seed_ull                    = (unsigned long long *)sd_seed;
    static unsigned long long sd_random_generator_offset_last= *sd_random_generator_offset;
    if (generator == NULL){
      curandCall(curandCreateGenerator(&generator, CURAND_RNG_PSEUDO_DEFAULT));
      curandCall(curandSetPseudoRandomGeneratorSeed(generator, *sd_seed_ull));
      curandCall(curandSetGeneratorOrdering( generator, CURAND_ORDERING_PSEUDO_BEST));
      curandCall(curandSetGeneratorOffset( generator, *sd_random_generator_offset));
    }
    if (* sd_random_generator_offset != sd_random_generator_offset_last){
      curandCall(curandSetGeneratorOffset( generator, *sd_random_generator_offset));
    }
    //#ifdef FLATNOISE
    // this does not work yet:
    //curandCall(curandGenerateUniformReal(generator, gaussian_d, num_of_rands, 0, sqrt(24.*temperature/time_step)));
    //#else
    sd_random_generator_offset[0]+=num_of_rands/2;
    curandCall(curandGenerateNormalReal(generator, gaussian, num_of_rands, 0, sqrt(4*temperature/time_step)));
    sd_random_generator_offset_last=sd_random_generator_offset[0];
    //#endif
    if (myInfo_h[2]){
      sd_set_zero<<<64,192>>>(brownian_force_nf,3*N);
      sd_compute_brownian_force_nearfield<<<numBlocks, numThreadsPerBlock>>>(r_d, gaussian_nf, N, L_d, radius,
									     1./(6.*M_PI*viscosity*radius), brownian_force_nf);
    }// end of near field
    int size=N*3;
    real alpha=1/sqrt(2.);
    // in the case of FF_ONLY this computes directly the displacement
    real E_cheby=sd_compute_brownian_force_farfield(mobility, gaussian_ff, ran_prec, brownian_force_ff);
    //real bff_nrm;
    //cublasCall(cublasRnrm2(cublas, size, brownian_force_ff, 1, &bff_nrm));
    if (E_cheby > 10*ran_prec){
      E_cheby=sd_compute_brownian_force_farfield(mobility, gaussian_ff, ran_prec, brownian_force_ff);
      //cublasCall(cublasRnrm2(cublas, size, brownian_force_ff, 1, &bff_nrm));
      if (warnings > 1) fprintf(stderr, "Recalculating the Chebyshev-polynome\n");
    }
    if ((E_cheby>100*ran_prec && warnings) || (E_cheby > 10*ran_prec &&  warnings > 1) ){
      fprintf(stderr, "The error of the Chebyshev-approximation was %7.3f%%\n",E_cheby*100);
      //printVectorDev(mobility.wavespace->sines,mobility.wavespace->max*mobility.ldd_short,"sines");
      //printVectorDev(mobility.wavespace->matrices,mobility.wavespace->max*6,"matrices");
    }
    //printVectorDev(gaussian_ff,size, "gaussian");
    //printVectorDev(brownian_force_ff, size, "brownian force");
    //printVectorDev(mobility_d, 1, "mobility");
    //printf("E_Cheby is %e   and bff is %e\n",E_cheby,bff_nrm);
    //printVectorDev(brownian_force_ff,size,"FBff: ");
    //printVectorDev(brownian_force_nf,size,"FBnf: ");
    real * brownian_disp;
#ifndef SD_FF_ONLY
    cuda_safe_mem(cudaMalloc( (void**)&brownian_disp, size*sizeof(real) ));    assert(brownian_disp != NULL);
    sd_set_zero<<<32,32>>>(brownian_disp,size);
    const real one=1;
    if (myInfo_h[2]){
      cublasCall(cublasRaxpy(cublas, size, &one, brownian_force_nf, 1, brownian_force_ff, 1 ));
    }
    sd_iterative_solver(mobility, resistance, brownian_force_ff, brownian_disp, ran_prec,0, false);//(E_cheby/2+err_force/2)*1e-3
#else
    //real one=1, zero=0;
    //cuda_safe_mem(cudaMalloc( (void**)&brownian_disp, size*sizeof(real) ));    assert(brownian_disp != NULL);
    //cublasCall(cublasRgemv(cublas, CUBLAS_OP_N, size, size, &one, mobility_d, lda, brownian_force_ff, 1, &zero, brownian_disp, 1));
    brownian_disp = brownian_force_ff;
#endif
    //real tmp;
    //printVectorDev(brownian_disp, size, "brownian disp");
    //printVectorDev(disp_d, size, "det disp");
    //cublasCall(cublasRnrm2(cublas, size, brownian_disp, 1, &tmp));
    //printf("Brownian disp is: %e   ",tmp);
    //cublasCall(cublasRnrm2(cublas, size, disp_d, 1, &tmp));
    //printf("Deterministic disp is: %e   \n",tmp);
    
    cublasCall(cublasRaxpy( cublas, size, &alpha, brownian_disp, 1, disp_d, 1 ));
#ifndef SD_FF_ONLY
    cuda_safe_mem(cudaFree((void*)brownian_disp));
#endif
    cuda_safe_mem(cudaFree((void*)brownian_force_ff));
    if (myInfo_h[2]){
      cuda_safe_mem(cudaFree((void*)brownian_force_nf));
    }
    cuda_safe_mem(cudaFree((void*)gaussian));
  }// end of brownian motion
  cudaCheckError("brownian motion error");
  
  // free everything
#ifndef SD_FF_ONLY
  if (use_buckets) {
    cuda_safe_mem(cudaFree((void*)bucket_size));
    cuda_safe_mem(cudaFree((void*)bucket_num));
    cuda_safe_mem(cudaFree((void*)particle_count));
    cuda_safe_mem(cudaFree((void*)particle_to_bucket_list));
    cuda_safe_mem(cudaFree((void*)particle_sorted_list));
  }
#endif
  cudaCheckError("in mobility");
}

/// this calls an iterative solver to solve the problem: 
/// disp * (1+resistance*mobility) = mobility_d *  force_d 
/// and returnes \param disp
/// \param mobility and \param resistance are square matrizes with \param size <size> and ldd <((size+31)/32)*32>
/// \param force and \param disp are vectors of \param size <size>
real sd_iterative_solver(const matrix& mobility, const matrix& resistance, const real * force, real * disp, real rel_tol, real abs_tol, bool recalc)
{
  int size=mobility.size;
  if (abs(abs_tol) > 1e-8){
    if (warnings > 1) fprintf(stderr,"Solving for brownian forces.\n");
  } else {
    if (warnings > 1) fprintf(stderr,"Solving for interparticle forces.\n");
  }
#ifdef SD_DEBUG
  assert(!hasAnyNanDev(mobility.data,mobility.size*mobility.ldd));
  assert(!hasAnyNanDev(resistance.data,resistance.size*mobility.ldd));
  assert(!hasAnyNanDev(force,size));
  assert(!hasAnyNanDev(disp,size));
#endif
#ifdef SD_PC
  static real * mat_a = NULL;
  if (mat_a==NULL){
    cuda_safe_mem(cudaMalloc( (void**)&mat_a, lda*size*sizeof(real) ));
	recalc=true;
  }       assert(mat_a != NULL);
  if (recalc){
    sd_set_zero_matrix<<<32,192>>>(mat_a,size,lda);
  }
#endif
  real * mob_force=NULL;
  cuda_safe_mem(cudaMalloc( (void**)&mob_force, size*sizeof(real) ));       assert(mob_force !=NULL);
  real * result_checker=NULL;
  cuda_safe_mem(cudaMalloc( (void**)&result_checker, size*sizeof(real) ));  assert(result_checker !=NULL);
  // vars for cuBLAS calls
  real alpha=1;
#ifdef SD_PC
  // mat_a = (1+resistance*mobility)
  if (recalc){
    assert(!resistance.is_sparse);
    assert(!mobility.is_sparse);
    assert(mobility.wavepart==NULL);
    cublasCall(cublasRgemm(cublas,CUBLAS_OP_N,CUBLAS_OP_N, size , size ,size, &alpha, mobility.data, lda,resistance.data, lda, &beta,mat_a, lda));
    sd_add_identity_matrix<<<64,192>>>(mat_a,size,lda);
  }
#endif
  // mob_force = mobility * force
  //cublasCall(cublasRgemv( cublas, CUBLAS_OP_N, size, size, &alpha, mobility.data, lda, force, 1, &beta, mob_force, 1));
  sd_Rgemv(&alpha, mobility, force,mob_force);
#ifdef SD_DEBUG
  assert(!hasAnyNanDev(mob_force,size));
#  ifdef SD_PC
  assert(!hasAnyNanDev(mat_a,size*lda));
#  endif
#endif
  int info;
  real res;
  //printVectorDev((real *)force,size,"Kraft");
  //printMatrixDev((real *)mobility, lda, size, "mobility");
  //printMatrixDev((real* )mat_a,lda,size,"A");
  //printVectorDev(disp,6,"before");
  static int last_solv_info=4;
  if (last_solv_info>1){
    sd_set_zero<<<16,192>>>(disp,size);
  }
  //info = sd_bicgstab_solver(cublas ,size, mat_a,lda, mob_force, rel_tol, abs_tol, 10*size+100, disp, &res);
#ifdef SD_PC
  info = sd_gmres_solver(cublas ,size, mat_a,                lda, mob_force, rel_tol, abs_tol, size, disp, &res);
  //info = sd_bicgstab_solver(cublas ,size, mat_a,                lda, mob_force, rel_tol, abs_tol, size, disp, &res);
#else
  info = sd_gmres_solver(mobility, resistance, mob_force, rel_tol, abs_tol, size, disp, &res);
  //info = sd_bicgstab_solver(cublas ,size, mobility, resistance, lda, mob_force, rel_tol, abs_tol, size, disp, &res);
#endif
  last_solv_info=info;
  //printVectorDev(disp,6,"after");
  // compary to expected result
  //cuda_safe_mem(cudaMemcpy(mat_a, mat_a_bak, lda*size*sizeof(real),cudaMemcpyDeviceToDevice));
  
  if (info != 0){
    if (info == 1){
      if (warnings>1) fprintf(stderr, "Iterative solver did not fully converge ... the residuum was %6e\nWe will continue anyway ...\n",res);
    }
    else{ // info == 2 || info == 4
      if (info == 1){
	if (warnings>1) fprintf(stderr, "Iterative solver did not fully converge ... the residuum was %6e\nWe will continue anyway ...\n",res);
      }
      else if (info == 2){
	if(warnings) fprintf(stdout, "Iterative solver failed ... the residuum was %6e\nWe will continue but the results may be problematic ...\n",res);
      }
    }
    // dgetrs is not better - the contrary: results are worse ...
    /*int ipiv[size];
      magma_dgetrf_gpu( size, size,mat_a, lda, ipiv, &info);
      assert(info==0);
      magma_dgetrs_gpu('N', size, 1,
      mat_a, lda, ipiv,
      disp, size, &info);
      assert(info==0);
      // compary to expected result
      cuda_safe_mem(cudaMemcpy(mat_a, mat_a_bak, lda*size*sizeof(real),cudaMemcpyDeviceToDevice));
      cublasCall(cublasRgemv( cublas, CUBLAS_OP_N, size, size, &alpha, mat_a, lda, disp, 1, &beta, result_checker, 1));
      alpha=-1;
      cublasCall(cublasRaxpy( cublas, size, &alpha, mob_force, 1, result_checker, 1));
      alpha=1;
      cublasCall(cublasRdot( cublas, size, result_checker, 1, result_checker, 1,&res));
      if (res > 1e-1){
      fprintf(stderr, "All methods failed :(. The residuum from getrs was %e\n",res);
      //cublasCall(cublasRgemv( cublas, CUBLAS_OP_N, size, size, &alpha, mat_a, lda, disp, 1, &beta, result_checker, 1));
      //printVectorDev(mob_force, size, "mob_force");
      //printVectorDev(result_checker, size, "result_checker");
      //printVectorDev(disp, size, "disp");
      //printMatrixDev((real *)mobility,lda,size,"mobility");
      //printMatrixDev((real *)resistance,lda,size,"res");
      //printMatrixDev((real *)mat_a,lda,size,"mat_a");
      }*/
    //magma_int_t magma_dgetrs_gpu( magma_trans_t trans, magma_int_t n, magma_int_t nrhs,
    //				  double *dA, magma_int_t ldda, magma_int_t *ipiv,
    //				  double *dB, magma_int_t lddb, magma_int_t *info);
  }
#ifdef SD_DEBUG
  assert(!hasAnyNanDev(disp,size));
#endif
  //assert(info==0);
  //cuda_safe_mem(cudaFree((void*)mat_a));
  //cuda_safe_mem(cudaFree((void*)mat_a_bak));
  cuda_safe_mem(cudaFree((void*)mob_force));
  cuda_safe_mem(cudaFree((void*)result_checker));
  return res;
}

// BICGSTAB-Solver
// GMRes works better, therefore no support for BiCGStab
// implimented as given in Numerik linearer Gleichungssysteme by Prof. Dr. Andreas Meister
// this solves A*x=b
// cublas a handle for cublas
// size   the size n of the matrix
// A      the given n*n matrix = A1*A2+1 (in)
// A1     the two matrices
// A2     
// lda    the leading demension of A(1,2)
// b      the given solution vector (in)
// tol    requested tolerance of the solution
// maxit  maximum number of iterations
// x      the requested solution with an initial guess (in/out)
// returns 0 on success, else error code
/*#define sd_bicgstab_solver_exit()  cuda_safe_mem(cudaFree((void*)r0)); \
  cuda_safe_mem(cudaFree((void*)r));cuda_safe_mem(cudaFree((void*)p));	\
  cuda_safe_mem(cudaFree((void*)v));cuda_safe_mem(cudaFree((void*)t));	\
  cuda_safe_mem(cudaFree((void*)tmp));					\
  res[0] = normr;							\
  if (warnings > 1) printf("BICSTAB used %d iterations.\n", iteration);	\
  if (tolb*1.01 > normr)     { return 0;}				\
  if (tolb*100 > normr)      { return 1;}				\
  if (initnorm*1e10 > normr) { return 2;}				\
  else                       { return 4;}
#ifdef SD_PC
int sd_bicgstab_solver(cublasHandle_t cublas ,int size, const real * A, int lda,
		       real * b, real tol,real abs_tol, int maxit, real * x, real * res){
#else
int sd_bicgstab_solver(cublasHandle_t cublas ,int size, const real * A1, const real * A2, int lda,
		       real * b, real tol,real abs_tol, int maxit, real * x, real * res){
#endif
  // vector malloc
  real * r0=NULL;
  cuda_safe_mem(cudaMalloc( (void**)&r0, size*sizeof(real) ));       assert(r0 != NULL);
  real * r=NULL;
  cuda_safe_mem(cudaMalloc( (void**)&r, size*sizeof(real) ));        assert(r != NULL);
  real * p=NULL;
  cuda_safe_mem(cudaMalloc( (void**)&p, size*sizeof(real) ));        assert(p != NULL);
  real * v=NULL;
  cuda_safe_mem(cudaMalloc( (void**)&v, size*sizeof(real) ));        assert(v != NULL);
  real * t=NULL;
  cuda_safe_mem(cudaMalloc( (void**)&t, size*sizeof(real) ));        assert(t != NULL);
  real * tmp=NULL;
  cuda_safe_mem(cudaMalloc( (void**)&tmp, 2*lda*sizeof(real) ));      assert(tmp != NULL);
  //printMatrixDev(A,lda,size,"mat_a");
  //printVectorDev(b,size,"force");
#ifdef SD_PC
  cublasCall(cublasRtrsv(cublas, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, size, A, lda, b, 1)); 
#endif
  // constants
  real eps;
  if (sizeof(real) == sizeof(double)){
    eps = 1e-15;
  } else {
    eps = 1e-7;
  }
  eps = min(eps,tol*1e-2);
  // other variables
  real alpha=1;
  real beta=0;
  real tolb;
  // compute the norm of b
  real normb;
  cublasCall(cublasRnrm2( cublas, size, b, 1, &normb));
  // compute the tolerance we want to achieve
  tolb=tol*normb+abs_tol;
  // r0 = b-A*x
  alpha = -1;
#ifdef SD_PC
  cublasCall(cublasRgemv(cublas, CUBLAS_OP_N, size, size, &alpha, A, lda, x, 1, &beta, r0, 1));
  cublasCall(cublasRtrsv(cublas, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, size, A, lda, r0, 1)); 
#else
  cublasCall(cublasRgemv(cublas, CUBLAS_OP_N, size, size, &alpha, A2, lda, x, 1, &beta, tmp, 1));
  alpha=1;
  cublasCall(cublasRgemv(cublas, CUBLAS_OP_N, size, size, &alpha, A1, lda, tmp, 1, &beta, r0, 1));
  alpha=-1;
  cublasCall(cublasRaxpy(cublas, size, &alpha, x, 1, r0, 1));
#endif
  alpha = 1;
  cublasCall(cublasRaxpy(cublas, size, &alpha, b, 1, r0, 1));
  // r = r0
  cublasCall(cublasRcopy(cublas, size,r0,1,r, 1));
  // rr0 = r*r0
  real rr0;
  cublasCall(cublasRdot( cublas, size, r0, 1, r0, 1, &rr0));
  // p =r
  cublasCall(cublasRcopy(cublas, size,r0,1,p, 1));
  // normr=norm(r)
  real normr=sqrt(rr0);
  int iteration=0;
  real lastnorm=normr;
  real initnorm=normr;
  // check for conversion or max iterations
  while (iteration < maxit && normr >= tolb){
    // v=A*p
#ifdef SD_PC
    cublasCall(cublasRgemv(cublas, CUBLAS_OP_N, size, size, &alpha, A, lda, p, 1, &beta, v, 1));
    cublasCall(cublasRtrsv(cublas, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, size, A, lda, v, 1)); 
#else
    cublasCall(cublasRgemv(cublas, CUBLAS_OP_N, size, size, &alpha, A1, lda, p, 1, &beta, tmp, 1));
    cublasCall(cublasRgemv(cublas, CUBLAS_OP_N, size, size, &alpha, A1, lda, tmp, 1, &beta, v, 1));
    cublasCall(cublasRaxpy(cublas, size, &alpha, p, 1, v, 1));
#endif
    // vr0 = v*r0
    real vr0;
    cublasCall(cublasRdot( cublas, size, v, 1, r0, 1, &vr0));
    if (fabs(vr0) < eps || rr0 == 0){
      if (fabs(vr0) < eps){
	if (warnings > 1) fprintf(stderr, "BICGSTAB break-down.\n");
      }else{
	if (warnings > 1) fprintf(stderr, "BICGSTAB solution stagnates.\n");
      }
      sd_bicgstab_solver_exit();
    }
    // alpha = rr0/vr0
    real myAlpha=rr0/vr0;
    real minusMyAlpha = -myAlpha;
    // s = r - alpha v
    //cublasCall(cublasRcopy(cublas, size,r,1,s, 1));
    cublasCall(cublasRaxpy(cublas, size, &minusMyAlpha, v, 1, r, 1)); //s->r
    // t = A * s
#ifdef SD_PC
    cublasCall(cublasRgemv(cublas, CUBLAS_OP_N, size, size, &alpha, A, lda, r, 1, &beta, t, 1));// s->r
    cublasCall(cublasRtrsv(cublas, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, size, A, lda, t, 1)); 
#else
    cublasCall(cublasRgemv(cublas, CUBLAS_OP_N, size, size, &alpha, A2, lda, r, 1, &beta, tmp, 1));// s->r
    cublasCall(cublasRgemv(cublas, CUBLAS_OP_N, size, size, &alpha, A1, lda, tmp, 1, &beta, t, 1));// s->r
    cublasCall(cublasRaxpy(cublas, size, &alpha, r, 1, t, 1));
#endif
    // ts = s * t
    real ts;
    cublasCall(cublasRdot( cublas, size, t, 1, r, 1, &ts));// s->r
    // tt = t * t
    real tt;
    cublasCall(cublasRdot( cublas, size, t, 1, t, 1, &tt));
    if (abs(tt)<eps || ts == 0){
      if (warnings > 1) fprintf(stderr, "Exit: abs(tt)<eps || ts == 0\n");
      if (warnings > 1) fprintf(stderr, "BICGSTAB break-down.\n");
      sd_bicgstab_solver_exit();
    }
    // omega = ts/tt
    real myOmega=ts/tt;
    // x = x + alpha p + omega s
    cublasCall(cublasRaxpy(cublas, size, &myAlpha, p, 1, x, 1));
    cublasCall(cublasRaxpy(cublas, size, &myOmega, r, 1, x, 1));
    // copyback of s to r
    // r = s - omega t
    real minusMyOmega=-1*myOmega;
    cublasCall(cublasRaxpy(cublas, size, &minusMyOmega, t, 1, r, 1));
    //myOmega*=-1;
    // r1r0 = r * r0
    real r1r0;
    cublasCall(cublasRdot( cublas, size, r, 1, r0, 1, &r1r0));
    // beta = (alpha * r1r0 ) / (omega rr0)
    real myBeta = (myAlpha*r1r0)/(myOmega*rr0);
    if (abs(myBeta)>1/eps){
      if (warnings > 1) fprintf(stderr,"Exit: abs(myBeta)<1/eps\n");
      sd_bicgstab_solver_exit()
    }
    // p = r + beta ( p - omega v)= beta p + r - beta omega v
    cublasCall(cublasRscal(cublas, size, &myBeta, p, 1));
    cublasCall(cublasRaxpy(cublas, size, &alpha, r, 1, p, 1));
    alpha=-myBeta*myOmega;
    cublasCall(cublasRaxpy(cublas, size, &alpha, v, 1, p, 1));
    alpha=1;
    rr0=r1r0;
    cublasCall(cublasRnrm2( cublas, size, r, 1, &normr));
    if (lastnorm*sqrt(eps) > normr){ // restart
      cublasCall(cublasRcopy(cublas, size,b,1,r, 1));
      alpha=-1;beta=1; // r <-- r-A*x
#ifdef SD_PC
      cublasCall(cublasRgemv(cublas, CUBLAS_OP_N, size, size, &alpha, A, lda, x, 1, &beta, r, 1));
      cublasCall(cublasRtrsv(cublas, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, size, A, lda, r, 1)); 
#else // TODO
      alpha=1, beta=0;
      cublasCall(cublasRgemv(cublas, CUBLAS_OP_N, size, size, &alpha, A2, lda, x, 1, &beta, tmp+lda, 1));// s->r
      cublasCall(cublasRgemv(cublas, CUBLAS_OP_N, size, size, &alpha, A1, lda, tmp+lda, 1, &beta, tmp, 1));// s->r
      cublasCall(cublasRaxpy(cublas, size, &alpha, x, 1, tmp, 1));
      alpha=-1;
      cublasCall(cublasRaxpy(cublas, size, &alpha, tmp, 1, r, 1));
#endif
      alpha= 1;beta=0;
      cublasCall(cublasRdot( cublas, size, r, 1, r, 1, &rr0));
      normr=sqrt(rr0);
      lastnorm = normr;
      // r = r0
      cublasCall(cublasRcopy(cublas, size,r,1,r0, 1));
      // p =r
      cublasCall(cublasRcopy(cublas, size,r,1,p, 1));
    }
    if (iteration%50000 == 0){ // enable debugging by setting this to a lower value
      real realnorm;
      {// recalculate normr
        alpha=-1;beta=0;
#ifdef SD_PC
	cublasCall(cublasRgemv(cublas, CUBLAS_OP_N, size, size, &alpha, A, lda, x, 1, &beta, tmp, 1));
	cublasCall(cublasRtrsv(cublas, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, size, A, lda, tmp, 1)); 
#else
	cublasCall(cublasRgemv(cublas, CUBLAS_OP_N, size, size, &alpha, A2, lda, x, 1, &beta, tmp+lda, 1));// s->r
	alpha=1;
	cublasCall(cublasRgemv(cublas, CUBLAS_OP_N, size, size, &alpha, A1, lda, tmp+lda, 1, &beta, tmp, 1));// s->r
	alpha=-1;
	cublasCall(cublasRaxpy(cublas, size, &alpha, r, 1, tmp, 1));
#endif
	alpha=1;
	cublasCall(cublasRaxpy(cublas, size, &alpha, b,1,tmp, 1));
	alpha= 1;beta=0;
	cublasCall(cublasRnrm2(cublas, size, tmp, 1, &realnorm));
      }
      fprintf(stderr,"  Iteration: %6d Residuum: %12f RealResiduum: %12f\n",iteration, normr, realnorm);
    }
    iteration++;
    if (initnorm*1e10 < normr){ // somehow our solution explodes ...
      if (warnings) fprintf(stderr, "BICGSTAB did not converge, residuum exploded. Aborting.\n");
      sd_bicgstab_solver_exit();
    }
  }
  res[0]=normr;
  if (normr > tolb*1.01){
    if (warnings) fprintf(stderr, "BICGSTAB solution did not converge after %d iterations. Error was %e1 %% to high.\n",iteration,(normr/tolb-1)*100);
  }
  //fprintf(stderr, "BICGSTAB solution did converge after %d iterations.\n",iteration);
  sd_bicgstab_solver_exit();
}
*/

// GMRes-Solver
// implimented as given in Numerik linearer Gleichungssysteme by Prof. Dr. Andreas Meister
// this solves A*x=b
// cublas a handle for cublas
// size   the size n of the matrix
// A      the given n*n matrix (in)
// lda    the leading demension of A
// b      the given solution vector (in)
// tol    requested tolerance of the solution
// maxit  maximum number of iterations
// x      the requested solution with an initial guess (in/out)
// returns 0 on success, else error code
#define sd_gmres_solver_exit()  cuda_safe_mem(cudaFree((void*)r0));	\
  cuda_safe_mem(cudaFree((void*)v));cuda_safe_mem(cudaFree((void*)w));	\
  cuda_safe_mem(cudaFree((void*)tmp)); res[0] = normr;			\
  if (warnings > 1) printf("GMRes used %d iterations.\n", iteration);	\
  if (tolb*1.01 >= normr)     { return 0;}				\
  printf("tolb: %e, normr: %e, eps: %e \n",tolb, normr, eps);		\
  if (tolb*100 > normr)       { return 1;}				\
  if (initnorm*1e10 >= normr) { return 2;}				\
  else                        { return 4;}
#ifdef SD_PC
int sd_gmres_solver(int size, const real * A,int lda, const real * b, real tol,real abs_tol, int maxit, real * x, real * res){
#else
int sd_gmres_solver(const matrix & A1, const matrix & A2 , const real * b, real tol,real abs_tol, int maxit, real * x, real * res){
#endif
  int size=A1.size;
  int lda =A1.ldd;
  const int m=60;
  // vector malloc device
  real * r0=NULL;
  cuda_safe_mem(cudaMalloc( (void**)&r0,  size*sizeof(real) ));      assert(r0  != NULL);
  real * v=NULL;
  cuda_safe_mem(cudaMalloc( (void**)&v,   m*lda*sizeof(real) ));     assert(v   != NULL);
  real * w=NULL;
  cuda_safe_mem(cudaMalloc( (void**)&w,   size*sizeof(real) ));      assert(w   != NULL);
  real * tmp=NULL;
  cuda_safe_mem(cudaMalloc( (void**)&tmp, size*sizeof(real) ));      assert(tmp != NULL);
  // malloc host
  const real minusone=-1;
  const real one=1;
  real gamma[m+1];
  real alpha[m];
  real h[(m+1)*(m+1)];
  real s[m];
  real c[m];
#ifdef SD_DEBUG
  for (int i=0;i<m;i++){
    s[i]=0;
    c[i]=0;
    alpha[i]=0;
    gamma[i]=0;
  }
  gamma[m]=0;
  for (int i=0;i<m*(m+1);i++){
    h[i]=0;
  }
#endif
#ifdef SD_PC
  // we use left preconditioner, so we dont need to do triang-solv on exit
  // on a random test upper reduced the conditional number more
  cublasCall(cublasRtrsv(cublas, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, size, A, lda, b, 1)); 
#endif
  // constants
  real eps;
  if (sizeof(real) == sizeof(double)){
    eps = 1e-15;
  } else {
    eps = 1e-7;
  }
  eps = min(eps,tol*1e-2);
  // other variables
  real tolb;
  // compute the norm of b
  real normb;
  cublasCall(cublasRnrm2( cublas, size, b, 1, &normb));
  // compute the tolerance we want to achieve
  tolb=tol*normb+abs_tol;
  // r0 = b-A*x
#ifdef SD_PC
  cublasCall(cublasRgemv(cublas, CUBLAS_OP_N, size, size, &minusone, A, lda, x, 1, &zero, r0, 1));
  cublasCall(cublasRtrsv(cublas, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, size, A, lda, r0, 1)); 
#else
  sd_Rgemv(&minusone, A2, x,  tmp);
  sd_Rgemv(&one,      A1, tmp, r0);
  //cublasCall(cublasRgemv(cublas, CUBLAS_OP_N, size, size, &one,      A1, lda, tmp, 1, &zero, r0, 1));
  cublasCall(cublasRaxpy(cublas, size, &minusone, x, 1, r0, 1));
#endif
  cublasCall(cublasRaxpy(cublas, size, &one, b, 1, r0, 1));
  real normr;
  cublasCall(cublasRnrm2(cublas, size, r0, 1, &normr));
  real initnorm=normr;
  if (warnings > 1) printf("inital normr: %e\n",normr);
  int iteration=0;
  if (normr==0){
    sd_gmres_solver_exit();
  }
  if (normb==0){
    sd_set_zero<<<192,192>>>(x,size);
    normr=0;
    sd_gmres_solver_exit();
  }
  do {
    cublasCall(cublasRcopy(cublas, size, r0, 1, v  ,1));
    gamma[0]=normr;
    real igamma=1/normr;
    cublasCall(cublasRscal(cublas, size, &igamma, v, 1));
    int j;
    for (j=0;j<m && normr > tolb ;j++){
#ifdef SD_PC
      cublasCall(cublasRgemv(cublas, CUBLAS_OP_N, size, size, &one, A, lda, v+lda*j, 1, &zero, w, 1));
      cublasCall(cublasRtrsv(cublas, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, size, A, lda, w, 1)); 
#else
      sd_Rgemv(&one, A2, v+lda*j, tmp);
      sd_Rgemv(&one, A1, tmp, w);
      //cublasCall(cublasRgemv(cublas, CUBLAS_OP_N, size, size, &one, A2, lda, v+lda*j, 1, &zero, tmp, 1));
      //cublasCall(cublasRgemv(cublas, CUBLAS_OP_N, size, size, &one, A1, lda, tmp, 1, &zero, w, 1));
      cublasCall(cublasRaxpy(cublas, size, &one, v+lda*j, 1, w, 1));
#endif
      for (int i=0;i<j+1;i++){
	cublasCall(cublasRdot(cublas, size, v+i*lda, 1, w, 1, h+j*(m+1)+i));
	real tmp=-h[(m+1)*j+i];
	cublasCall(cublasRaxpy(cublas, size, &tmp , v+lda*i, 1, w, 1));
      }
      cublasCall(cublasRnrm2(cublas, size, w, 1, h+j*(m+1)+j+1));
      for (int i=0;i<j;i++){
	real hij =h[j*(m+1)+i];
	real hipj=h[j*(m+1)+i+1];
	h[j*(m+1)+i]  = c[i]*hij +s[i]*hipj;
	h[j*(m+1)+i+1]=-s[i]*hij +c[i]*hipj;
      }
      real hjj =h[j*(m+1)+j];
      real hjjp=h[j*(m+1)+j+1];
      real beta=sqrt(SQR(hjj)+SQR(hjjp));
      s[j]     =hjjp/beta;
      c[j]     =hjj /beta;
      
      
      gamma[j+1]  =-s[j]*gamma[j];
      gamma[j]    = c[j]*gamma[j];
      
      normr=abs(gamma[j+1]);
      if (j+1<m){
	// TODO: optimize: set w = v+j+1*lda
	cublasCall(cublasRcopy(cublas, size, w, 1, v+(j+1)*lda, 1));
	real tmp =1/h[j*(m+1)+j+1];
	cublasCall(cublasRscal(cublas, size, &tmp, v+(j+1)*lda, 1));
      }
      h[j*(m+1)+j]=beta;
      h[j*(m+1)+j+1]=0;
    }
    for (int i=j-1;i>=0;i--){
      real tmp=gamma[i];
      for (int k=i+1;k<j;k++){
	tmp-=h[k*(m+1)+i]*alpha[k];
      }
      alpha[i]=tmp/(h[i*(m+1)+i]);
    }
    for (int i=0;i<j;i++){
      cublasCall(cublasRaxpy(cublas, size, alpha+i, v+i*lda, 1, x, 1));
    }
    if (warnings > 1) printf("Iteration: %d, j was %d, residiuum: %e\n", iteration, j, normr );
    iteration+=j;
    {
      real newnormr;
      // r0 = b-A*x
#ifdef SD_PC
      cublasCall(cublasRgemv(cublas, CUBLAS_OP_N, size, size, &minusone, A, lda, x, 1, &zero, r0, 1));
      cublasCall(cublasRtrsv(cublas, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, size, A, lda, r0, 1));
#else
      sd_Rgemv(&minusone, A2, x,  tmp);
      sd_Rgemv(&one,      A1, tmp, r0);
      //cublasCall(cublasRgemv(cublas, CUBLAS_OP_N, size, size, &minusone, A2, lda, x,  1, &zero, tmp, 1));
      //cublasCall(cublasRgemv(cublas, CUBLAS_OP_N, size, size, &one,      A1, lda, tmp, 1, &zero, r0, 1));
      cublasCall(cublasRaxpy(cublas, size, &minusone, x, 1, r0, 1));
#endif
      cublasCall(cublasRaxpy(cublas, size, &one, b, 1, r0, 1));
      cublasCall(cublasRnrm2(cublas, size, r0, 1, &newnormr));
      if ((newnormr/normr > 1.1 || newnormr/normr < 0.9) && newnormr > sqrt(eps)){
        if (warnings) printf("Our norm changed strangely: Expected value: %e, real value: %e - eps %e- Bug in GMRes Solver?\n",normr,newnormr, eps);
      }
      normr=newnormr;
    }
  } while (iteration < maxit && normr > tolb);
  sd_gmres_solver_exit();
}

// calculates the largest and snalles eigenvalue of the matrix
// size        : size of the eigenvector / the matrix           (IN)
// mobility_d  : handle of the mobility matrix (on the device)  (IN)
// lambda_min  : smalles eigenvalue                            (OUT)
// lambda_max  : largest eigenvalue                            (OUT)
void calculate_maxmin_eigenvalues(const matrix & mobility, real * lambda_min, real * lambda_max, real tol){
  int size=mobility.size;
  int lda = ((size+31)/32)*32;
  int maxit=max(500,size);
  int IDO;
  char BMAT='I'; // standard eigenvalue problem
  char WHICH[]="LR"; // start with largest eigenvalue
  int NEV = 1; // only one eigenvalue
  // TUNING: these could be adjusted?
  real TOL=tol;
  if (sizeof(double) == sizeof(real)){
    TOL=max(1e-12,TOL);
  } else {
    TOL=max(1e-6,TOL);
  }
  // TUNING: make some tests to find good value ...
  int NCV=min(size, 6); // must be at least 3, but a bit bigger should be better ...
  int LDV=lda;
  int LWORKL=3*NCV*(NCV + 2);
  int mode=1;
  int IPARAM[11] = {1,0,maxit,1,0,0,mode,0,0,0,0};
  int IPNTR[14];
  int INFO=0;
  real RESID[size];
  real V[LDV*NCV];
  real WORKD[3*size];
  real WORKL[LWORKL];
  real * vec_in_d;
  real * vec_out_d;
  cuda_safe_mem(cudaMalloc((void**)&vec_in_d , lda*sizeof(real)));
  cuda_safe_mem(cudaMalloc((void**)&vec_out_d, lda*sizeof(real)));
  for (int minmax=0;minmax<2;minmax++){
    IDO=0;
    if (minmax){
      sprintf(WHICH,"SR");
      INFO=0;
      IPARAM[2]=maxit;
      TOL=sqrt(tol);
    }
    while (IDO != 99){
      //dnaupd_(&IDO,&BMAT,&N,WHICH,&NEV,&TOL,RESID.memptr(),&NCV,V.memptr(),&LDV,IPARAM,IPNTR,WORKD,WORKL,&LWORKL,&INFONAUP);
      rnaupd(&IDO,&BMAT,&size, WHICH, &NEV, &TOL, RESID, &NCV, V, &LDV, IPARAM, IPNTR, WORKD, WORKL, &LWORKL, &INFO);
      switch (IDO){
      case 1:
	cuda_safe_mem(cudaMemcpy(vec_in_d,WORKD+IPNTR[0]-1,size*sizeof(real),cudaMemcpyHostToDevice));
	{
	  const real one=1;
	  //cublasCall(cublasRgemv( cublas, CUBLAS_OP_N, size, size, &alpha, mobility_d, lda, vec_in_d, 1, &zero, vec_out_d, 1));
	  sd_Rgemv(&one, mobility, vec_in_d, vec_out_d);
	}
	cuda_safe_mem(cudaMemcpy(WORKD+IPNTR[1]-1,vec_out_d,size*sizeof(real),cudaMemcpyDeviceToHost));
	break;
      case -1:
      case 2:
      case 3:
      case 4:
	fprintf(stderr,"Error in %s l. %d: unexpected work from rnaupd: %d: Not Implemented!\n",__FILE__,__LINE__,IDO);
	break;
      case 99: //we are done
	break;
      default:
	fprintf(stderr,"Error in %s l. %d: unexpected work from rnaupd: %d: Not Understood!\n",__FILE__,__LINE__,IDO);
	break;
      }
    }
    if (warnings > 1) fprintf(stderr,"calculating eigenvalue needed %d iterations and %d gemv operations (tolerance is %e, EW is %e).\n"
	    ,IPARAM[2], IPARAM[8], TOL,WORKL[IPNTR[5]-1]);
    if (INFO){
      if (INFO == 1 && warnings > 1){
	fprintf(stderr,"Maximum iterations taken to find Eigenvalues.\n");
      } else if( INFO > 1 && warnings){
	fprintf(stderr,"Unexpected return value in %s l. %d from rnaupd_: %d (debug info: %d %e)\n",__FILE__,__LINE__,INFO, minmax, V[0]);
      }
      /*
	c  INFO    Integer.  (INPUT/OUTPUT)
	c          If INFO .EQ. 0, a randomly initial residual vector is used.
	c          If INFO .NE. 0, RESID contains the initial residual vector,
	c                          possibly from a previous run.
	c          Error flag on output.
	c          =  0: Normal exit.
	c          =  1: Maximum number of iterations taken.
	c                All possible eigenvalues of OP has been found. IPARAM(5)
	c                returns the number of wanted converged Ritz values.
	c          =  2: No longer an informational error. Deprecated starting
	c                with release 2 of ARPACK.
	c          =  3: No shifts could be applied during a cycle of the
	c                Implicitly restarted Arnoldi iteration. One possibility
	c                is to increase the size of NCV relative to NEV.
	c                See remark 4 below.
	c          = -1: N must be positive.
	c          = -2: NEV must be positive.
	c          = -3: NCV-NEV >= 2 and less than or equal to N.
	c          = -4: The maximum number of Arnoldi update iteration
	c                must be greater than zero.
	c          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
	c          = -6: BMAT must be one of 'I' or 'G'.
	c          = -7: Length of private work array is not sufficient.
	c          = -8: Error return from LAPACK eigenvalue calculation;
	c          = -9: Starting vector is zero.
	c          = -10: IPARAM(7) must be 1,2,3,4.
	c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatable.
	c          = -12: IPARAM(1) must be equal to 0 or 1.
	c          = -9999: Could not build an Arnoldi factorization.
	c                   IPARAM(5) returns the size of the current Arnoldi
	c                   factorization.
      */
    }
    if (WORKL[IPNTR[5]-1]<0 && TOL > 1e-3){
      minmax--;
      TOL=min(tol,1e-4);
      tol=min(tol*tol, 1e-8);
      INFO=1;
      IPARAM[2]=maxit*100;
    }
    if (minmax){ // make them a bit larger/smaller to be sure that we are in the interval of interrest ...
      *lambda_min=WORKL[IPNTR[5]-1]*(1+TOL);
    } else {
      *lambda_max=WORKL[IPNTR[5]-1]*(1-TOL);
    }
  }
  /* FORTRAN Comments of the dnaupd
     c          IPNTR(6): pointer to the real part of the ritz value array     
     c                    RITZR in WORKL.                                          
     c          IPNTR(7): pointer to the imaginary part of the ritz value array    
     c                    RITZI in WORKL.                                          
     c          IPNTR(8): pointer to the Ritz estimates in array WORKL associated
     c                    with the Ritz values located in RITZR and RITZI in WORK
  */
  if (warnings > 1) {
    printf("Range of Farfield-Mobility: [%e, %e]\n", *lambda_min, *lambda_max);
  }
  cuda_safe_mem(cudaFree((void*)vec_in_d));
  cuda_safe_mem(cudaFree((void*)vec_out_d));
}

/// \brief caclulate the chebyshev coefficents for the squareroot of the matrix
/// It needs the largest and smalles eigenvalue which have to be computed before.
/// Using chebyshev-gausquadrature \link https://en.wikipedia.org/wiki/Chebyshev%E2%80%93Gauss_quadrature \endlink.
/// 
/// \param lambda_min   : the lower boundery
/// \param lambda_max   : the upper boundery of the interval
/// \param tol          : the given tollerance which should be achieved
/// \param coefficents  : the pointer where the data will be stored
real calculate_chebyshev_coefficents(real lambda_min, real lambda_max, real tol,real ** coefficents){
  const int steps=1024*128; // with 1024 it should fit in L1
  unsigned int N=1024; // guess the number of coefficents we need, if more are needed -> realocate
  if (*coefficents==NULL){
    *coefficents = (real *)Utils::malloc(N*sizeof(real));
  } else {
    *coefficents = (real *)Utils::realloc(*coefficents,N*sizeof(real));
  }
  real * current_polynome=NULL;
  real * last_polynome=NULL;
  last_polynome    = (real *) Utils::malloc(steps * sizeof(real));
  current_polynome = (real *) Utils::malloc(steps * sizeof(real));
  real x[steps];
  real weight_and_func[steps];
  real lambda_m=lambda_max-lambda_min;
  real lambda_p=lambda_max+lambda_min;
  //fprintf(stderr,"lambda_min: %e lambda_max: %e  ", lambda_min, lambda_max);
  //fprintf(stderr,"lambda_minus: %e lambda_plusminus: %e\n", lambda_m, lambda_pm);
  // init
  real fac=2./steps;
  //fprintf(stderr,"fac: %e\n", fac);
  double ai=0;
  double a1=0;
  for (int i=0;i<steps;i++){
    last_polynome[i]=1; //T0(x)=1
    //real x = -1.+(i*2.+1.)/steps;
    x[i]=cos((2.*i+1)/2./steps*M_PI);
    current_polynome[i]= x[i];  //T1(x)=x
    #ifdef SD_FF_ONLY
    weight_and_func[i]=fac*sqrt(x[i]*lambda_m/2.+lambda_p/2.);
    #else
    weight_and_func[i]=fac*1./sqrt(x[i]*lambda_m/2.+lambda_p/2.);// /sqrt(1-x*x);// this could be big, but should not be inf
    #endif
    ai+=weight_and_func[i];
    a1+=weight_and_func[i]*x[i];
  }
  real error;
  int loop=0;
  //fprintf(stderr,"%s l. %d: a[%d]: %e\n",__FILE__,__LINE__,0,ai);
  //fprintf(stderr,"%s l. %d: a[%d]: %e\n",__FILE__,__LINE__,1,a1);
  (*coefficents)[loop]=ai/2.;
  loop++;
  (*coefficents)[loop]=a1;
  //double sumfacmax=lambda_max;
  //double totalsum=abs(ai)+abs(a1)*sumfacmax;
  double totalsum=abs(ai)+abs(a1);
  //sumfacmax*=lambda_max;
  const int miniloop=10;
  do{
    error=0;
    do {
      std::swap(last_polynome,current_polynome);
      //{	real * tmp       = last_polynome; last_polynome    = current_polynome; current_polynome = tmp; }
      //printf("addresses: 0x%08x  0x%08x ",last_polynome,current_polynome);
      ai=0;
      for (int i=0;i<steps;i++){
	//real x = -1.+(i*2.+1.)/steps;
	current_polynome[i]=-1.*current_polynome[i]+((real)2)*x[i]*last_polynome[i];
	//printf ("%e %e %e\n", x, weight_and_func[i], current_polynome[i]);
	ai+=current_polynome[i]*weight_and_func[i];
      }
      //fprintf(stderr,"%s l. %d: a[%d]: %e\n",__FILE__,__LINE__,loop,ai);
      //printf("\n\n");
      //sumfacmax*=lambda_max;
      error+=abs(ai);//*sumfacmax;
      totalsum+=abs(ai);//*sumfacmax;
      loop++;
      (*coefficents)[loop]=ai;
    } while (loop%miniloop);
    if (loop+miniloop > N){
      N*=2;
      *coefficents=(real *)Utils::realloc(*coefficents,N*sizeof(real));
    }
  } while ((error > tol*totalsum/10 || loop < 20 ) && loop < sqrt(steps));
  if (loop >=steps/10 -1 ){
    fprintf(stderr,"to few steps to get sufficent results in %s l. %d\n",__FILE__,__LINE__);
  }
  error=0;
  while (error < tol*totalsum){ // approximate error
    loop--;
    error+=abs((*coefficents)[loop]);//*sumfacmax;
    //sumfacmax/=lambda_max;
  }
  //fprintf(stderr,"sum: %e   error: %e",totalsum,error);
  loop++;
  free(last_polynome);
  free(current_polynome);
  if (warnings > 1) fprintf(stderr,"chebyshev needs %d gemv\n",loop);
  return loop;
}
/// \brief Compute the brownian force resulting from the farfield
/// To compute the farfield contribution of the brownian motion the
/// squareroot of the inverse of the mobility is required.
/// Instead of calculating it directly, its action on the vector
/// containing the randomnumbers is computed via the chebyshev-polynomials.
/// In the case of FF_ONLY this computes directly the displacement!
real sd_compute_brownian_force_farfield(const matrix & mobility, const real * gaussian_ff,
					real tol, real * brownian_force_ff ){

  cudaCheckError("Unknown Error");
  int size=mobility.size;
  static real * cheby_coefficents=NULL;
  static int N_chebyshev;
  static bool recalc_ew = true;
  static real lambda_min, lambda_max;
  static int count_recalc=0;
  static int count_total=0;
  count_total++;
  cudaCheckError("Unknown Error");
  if (recalc_ew){
    count_recalc++;
    calculate_maxmin_eigenvalues(mobility,&lambda_min, &lambda_max, tol);
    //printf("lambda min: %e, lambda_max: %e",lambda_min, lambda_max);
    N_chebyshev = calculate_chebyshev_coefficents(lambda_min, lambda_max,tol*tol,&cheby_coefficents);
    //printf("\nWe need %d iteration for Chebyshev.\n", N_chebyshev);
    recalc_ew=false;
    if (warnings>1){printf("recalc ratio: %6f\n",(((double)count_recalc)/((double)count_total)));}
  }
  if (lambda_min < 0){
    std::cerr << "Mobility has negative eigenvalues!\n" << std::endl;
    std::cout << "Mobility has negative eigenvalues!\n" << std::endl;
    if (mobility.wavespace != NULL){
      //printMatrixDev(mobility);
    }
    if (warnings > 2){
      //printMatrixDev(mobility);
    }
    //printMatrixDev(mobility.data,mobility.ldd,mobility.size,"Mobility has negative eigenvalues!\n");
    errexit();
  }
  real * chebyshev_vec_curr, * chebyshev_vec_last, * chebyshev_vec_next;
  sd_set_zero<<<64,192>>>(brownian_force_ff,size);
  cudaThreadSynchronize();
  cudaCheckError("set zero");
  chebyshev_vec_curr=NULL;
  cuda_safe_mem(cudaMalloc( (void**)&chebyshev_vec_curr, size*sizeof(real) ));    assert(chebyshev_vec_curr != NULL);
#ifdef SD_FF_ONLY
  const real one  =1;
  real gMg;
  //cublasCall(cublasRgemv(cublas, CUBLAS_OP_N, size, size, &one, mobility_d, lda, gaussian_ff, 1, &zero,  chebyshev_vec_curr, 1));
  sd_Rgemv(&one, mobility, gaussian_ff, chebyshev_vec_curr);
  cublasCall(cublasRdot(cublas, size, chebyshev_vec_curr, 1, gaussian_ff, 1, &gMg));
#else
  real gaussian_ff_norm;
  cublasCall(cublasRnrm2(cublas, size, gaussian_ff, 1, &gaussian_ff_norm));
#endif
  cublasCall(cublasRcopy( cublas, size, gaussian_ff, 1, chebyshev_vec_curr, 1));
  cublasCall(cublasRaxpy( cublas, size, cheby_coefficents+0, chebyshev_vec_curr, 1, brownian_force_ff, 1 ));
  //chebyshev_vec_last=chebyshev_vec_curr;
  chebyshev_vec_last=NULL;
  cuda_safe_mem(cudaMalloc( (void**)&chebyshev_vec_last, size*sizeof(real) ));    assert(chebyshev_vec_last != NULL);
  chebyshev_vec_next=NULL;
  cuda_safe_mem(cudaMalloc( (void**)&chebyshev_vec_next, size*sizeof(real) ));    assert(chebyshev_vec_next != NULL);
  //sd_set_zero<<<64,192>>>(chebyshev_vec_????,size);
  real lambda_minus=lambda_max-lambda_min;
  real alpha=2./lambda_minus;
  sd_Rgemv(&alpha, mobility,gaussian_ff, chebyshev_vec_next);
  //cublasCall(cublasRgemv( cublas, CUBLAS_OP_N, size, size, &alpha, mobility_d, lda, gaussian_ff, 1, &zero, chebyshev_vec_next , 1));
  alpha=-(lambda_min+lambda_max)/lambda_minus;
  cublasCall(cublasRaxpy( cublas, size, &alpha, chebyshev_vec_curr, 1, chebyshev_vec_next, 1 ));
  std::swap(chebyshev_vec_curr,chebyshev_vec_next);
  std::swap(chebyshev_vec_last,chebyshev_vec_next);
  cublasCall(cublasRaxpy( cublas, size, cheby_coefficents+1, chebyshev_vec_curr, 1, brownian_force_ff, 1 ));
  for (int i=2;i<=N_chebyshev;i++){
    alpha=4./lambda_minus;
    sd_Rgemv(&alpha, mobility, chebyshev_vec_curr, chebyshev_vec_next);
    //cublasCall(cublasRgemv( cublas, CUBLAS_OP_N, size, size, &alpha, mobility_d, lda, chebyshev_vec_curr, 1, &zero, chebyshev_vec_next , 1));
    alpha=-2*(lambda_min+lambda_max)/lambda_minus;
    cublasCall(cublasRaxpy( cublas, size, &alpha, chebyshev_vec_curr, 1, chebyshev_vec_next, 1 ));
    alpha=-1;
    cublasCall(cublasRaxpy( cublas, size, &alpha, chebyshev_vec_last, 1, chebyshev_vec_next, 1 ));
    std::swap(chebyshev_vec_curr,chebyshev_vec_next);
    std::swap(chebyshev_vec_last,chebyshev_vec_next);
    cublasCall(cublasRaxpy( cublas, size, cheby_coefficents+i, chebyshev_vec_curr, 1, brownian_force_ff, 1 ));
  }
#ifdef SD_DEBUG
  assert(isSymmetricDev(mobility.data,mobility.ldd,mobility.size));
#endif
  // errorcheck of chebyshev polynomial
#ifdef SD_FF_ONLY
  real dispdisp;
  cublasCall(cublasRdot(cublas, size, brownian_force_ff, 1,brownian_force_ff, 1, &dispdisp));
  real E_cheby = sqrt(abs(dispdisp-gMg)/gMg);
#else
  real zMz;
  alpha = 1;
  sd_Rgemv(&alpha, mobility, brownian_force_ff, chebyshev_vec_last);
  //cublasCall(cublasRgemv(cublas, CUBLAS_OP_N, size, size, &alpha, mobility_d, lda, brownian_force_ff, 1, &zero,  chebyshev_vec_last, 1));
  cublasCall(cublasRdot(cublas, size, chebyshev_vec_last, 1, brownian_force_ff, 1, &zMz));
  real E_cheby = sqrt(abs(zMz-gaussian_ff_norm*gaussian_ff_norm))/gaussian_ff_norm;
#endif
  if (E_cheby > tol){
    recalc_ew=true;
  }
  cuda_safe_mem(cudaFree((void*)chebyshev_vec_last));
  cuda_safe_mem(cudaFree((void*)chebyshev_vec_curr));
  cuda_safe_mem(cudaFree((void*)chebyshev_vec_next));
  return E_cheby;
}

void sd_Rgemv(const real * factor, const matrix & A, const real * in, real * out){
  A.assert_all();
  if(A.is_sparse){
    int threads=A.size;
    const int tpB = numThreadsPerBlock;
    int blocks = ((threads+tpB-1)/tpB);
    sd_multiply_sparse_matrix_vector<<<blocks, tpB>>>(A.size, *factor, A.data, A.ldd, A.ldd_short, A.col_idx, A.row_l, in, out );
    cudaThreadSynchronize(); 
    real tmp;
    cublasCall(cublasRnrm2(cublas, A.size, out,1, &tmp));
    cudaCheckError("Sparse Matrix Vector Multiplication");
  }
  else{
    assert(A.ldd);
    const real zero=0;
    #ifdef SD_USE_FLOAT
    cublasCall(cublasSgemv(cublas, CUBLAS_OP_N, A.size, A.size, factor, A.data, A.ldd, in,  1, &zero, out, 1));
    #else
    cublasCall(cublasDgemv(cublas, CUBLAS_OP_N, A.size, A.size, factor, A.data, A.ldd, in,  1, &zero, out, 1));
    #endif
  }
  if(A.wavespace != NULL){
    /*real insize;
    cublasCall(cublasRnrm2(cublas, A.size, in, 1, &insize));
    real before;
    cublasCall(cublasRnrm2(cublas, A.size, out, 1, &before));*/
    wavepart  & wave=* A.wavespace;
    int buf_size = A.wavespace->max*3;
    real * sin_sum=NULL;
    cuda_safe_mem(cudaMalloc( (void**)& sin_sum, buf_size*sizeof(real)));
    real max;
    sd_nrm1<<<1,numThreadsPerBlock>>>(A.size,in,sin_sum);
    cuda_safe_mem(cudaMemcpy(&max, sin_sum, sizeof(real), cudaMemcpyDeviceToHost));
    sd_set_zero<<<64,192>>>(sin_sum,buf_size);
    //fprintf(stderr,"buf_size is %d\n",buf_size);
    real * cos_sum=NULL;
    cuda_safe_mem(cudaMalloc( (void**)& cos_sum, buf_size*sizeof(real)));
    sd_set_zero<<<64,192>>>(cos_sum,buf_size);
    cudaThreadSynchronize();
    cudaCheckError("");
    // thread number has to be a power of two!
    //sd_wavepart_sum<<<wave.num,64,64*3*sizeof(real)>>>(in, wave.vecs, wave.num, wave.matrices, wave.sines, wave.cosines, \
    //fprintf(stderr,"calling sd_wavepart_sum<<<%d,%d,%d>>>\n",wave.num,8,8*3*sizeof(real));
    //                         ,------ has to be power of two and at least 8 (and there might be still a bug if larger than 32)
    //                        \|/
    //                         v   
    sd_wavepart_sum<<<wave.num,32,32*3*sizeof(real)>>>(in, wave.vecs, wave.num, wave.matrices, wave.sines, wave.cosines, \
						       A.ldd_short, A.size/3, sin_sum, cos_sum, max);
    cudaThreadSynchronize();
    cudaCheckError("");
    int tpB = 32;
    int blocks = (A.ldd_short)/tpB;
    sd_wavepart_assemble<<<blocks,tpB,tpB*3*sizeof(real)>>>(wave.num, wave.sines, wave.cosines, sin_sum, cos_sum,
    							    A.ldd_short, out, max, A.size/3, *factor);
    //sd_wavepart_assemble_block<<<A.size/3,tpB,tpB*3*sizeof(real)>>>(wave.num, wave.sines, wave.cosines, sin_sum, cos_sum,
    //							              A.ldd_short, out, max, A.size/3, *factor);
    cudaThreadSynchronize();
    cudaCheckError("");
    /*#ifdef SD_DEBUG
    if (A.dense){
      printVectorDev(sin_sum,buf_size,"insin_sum");
      printVectorDev(cos_sum,buf_size,"incos_sum");
      real * correct;
      cuda_safe_mem(cudaMalloc( (void**)&correct,    A.ldd*sizeof(real) ));
      const real zero=0;
      #ifdef SD_USE_FLOAT
      cublasCall(cublasSgemv(cublas, CUBLAS_OP_N, A.size, A.size, factor, A.dense, A.ldd, in,  1, &zero, correct, 1));
      #else
      cublasCall(cublasDgemv(cublas, CUBLAS_OP_N, A.size, A.size, factor, A.dense, A.ldd, in,  1, &zero, correct, 1));
      #endif
      const real minusone = -1;
      cublasCall(cublasRaxpy(cublas, A.size, &minusone, out, 1 ,correct, 1));
      real erg;
      cublasCall(cublasRnrm2(cublas, A.size, correct, 1 , &erg));
      if (erg < 1e-3){
	std::cout << "s";
	printVectorDev(in     ,A.size,"\nin, worked  ");
	printVectorDev(correct,A.size,"correct     ");
	printVectorDev(out    ,A.size,"out         ");
      } else {
	A.print();
	cublasCall(cublasRaxpy(cublas, A.size, &minusone, out, 1 ,correct, 1));
	printVectorDev(in     ,A.size,"in, failed  ");
	printVectorDev(correct,A.size,"correct     ");
	printVectorDev(out    ,A.size,"out         ");
	
      }
      cuda_safe_mem(cudaFree((void*) correct));
    }
    if (warnings > 20){
      printVectorDev(sin_sum,buf_size,"sin_sum");
      printVectorDev(cos_sum,buf_size,"cos_sum");
    }
    #endif*/
    cuda_safe_mem(cudaFree((void*) sin_sum));
    cuda_safe_mem(cudaFree((void*) cos_sum));
  }
}


void sd_compute_mobility_matrix_wave_cpu(int kmax, real kmax2, matrix & mobility, real a,const real * L_h, real xi, real selfmobility, real ew_prec){
#ifdef SD_DEBUG
  //kmax=1;
  //kmax2=kmax*kmax;
#endif
  cudaCheckError("");
  real fak=selfmobility;
  for(int k=0;k<3;k++){
    fak/=L_h[k];
  }
  //std::cerr <<fak<<std::endl;
  int cki[3];
  real k[3];
  real ki[3];
  real xi2=1/xi/xi;
  real kvol=1;
  for (int i =0 ; i < 3; i++){
    ki[i]=2*M_PI/L_h[i];
    kvol*=ki[i];
  }
  // count how many wave vectors we need
  int vc = 0;
  for (cki[0]=-kmax;cki[0]<kmax+1;cki[0]+=1){
    for (cki[1]=-kmax;cki[1]<kmax+1;cki[1]+=1){
      for (cki[2]=-kmax;cki[2]<kmax+1;cki[2]+=1){
	int cki2=cki[0]*cki[0]+ cki[1]*cki[1] + cki[2]*cki[2];
	if (cki2 <= kmax2 && cki2 > 0 ){
	  vc++;
	}
      }
    }
  }
  if (!vc){ // nothing to do
    return;
  }
  if (mobility.wavespace == NULL){
    mobility.wavespace = new wavepart(mobility.ldd_short);
    //mobility.wavespace = (wavepart *) Utils::malloc(sizeof(wavepart));
    //mobility.wavespace->wavepart();
  }
  mobility.wavespace->num      = vc;
  int max = vc;
  if (max !=  mobility.wavespace->num){
    mobility.wavespace->num=max;
  }
  assert(max);
  if (max > mobility.wavespace->max){
    mobility.wavespace->max = ((max+31)/32)*32;
    max=mobility.wavespace->max;
    //fprintf(stderr,"there are %d wavespace contribution taking %d places",vc, max);
    cuda_safe_mem(cudaMalloc( (void**)& (mobility.wavespace->vecs)     , max*3*sizeof(real) ));
    cuda_safe_mem(cudaMalloc( (void**)& (mobility.wavespace->matrices) , max*6*sizeof(real) ));
    cuda_safe_mem(cudaMalloc( (void**)& (mobility.wavespace->sines)    , max*mobility.ldd_short*sizeof(real) ));
    cuda_safe_mem(cudaMalloc( (void**)& (mobility.wavespace->cosines)  , max*mobility.ldd_short*sizeof(real) ));
    cudaCheckError("");
    sd_set_zero<<<640,64>>>(mobility.wavespace->sines   , max*mobility.ldd_short);
    sd_set_zero<<<640,64>>>(mobility.wavespace->cosines , max*mobility.ldd_short);
    cudaCheckError("");
  }
  int ldd = mobility.wavespace->max;
  real vecs_h[3 * ldd];
  real matrices_h[6 * ldd];
  vc=0;
  for (cki[0]=-kmax;cki[0]<kmax+1;cki[0]+=1){
    k[0]=cki[0]*ki[0];
    for (cki[1]=-kmax;cki[1]<kmax+1;cki[1]+=1){
      k[1]=cki[1]*ki[1];
      for (cki[2]=-kmax;cki[2]<kmax+1;cki[2]+=1){
	int cki2=cki[0]*cki[0]+ cki[1]*cki[1] + cki[2]*cki[2];
	if (cki2 <= kmax2 && cki2 > 0){
	  k[2]=cki[2]*ki[2];
	  real k2=k[0]*k[0]+ k[1]*k[1] + k[2]*k[2];
	  for (int i = 0;i<3;i++){
	    vecs_h[vc*3+i]=k[i];
	  }
	  real scal = (a-1./3.*a*a*a*k2)*(1+0.25*xi2*k2+0.125*xi2*xi2*k2*k2)*6.*M_PI/k2*exp(-0.25*xi2*k2)*fak;
	  //if (cki[0]==-kmax){
	  //  printf("scal is %e, exp is %e\n",scal/fak,exp(-0.25*xi2*k2));
	  //}
	  //if (scal/kvol > 1e-6*ew_prec){
	    for (int i = 0 ; i < 3 ; i++ ){
	      matrices_h[vc*6+i] = (1-(k[i]*k[i]/k2))*scal;
	    }
	    matrices_h[vc*6+3] = (-k[0]*k[1]/k2)*scal;
	    matrices_h[vc*6+4] = (-k[0]*k[2]/k2)*scal;
	    matrices_h[vc*6+5] = (-k[1]*k[2]/k2)*scal;
	    /*if (scal/fak > 1e-5){
	      printf("scal is %e, k is %e",scal/fak,sqrt(k2));
	      printf("\tk=[%4.0f, %4.0f, %4.0f]  ",k[0]*L_h[0]/2/M_PI,k[1]*L_h[1]/2/M_PI,k[2]*L_h[2]/2/M_PI);
	      printf("  m=[%12.4e, %12.4e, %12.4e,%12.4e, %12.4e, %12.4e]  \n",matrices_h[vc*6],matrices_h[vc*6+1],matrices_h[vc*6+2],matrices_h[vc*6+3],matrices_h[vc*6+4],matrices_h[vc*6+5]);
	      }*/
	    vc++;
	    //}
	}
      }
    }
  }
  mobility.wavespace->num      = vc;
  {
    static int printed=0;
    if ((printed%reprint)==1){
      printf("\nwe need %d k vectors\n",vc);
    }
    printed++;
  }
  cudaCheckError("");
  //fprintf(stderr,"h: 0x%x d: 0x%x, size: %d",vecs_h,mobility.wavespace->vecs,  ldd*3*sizeof(real));
  cuda_safe_mem(cudaMemcpy(mobility.wavespace->vecs,    vecs_h,    ldd*3*sizeof(real),cudaMemcpyHostToDevice));
  cuda_safe_mem(cudaMemcpy(mobility.wavespace->matrices,matrices_h,ldd*6*sizeof(real),cudaMemcpyHostToDevice));
  //mobility.wavespace->hash_vecs();
}

void _cuda_check_error(char *file, unsigned int line){
  cudaError_t err=cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf(stderr, "Error found during memory operation. Possibly however from an failed operation before. %s:%u.\n", file, line);
    printf("CUDA error: %s\n", cudaGetErrorString(err));
    if ( err == cudaErrorInvalidValue )
      fprintf(stderr, "You may have tried to allocate zero memory before %s:%u.\n", file, line);
    exit(EXIT_FAILURE);
  }
}



int sd_find_max(int * data, int length){
  int tlength=length;
  int * tdata=data;
  if (tlength >= 512){
    int bs=32;
    if (tlength > 4*1024){
      bs=128;
    } else if (tlength > 1024){
      bs=64;
    }
    tlength = (length+bs-1)/bs;
    cuda_safe_mem(cudaMalloc((void **) &tdata,tlength*sizeof(int)));
    sd_nrm_inf<<<tlength,bs,bs*sizeof(int)>>>(length,data,tdata);
  }
  //if (length < 64){ // copy to host and compare
  int tmp[tlength];
  cuda_safe_mem(cudaMemcpy(tmp,tdata,tlength*sizeof(int),cudaMemcpyDeviceToHost));
  int result=tmp[0];
  for (int i=1;i<tlength;i++){
    if (tmp[i] > result){
      result=tmp[i];
    }
  }
  if (length>=512) {
    cuda_safe_mem(cudaFree((void*)tdata));
  }
  return result;
}



#endif
