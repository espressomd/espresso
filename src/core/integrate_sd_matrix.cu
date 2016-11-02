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

#include "config.hpp"

#ifdef SD

#include <stdio.h>
#include <assert.h>
#include "cuda_runtime.h"

#include "integrate_sd.hpp"
#include "integrate_sd_cuda.hpp"
#include "integrate_sd_matrix.hpp"

#include <openssl/md5.h>

/// constructor for wave part
/// _ldd_short : ldd_short of the matrix
wavepart::wavepart(int _ldd_short):ldd_short(_ldd_short){
  vecs          = NULL;
  matrices      = NULL;
  cosines       = NULL;
  sines         = NULL;
  num           = 0;
  max           = 0;
#ifdef SD_DEBUG
  vecs_hash     = NULL;
  matrices_hash = NULL;
  cosines_hash  = NULL;
  sines_hash    = NULL;
#endif
}

/// destructor of the wave part
wavepart::~wavepart(){
  if (vecs){
    cuda_safe_mem(cudaFree((void*)vecs));
    vecs=NULL;
  }
  if (matrices){
    cuda_safe_mem(cudaFree((void*)matrices));
    matrices=NULL;
  }
  if (cosines){
    cuda_safe_mem(cudaFree((void*)cosines));
    cosines=NULL;
  }
  if (sines){
    cuda_safe_mem(cudaFree((void*)sines));
    sines=NULL;
  }
#ifdef SD_DEBUG
  if (vecs_hash){
    free(vecs_hash);
    vecs_hash=NULL;
  }
  if (matrices_hash){
    free(matrices_hash);
    matrices_hash=NULL;
  }
  if (sines_hash){
    free(sines_hash);
    sines_hash=NULL;
  }
  if (cosines_hash){
    free(cosines_hash);
    cosines_hash=NULL;
  }
#endif
  max=0;
}


/// free the wave part of the matrix
void matrix::_free_wavespace(){
  if (wavespace){
    //wavespace->~wavepart();
    //free((void*)wavespace);
    delete wavespace;
    wavespace=NULL;
  }
}

/// constructor of the matrix
matrix::matrix(){
  data=NULL;
  col_idx=NULL;
  row_l=NULL;
  wavespace=NULL;
  is_sparse=false;
  size=0;
  ldd=0;
  ldd_short=0;
#ifdef SD_DEBUG
  dense=NULL;
  data_hash=NULL;
  dense_hash=NULL;
#endif
}

/// destructor of the matrix
matrix::~matrix(){
  if (data){
    cuda_safe_mem(cudaFree((void*)data));
    data=NULL;
  }
  if (col_idx){
    cuda_safe_mem(cudaFree((void*)col_idx));
    col_idx=NULL;
  }
  if (row_l){
    cuda_safe_mem(cudaFree((void*)row_l));
    row_l=NULL;
  }
#ifdef SD_DEBUG
  if (dense){
    cuda_safe_mem(cudaFree((void*)dense));
    dense=NULL;
  }
  if (data_hash){
    free(data_hash);
    data_hash=NULL;
  }
  if (dense_hash){
    free(dense_hash);
    dense_hash=NULL;
  }
#endif
  _free_wavespace();
}


/// print some basic properties of the matrix, if SD_DEBUG is enabled
void matrix::printStats() const{
#ifdef SD_DEBUG
  printf("addr:%p\tsize:%d\tldd:%d\tsparse:%d",this,size,ldd,is_sparse);
#endif
}

void getAndHash(const real * data, const int size, unsigned char result[MD5_DIGEST_LENGTH]){
#ifdef SD_DEBUG
  real host[size];
  cuda_safe_mem(cudaMemcpy( host, data, size*sizeof(real), cudaMemcpyDeviceToHost ));
  MD5((unsigned char*) host, size, result);
#endif
}

void matrix::hash_data(){
#ifdef SD_DEBUG
  if (data_hash == NULL){
    data_hash=(unsigned char *) Utils::malloc(MD5_DIGEST_LENGTH);
  }
  getAndHash(data, size*ldd, data_hash);
#endif
}

void matrix::hash_dense(){
#ifdef SD_DEBUG
  if (dense_hash == NULL){
    dense_hash=(unsigned char *) Utils::malloc(MD5_DIGEST_LENGTH);
  }
  getAndHash(dense, size*ldd, dense_hash);
#endif
}

void matrix::assert_all() const{
#ifdef SD_DEBUG
  unsigned char tmp[MD5_DIGEST_LENGTH];
  if (data_hash){
    getAndHash(data, size*ldd, tmp);
    for (int i=0; i < MD5_DIGEST_LENGTH; i++){
      assert (tmp[i]==data_hash[i]);
      printf("%02x",tmp[i]);
    }
  }
  if (dense_hash){
    getAndHash(dense, size*ldd, tmp);
    for (int i=0; i < MD5_DIGEST_LENGTH; i++){
      assert (tmp[i]==dense_hash[i]);
      printf("%02x",tmp[i]);
    }
  }
  if (wavespace){
    wavespace->assert_all();
  }
#endif
}


void wavepart::hash_vecs(){
#ifdef SD_DEBUG
  fprintf(stderr,"going to allocate %d bytes, vecs_hash is %p\n",MD5_DIGEST_LENGTH, vecs_hash);
  if (vecs_hash == NULL){
    vecs_hash=(unsigned char *) Utils::malloc(MD5_DIGEST_LENGTH);
    fprintf(stderr,"allocated %d bytes\n",MD5_DIGEST_LENGTH);
  }
  getAndHash(vecs, num*3 , vecs_hash);
#endif
}
void wavepart::hash_matrices(){
#ifdef SD_DEBUG
  if (matrices_hash == NULL){
    matrices_hash=(unsigned char *) Utils::malloc(MD5_DIGEST_LENGTH);
  }
  getAndHash(matrices, num*6, matrices_hash);
#endif
}
void wavepart::hash_cosines(){
#ifdef SD_DEBUG
  if (cosines_hash == NULL){
    cosines_hash=(unsigned char *) Utils::malloc(MD5_DIGEST_LENGTH);
  }
  getAndHash(cosines, num*ldd_short, cosines_hash);
#endif
}
void wavepart::hash_sines(){
#ifdef SD_DEBUG
  if (sines_hash == NULL){
    sines_hash=(unsigned char *) Utils::malloc(MD5_DIGEST_LENGTH);
  }
  getAndHash(sines, num*ldd_short, sines_hash);
#endif
}

void wavepart::assert_all() const{
#ifdef SD_DEBUG
  unsigned char tmp[MD5_DIGEST_LENGTH];
  if (vecs_hash){
    getAndHash(vecs, 3*num, tmp);
    for (int i=0; i < MD5_DIGEST_LENGTH; i++){
      assert (tmp[i]==vecs_hash[i]);
      printf("%02x",tmp[i]);
    }
  }
  if (matrices_hash){
    getAndHash(matrices, 6*num, tmp);
    for (int i=0; i < MD5_DIGEST_LENGTH; i++){
      assert (tmp[i]==matrices_hash[i]);
      printf("%02x",tmp[i]);
    }
  }
  if (sines_hash){
    getAndHash(sines, num*ldd_short, tmp);
    for (int i=0; i < MD5_DIGEST_LENGTH; i++){
      assert (tmp[i]==sines_hash[i]);
      printf("%02x",tmp[i]);
    }
  }
  if (cosines_hash){
    getAndHash(cosines, num*ldd_short, tmp);
    for (int i=0; i < MD5_DIGEST_LENGTH; i++){
      assert (tmp[i]==cosines_hash[i]);
      printf("%02x",tmp[i]);
    }
  }
#endif
}

#endif /* SD */

