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

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <cublas_v2.h>
#include <curand.h>

#include "integrate_sd_cuda.hpp"
#include "integrate_sd_cuda_debug.hpp"
#include "integrate_sd_cuda_kernel.hpp"
//#include "errorhandling.hpp"
#include "cuda_utils.hpp"


/* *************************************************************************************************************** *
 * ********************************************    DEBUGING & TESTs   ******************************************** *
 * *************************************************************************************************************** */
#ifdef SD
 
int sd_test(int size, int type){
  assert(size%3 == 0);
  int lda=((size+31)/32)*32;
  matrix mat;
  mat.ldd=lda;
  mat.size=size;
  real * x1=NULL, * x2=NULL, * b =NULL;
  cuda_safe_mem(cudaMalloc((void**) &mat.data   ,size*lda*sizeof(real)));
  cuda_safe_mem(cudaMalloc((void**) &mat.col_idx,size*lda*sizeof(int)));
  cuda_safe_mem(cudaMalloc((void**) &mat.row_l  ,size*sizeof(int)));
  cuda_safe_mem(cudaMalloc((void**) &x1         ,size*sizeof(real)));
  cuda_safe_mem(cudaMalloc((void**) &x2         ,size*sizeof(real)));
  cuda_safe_mem(cudaMalloc((void**) &b          ,size*sizeof(real))); 
  sd_set_zero_matrix<<<32,32>>>(mat.data,mat.size,mat.ldd);
  sd_set_zero<<<32,32>>>(b,size);
  for (int i=0;i<size/3;i++){
    sd_set_value<<<32,32>>>(mat.col_idx+i*mat.ldd,lda,i);
  }
  sd_set_value<<<32,32>>>(mat.row_l,mat.ldd, mat.size/3);
  //cublasHandle_t cublas=NULL; // use global
  cublasCall(cublasCreate(&cublas));
  curandGenerator_t generator = NULL;
  curandCall(curandCreateGenerator(&generator, CURAND_RNG_PSEUDO_DEFAULT));
  curandCall(curandSetPseudoRandomGeneratorSeed(generator, (unsigned long long)('E'+'S'+'P'+'R'+'e'+'s'+'S'+'o')));
  curandCall(curandSetGeneratorOrdering( generator, CURAND_ORDERING_PSEUDO_BEST));
  curandCall(curandSetGeneratorOffset( generator, 0));
  curandCall(curandGenerateNormalReal( generator, b,  size, 0, sqrt(2.)));
  curandCall(curandGenerateNormalReal( generator, x1, size, 0, sqrt(2.)));
  curandCall(curandGenerateNormalReal( generator, x2, size, 0, sqrt(2.)));
  if (type==0){
    ;
  }
  else{
    curandCall(curandGenerateNormalReal( generator, mat.data, mat.size*mat.ldd, 0, sqrt(2.)));
  }
  mat.is_sparse=false;
  const real one=1;
  const real minusone=-1;
  sd_Rgemv(&one, mat, b,x1);
  mat.is_sparse=true;
  sd_Rgemv(&one, mat, b,x2);
  real tmp;
  cublasCall(cublasRnrm2(cublas, size, x1, 1, &tmp));
  std::cout << "x1 norm: " << tmp;
  cublasCall(cublasRnrm2(cublas, size, x2, 1, &tmp));
  std::cout << "   x2 norm: " << tmp;
  cublasCall(cublasRaxpy(cublas, size, &minusone, x1,1,x2,1));
  cublasCall(cublasRnrm2(cublas, size, x2, 1, &tmp));
  std::cout << "   diff norm: " << tmp;
  return 0;
}


void printVectorDev( const real * data, const int m, const char * msg){
  std::cout << msg;
  printVectorDev(data,m);
}
void printVectorDev( const real * data, const int m){
  real host[m];
  assert(m>0);
  cuda_safe_mem(cudaMemcpy( host, data, m*sizeof(*data), cudaMemcpyDeviceToHost ));
  int count = 1;
  char komma[2]={0,0};
  std::cout << " = [";
  for (int j=0; j < m;j++){
    bool print=true;
    //std::cout << ",\t" << std::setw(20) << std::setprecision(15) << host[j];
    if (j < m-1){
      if (host[j]==host[j+1]){
	count++;
	print=false;
      }
    }
    if (print){
      if (count > 3){
	printf("%s\t%12.4e <repeats %d times>",komma,host[j],count);
	komma[0]=',';
	count = 1;
      } else {
	for (; count > 0 ; count--){
	  printf("%s\t%12.4e",komma,host[j]);
	  komma[0]=',';
	}
	count = 1;
      }
    }
  }
  std::cout << "];\n";
}


void printVectorHost( const real * host, const int m, const char * msg){
  std::cout << msg;
  printVectorHost(host,m);
}

void printVectorHost( const real * host, const int m){
  for (int j=0; j < m;j++){
    std::cout << ",\t" << std::setw(20) << std::setprecision(15) << host[j];
  }
  std::cout << "\n";
}

void printPosDev( real * data, int m, const char * msg){
  std::cout << msg;
  real host[m*3];
  cuda_safe_mem(cudaMemcpy( host, data,3*m*sizeof(*data), cudaMemcpyDeviceToHost ));
  for (int l=0; l<m;l++){
    printf("%5d : [ %e, %e, %e]\n",l,host[3*l],host[3*l+1],host[3*l+2]);
  }
}

/*void _cudaCheckError(const char *msg, const char * file, const int line)
{
  cudaError_t err = cudaGetLastError();
  if( cudaSuccess != err)
    {
      std::cerr <<  "Cuda error:" <<  msg << ": " <<  cudaGetErrorString( err) << " in "<<file << "l. "<<line<<"\n";
      errexit();
    }
    }*/

 
// check whether there was any cuda error so far.
// do not use this function directly but use the macro cudaCheckError(const char *msg);
// which requires only the first paramter
// PARAMTERS:
// msg   : the message which should be printed in case of an error
// file  : the file in which the function is called
// line  : the line in which the function is called
void _cudaCheckError(const char *msg, const char * file, const int line)
{
  cudaError_t err = cudaGetLastError();
  if( cudaSuccess != err)
    {
      std::cerr <<  "Cuda error:" <<  msg << ": '" <<  cudaGetErrorString( err) << "' in "<<file << " l. "<<line<<"\n";
      errexit();
    }
}



myTimer::myTimer(unsigned int _numCounters):numCounters(_numCounters),counters(NULL){
  counters=(unsigned long long *)Utils::malloc(sizeof(unsigned long long)*numCounters);
  assert(counters!=NULL);
  for (unsigned int i=0;i<numCounters;i++){
    counters[i]=0;
  }
  clock_gettime(CLOCK_MONOTONIC_RAW,&last);
}

myTimer::~myTimer(){
  free(counters);
}

void myTimer::add(unsigned int _counter){
  struct timespec current;
  clock_gettime(CLOCK_MONOTONIC_RAW,&current);
  // use usec
  counters[_counter]+=(current.tv_sec-last.tv_sec)*1e6+(current.tv_nsec-last.tv_nsec+500)*1e-3;
  clock_gettime(CLOCK_MONOTONIC_RAW,&last);
}

void myTimer::print(unsigned int _counter,const char * msg){
  std::cout << msg << " took "<< counters[_counter]*1e-3<<" ms\n";
}

void myTimer::printAll(const char *msg){
  std::cout << "\n" << msg;
  for (unsigned int i=0;i<numCounters;i++){
    std::cout << "\t" << counters[i]*1e-3;
  }
  std::cout << "\n";
}



real getAbsMax(const real * data, int length){
  real max=0;
  for (int i=0;i<length;i++){
    if (abs(data[i])>max){
      max=abs(data[i]);
    }
  }
  return max;
}


real getAbsMaxDev(const real * data_d, int length){
  real host[length];
  cuda_safe_mem(cudaMemcpy( host, data_d, length*sizeof(*data_d), cudaMemcpyDeviceToHost ));
  return getAbsMax(host,length);
}


int countMinusDev(const real * data_d, int length){
  real host[length];
  cuda_safe_mem(cudaMemcpy( host, data_d, length*sizeof(*data_d), cudaMemcpyDeviceToHost ));
  int c=0;
  for (int i=0;i<length;i++){
    if (host[i]<0)
      c++;
  }
  return c;
}

bool hasAnyNanDev(const real * data, int length){
  real host[length];
  cuda_safe_mem(cudaMemcpy( host, data, length*sizeof(*data), cudaMemcpyDeviceToHost ));
  for (int c=0;c<length;c++){
    if (isnan(host[c])){
      return true;
    }
  }
  return false;
}

bool isSymmetricDev(const matrix & mat){
  return isSymmetricDev(mat.data,mat.ldd,mat.size);
}

bool  isSymmetricDev(const real * data, int lda, int size){
  #ifdef SD_DEBUG
  real host[lda*size];
  cuda_safe_mem(cudaMemcpy( host, data, size*lda*sizeof(real), cudaMemcpyDeviceToHost ));
  for (int i=0;i<size;i++){
    for (int j=0;j<size;j++){
      if (fabs(host[i*lda+j]-host[j*lda+i])/(fabs(host[i*lda+j])+fabs(host[j*lda+i])) > 1e-8){
	fprintf(stderr,"not symmetric: elements (%d,%d) and (%d,%d)\n",i,j,j,i);
	return false;
      }
    }
  }
  #endif
  return true;
}



void printMatrixDev( const matrix & mat){
  /*if (mat.wavespace){
    printf("wavespace:: max entries: %d, actual entries: %d\n",mat.wavespace->max,mat.wavespace->num);
    printVectorDev(mat.wavespace->matrices,mat.wavespace->max*6, "matrices");
    printVectorDev(mat.wavespace->vecs    ,mat.wavespace->max*3, "vecs");
  }*/
  mat.printWavespace();
  printMatrixDev(mat.data,mat.ldd,mat.size);
}


void  matrix::printWavespace() const{
  if (wavespace){
    wavespace->print(size/3,ldd_short);
  }else {
    printf("No wavespace available!\n");
  }
}
void matrix::print() const{
  printMatrixDev(*this);
}

void wavepart::print(int N, int ldd) const{
  printVectorDev(sines,max*ldd);
  printVectorDev(cosines,max*ldd);
  assert(max > 0);
  real vec_h[max*3];
  real mat_h[max*6];
  real sin_h[max*ldd];
  real cos_h[max*ldd];
  //real box_l[]={8,8,8};
  extern double box_l[3];
  cuda_safe_mem(cudaMemcpy( vec_h, vecs,     max*3*sizeof(real),   cudaMemcpyDeviceToHost ));
  cuda_safe_mem(cudaMemcpy( mat_h, matrices, max*6*sizeof(real),   cudaMemcpyDeviceToHost ));
  cuda_safe_mem(cudaMemcpy( sin_h, sines,    max*ldd*sizeof(real), cudaMemcpyDeviceToHost ));
  cuda_safe_mem(cudaMemcpy( cos_h, cosines,  max*ldd*sizeof(real), cudaMemcpyDeviceToHost ));
  printf("printing wavespace:\n");
  for (int i=0; i<num; i++){
    printf("\tk=[%4.0f, %4.0f, %4.0f]  ",vec_h[i*3]*box_l[0]/2/M_PI,vec_h[i*3+1]*box_l[1]/2/M_PI,vec_h[i*3+2]*box_l[2]/2/M_PI);
    printf("  m=[%12.4e, %12.4e, %12.4e,%12.4e, %12.4e, %12.4e]  ",mat_h[i*6],mat_h[i*6+1],mat_h[i*6+2],mat_h[i*6+3],mat_h[i*6+4],mat_h[i*6+5]);
    printf("\tsines:  \t");
    for (int p=0;p<N;p++){
      printf("%12.4e,  ",sin_h[p+ldd*i]);
    }
    printf("\tcosines:\t");
    for (int p=0;p<N;p++){
      printf("%12.4e,  ",cos_h[p+ldd*i]);
    }
    printf("\n");
  }
}

#endif /* SD */
