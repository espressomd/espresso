#include "config.hpp"

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <cublas_v2.h>
#include <curand.h>

#include "integrate_sd_cuda.cuh"
#include "integrate_sd_cuda_debug.cuh"
#include "integrate_sd_cuda_kernel.cuh"
//#include "errorhandling.hpp"
#include "cuda_utils.hpp"


/* *************************************************************************************************************** *
 * ********************************************    DEBUGING & TESTs   ******************************************** *
 * *************************************************************************************************************** */

 
int sd_test(int size, int type){
  assert(size%3 == 0);
  int lda=((size+31)/32)*32;
  matrix mat={NULL,NULL,NULL};
  mat.ldd=lda;
  mat.size=size;
  real * x1=NULL, * x2=NULL, * b =NULL;
  cuda_safe_mem(cudaMalloc((void**) &mat.data   ,size*lda*sizeof(real)));
  cuda_safe_mem(cudaMalloc((void**) &mat.col_idx,size*lda*sizeof(int)));
  cuda_safe_mem(cudaMalloc((void**) &mat.row_l  ,size*sizeof(int)));
  cuda_safe_mem(cudaMalloc((void**) &x1         ,size*sizeof(real)));
  cuda_safe_mem(cudaMalloc((void**) &x2         ,size*sizeof(real)));
  cuda_safe_mem(cudaMalloc((void**) &b          ,size*sizeof(real))); 
  sd_set_zero_matrix<<<32,32>>>(mat.data,size);
  sd_set_zero<<<32,32>>>(b,size);
  for (int i=0;i<size/3;i++){
    sd_set_int<<<32,32>>>(mat.col_idx+i*mat.ldd,lda,i);
  }
  sd_set_int<<<32,32>>>(mat.row_l,mat.ldd, mat.size/3);
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


void printVectorDev( real * data, int m, const char * msg){
  std::cout << msg;
  printVectorDev(data,m);
}
void printVectorDev( real * data, int m){
  real host[m];
  assert(m>0);
  cuda_safe_mem(cudaMemcpy( host, data, m*sizeof(*data), cudaMemcpyDeviceToHost ));
  for (int j=0; j < m;j++){
    //std::cout << ",\t" << std::setw(20) << std::setprecision(15) << host[j];
    printf(",\t%20.12e",host[j]);
  }
  std::cout << "\n";
}


void printVectorHost( real * host, int m, const char * msg){
  std::cout << msg;
  printVectorHost(host,m);
}

void printVectorHost( real * host, int m){
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
  counters=(unsigned long long *)malloc(sizeof(unsigned long long)*numCounters);
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



real getAbsMax(real * data, int length){
  real max=0;
  for (int i=0;i<length;i++){
    if (abs(data[i])>max){
      max=abs(data[i]);
    }
  }
  return max;
}


real getAbsMaxDev(real * data_d, int length){
  real host[length];
  cuda_safe_mem(cudaMemcpy( host, data_d, length*sizeof(*data_d), cudaMemcpyDeviceToHost ));
  return getAbsMax(host,length);
}


int countMinusDev(real * data_d, int length){
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

bool  isSymmetricDev(const real * data, int lda, int size){
  real host[lda*size];
  cuda_safe_mem(cudaMemcpy( host, data, size*lda*sizeof(*data), cudaMemcpyDeviceToHost ));
  for (int i=0;i<size;i++){
    for (int j=0;j<size;j++){
      if (fabs(host[i*lda+j]-host[j*lda+i]) > 1e-8){
	return false;
      }
    }
  }
  return true;
}
