#include <iostream>
#include <iomanip>
#include <stdio.h>

#include "integrate_sd_cuda_debug.cuh"
//#include "errorhandling.hpp"
#include "cuda_utils.hpp"

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
