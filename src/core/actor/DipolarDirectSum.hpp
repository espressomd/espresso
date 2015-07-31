#include "config.hpp"

#ifdef DIPOLAR_DIRECT_SUM

#include "SystemInterface.hpp"
#include <iostream>
#include "Actor.hpp"
#include "DipolarDirectSum_cuda.hpp"
#include "grid.hpp" 
#include "cuda_interface.hpp"

#ifndef ACTOR_DIPOLARDIRECTSUM_HPP
#define ACTOR_DIPOLARDIRECTSUM_HPP

//This needs to be done in the .cu file too!!!!!
typedef float dds_float;


void DipolarDirectSum_kernel_wrapper_energy(dds_float k, int n, float *pos, float *dip, dds_float box_l[3],int periodic[3],float* E); 
void DipolarDirectSum_kernel_wrapper_force(dds_float k, int n, float *pos, float *dip, float* f, float* torque, dds_float box_l[3],int periodic[3]); 



class DipolarDirectSum : public Actor {
public:
  DipolarDirectSum(SystemInterface &s) {
    if(!s.requestFGpu())
      std::cerr << "DipolarDirectSum needs access to forces on GPU!" << std::endl;

    if(!s.requestRGpu())
      std::cerr << "DipolarDirectSum needs access to positions on GPU!" << std::endl;
    
    if(!s.requestDipGpu())
      std::cerr << "DipolarDirectSum needs access to dipoles on GPU!" << std::endl;

  }; 
  void computeForces(SystemInterface &s) {
    dds_float box[3];
    int per[3];
    for (int i=0;i<3;i++)
    {
     box[i]=s.box()[i];
     per[i] = (PERIODIC(i));
    }
    DipolarDirectSum_kernel_wrapper_force(k,s.npart_gpu(),
					 s.rGpuBegin(), s.dipGpuBegin(), s.fGpuBegin(),s.torqueGpuBegin(),box,per);
  };
  void computeEnergy(SystemInterface &s) {
    dds_float box[3];
    int per[3];
    for (int i=0;i<3;i++)
    {
     box[i]=s.box()[i];
     per[i] = (PERIODIC(i));
    }
    DipolarDirectSum_kernel_wrapper_energy(k,s.npart_gpu(),
					 s.rGpuBegin(), s.dipGpuBegin(), box,per,(&(((CUDA_energy*)s.eGpu())->dipolar)));
 };
protected:
  float k;
};

void activate_dipolar_direct_sum_gpu();
void deactivate_dipolar_direct_sum_gpu();

extern DipolarDirectSum *dipolarDirectSum;

#endif

#endif
