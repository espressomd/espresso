#include "config.hpp"

#ifdef DIPOLAR_BARNES_HUT

#include "SystemInterface.hpp"
#include "EspressoSystemInterface.hpp"
#include <iostream>
#include "Actor.hpp"
#include "DipolarBarnesHut_cuda.cuh"
#include "grid.hpp" 
#include "cuda_interface.hpp"
#include "interaction_data.hpp"

#ifndef ACTOR_DIPOLARBARNESHUT_HPP
#define ACTOR_DIPOLARBARNESHUT_HPP

//This needs to be done in the .cu file too
typedef float dds_float;

class DipolarBarnesHut : public Actor {
public:
  DipolarBarnesHut(SystemInterface &s, float epssq, float itolsq)
  {
	k = coulomb.Dprefactor;
	m_epssq = epssq;
	m_itolsq = itolsq;
	setBHPrecision(m_epssq,m_itolsq);
	if(!s.requestFGpu())
      std::cerr << "DipolarBarnesHut needs access to forces on GPU!" << std::endl;

    if(!s.requestRGpu())
      std::cerr << "DipolarBarnesHut needs access to positions on GPU!" << std::endl;

    if(!s.requestDipGpu())
      std::cerr << "DipolarBarnesHut needs access to dipoles on GPU!" << std::endl;
  };

  void computeForces(SystemInterface &s) {
    dds_float box[3];
    int per[3];

    for (int i=0;i<3;i++)
    {
     box[i]=s.box()[i];
     per[i] = (PERIODIC(i));
    }

    fillConstantPointers(s.rGpuBegin(), s.dipGpuBegin(),
                s.npart_gpu(), s.bhnnodes(), s.BHarrl(), s.BHboxl(), s.massGpuBegin());
    initBHgpu(s.blocksGpu());

	buildBoxBH(s.blocksGpu());
	buildTreeBH(s.blocksGpu());
	summarizeBH(s.blocksGpu());
	sortBH(s.blocksGpu());
	forceBH(s.blocksGpu(),k,s.fGpuBegin(),s.torqueGpuBegin(),box,per);
  };
  void computeEnergy(SystemInterface &s) {
    dds_float box[3];
    int per[3];

    for (int i=0;i<3;i++)
    {
     box[i]=s.box()[i];
     per[i] = (PERIODIC(i));
    }

    fillConstantPointers(s.rGpuBegin(), s.dipGpuBegin(),
                    s.npart_gpu(), s.bhnnodes(), s.BHarrl(), s.BHboxl(), s.massGpuBegin());
    initBHgpu(s.blocksGpu());
    buildBoxBH(s.blocksGpu());
    buildTreeBH(s.blocksGpu());
    summarizeBH(s.blocksGpu());
    sortBH(s.blocksGpu());
    energyBH(s.blocksGpu(),k,box,per,(&(((CUDA_energy*)s.eGpu())->dipolar)));
 };

protected:
  float k;
  float m_epssq;
  float m_itolsq;
};

void activate_dipolar_barnes_hut(float epssq, float itolsq);
void deactivate_dipolar_barnes_hut();

extern DipolarBarnesHut *dipolarBarnesHut;

#endif

#endif
