#include "config.hpp"

#ifdef BARNES_HUT

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

//This needs to be done in the .cu file too!!!!!
typedef float dds_float;

class DipolarBarnesHut : public Actor {
public:
  //DipolarBarnesHut(SystemInterface &s, float epssq, float itolsq)
  DipolarBarnesHut(SystemInterface &s)
  {

	k = coulomb.Dprefactor;
	//m_epssq = epssq;
	//m_itolsq = itolsq;
	//setBHPrecision(m_epssq,m_itolsq);
	if(!s.requestFGpu())
      std::cerr << "DipolarBarnesHut needs access to forces on GPU!" << std::endl;

    if(!s.requestRGpu())
      std::cerr << "DipolarBarnesHut needs access to positions on GPU!" << std::endl;

    if(!s.requestDipGpu())
      std::cerr << "DipolarBarnesHut needs access to dipoles on GPU!" << std::endl;

    //std::cout << "Trace DipolarBarnesHut init" << std::endl;

	/*fillConstantPointers(s.rxGpuBegin(), s.ryGpuBegin(), s.rzGpuBegin(),
			s.dipxGpuBegin(), s.dipyGpuBegin(), s.dipzGpuBegin(),
			s.npart_gpu(), s.bhnnodes(), s.arrl(), s.boxl(), s.massGpuBegin());
	initBH(s.blocksGpu());*/

  };

//virtual ~DipolarBarnesHut() {}

  void computeForces(SystemInterface &s) {
    dds_float box[3];
    int per[3];

    //std::cout << "Trace computeForces 1" << std::endl;

    for (int i=0;i<3;i++)
    {
     box[i]=s.box()[i];
     per[i] = (PERIODIC(i));
    }

	//cudaThreadSynchronize();

    //printf("\n blocks = %d",s.blocksGpu());
	buildBoxBH(s.blocksGpu());
	//cudaThreadSynchronize();

	buildTreeBH(s.blocksGpu());
	//cudaThreadSynchronize();

	summarizeBH(s.blocksGpu());
	//cudaThreadSynchronize();

	sortBH(s.blocksGpu());
	//cudaThreadSynchronize();

	//std::cout << "Trace computeForces 2" << std::endl;
	//printf("Trace computeForces 2");
	forceBH(s.blocksGpu(),k,s.fGpuBegin(),s.torqueGpuBegin(),box,per);
	//cudaThreadSynchronize();
  };
  void computeEnergy(SystemInterface &s) {
    dds_float box[3];
    int per[3];

    //std::cout << "Trace computeEnergy 1" << std::endl;

    for (int i=0;i<3;i++)
    {
     box[i]=s.box()[i];
     per[i] = (PERIODIC(i));
    }
    //cudaThreadSynchronize();

    //printf("\n blocks = %d",s.blocksGpu());
    buildBoxBH(s.blocksGpu());
    //cudaThreadSynchronize();

    buildTreeBH(s.blocksGpu());
    //cudaThreadSynchronize();

    summarizeBH(s.blocksGpu());
    //cudaThreadSynchronize();

    sortBH(s.blocksGpu());
    //cudaThreadSynchronize();

    //std::cout << "Trace computeEnergy 2" << std::endl;
    energyBH(s.blocksGpu(),k,box,per,(&(((CUDA_energy*)s.eGpu())->dipolar)));
    //cudaThreadSynchronize();
 };

protected:
  float k;
  float m_epssq;
  float m_itolsq;
};

//void activate_dipolar_barnes_hut(float epssq, float itolsq);
void activate_dipolar_barnes_hut();
void deactivate_dipolar_barnes_hut();

extern DipolarBarnesHut *dipolarBarnesHut;

#endif

#endif
