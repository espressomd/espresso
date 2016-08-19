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
  DipolarBarnesHut(EspressoSystemInterface &s) {

	k = coulomb.Dbjerrum;

	if(!s.requestFGpu())
      std::cerr << "DipolarDirectSum needs access to forces on GPU!" << std::endl;

    if(!s.requestRGpu())
      std::cerr << "DipolarDirectSum needs access to positions on GPU!" << std::endl;

    if(!s.requestDipGpu())
      std::cerr << "DipolarDirectSum needs access to dipoles on GPU!" << std::endl;

	/*fillConstantPointers(s.rxGpuBegin(), s.ryGpuBegin(), s.rzGpuBegin(),
			s.dipxGpuBegin(), s.dipyGpuBegin(), s.dipzGpuBegin(),
			s.npart_gpu(), s.bhnnodes(), s.arrl(), s.boxl(), s.massGpuBegin());
	initBH(s.blocksGpu());*/

  };

  virtual ~DipolarBarnesHut() {} // TODO: any memory cleanup?

  void computeForces(EspressoSystemInterface &s) {
    dds_float box[3];
    int per[3];
    for (int i=0;i<3;i++)
    {
     box[i]=s.box()[i];
     per[i] = (PERIODIC(i));
    }

	//cudaThreadSynchronize();

	buildBoxBH(s.blocksGpu());
	//cudaThreadSynchronize();

	buildTreeBH(s.blocksGpu());
	//cudaThreadSynchronize();

	summarizeBH(s.blocksGpu());
	//cudaThreadSynchronize();

	sortBH(s.blocksGpu());
	//cudaThreadSynchronize();

	forceBH(s.blocksGpu(),k, s.fGpuBegin(),s.torqueGpuBegin(),box,per);
	//cudaThreadSynchronize();
  };
  void computeEnergy(EspressoSystemInterface &s) {
    dds_float box[3];
    int per[3];
    for (int i=0;i<3;i++)
    {
     box[i]=s.box()[i];
     per[i] = (PERIODIC(i));
    }
    //cudaThreadSynchronize();

    buildBoxBH(s.blocksGpu());
    //cudaThreadSynchronize();

    buildTreeBH(s.blocksGpu());
    //cudaThreadSynchronize();

    summarizeBH(s.blocksGpu());
    //cudaThreadSynchronize();

    sortBH(s.blocksGpu());
    //cudaThreadSynchronize();

    energyBH(s.blocksGpu(),k,box,per,(&(((CUDA_energy*)s.eGpu())->dipolar)));
    //cudaThreadSynchronize();
 };

protected:
  float k;
};

void activate_dipolar_barnes_hut();
void deactivate_dipolar_barnes_hut();

extern DipolarBarnesHut *dipolarBarnesHut;

#endif

#endif
