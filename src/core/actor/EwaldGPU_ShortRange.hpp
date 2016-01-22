#ifndef EWALDGPUFORCESHORTRANGE_HPP_
#define EWALDGPUFORCESHORTRANGE_HPP_

#ifdef EWALD_GPU
#include "actor/EwaldGPU.hpp"

#define SQR(A) ((A)*(A))

//Add energy
inline double ewaldgpu_coulomb_pair_energy(double chgfac, double *d,double dist2,double dist)
{
  if (dist < ewaldgpu_params.rcut)
  {
     return coulomb.prefactor*chgfac*erfc(ewaldgpu_params.alpha*dist)/dist;
  }
  return 0.0;
}

//Add forces
inline void add_ewald_gpu_coulomb_pair_force(Particle *p1, Particle *p2, double d[3], double dist, double force[3])
{
  const double rcut=ewaldgpu_params.rcut;
  if(dist < rcut)
  {
    int j;
    double fac;

    fac=coulomb.prefactor * p1->p.q * p2->p.q * (  2*ewaldgpu_params.alpha*wupii*exp(-SQR(ewaldgpu_params.alpha *dist)) + erfc(ewaldgpu_params.alpha*dist)/dist )/SQR(dist);

    for(j=0;j<3;j++)
      force[j] += fac * d[j];

    ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: EWALD_GPU   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist,fac));
    ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: EWALD_GPU   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist,fac));
  }
}

#endif

#endif /* EWALDGPUFORCESHORTRANGE_HPP_ */
