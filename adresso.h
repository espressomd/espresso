// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.

#ifndef ADRESSO_H
#define ADRESSO_H
/** \file adresso.h
    This is the place for adaptive resolution scheme (adress)
    Implementation of adresso.h
    <b>Responsible:</b>
    <a href="mailto:junghans@mpip-mainz.mpg.de">Axel</a>
*/

#include <tcl.h>
#include "particle_data.h"
#include "virtual_sites.h"

/** \name Exported Variables */
/************************************************************/
/*@{*/
extern double adress_vars[7];
/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/
/** Implements the Tcl command \ref tcl_adress. This allows for seetings for adress
*/
int adress_tcl(ClientData data, Tcl_Interp *interp, int argc, char **argv);

#ifdef ADRESS
/** Calc adress weight function of a vector
    @param x[3] vector
    @return weight of the vector
*/
double adress_wf_vector(double x[3]);

/** Calc adress weight function of a particle
    @param x[3] vector
    @return weight of the particle
*/
MDINLINE double adress_wf_particle(Particle *p){
   if (p==NULL) return 0.0;
   if (ifParticleIsVirtual(p)){
      return p->p.adress_weight;
   }
   else{
      return adress_wf_particle(get_mol_com_particle(p));
   }
}

/** Update adress weight of all particles
*/
void adress_update_weights();

MDINLINE double adress_non_bonded_force_weight(Particle *p1,Particle *p2){
  double adress_weight_1,adress_weight_2,force_weight;
  int virtual_1,virtual_2;

  //NOTE this is in order of probability to appear
  adress_weight_1=adress_wf_particle(p1);
  virtual_1=ifParticleIsVirtual(p1);

   //if particles 1 is ex, but in the cg regime
  if ( (adress_weight_1<ROUND_ERROR_PREC) && (virtual_1==0) ) return 0.0;

  adress_weight_2=adress_wf_particle(p2);
  virtual_2=ifParticleIsVirtual(p2);

  //if particles 2 is ex, but in the cg regime
  if ( (adress_weight_2<ROUND_ERROR_PREC) && (virtual_2==0) ) return 0.0;

  //mixed case is captured by cg-cg interation
  if ((virtual_1+virtual_2)==1) return 0.0;

  force_weight=adress_weight_1*adress_weight_2;

  //both are cg
  if ((virtual_1+virtual_2)==2) {
     //both are in ex regime
     if (force_weight>1-ROUND_ERROR_PREC) return 0.0;
     force_weight=1-force_weight;
  }
  //both are ex -> force_weight is already set
  //if ((virtual_1+virtual_2)==0) force_weight=force_weight;

  return force_weight;
}

MDINLINE double adress_bonded_force_weight(Particle *p1){
  double adress_weight_1=adress_wf_particle(p1);
  double force_weight;
  //NOTE only ex particles have bonded interations and bonded interactions are only inside a molecule

  //particle is cg
  if (adress_weight_1<ROUND_ERROR_PREC) return 0.0;

  //both
  force_weight=adress_weight_1*adress_weight_1;
  return force_weight;
}
#endif
/*@}*/

#endif
