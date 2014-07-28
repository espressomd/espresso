/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
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
#ifndef REACTION_FIELD_H
#define REACTION_FIELD_H
/** \file reaction_field.hpp
 *  Routines to calculate the Reaction Field Energy or/and force 
 *  for a particle pair.
 *  M. Neumann, J. Chem. Phys 82, 5663 (1985)
 *  \ref forces.cpp
 *
*/

#include "utils.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "mol_cut.hpp"

#ifdef ELECTROSTATICS

/** Structure to hold Reaction Field Parameters. */
typedef struct {
  /** ionic strength . */
  double kappa;
  /** epsilon1 (continuum dielectric constant inside) . */
  double epsilon1;
  /** epsilon2 (continuum dielectric constant outside) . */
  double epsilon2;
  /** Cutoff for Reaction Field interaction. */
  double r_cut;
  /** B important prefactor . */
  double B;
} Reaction_field_params;

/** Structure containing the Reaction Field parameters. */
extern Reaction_field_params rf_params;

/** \name Functions */
/************************************************************/
/*@{*/

///
int rf_set_params(double kappa,double epsilon1,double epsilon2, double r_cut);

inline void add_rf_coulomb_pair_force_no_cutoff(Particle *p1, Particle *p2, double d[3], double dist, double force[3])
{
   int j;
   double fac;
   fac = 1.0 / (dist*dist*dist)  +  rf_params.B / (rf_params.r_cut*rf_params.r_cut*rf_params.r_cut);
   fac *= coulomb.prefactor * p1->p.q * p2->p.q;
   
   for (j=0;j<3;j++)
         force[j] += fac * d[j];
   
   ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: RF   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist,fac));
   ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: RF   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist,fac));
}

/** Computes the Reaction Field pair force and adds this
    force to the particle forces (see \ref tclcommand_inter). 
    @param p1        Pointer to first particle.
    @param p2        Pointer to second/middle particle.
    @param d         Vector pointing from p1 to p2.
    @param dist      Distance between p1 and p2.
    @param force     returns the force on particle 1.
*/
inline void add_rf_coulomb_pair_force(Particle *p1, Particle *p2, double d[3], double dist, double force[3])
{
   if (dist < rf_params.r_cut)
   {
      add_rf_coulomb_pair_force_no_cutoff(p1,p2,d,dist,force);
   }
}

inline double rf_coulomb_pair_energy_no_cutoff(Particle *p1, Particle *p2, double dist)
{
   double  fac;
   fac = 1.0 / dist  -  (rf_params.B*dist*dist) / (2*rf_params.r_cut*rf_params.r_cut*rf_params.r_cut);
   //cut off part
   fac -= (1-rf_params.B/2)  / rf_params.r_cut;
   fac *= coulomb.prefactor * p1->p.q * p2->p.q;
   return fac;
}

inline double rf_coulomb_pair_energy(Particle *p1, Particle *p2, double dist)
{
  if (dist < rf_params.r_cut)
  {
     return rf_coulomb_pair_energy_no_cutoff(p1,p2,dist);
  }
  return 0.0;
}

/*from I. G. Tironi et al., J. Chem. Phys. 102, 5451 (1995)*/
#ifdef INTER_RF

///
int interrf_set_params(int part_type_a, int part_type_b,int rf_on);

inline void add_interrf_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist, double force[3])
{
#ifdef ONEPART_DEBUG
  double fac=0.0 ; /* TODO: this  variable was not declared. Now the code compiles, but I have no idea of what value to assign to it (MS) */
#endif
  if ((ia_params->rf_on ==1) && CUTOFF_CHECK(dist < rf_params.r_cut)) {
     add_rf_coulomb_pair_force_no_cutoff(p1,p2,d, dist,force);
  }

  ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: INTER_RF   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist,fac));
  ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: INTER_RF   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist,fac));
}

inline double interrf_pair_energy(Particle *p1, Particle *p2,IA_parameters *ia_params, double dist)
{
  double val;
  if ((ia_params->rf_on==1) && CUTOFF_CHECK(dist < rf_params.r_cut)) {
     val=rf_coulomb_pair_energy_no_cutoff(p1,p2,dist);
     return val;
  }
  return 0.0;
}

#endif

/*@}*/
#endif

#endif
