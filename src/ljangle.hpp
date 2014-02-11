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
#ifndef _LJANGLE_H
#define _LJANGLE_H

/** \file ljangle.hpp
 *  Routines to calculate the lennard-jones 12-10 with angular dependance.
 *  The potential is a product of a 12-10 LJ potential with two cos^2.
 *  The potential actually relies on 6 particles: the 2 primary beads
 *  and each bead needs two other particles to define an orientation.
 *  We calculate the interaction explicitly *without* the use of the ROTATION feature.
 *
 *  Optional: simulate two different environments in which the interaction
 *  strengths differ. For example: simulate hydrogen-bonds both in water and
 *  inside a bilayer. The two environments are distinguished by their
 *  z-position. Input: the midplane of the 2nd environment, its total
 *  thickness, the thickness of the interface, and the value of the
 *  interaction strength in this medium. The interaction strengh of the second
 *  environment must be *stronger* than of the first one.
 *
 *  \ref forces.cpp
 */

#include "utils.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "mol_cut.hpp"
#include "grid.hpp"
#include "forcecap.hpp"


#ifdef LJ_ANGLE
#include <cmath>

/** For the warmup you can cap the singularity of the directionnal LJ
    potential at r=0. look into the warmup documentation for more
    details (who wants to write that?).*/
    /* replaced by global force_cap -- to be removed */

int ljangle_set_params(int part_type_a, int part_type_b,
				double eps, double sig, double cut,
				int b1p, int b1n, int b2p, int b2n,
				double cap_radius, double z0, double dz, 
				double kappa, double epsprime);


/** Calculate lj-angle force between particle p1 and p2 
    Involves 6 particles total */
inline void add_ljangle_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
				     double d[3], double dist)
{
  if(!CUTOFF_CHECK(dist < ia_params->LJANGLE_cut)) 
    return;
  int j;
  double frac2=0.0, frac10=0.0, rad=0.0, radprime=0.0;
  double r31[3], r41[3], r52[3], r62[3], rij[3], rik[3], rkn[3];
  double l_rij, l_rik, l_rkn, l_rij2, l_rik2, l_rkn2;
  double cos_jik, cos_ikn;
  double angular_jik, angular_ikn, angular_jik_prime, angular_ikn_prime;
  /* Recreate angular dependence of potential by including 6 particles instead of 2. */    
  Particle *p3=NULL, *p4=NULL, *p5=NULL, *p6=NULL;
  int part1p, part1n, part2p, part2n;
  /* Optional 2nd environment */
  double effective_eps=ia_params->LJANGLE_eps, z1, z2, z_middle, z_ref, z_ref5,
    localz0=ia_params->LJANGLE_z0, localdz2=ia_params->LJANGLE_dz/2.,
    localkappa6=pow(1/ia_params->LJANGLE_kappa,6);


	  
    /* Retrieve the bonded partners from parsing */
    if (ia_params->LJANGLE_bonded1type == p1->p.type) {
      part1p = p1->p.identity + ia_params->LJANGLE_bonded1pos;
      part1n = p1->p.identity + ia_params->LJANGLE_bonded1neg;
      part2p = p2->p.identity + ia_params->LJANGLE_bonded2pos;
      part2n = p2->p.identity + ia_params->LJANGLE_bonded2neg;
    } else {
      part1p = p1->p.identity + ia_params->LJANGLE_bonded2pos;
      part1n = p1->p.identity + ia_params->LJANGLE_bonded2neg;
      part2p = p2->p.identity + ia_params->LJANGLE_bonded1pos;
      part2n = p2->p.identity + ia_params->LJANGLE_bonded1neg;
    }

    if (part1p >= 0 && part1p < n_part &&
	part1n >= 0 && part1n < n_part &&
	part2p >= 0 && part2p < n_part &&
	part2n >= 0 && part2n < n_part     ) {
	
      p3 = local_particles[part1p];
      p4 = local_particles[part1n];
      p5 = local_particles[part2p];
      p6 = local_particles[part2n];

      /* Check whether pointers have been allocated.
       * Otherwise, there's a communication error (verlet skin too small).
       */			  
      if (p3==NULL ||p4==NULL ||p5==NULL ||p6==NULL) 
	fprintf(stderr, "LJANGLE - Communication error, all particles cannot be reached locally.\n");
			

			
      get_mi_vector(r31, p1->r.p, p3->r.p);
      get_mi_vector(r41, p1->r.p, p4->r.p);
      get_mi_vector(r52, p2->r.p, p5->r.p);
      get_mi_vector(r62, p2->r.p, p6->r.p);
	
		  
      /* Bead i represents the central particle of monomer 1.
       * Bead j is the virtual particle that gives the orientation of monomer 1.
       * Bead k represents the central particle of monomer 2.
       * Bead n is the virtual particle that gives the orientation of monomer 2.
       */
      for(j=0;j<3;++j){
	rij[j] = r31[j] + r41[j];
	rik[j] = -d[j]; /* At this point, rik[3] has length dist */
	rkn[j] = r52[j] + r62[j];
      }
	  
      l_rij2 = sqrlen(rij);
      l_rik2 = dist*dist;
      l_rkn2 = sqrlen(rkn);
      l_rij  = sqrt(l_rij2);
      l_rik  = dist;
      l_rkn  = sqrt(l_rkn2);
	  
      cos_jik = scalar(rij,rik)/(l_rij*l_rik);
      cos_ikn = -scalar(rik,rkn)/(l_rik*l_rkn);
	  
      if(cos_jik>0. && cos_ikn>0.) { 
	angular_jik       = pow(cos_jik,2);
	angular_ikn       = pow(cos_ikn,2);
	angular_jik_prime = -2*cos_jik;
	angular_ikn_prime = -2*cos_ikn;
	      
	/* Optional 2nd environment */
	if (localdz2 > 0.) {
	  /* calculate center position of the interaction (calculate minimal distance) */
	  z1  = p1->r.p[2];
	  z2  = p2->r.p[2];
	  z_middle  = z1 + z2;
	  z_middle /= 2.;
#ifdef PARTIAL_PERIODIC
	  if (PERIODIC(2))
#endif
	    if (z_middle > box_l[2] || z_middle < 0.) 
	      z_middle -= dround(z_middle *box_l_i[2])*box_l[2];

	  /* If we're in the environment #2 region, calculate the new interaction strength */
	  if (z_middle > localz0 - localdz2 && z_middle <= localz0) {
	    z_ref = z_middle - localz0 + localdz2;
	    z_ref5 = pow(z_ref,5);
	    effective_eps += (ia_params->LJANGLE_epsprime - ia_params->LJANGLE_eps) *
	      localkappa6*z_ref5*fabs(z_ref) / 
	      (1 + localkappa6*z_ref5*z_ref);
	  } else if (z_middle > localz0 && z_middle < localz0 + localdz2) {
	    z_ref = z_middle - localz0 - localdz2;
	    z_ref5 = pow(z_ref,5);
	    effective_eps -= (ia_params->LJANGLE_epsprime - ia_params->LJANGLE_eps) *
	      localkappa6*z_ref5*fabs(z_ref) / 
	      (1 + localkappa6*z_ref5*z_ref);
	  }
	}  

	/* normal case: resulting force/energy smaller than capping. */
	if(dist > ia_params->LJANGLE_capradius) {
	  frac2 = SQR(ia_params->LJANGLE_sig/dist);
	  frac10 = frac2*frac2*frac2*frac2*frac2;
	  rad        = effective_eps * frac10*(5.0 * frac2 - 6.0);
	  radprime   = 60.0 * effective_eps * frac10*(1.0 - frac2) / dist;
		  
#ifdef LJ_WARN_WHEN_CLOSE
	  if(radprime > 1000) fprintf(stderr,"%d: LJANGLE-Warning: Pair (%d-%d) force=%f dist=%f\n",
				      this_node,p1->p.identity,p2->p.identity,radprime,dist);
#endif
	}
	/* capped part of lj-angle potential. */
	else if(dist > 0.0) {
	  /* set the length of rik to capradius */
	  for (j=0;j<3;++j) 
	    rik[j] *= ia_params->LJANGLE_capradius/dist;
		  
	  l_rik = ia_params->LJANGLE_capradius;
		  
	  frac2 = SQR(ia_params->LJANGLE_sig/ia_params->LJANGLE_capradius);
	  frac10 = frac2*frac2*frac2*frac2*frac2;
		  
	  rad        = effective_eps * frac10*(5.0 * frac2 - 6.0);
	  radprime   = 60.0 * effective_eps * frac10*(1.0 - frac2) / (ia_params->LJANGLE_capradius);
	}
	      
	/* Regroup the last two cases in one */
	/* Propagate all forces in this function rather than in the forces.hpp file */
	if (dist > 0.0){
	  for(j=0;j<3;++j){
	    p1->f.f[j] += angular_jik * angular_ikn * rik[j]/l_rik *radprime 
	      + angular_ikn * rad * angular_jik_prime 
	      * ( ( 2*rik[j]-rij[j] )/( l_rij*l_rik )
		  - cos_jik*( 2*rij[j]/l_rij2 - rik[j]/l_rik2 ) )
	      + angular_jik * rad * angular_ikn_prime
	      * ( rkn[j]/( l_rik*l_rkn )
		  + cos_ikn*rik[j]/l_rik2 );
	    p2->f.f[j] += -angular_jik * angular_ikn * rik[j]/l_rik *radprime
	      + angular_ikn * rad * angular_jik_prime 
	      * ( rij[j]/( l_rij*l_rik )
		  - cos_jik*rik[j]/l_rik2 )
	      + angular_jik * rad * angular_ikn_prime
	      * ( -(2*rik[j]+rkn[j])/(l_rik*l_rkn)
		  - cos_ikn*( 2*rkn[j]/l_rkn2 + rik[j]/l_rik2 ) );
	    p3->f.f[j] += angular_ikn * rad * angular_jik_prime 
	      * ( -rik[j]/(l_rij*l_rik) 
		  + cos_jik*rij[j]/l_rij2 );
	    p4->f.f[j] += angular_ikn * rad * angular_jik_prime 
	      * ( -rik[j]/(l_rij*l_rik) 
		  + cos_jik*rij[j]/l_rij2 );
	    p5->f.f[j] += angular_jik * rad * angular_ikn_prime
	      * ( rik[j]/(l_rik*l_rkn)
		  + cos_ikn*rkn[j]/l_rkn2 );
	    p6->f.f[j] += angular_jik * rad * angular_ikn_prime
	      * ( rik[j]/(l_rik*l_rkn)
		  + cos_ikn*rkn[j]/l_rkn2 );
	  }
	}
	      
	      
	/* this should not happen! In this case consider only radial potential*/
	else {
	  LJ_TRACE(fprintf(stderr, "%d: LJ-angle warning: Particles id1=%d id2=%d exactly on top of each other\n",this_node,p1->p.identity,p2->p.identity));
		  
	  frac2 = SQR(ia_params->LJANGLE_sig/ia_params->LJANGLE_capradius);
	  frac10 = frac2*frac2*frac2*frac2*frac2;
	  rad   = effective_eps * frac10*(5.0 * frac2 - 6.0) / ia_params->LJANGLE_capradius;
		  
	  p1->f.f[0] += rad * ia_params->LJANGLE_capradius;
	}
	      
	      
	ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: LJANGLE   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist,radprime));
	ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: LJANGLE   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist,radprime));
	      
	LJ_TRACE(fprintf(stderr,"%d: LJANGLE: Pair (%d-%d) dist=%.3f: force+-: (%.3e,%.3e,%.3e)\n",
			 this_node,p1->p.identity,p2->p.identity,dist,rad*d[0],rad*d[1],rad*d[2]));
      }
    }
}

/** calculate Lennard jones energy between particle p1 and p2. */
inline double ljangle_pair_energy(Particle *p1, Particle *p2, IA_parameters *ia_params,
				    double d[3], double dist)
{
  if(!CUTOFF_CHECK(dist < ia_params->LJANGLE_cut)) 
    return 0.0;

  int j;
  double frac2, frac10;
  double r31[3], r41[3], r52[3], r62[3], rij[3], rik[3], rkn[3];
  double l_rij, l_rik, l_rkn;
  double cos_jik, cos_ikn;
  double angular_jik, angular_ikn;
  /* Recreate angular dependence of potential by including 6 particles instead of 2. */    
  Particle *p3=NULL, *p4=NULL, *p5=NULL, *p6=NULL;
  int part1p, part1n, part2p, part2n;
  /* Optional 2nd environment */
  double effective_eps=ia_params->LJANGLE_eps, z1, z2, z_middle, z_ref, z_ref5,
    localz0=ia_params->LJANGLE_z0, localdz2=ia_params->LJANGLE_dz/2.,
    localkappa6=pow(1/ia_params->LJANGLE_kappa,6);

	
    /* Retrieve the bonded partners from parsing */
    if (ia_params->LJANGLE_bonded1type == p1->p.type) {
      part1p = p1->p.identity + ia_params->LJANGLE_bonded1pos;
      part1n = p1->p.identity + ia_params->LJANGLE_bonded1neg;
      part2p = p2->p.identity + ia_params->LJANGLE_bonded2pos;
      part2n = p2->p.identity + ia_params->LJANGLE_bonded2neg;
    } else {
      part1p = p1->p.identity + ia_params->LJANGLE_bonded2pos;
      part1n = p1->p.identity + ia_params->LJANGLE_bonded2neg;
      part2p = p2->p.identity + ia_params->LJANGLE_bonded1pos;
      part2n = p2->p.identity + ia_params->LJANGLE_bonded1neg;
    }
	
    if (part1p >= 0 && part1p < n_part &&
	part1n >= 0 && part1n < n_part &&
	part2p >= 0 && part2p < n_part &&
	part2n >= 0 && part2n < n_part     ) {
	    
      p3 = local_particles[part1p];
      p4 = local_particles[part1n];
      p5 = local_particles[part2p];
      p6 = local_particles[part2n];
	    
      get_mi_vector(r31, p1->r.p, p3->r.p);
      get_mi_vector(r41, p1->r.p, p4->r.p);
      get_mi_vector(r52, p2->r.p, p5->r.p);
      get_mi_vector(r62, p2->r.p, p6->r.p);
	    
      for(j=0;j<3;++j){
	rij[j] = r31[j] + r41[j];
	rik[j] = -d[j]; /* At this point, rik[3] has length dist */
	rkn[j] = r52[j] + r62[j];
      }
	    
      l_rij = sqrt(sqrlen(rij));
      l_rik = dist;
      l_rkn = sqrt(sqrlen(rkn));
	    
      cos_jik =  scalar(rij,rik)/(dist*l_rij);
      cos_ikn = -scalar(rik,rkn)/(l_rik*l_rkn);
	    
      angular_jik = pow(cos_jik,2);
      angular_ikn = pow(cos_ikn,2);
	    
      if(cos_jik>0. && cos_ikn>0.) {
	/* Optional 2nd environment */
	if (localdz2 > 0.) {
	  /* calculate center position of the interaction (calculate minimal distance) */
	  z1  = p1->r.p[2];
	  z2  = p2->r.p[2];
	  z_middle  = z1 + z2;
	  z_middle /= 2.;
#ifdef PARTIAL_PERIODIC
	  if (PERIODIC(2))
#endif
	    if (z_middle > box_l[2] || z_middle < 0.) 
	      z_middle -= dround(z_middle *box_l_i[2])*box_l[2];

	  /* If we're in the environment #2 region, calculate the new interaction strength */
	  if (z_middle > localz0-localdz2 && z_middle <= localz0) {
	    z_ref = z_middle - localz0 + localdz2;
	    z_ref5 = pow(z_ref,5);
	    effective_eps += (ia_params->LJANGLE_epsprime - ia_params->LJANGLE_eps) *
	      localkappa6 * z_ref5 * fabs(z_ref) / 
	      (1 + localkappa6 * z_ref5 * z_ref);
	  } else if (z_middle > localz0 && z_middle < localz0+localdz2) {
	    z_ref = z_middle - localz0 - localdz2;
	    z_ref5 = pow(z_ref,5);
	    effective_eps -= (ia_params->LJANGLE_epsprime - ia_params->LJANGLE_eps) *
	      localkappa6 * z_ref5 * fabs(z_ref) / 
	      (1 + localkappa6 * z_ref5 * z_ref);
	  }
	}  

	/* normal case: resulting force/energy smaller than capping. */
	if(dist > ia_params->LJANGLE_capradius) {
	  frac2  = SQR(ia_params->LJANGLE_sig/dist);
	  frac10 = frac2*frac2*frac2*frac2*frac2;
	  return effective_eps * frac10*(5.0 * frac2 - 6.0) * angular_jik * angular_ikn;
	}
	/* capped part of lj potential. */
	else if(dist > 0.0) {
	  frac2  = SQR(ia_params->LJANGLE_sig/ia_params->LJANGLE_capradius);
	  frac10 = frac2*frac2*frac2*frac2*frac2;
	  return effective_eps * frac10*(5.0 * frac2 - 6.0) * angular_jik * angular_ikn;
	}
	/* this should not happen!  In this case consider only radial potential*/
	else {
	  frac2  = SQR(ia_params->LJANGLE_sig/ia_params->LJANGLE_capradius);
	  frac10 = frac2*frac2*frac2*frac2*frac2;
	  return effective_eps * frac10*(5.0 * frac2 - 6.0);
	}
      }
    }
  
  return 0.0;
}


/** calculate ljangle_capradius from force_cap */
/* This routine does not take the optional 2nd environment into account. */
void calc_ljangle_cap_radii();


#endif /* ifdef LJ_ANGLE */
/* LJANGLE_H */
#endif 

