/*
  Copyright (C) 2010,2012,2013 The ESPResSo project
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
/** \file pressure.hpp
    Pressure calculation. Really similar to \ref energy.hpp "energy.h".
*/

#ifndef _PRESSURE_H
#define _PRESSURE_H

#include "utils.hpp"
#include "integrate.hpp"
#include "statistics.hpp"
#include "thermostat.hpp"
#include "forces.hpp"
#include "npt.hpp"

/** \name Exported Variables */
/************************************************************/
/*@{*/
///
extern Observable_stat virials, total_pressure, p_tensor, total_p_tensor;
///
extern Observable_stat_non_bonded virials_non_bonded, total_pressure_non_bonded, p_tensor_non_bonded, total_p_tensor_non_bonded;
/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/
void init_virials(Observable_stat *stat);
void init_virials_non_bonded(Observable_stat_non_bonded *stat_nb);
void init_p_tensor_non_bonded(Observable_stat_non_bonded *stat_nb);
void init_p_tensor(Observable_stat *stat);
void master_pressure_calc(int v_comp);


/** Calculates the pressure in the system from a virial expansion using the terms from \ref calculate_verlet_virials or \ref nsq_calculate_virials dependeing on the used cell system.<BR>
    @param result here the data about the scalar pressure are stored
    @param result_t here the data about the stress tensor are stored
    @param result_nb here the data about the intra- and inter- molecular nonbonded contributions to scalar pressure are stored
    @param result_t_nb here the data about the intra- and inter- molecular nonbonded contributions to stress tensor are stored
    @param v_comp flag which enables (1) compensation of the velocities required
		  for deriving a pressure reflecting \ref nptiso_struct::p_inst
		  (hence it only works with domain decomposition); naturally it
		  therefore doesn't make sense to use it without NpT.
*/
void pressure_calc(double *result, double *result_t, double *result_nb, double *result_t_nb, int v_comp);

/** Calculate non bonded energies between a pair of particles.
    @param p1        pointer to particle 1.
    @param p2        pointer to particle 2.
    @param d         vector between p1 and p2.
    @param dist      distance between p1 and p2.
    @param dist2     distance squared between p1 and p2. */
inline void add_non_bonded_pair_virials(Particle *p1, Particle *p2, double d[3],
					  double dist, double dist2)
{
  int p1molid, p2molid, k, l;
  double force[3] = {0, 0, 0};

  calc_non_bonded_pair_force(p1, p2,d, dist, dist2, force);

  *obsstat_nonbonded(&virials, p1->p.type, p2->p.type) += d[0]*force[0] + d[1]*force[1] + d[2]*force[2];

 /* stress tensor part */
  for(k=0;k<3;k++)
    for(l=0;l<3;l++)
      obsstat_nonbonded(&p_tensor, p1->p.type, p2->p.type)[k*3 + l] += force[k]*d[l];

  p1molid = p1->p.mol_id;
  p2molid = p2->p.mol_id;
  if ( p1molid == p2molid ) {
    *obsstat_nonbonded_intra(&virials_non_bonded, p1->p.type, p2->p.type) += d[0]*force[0] + d[1]*force[1] + d[2]*force[2];
    
    for(k=0;k<3;k++)
      for(l=0;l<3;l++)
        obsstat_nonbonded_intra(&p_tensor_non_bonded, p1->p.type, p2->p.type)[k*3 + l] += force[k]*d[l];
  } 
  if ( p1molid != p2molid ) {
    *obsstat_nonbonded_inter(&virials_non_bonded, p1->p.type, p2->p.type) += d[0]*force[0] + d[1]*force[1] + d[2]*force[2];
    
    for(k=0;k<3;k++)
      for(l=0;l<3;l++)
        obsstat_nonbonded_inter(&p_tensor_non_bonded, p1->p.type, p2->p.type)[k*3 + l] += force[k]*d[l];
  }
  
#ifdef ELECTROSTATICS
  /* real space coulomb */
  if (coulomb.method != COULOMB_NONE) {
    switch (coulomb.method) {
#ifdef P3M
    case COULOMB_P3M_GPU:
    case COULOMB_P3M:
      virials.coulomb[0] += p3m_pair_energy(p1->p.q*p2->p.q,d,dist2,dist);
      break;
#endif

    /* short range potentials, where we use the virial */
    /***************************************************/
    case COULOMB_DH: {
      double force[3] = {0, 0, 0};
    
      add_dh_coulomb_pair_force(p1,p2,d,dist, force);
      for(k=0;k<3;k++)
	for(l=0;l<3;l++)
	  p_tensor.coulomb[k*3 + l] += force[k]*d[l];
      virials.coulomb[0] += force[0]*d[0] + force[1]*d[1] + force[2]*d[2];
      break;
    }
    case COULOMB_RF: {
      double force[3] = {0, 0, 0};
    
      add_rf_coulomb_pair_force(p1,p2,d,dist, force);
      for(k=0;k<3;k++)
	for(l=0;l<3;l++)
	  p_tensor.coulomb[k*3 + l] += force[k]*d[l];
      virials.coulomb[0] += force[0]*d[0] + force[1]*d[1] + force[2]*d[2];
      break;
    }
    case COULOMB_INTER_RF:
      // this is done together with the other short range interactions
      break;
    default:
      fprintf(stderr,"calculating pressure for electrostatics method that doesn't have it implemented\n");
      break;
    }
  }
#endif /*ifdef ELECTROSTATICS */

#ifdef DIPOLES
  /* real space magnetic dipole-dipole */
  if (coulomb.Dmethod != DIPOLAR_NONE) {
    fprintf(stderr,"calculating pressure for magnetostatics which doesn't have it implemented\n");
  }
#endif /*ifdef DIPOLES */
}

inline void calc_bonded_force(Particle *p1, Particle *p2, Bonded_ia_parameters *iaparams, int *i, double dx[3], double force[3]) {
#ifdef TABULATED
  char* errtxt;
#endif

  /* Calculates the bonded force between two particles */
    switch(iaparams->type) {
    case BONDED_IA_FENE:
      calc_fene_pair_force(p1,p2,iaparams,dx,force);
      break;
    case BONDED_IA_HARMONIC:
      calc_harmonic_pair_force(p1,p2,iaparams,dx,force);
      break;
#ifdef LENNARD_JONES
    case BONDED_IA_SUBT_LJ:
      calc_subt_lj_pair_force(p1,p2,iaparams,dx,force);
      break;
#endif
      /* since it is not clear at the moment how to handle a many body interaction here, I skip it */
    case BONDED_IA_ANGLE_HARMONIC:
      (*i)++; force[0] = force[1] = force[2] = 0; break;
    case BONDED_IA_ANGLE_COSINE:
      (*i)++; force[0] = force[1] = force[2] = 0; break;
    case BONDED_IA_ANGLE_COSSQUARE:
      (*i)++; force[0] = force[1] = force[2] = 0; break;
    case BONDED_IA_ANGLEDIST:
      (*i)++; force[0] = force[1] = force[2] = 0; break;
    case BONDED_IA_DIHEDRAL:
      (*i)+=2; force[0] = force[1] = force[2] = 0; break;

#ifdef TABULATED
    case BONDED_IA_TABULATED:
      // printf("BONDED TAB, Particle: %d, P2: %d TYPE_TAB: %d\n",p1->p.identity,p2->p.identity,iparams->p.tab.type);
      switch(iaparams->p.tab.type) {
        case TAB_BOND_LENGTH:
	  calc_tab_bond_force(p1, p2, iaparams, dx, force); break;
        case TAB_BOND_ANGLE:
          (*i)++; force[0] = force[1] = force[2] = 0; break;
        case TAB_BOND_DIHEDRAL:
          (*i)+=2; force[0] = force[1] = force[2] = 0; break;
        default:
	  errtxt = runtime_error(128 + ES_INTEGER_SPACE);
	  ERROR_SPRINTF(errtxt,"{081 calc_bonded_force: tabulated bond type of atom %d unknown\n", p1->p.identity);
	  return;
      }
      break;
#endif
#ifdef OVERLAPPED
    case BONDED_IA_OVERLAPPED:
      // printf("BONDED OVERLAP, Particle: %d, P2: %d TYPE_OVERLAP: %d\n",p1->p.identity,p2->p.identity,iparams->p.tab.type);
      char *errtxt;
      switch(iaparams->p.overlap.type) {
        case OVERLAP_BOND_LENGTH:
          calc_overlap_bond_force(p1, p2, iaparams, dx, force); break;
        case OVERLAP_BOND_ANGLE:
          (*i)++; force[0] = force[1] = force[2] = 0; break;
        case OVERLAP_BOND_DIHEDRAL:
          (*i)+=2; force[0] = force[1] = force[2] = 0; break;
        default:
          errtxt = runtime_error(128 + ES_INTEGER_SPACE);
          ERROR_SPRINTF(errtxt,"{081 calc_bonded_force: overlapped bond type of atom %d unknown\n", p1->p.identity);
          return;
      }
      break;
#endif
#ifdef BOND_CONSTRAINT
    case BONDED_IA_RIGID_BOND:
      force[0] = force[1] = force[2] = 0; break;
#endif
#ifdef BOND_VIRTUAL
    case BONDED_IA_VIRTUAL_BOND:
      force[0] = force[1] = force[2] = 0; break;
#endif
    default :
      //      fprintf(stderr,"add_bonded_virials: WARNING: Bond type %d of atom %d unhandled\n",bonded_ia_params[type_num].type,p1->p.identity);
      fprintf(stderr,"add_bonded_virials: WARNING: Bond type %d , atom %d unhandled, Atom 2: %d\n",iaparams->type,p1->p.identity,p2->p.identity);
      force[0] = force[1] = force[2] = 0;
      break;
    }
}


/* calc_three_body_bonded_forces is called by add_three_body_bonded_stress. This
   routine is only entered for angular potentials. */
inline void calc_three_body_bonded_forces(Particle *p1, Particle *p2, Particle *p3,
              Bonded_ia_parameters *iaparams, double force1[3], double force2[3], double force3[3]) {

#ifdef TABULATED
  char* errtxt;
#endif

  switch(iaparams->type) {
#ifdef BOND_ANGLE_OLD
  case BONDED_IA_ANGLE_OLD:
    // p1 is *p_mid, p2 is *p_left, p3 is *p_right
    calc_angle_3body_forces(p1, p2, p3, iaparams, force1, force2, force3);
    break;
#endif
#ifdef BOND_ANGLE
  case BONDED_IA_ANGLE_HARMONIC:
    // p1 is *p_mid, p2 is *p_left, p3 is *p_right
    calc_angle_harmonic_3body_forces(p1, p2, p3, iaparams, force1, force2, force3);
    break;
 case BONDED_IA_ANGLE_COSINE:
    // p1 is *p_mid, p2 is *p_left, p3 is *p_right
    calc_angle_cosine_3body_forces(p1, p2, p3, iaparams, force1, force2, force3);
    break;
  case BONDED_IA_ANGLE_COSSQUARE:
    // p1 is *p_mid, p2 is *p_left, p3 is *p_right
    calc_angle_cossquare_3body_forces(p1, p2, p3, iaparams, force1, force2, force3);
    break;
#endif
#ifdef TABULATED
  case BONDED_IA_TABULATED:
    switch(iaparams->p.tab.type) {
      case TAB_BOND_ANGLE:
        // p1 is *p_mid, p2 is *p_left, p3 is *p_right
        calc_angle_3body_tabulated_forces(p1, p2, p3, iaparams, force1, force2, force3);
        break;
      default:
        errtxt = runtime_error(128 + ES_INTEGER_SPACE);
        ERROR_SPRINTF(errtxt,"{081 calc_bonded_force: tabulated bond type of atom %d unknown\n", p1->p.identity);
        return;
      }
      break;
#endif
  default :
    fprintf(stderr,"calc_three_body_bonded_forces: \
            WARNING: Bond type %d , atom %d unhandled, Atom 2: %d\n", \
            iaparams->type, p1->p.identity, p2->p.identity);
    break;
  }
}


/** Calculate bonded virials for one particle.
    For performance reasons the force routines add their values directly to the particles.
    So here we do some tricks to get the value out without changing the forces.
    @param p1 particle for which to calculate virials
*/
inline void add_bonded_virials(Particle *p1)
{
  double dx[3], force[3] = {0,0,0};
  char *errtxt;
  Particle *p2;
  Bonded_ia_parameters *iaparams;

  int i, k, l;
  int type_num;

  i = 0;
  while(i<p1->bl.n) {
    type_num = p1->bl.e[i++];
    iaparams = &bonded_ia_params[type_num];

    /* fetch particle 2 */
    p2 = local_particles[p1->bl.e[i++]];
    if ( ! p2 ) {
      // for harmonic spring:
      // if cutoff was defined and p2 is not there it is anyway outside the cutoff, see calc_maximal_cutoff()
      if ((type_num==BONDED_IA_HARMONIC)&&(iaparams->p.harmonic.r_cut>0)) return;
      errtxt = runtime_error(128 + 2*ES_INTEGER_SPACE);
      ERROR_SPRINTF(errtxt,"{088 bond broken between particles %d and %d (particles not stored on the same node)} ",
		    p1->p.identity, p1->bl.e[i-1]);
      return;
    }

    get_mi_vector(dx, p1->r.p, p2->r.p);
    calc_bonded_force(p1,p2,iaparams,&i,dx,force);
    *obsstat_bonded(&virials, type_num) += dx[0]*force[0] + dx[1]*force[1] + dx[2]*force[2];

 /* stress tensor part */
    for(k=0;k<3;k++)
      for(l=0;l<3;l++)
	obsstat_bonded(&p_tensor, type_num)[k*3 + l] += force[k]*dx[l];

  }
}

/** Calculate the contribution to the stress tensor from angular potentials.
    The central particle of the three-particle interaction is responsible
    for the contribution of the entire interaction - this is the coding
    not the physics.
*/
inline void add_three_body_bonded_stress(Particle *p1) {
  double dx12[3]; // espresso notation
  double dx21[3];
  double dx31[3];
  double force1[3];
  double force2[3];
  double force3[3];

  char *errtxt;
  Particle *p2;
  Particle *p3;
  Bonded_ia_parameters *iaparams;

  int i, k, j, l;
  int type_num;
  int type;

  i = 0;
  while(i < p1->bl.n) {
    /* scan bond list for angular interactions */
    type_num = p1->bl.e[i];
    iaparams = &bonded_ia_params[type_num];
    type = iaparams->type;

    if(type == BONDED_IA_ANGLE_HARMONIC || type == BONDED_IA_ANGLE_COSINE || type == BONDED_IA_ANGLE_COSSQUARE){
      p2 = local_particles[p1->bl.e[++i]];
      p3 = local_particles[p1->bl.e[++i]];

      get_mi_vector(dx12, p1->r.p, p2->r.p);
      for(j = 0; j < 3; j++)
        dx21[j] = -dx12[j];

      get_mi_vector(dx31, p3->r.p, p1->r.p);

      for(j = 0; j < 3; j++) {
        force1[j] = 0.0;
        force2[j] = 0.0;
        force3[j] = 0.0;
      }

      calc_three_body_bonded_forces(p1, p2, p3, iaparams, force1, force2, force3);

      /* uncomment the next line to see that the virial is indeed zero */
      //printf("W = %g\n", scalar(force2, dx21) + scalar(force3, dx31));

      /* three-body bonded interactions contribute to the stress but not the scalar pressure */
      for(k = 0; k < 3; k++) {
        for(l = 0; l < 3; l++) {
          obsstat_bonded(&p_tensor, type_num)[3 * k + l] += force2[k] * dx21[l] + force3[k] * dx31[l];
        }
      }
      i = i + 1;
    }
    // skip over non-angular interactions
    else if(type == BONDED_IA_FENE) {
      i = i + 2;
    }
    else if(type == BONDED_IA_STRETCHING_FORCE) {
      i = i + 2;
    }
    else if(type == BONDED_IA_STRETCHLIN_FORCE) {
      i = i + 2;
    }
    else if(type == BONDED_IA_AREA_FORCE_LOCAL) {
      i = i + 3;
    }
    else if(type == BONDED_IA_AREA_FORCE_GLOBAL) {
      i = i + 3;
    }
    else if(type == BONDED_IA_BENDING_FORCE) {
      i = i + 4;
    }
    else if(type == BONDED_IA_VOLUME_FORCE) {
      i = i + 3;
    }
    else if(type == BONDED_IA_HARMONIC) {
      i = i + 2;
    }
#ifdef LENNARD_JONES
    else if(type == BONDED_IA_SUBT_LJ) {
      i = i + 2;
    }
#endif
    else if(type == BONDED_IA_ANGLEDIST) {
      i = i + 3;
    }
    else if(type == BONDED_IA_DIHEDRAL) {
      i = i + 4;
    }
#ifdef TABULATED
    else if(type == BONDED_IA_TABULATED) {
      if(iaparams->p.tab.type == TAB_BOND_LENGTH) {
        i = i + 2;
      }
      else if(iaparams->p.tab.type == TAB_BOND_ANGLE) {
        p2 = local_particles[p1->bl.e[++i]];
        p3 = local_particles[p1->bl.e[++i]];

        get_mi_vector(dx12, p1->r.p, p2->r.p);
        for(j = 0; j < 3; j++)
          dx21[j] = -dx12[j];

        get_mi_vector(dx31, p3->r.p, p1->r.p);

        for(j = 0; j < 3; j++) {
          force1[j] = 0.0;
          force2[j] = 0.0;
          force3[j] = 0.0;
        }

        calc_three_body_bonded_forces(p1, p2, p3, iaparams, force1, force2, force3);

        /* uncomment the next line to see that the virial is indeed zero */
        //printf("W = %g\n", scalar(force2, dx21) + scalar(force3, dx31));

        /* three-body bonded interactions contribute to the stress but not the scalar pressure */
        for(k = 0; k < 3; k++) {
          for(l = 0; l < 3; l++) {
            obsstat_bonded(&p_tensor, type_num)[3 * k + l] += force2[k] * dx21[l] + force3[k] * dx31[l];
          }
        }
        i = i + 1;
      }
      else if (iaparams->p.tab.type == TAB_BOND_DIHEDRAL) {
        i = i + 4;
      }
      else {
        errtxt = runtime_error(128 + ES_INTEGER_SPACE);
        ERROR_SPRINTF(errtxt,"add_three_body_bonded_stress: match not found for particle %d.\n", p1->p.identity);
      }
    }
#endif
#ifdef BOND_CONSTRAINT
    else if(type == BONDED_IA_RIGID_BOND) {
      i = i + 2;
    }
#endif
#ifdef BOND_VIRTUAL
    else if(type == BONDED_IA_VIRTUAL_BOND) {
      i = i + 2;
    }
#endif
    else {
      errtxt = runtime_error(128 + ES_INTEGER_SPACE);
      ERROR_SPRINTF(errtxt,"add_three_body_bonded_stress: match not found for particle %d.\n", p1->p.identity);
    }
  } 
}
 
/** Calculate kinetic pressure (aka energy) for one particle.
    @param p1 particle for which to calculate pressure
    @param v_comp flag which enables (1) compensation of the velocities required
		  for deriving a pressure reflecting \ref nptiso_struct::p_inst
		  (hence it only works with domain decomposition); naturally it
		  therefore doesn't make sense to use it without NpT.
*/
inline void add_kinetic_virials(Particle *p1,int v_comp)
{
  int k, l;
  /* kinetic energy */
  if(v_comp)
    virials.data.e[0] += (SQR(p1->m.v[0] - p1->f.f[0]) + SQR(p1->m.v[1] - p1->f.f[1]) + SQR(p1->m.v[2] - p1->f.f[2]))*PMASS(*p1);
  else
    virials.data.e[0] += (SQR(p1->m.v[0]) + SQR(p1->m.v[1]) + SQR(p1->m.v[2]))*PMASS(*p1);


  /* ideal gas contribution (the rescaling of the velocities by '/=time_step' each will be done later) */
  for(k=0;k<3;k++)
    for(l=0;l<3;l++)
      p_tensor.data.e[k*3 + l] += (p1->m.v[k])*(p1->m.v[l])*PMASS(*p1);

}

/** implementation of 'analyse local_stress_tensor */
int local_stress_tensor_calc (DoubleList *TensorInBin, int bins[3], int periodic[3], double range_start[3], double range[3]);

/** function to calculate stress tensor for the observables */
int observable_compute_stress_tensor(int v_comp, double *A, unsigned int n_A);


/*@}*/

#endif
