// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
/** \file pressure.h
    Pressure calculation. Really similar to \ref energy.h "energy.h".
*/

#ifndef PRESSURE_H
#define PRESSURE_H


#include "utils.h"
#include "integrate.h"
#include "statistics.h"
#include "thermostat.h"
#include "communication.h"
#include "adresso.h"

/************************************************
 * data types
 ************************************************/

/** Structure to hold all variables related to the isotropic NpT-integration scheme. */
typedef struct { 
  /** mass of a virtual piston representing the shaken box */
  double piston;
  /** inverse of piston */
  double inv_piston;
  /** isotropic volume.  Note that we use the term volume throughout
      although for a 2d or 1d system we mean Area and Length
      respectively */
  double volume;

  /** desired pressure to which the algorithm strives to */
  double p_ext;
  /** instantaneous pressure the system currently has */
  double p_inst;
  /** instantaneous pressure averaged over current integration cycle */
  double p_inst_av;
  /** difference between \ref p_ext and \ref p_inst */
  double p_diff;
  /** virial (short-range) components of \ref p_inst */
  double p_vir[3];
  /** ideal gas components of \ref p_inst, derived from the velocities */
  double p_vel[3];
  /** flag which indicates if \ref p_vel may (0) or may not (1) be used
      in offline pressure calculations such as 'analyze p_inst' */ 
  int invalidate_p_vel;
  /** geometry information for the npt integrator.  Holds the vector
      <dir, dir ,dir> where a positive value for dir indicates that
      box movement is allowed in that direction. To check whether a
      given direction is turned on use bitwise comparison with \ref
      nptgeom_dir */
  int geometry;
  /** bitwise comparison values corresponding to different directions*/
  int nptgeom_dir[3];
  /** The number of dimensions in which npt boxlength motion is coupled to particles */
  int dimension;
  /** Set this flag if you want all box dimensions to be identical. Needed for electrostatics and magnetostatics.  
      If the value of dimension is less than 3 then box length motion in one or more
      directions will be decoupled from the particle motion */
  int cubic_box;
  /** An index to one of the non_constant dimensions. handy if you just want the variable box_l */
  int non_const_dim;
} nptiso_struct;
extern nptiso_struct nptiso;

/** Allowable values for nptiso.geometry*/
#define NPTGEOM_XDIR 1
#define NPTGEOM_YDIR 2
#define NPTGEOM_ZDIR 4

/* include the potential files */
#include "p3m.h"
#include "lj.h"
#include "ljgen.h"
#include "ljangle.h"
#include "steppot.h"
#include "bmhtf-nacl.h"
#include "buckingham.h"
#include "morse.h"
#include "soft_sphere.h"
#include "ljcos.h"
#include "ljcos2.h"
#include "tab.h"
#include "gb.h"
#include "fene.h"
#include "harmonic.h"
#include "subt_lj.h"
#include "angle.h"
#include "angledist.h"
#include "dihedral.h"
#include "debye_hueckel.h"
#include "endangledist.h"
#include "reaction_field.h"
#include "mmm1d.h"
#include "mol_cut.h"
#include "tunable_slip.h"

/** \name Exported Variables */
/************************************************************/
/*@{*/
///
extern Observable_stat virials, total_pressure;
///
extern Observable_stat p_tensor;
///
extern Observable_stat_non_bonded virials_non_bonded, total_pressure_non_bonded;
///
extern Observable_stat_non_bonded p_tensor_non_bonded;
/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/

/** Callback for setting \ref nptiso_struct::piston */
int piston_callback(Tcl_Interp *interp, void *_data);
/** Callback for setting \ref nptiso_struct::p_ext */
int p_ext_callback(Tcl_Interp *interp, void *_data);
/** Callback for setting \ref nptiso_struct::p_diff */
int p_diff_callback(Tcl_Interp *interp, void *_data);

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

MDINLINE void calc_non_bonded_pair_force_parts(Particle *p1, Particle *p2, IA_parameters *ia_params,double d[3],
					 double dist, double dist2, double force[3],double torgue1[3],double torgue2[3])
{
#ifdef NO_INTRA_NB
  if (p1->p.mol_id==p2->p.mol_id) return;
#endif
  /* lennard jones */
#ifdef LENNARD_JONES
  add_lj_pair_force(p1,p2,ia_params,d,dist, force);
#endif
  /* lennard jones generic */
#ifdef LENNARD_JONES_GENERIC
  add_ljgen_pair_force(p1,p2,ia_params,d,dist, force);
#endif
  /* Directional LJ */
#ifdef LJ_ANGLE
  /* The forces are propagated within the function */
  add_ljangle_pair_force(p1,p2,ia_params,d,dist);
#endif
  /* smooth step */
#ifdef SMOOTH_STEP
  add_SmSt_pair_force(p1,p2,ia_params,d,dist,dist2, force);
#endif
  /* BMHTF NaCl */
#ifdef BMHTF_NACL
  add_BMHTF_pair_force(p1,p2,ia_params,d,dist,dist2, force);
#endif
  /* buckingham*/
#ifdef BUCKINGHAM
  add_buck_pair_force(p1,p2,ia_params,d,dist,force);
#endif
  /* morse*/
#ifdef MORSE
  add_morse_pair_force(p1,p2,ia_params,d,dist,force);
#endif
 /*soft-sphere potential*/
#ifdef SOFT_SPHERE
  add_soft_pair_force(p1,p2,ia_params,d,dist,force);
#endif
  /* lennard jones cosine */
#ifdef LJCOS
  add_ljcos_pair_force(p1,p2,ia_params,d,dist,force);
#endif
  /* lennard jones cosine */
#ifdef LJCOS2
  add_ljcos2_pair_force(p1,p2,ia_params,d,dist,force);
#endif
  /* tabulated */
#ifdef TABULATED
  add_tabulated_pair_force(p1,p2,ia_params,d,dist,force);
#endif
  /* Gay-Berne */
#ifdef ROTATION
  add_gb_pair_force(p1,p2,ia_params,d,dist,force,torgue1,torgue2);
#endif
#ifdef INTER_RF
  add_interrf_pair_force(p1,p2,ia_params,d,dist, force);
#endif
#ifdef ADRESS
#ifdef INTERFACE_CORRECTION
  add_adress_tab_pair_force(p1,p2,ia_params,d,dist,force);
#endif
#endif
}

MDINLINE void calc_non_bonded_pair_force(Particle *p1,Particle *p2,IA_parameters *ia_params,double d[3],double dist,double dist2,double force[3],double t1[3],double t2[3]){
#ifdef MOL_CUT
   //You may want to put a correction factor and correction term for smoothing function else then theta
   if (checkIfParticlesInteractViaMolCut(p1,p2,ia_params)==1)
#endif
   {
      calc_non_bonded_pair_force_parts(p1, p2, ia_params,d, dist, dist2,force,t1,t2);
   }
}

MDINLINE void calc_non_bonded_pair_force_simple(Particle *p1,Particle *p2,double d[3],double dist,double dist2,double force[3]){
   IA_parameters *ia_params = get_ia_param(p1->p.type,p2->p.type);
   double t1[3],t2[3];
#ifdef ADRESS
  int j;
  double force_weight=adress_non_bonded_force_weight(p1,p2);
  if (force_weight<ROUND_ERROR_PREC) return;
#endif
  calc_non_bonded_pair_force(p1,p2,ia_params,d,dist,dist2,force,t1,t2);
#ifdef ADRESS
   for (j=0;j<3;j++){
      force[j]*=force_weight;
   }
#endif
}

MDINLINE void calc_non_bonded_pair_force_from_partcfg(Particle *p1,Particle *p2,IA_parameters *ia_params,double d[3],double dist,double dist2,double force[3],double t1[3],double t2[3]){
#ifdef MOL_CUT
   //You may want to put a correction factor and correction term for smoothing function else then theta
   if (checkIfParticlesInteractViaMolCut_partcfg(p1,p2,ia_params)==1)
#endif
   {
     calc_non_bonded_pair_force_parts(p1, p2, ia_params,d, dist, dist2,force,t1,t2);
   }
}

MDINLINE void calc_non_bonded_pair_force_from_partcfg_simple(Particle *p1,Particle *p2,double d[3],double dist,double dist2,double force[3]){
   IA_parameters *ia_params = get_ia_param(p1->p.type,p2->p.type);
   double t1[3],t2[3];
   calc_non_bonded_pair_force_from_partcfg(p1,p2,ia_params,d,dist,dist2,force,t1,t2);
}

/** Calculate non bonded energies between a pair of particles.
    @param p1        pointer to particle 1.
    @param p2        pointer to particle 2.
    @param d         vector between p1 and p2.
    @param dist      distance between p1 and p2.
    @param dist2     distance squared between p1 and p2. */
MDINLINE void add_non_bonded_pair_virials(Particle *p1, Particle *p2, double d[3],
					  double dist, double dist2)
{
  int p1molid, p2molid, k, l;
  double force[3] = {0, 0, 0};
#if defined(ELECTROSTATICS)  || defined(MAGNETOSTATICS)
  double ret=0;
#endif
  calc_non_bonded_pair_force_simple(p1, p2,d, dist, dist2,force);

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
#ifdef ELP3M
    case COULOMB_P3M:
      ret = p3m_coulomb_pair_energy(p1->p.q*p2->p.q,d,dist2,dist);
      break;
#endif
    case COULOMB_DH:
      ret = dh_coulomb_pair_energy(p1,p2,dist);
      break;
    case COULOMB_RF:
      ret = rf_coulomb_pair_energy(p1,p2,dist);
      break;
    case COULOMB_INTER_RF:
      //this is done elsewhere
      ret = 0;
      break;
    case COULOMB_MMM1D:
      ret = mmm1d_coulomb_pair_energy(p1,p2,d, dist2,dist);
      break;
    default:
      ret = 0;
    }
    virials.coulomb[0] += ret;
  }
  /* stress tensor part */
  if (coulomb.method == COULOMB_DH) {
    int i;
    for (i = 0; i < 3; i++)
      force[i] = 0;
    
    add_dh_coulomb_pair_force(p1,p2,d,dist, force);
    for(k=0;k<3;k++)
      for(l=0;l<3;l++)
	p_tensor.coulomb[k*3 + l] += force[k]*d[l];
  }
  if (coulomb.method == COULOMB_RF) {
    int i;
    for (i = 0; i < 3; i++)
      force[i] = 0;
    
    add_rf_coulomb_pair_force(p1,p2,d,dist, force);
    for(k=0;k<3;k++)
      for(l=0;l<3;l++)
	p_tensor.coulomb[k*3 + l] += force[k]*d[l];
  }
  if (coulomb.method == COULOMB_INTER_RF) {
     //this is done elsewhere
  }
#endif /*ifdef ELECTROSTATICS */

#ifdef MAGNETOSTATICS
  /* real space magnetic dipole-dipole */
  if (coulomb.Dmethod != DIPOLAR_NONE) {
    switch (coulomb.Dmethod) {
#ifdef ELP3M
    case  DIPOLAR_P3M:
        /*ret = p3m_dipolar_pair_energy(p1->p.q*p2->p.q,d,dist2,dist);*/
	fprintf(stderr,"virials Not working for dipoles P3M .... pressure.h \n");
	ret=0;
        break; 
#endif
#ifdef DAWAANR
    case  DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA:
        /*ret = p3m_dipolar_pair_energy(p1->p.q*p2->p.q,d,dist2,dist);*/
	fprintf(stderr,"virials Not working for dipoles DAWAANR .... pressure.h \n");
	ret=0;
        break; 
#endif
#ifdef MAGNETIC_DIPOLAR_DIRECT_SUM
    case  DIPOLAR_DS:
        /*ret = p3m_dipolar_pair_energy(p1->p.q*p2->p.q,d,dist2,dist);*/
	fprintf(stderr,"virials Not working for dipoles MAGNETIC DIRECT SUM .... pressure.h \n");
	ret=0;
        break; 
#endif

      default:
      ret = 0;
    }
    virials.dipolar[0] += ret;
  }  
#endif /*ifdef MAGNETOSTATICS */
}

MDINLINE void calc_bonded_force(Particle *p1, Particle *p2, Bonded_ia_parameters *iaparams, int *i, double dx[3], double force[3]) {
#ifdef TABULATED
  char* errtxt;
#endif

#ifdef ADRESS
  int j;
  double force_weight=1;
  //adress_bonded_force_weight(p1);
  if (force_weight<ROUND_ERROR_PREC) return;
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
    case BONDED_IA_ANGLE:
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
	  errtxt = runtime_error(128 + TCL_INTEGER_SPACE);
	  ERROR_SPRINTF(errtxt,"{081 calc_bonded_force: tabulated bond type of atom %d unknown\n", p1->p.identity);
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
#ifdef ADRESS
    if((get_mol_com_particle(p1))->p.identity == (get_mol_com_particle(p2))->p.identity)
      force_weight = 1.0;
    else
      force_weight=adress_non_bonded_force_weight(p1,p2);
    for (j=0;j<3;j++){
      force[j]*=force_weight;
    }
#endif
}


/* calc_three_body_bonded_forces is called by add_three_body_bonded_stress. This
   routine is only entered for angular potentials. */
MDINLINE void calc_three_body_bonded_forces(Particle *p1, Particle *p2, Particle *p3,
              Bonded_ia_parameters *iaparams, double force1[3], double force2[3], double force3[3]) {

#ifdef TABULATED
  char* errtxt;
#endif

  switch(iaparams->type) {
  case BONDED_IA_ANGLE:
    // p1 is *p_mid, p2 is *p_left, p3 is *p_right
    calc_angle_3body_forces(p1, p2, p3, iaparams, force1, force2, force3);
    break;
#ifdef TABULATED
  case BONDED_IA_TABULATED:
    switch(iaparams->p.tab.type) {
      case TAB_BOND_ANGLE:
        // p1 is *p_mid, p2 is *p_left, p3 is *p_right
        calc_angle_3body_tabulated_forces(p1, p2, p3, iaparams, force1, force2, force3);
        break;
      default:
        errtxt = runtime_error(128 + TCL_INTEGER_SPACE);
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
MDINLINE void add_bonded_virials(Particle *p1)
{
  double dx[3], force[3];
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
      errtxt = runtime_error(128 + 2*TCL_INTEGER_SPACE);
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
MDINLINE void add_three_body_bonded_stress(Particle *p1) {
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

    if(type == BONDED_IA_ANGLE) {
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
        errtxt = runtime_error(128 + TCL_INTEGER_SPACE);
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
      errtxt = runtime_error(128 + TCL_INTEGER_SPACE);
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
MDINLINE void add_kinetic_virials(Particle *p1,int v_comp)
{
  int k, l;
  /* kinetic energy */
  if(v_comp)
    virials.data.e[0] += (SQR(p1->m.v[0] - p1->f.f[0]) + SQR(p1->m.v[1] - p1->f.f[1]) + SQR(p1->m.v[2] - p1->f.f[2]))*PMASS(*p1);
  else
    virials.data.e[0] += (SQR(p1->m.v[0]) + SQR(p1->m.v[1]) + SQR(p1->m.v[2]))*PMASS(*p1);

#ifdef ROTATION
  virials.data.e[0] += (SQR(p1->m.omega[0]) + SQR(p1->m.omega[1]) + SQR(p1->m.omega[2]))*SQR(time_step);
#endif

  /* ideal gas contribution (the rescaling of the velocities by '/=time_step' each will be done later) */
  for(k=0;k<3;k++)
    for(l=0;l<3;l++)
      p_tensor.data.e[k*3 + l] += (p1->m.v[k])*(p1->m.v[l])*PMASS(*p1);

}

/** implementation of 'analyze pressure'
    @param interp Tcl interpreter
    @param argc   arguments
    @param argv   arguments
    @param v_comp flag which enables (1) compensation of the velocities required
		  for deriving a pressure reflecting \ref nptiso_struct::p_inst
		  (hence it only works with domain decomposition); naturally it
		  therefore doesn't make sense to use it without NpT. */
int parse_and_print_pressure(Tcl_Interp *interp, int v_comp, int argc, char **argv);

/** Implementation of 'analyze bins' */
int parse_bins(Tcl_Interp *interp, int argc, char **argv);

/** implementation of 'analyze p_IK1' */
int parse_and_print_p_IK1(Tcl_Interp *interp, int argc, char **argv);

/** implementation of 'analyze stress_tensor' */
int parse_and_print_stress_tensor(Tcl_Interp *interp, int v_comp, int argc, char **argv);

/** implementation of 'analyse local_stress_tensor */
int local_stress_tensor_calc (DoubleList *TensorInBin, int bins[3], int periodic[3], double range_start[3], double range[3]);
int parse_local_stress_tensor(Tcl_Interp *interp, int argc, char **argv);

/*@}*/

#endif
