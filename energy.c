// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
#include "energy.h"

#include <mpi.h>
#include "parser.h"
#include "communication.h"
#include "grid.h"
#include "integrate.h"
#include "particle_data.h"
#include "interaction_data.h"
#include "cells.h"
#include "verlet.h"
#include "p3m.h"
#include "thermostat.h"
/* include the force files */
#include "lj.h"
#include "tab.h"
#ifdef LJCOS
  #include "ljcos.h"
#endif
#include "gb.h"
#include "fene.h"
#include "harmonic.h"
#include "subt_lj_harm.h"
#include "subt_lj_fene.h"
#include "angle.h"
#include "debye_hueckel.h"
#include "forces.h"
#include "constraint.h"

Observable_stat energy = {0, {NULL,0,0}, {NULL,0,0}, 0,0,0,0,0,0};

void calc_energy()
{
  Cell *cell;
  Particle *p, **pairs;
  Particle *p1, *p2;
  int i, j, k,  m, n, o, np, size;
  double d[3], dist2, dist;
  IA_parameters *ia_params;
  /* bonded interactions */
  int type_num;
  /* non bonded */ 
  int type1,type2;
  /* energy classification numbers */
  int s_bonded,s_non_bonded,s_coulomb,max;

  s_bonded = energy.n_pre;
  s_non_bonded = s_bonded + energy.n_bonded;
  s_coulomb = s_non_bonded + energy.n_non_bonded;
  max = s_coulomb + energy.n_coulomb;

  for(i=0;i<max;i++) {
    energy.node.e[i] = 0.0;
    energy.sum.e[i]  = 0.0;
  }

#ifdef ELECTROSTATICS
  switch (coulomb.method) {
  case COULOMB_NONE:
  case COULOMB_P3M:
  case COULOMB_DH: break;
  default:
    fprintf(stderr, "calc_energy: cannot calculate energy for coulomb method %d\n",
	    coulomb.method);
    errexit();
  }
#endif

  /* energy calculation loop. */
  INNER_CELLS_LOOP(m, n, o) {
    cell = CELL_PTR(m, n, o);
    p  = cell->pList.part;
    np = cell->pList.n;
    
    if(energy.ana_num < s_non_bonded ) {
      /* calculate bonded interactions and kinetic energy (loop local particles) */
      for(j = 0; j < np; j++) {
	p1 = &p[j];
	/* kinetic energy */
	energy.node.e[1] += SQR(p1->v[0]) + SQR(p1->v[1]) + SQR(p1->v[2]);
#ifdef ROTATION
        /* the rotational part is added to the total kinetic energy;
	   at the moment, we assume unit inertia tensor I=(1,1,1)  */
        energy.node.e[1] += (SQR(p1->omega[0]) + SQR(p1->omega[1]) + SQR(p1->omega[2]))*time_step*time_step;
#endif	
	/* bonded interaction energies */
	i=0;
	while(i<p1->bl.n) {
	  type_num = p1->bl.e[i];
	  switch(bonded_ia_params[type_num].type) {
	  case BONDED_IA_FENE:
	    energy.node.e[type_num + energy.n_pre] +=
	      fene_pair_energy(p1, checked_particle_ptr(p1->bl.e[i+1]), type_num);
	    i+=2; break;
	  case BONDED_IA_HARMONIC:
	    energy.node.e[type_num + energy.n_pre] +=
	      harmonic_pair_energy(p1, checked_particle_ptr(p1->bl.e[i+1]), type_num);
	    i+=2; break;
	  case BONDED_IA_SUBT_LJ_HARM:
	    energy.node.e[type_num + energy.n_pre] +=
	      subt_lj_harm_pair_energy(p1, checked_particle_ptr(p1->bl.e[i+1]), type_num);
	    i+=2; break; 
	  case BONDED_IA_SUBT_LJ_FENE:
	    energy.node.e[type_num + energy.n_pre] +=
	      subt_lj_fene_pair_energy(p1, checked_particle_ptr(p1->bl.e[i+1]), type_num);
	    i+=2; break; 
	  case BONDED_IA_ANGLE:
	    energy.node.e[type_num + energy.n_pre] +=
	      angle_energy(p1, checked_particle_ptr(p1->bl.e[i+1]),
			   checked_particle_ptr(p1->bl.e[i+2]), type_num);
	    i+=3; break;
	  default :
	    fprintf(stderr,"WARNING: Bonds of atom %d unknown\n",p1->r.identity);
	    i = p1->bl.n; 
	    break;
	  }
	}
#ifdef CONSTRAINTS
	/* constaint energies */
	for (i=0; i< n_constraints ; i++) {
    
	  type1 = p1->r.type;
	  type2 = (&constraints[i].part_rep)->r.type;
	  ia_params=get_ia_param(type1,type2);

	  if(ia_params->LJ_cut > 0. ) {
            type_num = s_non_bonded + ((2 * n_particle_types - 1 - type1) * type1) / 2  +  type2;
            energy.node.e[type_num] += add_constraints_energy(p1,i);
	  }
	}
#endif
      }
    }

    if(energy.ana_num == 0 || 
       (energy.ana_num >= s_non_bonded && energy.ana_num < max ) ) {
      /* calculate non bonded interactions (loop verlet lists of neighbors) */
      for (k = 0; k < cell->n_neighbors; k++) {
	pairs = cell->nList[k].vList.pair;  /* verlet list */
	np    = cell->nList[k].vList.n;     /* length of verlet list */
	
	/* verlet list loop */
	for(i=0; i<2*np; i+=2) {
	  p1 = pairs[i];                    /* pointer to particle 1 */
	  p2 = pairs[i+1];                  /* pointer to particle 2 */
	  ia_params = get_ia_param(p1->r.type,p2->r.type);
	  
	  if(p1->r.type > p2->r.type) { type1 = p2->r.type; type2 = p1->r.type; }
	  else { type2 = p2->r.type; type1 = p1->r.type; }
	  type_num = s_non_bonded + ((2 * n_particle_types - 1 - type1) * type1) / 2  +  type2 ;

	  /* distance calculation */
	  for(j=0; j<3; j++) d[j] = p1->r.p[j] - p2->r.p[j];
	  dist2 = SQR(d[0]) + SQR(d[1]) + SQR(d[2]);
	  dist  = sqrt(dist2);
	  
#ifdef TABULATED
	  /* Tabulated */
	  energy.node.e[type_num] += tabulated_pair_energy(p1,p2,ia_params,d,dist);
#endif
	  
	  /* lennard jones */
	  energy.node.e[type_num] += lj_pair_energy(p1,p2,ia_params,d,dist);
	  
#ifdef LJCOS
	  /* lennnard jones cosine */
	  energy.node.e[type_num] += ljcos_pair_energy(p1,p2,ia_params,d,dist);
#endif

#ifdef ROTATION	  
	  /* gay-berne */
	  energy.node.e[type_num] += gb_pair_energy(p1,p2,ia_params,d,dist);
#endif
 
#ifdef ELECTROSTATICS
	  /* real space coulomb */
	  if(coulomb.method==COULOMB_P3M) 
	    energy.node.e[s_coulomb+1] += p3m_coulomb_pair_energy(p1,p2,d,dist2,dist);
	  else if(coulomb.method==COULOMB_DH)
	    energy.node.e[s_coulomb] += dh_coulomb_pair_energy(p1,p2,dist);
#endif

	} 
      }
    }
  }

#ifdef ELECTROSTATICS
  /* calculate k-space part of electrostatic interaction. */ 
  if(coulomb.method==COULOMB_P3M && (energy.ana_num == 0 || energy.ana_num >= s_coulomb) ) {
    energy.node.e[s_coulomb+2] = P3M_calc_kspace_forces(0,1);
    energy.node.e[s_coulomb] = energy.node.e[s_coulomb+1]+energy.node.e[s_coulomb+2];
  }
#endif

  /* rescale kinetic energy */
  energy.node.e[1] /= (2.0*time_step*time_step);
  /* sum energies over nodes */
  if(energy.ana_num==0) size=energy.n;
  else {                
    size = 1;
    energy.node.e[0] = energy.node.e[energy.ana_num];
  }
  MPI_Reduce(energy.node.e, energy.sum.e, size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(energy.ana_num==0 && this_node==0) {
#ifdef ELECTROSTATICS
    if(coulomb.method==COULOMB_P3M) {
      for(i=1;i<energy.n-2;i++)
	energy.sum.e[0] += energy.sum.e[i];
    } else 
#endif
      for(i=1;i<energy.n;i++)
	energy.sum.e[0] += energy.sum.e[i];
  }
}

void init_energies()
{
  if (energy.init_status != 0 && ! interactions_changed)
    return;

  energy.n_pre        = 2;
  energy.n_bonded     = n_bonded_ia;
  energy.n_non_bonded = (n_particle_types*(n_particle_types+1))/2;
#ifdef ELECTROSTATICS
  if(coulomb.bjerrum > 0.0)         energy.n_coulomb =  1;
  if(coulomb.method==COULOMB_P3M) energy.n_coulomb += 2;
#endif
  energy.n = energy.n_pre+energy.n_bonded+energy.n_non_bonded+energy.n_coulomb;
  realloc_doublelist(&(energy.node),energy.n);
  realloc_doublelist(&(energy.sum),energy.n);

  energy.init_status=0;
}

void calc_energies()
{
  /* check integrator status */
  if (parameter_changed || interactions_changed || topology_changed || particle_changed) {
    mpi_integrate(0);
  }
  /* calculate energies (in parallel) */
  mpi_gather_stats(1, NULL);
}

/****************************************************************************************
 *                                 parser
 ****************************************************************************************/

int parse_and_print_energy(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze energy [{ fene <type_num> | harmonic <type_num> | subt_lj_harm <type_num> | subt_lj_fene <type_num> / lj <type1> <type2> | ljcos <type1> <type2> | gb <type1> <type2> | coulomb | kinetic | total }]' */
  char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE + 2];
  int i, j, p;
  int out_mode=0;

  if (n_total_particles == 0) {
    Tcl_AppendResult(interp, "(no particles)",
		     (char *)NULL);
    return (TCL_OK);
  }
 
  /* check energy initialization */
  init_energies();


  if(argc>0) {
    if      (ARG0_IS_S("kinetic")) energy.ana_num = 1;
    else if (ARG0_IS_S("bonded") || ARG0_IS_S("fene") ||
	     ARG0_IS_S("harmonic") || ARG0_IS_S("subt_lj_harm") ||
	     ARG0_IS_S("subt_lj_fene")) {
      if(argc<2) {
	Tcl_AppendResult(interp, "wrong # arguments for: analyze energy bonded <type_num>",
			 (char *)NULL);
	return (TCL_ERROR);
      }
      if(! ARG1_IS_I(i)) return (TCL_ERROR);
      if(i >= energy.n_bonded) {
	Tcl_AppendResult(interp, "bond type does not exist",
			 (char *)NULL);
	return (TCL_ERROR);
      }
      energy.ana_num = energy.n_pre+i;
    }
    else if (ARG0_IS_S("nonbonded") || ARG0_IS_S("lj") ||
	     ARG0_IS_S("lj-cos") || ARG0_IS_S("gb") || ARG0_IS_S("tabulated")){
      if(argc<3) {
	Tcl_AppendResult(interp, "wrong # arguments for: analyze energy nonbonded <type1> <type2>",
			 (char *)NULL);
	return (TCL_ERROR);
      }
      if(! ARG_IS_I(1, i)) return (TCL_ERROR);
      if(! ARG_IS_I(2, j)) return (TCL_ERROR);
      if(i >= n_particle_types || j >= n_particle_types) {
	Tcl_AppendResult(interp, "particle type does not exist",
			 (char *)NULL);
	return (TCL_ERROR);
      }
      energy.ana_num = energy.n_pre+energy.n_bonded + ((2 * n_particle_types - 1 - i) * i) / 2  +  j;
    }
    else if( ARG0_IS_S("coulomb")) {
#ifdef ELECTROSTATICS
      energy.ana_num = energy.n_pre+energy.n_bonded+energy.n_non_bonded;
#else
      Tcl_AppendResult(interp, "ELECTROSTATICS not compiled (see config.h)\n", (char *)NULL);
#endif
    }
    else if (ARG0_IS_S("total")) {
      energy.ana_num = 0; 
      out_mode = 1;
    }
    else {
      Tcl_AppendResult(interp, "unknown feature of: analyze energy",
		       (char *)NULL);
      return (TCL_ERROR);
    }
  }
  else
    energy.ana_num=0;

  calc_energies();

  /* print out results */
  if(energy.ana_num > 0) {
    Tcl_PrintDouble(interp, energy.sum.e[0], buffer);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
  }
  else {
    if(out_mode == 1) {
      Tcl_PrintDouble(interp, energy.sum.e[0], buffer);
      Tcl_AppendResult(interp, buffer, (char *)NULL);
    }
    else {
      Tcl_PrintDouble(interp, energy.sum.e[0], buffer);
      Tcl_AppendResult(interp, "{ energy ", buffer, " } ", (char *)NULL);
      Tcl_PrintDouble(interp, energy.sum.e[1], buffer);
      Tcl_AppendResult(interp, "{ kinetic ", buffer, " } ", (char *)NULL);
      for(i=0;i<n_bonded_ia;i++) {
	sprintf(buffer, "%d ", i);
	Tcl_AppendResult(interp, "{ ", buffer, (char *)NULL);
	Tcl_PrintDouble(interp, energy.sum.e[energy.n_pre+i], buffer);
	Tcl_AppendResult(interp,
			 get_name_of_bonded_ia(bonded_ia_params[i].type),
			 " ", buffer, " } ", (char *) NULL);
      }
      p = energy.n_pre+energy.n_bonded;
      for (i = 0; i < n_particle_types; i++)
	for (j = i; j < n_particle_types; j++) {
	  if (checkIfParticlesInteract(i, j)) {
	    sprintf(buffer, "%d ", i);
	    Tcl_AppendResult(interp, "{ ", buffer, (char *)NULL);
	    sprintf(buffer, "%d ", j);
	    Tcl_AppendResult(interp, " ", buffer, (char *)NULL);
	    Tcl_PrintDouble(interp, energy.sum.e[p], buffer);
	    Tcl_AppendResult(interp, "nonbonded ", buffer, " } ", (char *)NULL);	    
	  }
	  p++;
	}
#ifdef ELECTROSTATICS
      if(coulomb.bjerrum > 0.0) {
	Tcl_PrintDouble(interp, energy.sum.e[p], buffer);
	Tcl_AppendResult(interp, "{ coulomb ", buffer, (char *)NULL);
	if(coulomb.method == COULOMB_P3M) {
	  Tcl_PrintDouble(interp, energy.sum.e[p+1], buffer);
	  Tcl_AppendResult(interp, " ", buffer, (char *)NULL);
	  Tcl_PrintDouble(interp, energy.sum.e[p+2], buffer);
	  Tcl_AppendResult(interp, " ", buffer, (char *)NULL);
	}
	Tcl_AppendResult(interp, " }", (char *)NULL);
      }
#endif
    }
  }

  energy.init_status=1;
  return (TCL_OK);
}
