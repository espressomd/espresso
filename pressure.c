// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
#include "pressure.h"
#include "parser.h"
#include <string.h>
#include <mpi.h>
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
#include "ljcos.h"
#include "gb.h"
#include "fene.h"
#include "harmonic.h"
#include "angle.h"
#include "debye_hueckel.h"
#include "forces.h"
#include "constraint.h"


Observable_stat virials = {0, {NULL,0,0}, {NULL,0,0}, 0,0,0,0,0,0};
Observable_stat p_tensor= {0, {NULL,0,0}, {NULL,0,0}, 0,0,0,0,0,0};


/*******************/
/* Scalar Pressure */
/*******************/

int parse_and_print_pressure(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze pressure [{ <bonded> [<type_num>] | <nonbonded> [<type1> <type2>] | coulomb | ideal | total[s] }]' */
  /**************************************************************************************************************************/
  char buffer[3*TCL_DOUBLE_SPACE + 256];
  double *buf;
  int i, j, p, k;

  if (n_total_particles == 0) { Tcl_AppendResult(interp, "(no particles)",(char *)NULL); return (TCL_OK); }

  init_virials();

  if(argc == 0)
    virials.ana_num=0;
  else {
    if     (ARG0_IS_S("totals")) virials.ana_num=0;
    else if(ARG0_IS_S("ideal")) virials.ana_num=1;
    else if(ARG0_IS_S("bonded") || ARG0_IS_S("fene") ||
	    ARG0_IS_S("harmonic")) {
      if(argc<2) { virials.ana_num=0; }
      else {
	if(!ARG1_IS_I(i)) return (TCL_ERROR);
	if(i >= virials.n_bonded) { 
	  Tcl_AppendResult(interp,"bond type does not exist!",(char *)NULL);
	  return (TCL_ERROR);
	}
      }
      virials.ana_num = virials.n_pre+i;
    }
    else if(ARG0_IS_S("nonbonded") || ARG0_IS_S("lj") || ARG0_IS_S("lj-cos") || ARG0_IS_S("gb")) {
      if(argc<3) { virials.ana_num=0; }
      else {
	if(!ARG_IS_I(1, i)) return (TCL_ERROR);
	if(!ARG_IS_I(2, j)) return (TCL_ERROR);
	if(i >= n_particle_types || j >= n_particle_types) {
	  Tcl_AppendResult(interp, "particle type does not exist",
			   (char *)NULL);
	  return (TCL_ERROR);
	}
	virials.ana_num = virials.n_pre+virials.n_bonded + j -i;
	while(i>0) {
	  virials.ana_num += n_particle_types - (i-1); i--;
	}
      }
    }
    else if(ARG0_IS_S("coulomb")) {
#ifdef ELECTROSTATICS
      virials.ana_num = virials.n_pre+virials.n_bonded+virials.n_non_bonded; 
#else
      Tcl_AppendResult(interp, "ELECTROSTATICS not compiled (see config.h)\n", (char *)NULL);
#endif
    }
    else {
      Tcl_AppendResult(interp, "unknown feature of analyze pressure",(char *)NULL); return (TCL_ERROR);
    }
  }

  calc_pressure();

  buf = malloc(6*sizeof(double)); buf[0]=buf[1]=buf[2]=buf[3]=buf[4]=buf[5] = 0.0;
  if(argc > 0) {
    if(ARG0_IS_S("total")) {
      Tcl_PrintDouble(interp, virials.sum.e[virials.ana_num], buffer); 
      Tcl_AppendResult(interp, buffer, (char *)NULL); }
    else if(ARG0_IS_S("totals")) {
      sprintf(buffer,"%f %f",virials.sum.e[virials.ana_num],virials.node.e[virials.ana_num]); 
      Tcl_AppendResult(interp, buffer, (char *)NULL); }
    else if(ARG0_IS_S("ideal")) {
      Tcl_PrintDouble(interp, virials.sum.e[virials.ana_num], buffer); 
      Tcl_AppendResult(interp, buffer, (char *)NULL); }
    else if(virials.ana_num > 0) {  /* this covers case 'coulomb' as well */
      sprintf(buffer,"%f %f",virials.sum.e[virials.ana_num],virials.node.e[virials.ana_num]);
      Tcl_AppendResult(interp, buffer, (char *)NULL); }
    else if(ARG0_IS_S("bonded")) {
      for(i=0;i<n_bonded_ia;i++) {
	buf[0] += virials.sum.e[virials.n_pre+i];
	buf[1] += virials.node.e[virials.n_pre+i]; 
      }
      sprintf(buffer,"%f %f",buf[0],buf[1]); Tcl_AppendResult(interp,buffer,(char *)NULL); }
    else if(ARG0_IS_S("fene")) {
      for(i=0;i<n_bonded_ia;i++)
	if(bonded_ia_params[i].type == BONDED_IA_FENE) { buf[0] += virials.sum.e[virials.n_pre+i]; buf[1] += virials.node.e[virials.n_pre+i]; }
      sprintf(buffer,"%f %f",buf[0],buf[1]); Tcl_AppendResult(interp,buffer,(char *)NULL); }
    else if(ARG0_IS_S("harmonic")) {
      for(i=0;i<n_bonded_ia;i++)
	if(bonded_ia_params[i].type == BONDED_IA_HARMONIC) { buf[4] += virials.sum.e[virials.n_pre+i]; buf[5] += virials.node.e[virials.n_pre+i]; }
      sprintf(buffer,"%f %f",buf[4],buf[5]); Tcl_AppendResult(interp,buffer,(char *)NULL); }
    else if(ARG0_IS_S("nonbonded") || ARG0_IS_S("lj")) {
      p = virials.n_pre+virials.n_bonded;
      if(argc == 2) Tcl_GetInt(interp, argv[1], &k); else k=-1;
      for (i = 0; i < n_particle_types; i++) 
	if((k==-1) || (k==i))
	  for (j = i; j < n_particle_types; j++) {
	    if(checkIfParticlesInteract(i,j)) { buf[0] += virials.sum.e[p]; buf[1] += virials.node.e[p]; }
	    p++;
	  }
      sprintf(buffer,"%f %f",buf[0],buf[1]); Tcl_AppendResult(interp,buffer,(char *)NULL); }
    else { Tcl_AppendResult(interp, "unknown feature of analyze pressure",(char *)NULL); return (TCL_ERROR); } }
  else {
    int buf0, buf2, buf4;
    sprintf(buffer,"%f %f { ideal %f } { ",virials.sum.e[0],virials.node.e[0],virials.sum.e[1]); 
    Tcl_AppendResult(interp, buffer, (char *)NULL);
    buf0 = buf2 = buf4 = 0;
    for(i=0;i<n_bonded_ia;i++) {
      switch (bonded_ia_params[i].type) {
      case BONDED_IA_FENE:
	buf[0] += virials.sum.e[virials.n_pre+i]; buf[1] += virials.node.e[virials.n_pre+i]; buf0 = 1; break;
      case BONDED_IA_ANGLE:
	buf[2] += virials.sum.e[virials.n_pre+i]; buf[3] += virials.node.e[virials.n_pre+i]; buf2 = 1; break;
      case BONDED_IA_HARMONIC:
	buf[4] += virials.sum.e[virials.n_pre+i]; buf[5] += virials.node.e[virials.n_pre+i]; buf4 = 1; break;
      default: break; }
    }
    if(buf0 != 0) { sprintf(buffer,"{ FENE %f %f } ",buf[0],buf[1]); Tcl_AppendResult(interp,buffer,(char *)NULL); }
    if(buf2 != 0) { sprintf(buffer,"{ angle %f %f } ",buf[2],buf[3]); Tcl_AppendResult(interp,buffer,(char *)NULL); }
    if(buf4 != 0) { sprintf(buffer,"{ harmonic %f %f } ",buf[4],buf[5]); Tcl_AppendResult(interp,buffer,(char *)NULL); }
    buf[0] = buf[1] = 0.0; buf0 = 0; Tcl_AppendResult(interp, "}  { ", (char *)NULL);
    p = virials.n_pre+virials.n_bonded;
    for (i = 0; i < n_particle_types; i++) {
      for (j = i; j < n_particle_types; j++) {
	if (checkIfParticlesInteract(i, j)) { buf[0] += virials.sum.e[p];  buf[1] += virials.node.e[p];  buf0 = 1; }
	p++;
      }
    }
    if(buf0 != 0) { sprintf(buffer, "lj %f %f }  ",buf[0],buf[1]); Tcl_AppendResult(interp, buffer, (char *)NULL); }
    else Tcl_AppendResult(interp, "}  ", (char *)NULL);
#ifdef ELECTROSTATICS
    if(coulomb.bjerrum > 0.0) {
      sprintf(buffer, "{ coulomb %f %f } ",virials.sum.e[p],virials.node.e[p]);
      Tcl_AppendResult(interp, buffer,  (char *)NULL);
    }
#endif
  }
  virials.init_status=1;
  free(buf);
  return (TCL_OK);
}


/* Initialize the Virials */
/**************************/
void init_virials() {
  if (virials.init_status != 0 && ! interactions_changed)
    return;
  virials.n_pre        = 2;
  virials.n_bonded     = n_bonded_ia;
  virials.n_non_bonded = (n_particle_types*(n_particle_types+1))/2;

#ifdef ELECTROSTATICS
  switch (coulomb.method) {
  case COULOMB_NONE: virials.n_coulomb = 0; break;
  case COULOMB_P3M: virials.n_coulomb = 3; break;
  case COULOMB_DH:  virials.n_coulomb = 1; break;
  default:
    fprintf(stderr, "%d: init_virials: cannot calculate virials for coulomb method %d\n",
	    this_node, coulomb.method);
    errexit();
  }
#endif
  virials.n = virials.n_pre+virials.n_bonded+virials.n_non_bonded+virials.n_coulomb;
  realloc_doublelist(&(virials.node),virials.n);
  realloc_doublelist(&(virials.sum),virials.n);
  virials.init_status = 0;
}


/* Derive the Virials */
/**********************/
void calc_virials() {
  Cell *cell;
  Particle *p, **pairs, *p1,*p2,*p3;
  IA_parameters *ia_params;
  int i,j,k,  m,n,o, np, size;
  double d[3], dist,dist2, f1[3],f2[3],f3[3];
  /* bonded interactions */
  int type_num;
  /* non bonded */ 
  int type1,type2;
  /* virial classification numbers */
  int v_bonded,v_non_bonded,v_coulomb,max;

  v_bonded = virials.n_pre;
  v_non_bonded = v_bonded + virials.n_bonded;
  v_coulomb = v_non_bonded + virials.n_non_bonded;
  max = v_coulomb + virials.n_coulomb;
  for(i=0;i<max;i++) {
    virials.node.e[i] = 0.0;
    virials.sum.e[i]  = 0.0;
  }

  /* virials calculation loop. */
  INNER_CELLS_LOOP(m, n, o) {
    cell = CELL_PTR(m, n, o);
    p  = cell->pList.part;
    np = cell->pList.n;

    if(virials.ana_num < v_non_bonded ) {
      /* bonded interactions (loop local particles) */
      for(j = 0; j < np; j++) {
	p1 = &p[j]; i=0;
	/* kinetic energy */
	virials.node.e[1] += SQR(p1->v[0]) + SQR(p1->v[1]) + SQR(p1->v[2]);
	/* bonded interaction virials */
	while(i < p1->bl.n) {
	  p2 = checked_particle_ptr(p1->bl.e[i+1]);
	  f1[0] = p1->f[0]; f1[1] = p1->f[1]; f1[2] = p1->f[2];
	  f2[0] = p2->f[0]; f2[1] = p2->f[1]; f2[2] = p2->f[2];
	  get_mi_vector(d, p1->r.p, p2->r.p);
	  type_num = p1->bl.e[i];
	  switch(bonded_ia_params[type_num].type) {
	  case BONDED_IA_FENE:
	    add_fene_pair_force(p1,p2,type_num);
	    for(k=0;k<3;k++) { p1->f[k] -= (f1[k] = p1->f[k] - f1[k]); p2->f[k] -= (f2[k] = p2->f[k] - f2[k]); }
	    virials.node.e[type_num + virials.n_pre] += d[0]*f1[0] + d[1]*f1[1] + d[2]*f1[2];
	    i+=2; break;
	  case BONDED_IA_HARMONIC:
	    add_harmonic_pair_force(p1,p2,type_num);
	    for(k=0;k<3;k++) { p1->f[k] -= (f1[k] = p1->f[k] - f1[k]); p2->f[k] -= (f2[k] = p2->f[k] - f2[k]); }
	    virials.node.e[type_num + virials.n_pre] += d[0]*f1[0] + d[1]*f1[1] + d[2]*f1[2];
	    i+=2; break;
	  case BONDED_IA_ANGLE:
	    p3 = checked_particle_ptr(p1->bl.e[i+2]);
	    f3[0] = p3->f[0]; f3[1] = p3->f[1]; f3[2] = p3->f[2];
	    add_angle_force(p1,p2,p3,type_num);
	    for(k=0;k<3;k++) { p1->f[k] -= (f1[k] = p1->f[k] - f1[k]); p2->f[k] -= (f2[k] = p2->f[k] - f2[k]); }
	    virials.node.e[type_num + virials.n_pre] += -d[0]*f2[0] - d[1]*f2[1] - d[2]*f2[2];
	    for(k=0;k<3;k++) { p3->f[k] -= (f3[k] = p3->f[k] - f3[k]); }
	    get_mi_vector(d, p1->r.p, p3->r.p);
	    virials.node.e[type_num + virials.n_pre] += -d[0]*f3[0] - d[1]*f3[1] - d[2]*f3[2];
	    i+=3; break;
	  default :
	    fprintf(stderr,"WARNING: Bonds of atom %d unknown\n",p1->r.identity);
	    i = p1->bl.n; break;
	  }
	}
      }
    }

    if(virials.ana_num == 0 || (virials.ana_num >= v_non_bonded && virials.ana_num < max ) ) {
      /* calculate non bonded interactions (loop verlet lists of neighbors) */
      for (k = 0; k < cell->n_neighbors; k++) {
	pairs = cell->nList[k].vList.pair;  /* verlet list */
	np    = cell->nList[k].vList.n;     /* length of verlet list */
	
	/* verlet list loop */
	for(i=0; i<2*np; i+=2) {
	  p1 = pairs[i];                    /* pointer to particle 1 */
	  p2 = pairs[i+1];                  /* pointer to particle 2 */
	  ia_params = get_ia_param(p1->r.type,p2->r.type);

	  /* derive index 'type_num' */
	  if(p1->r.type > p2->r.type) { type1 = p2->r.type; type2 = p1->r.type; } else { type2 = p2->r.type; type1 = p1->r.type; }
	  type_num = v_non_bonded + ((2 * n_particle_types - 1 - type1) * type1) / 2  +  type2;

	  /* distance calculation */
	  for(j=0; j<3; j++) d[j] = p1->r.p[j] - p2->r.p[j];
	  dist2 = SQR(d[0]) + SQR(d[1]) + SQR(d[2]);
	  dist  = sqrt(dist2);
	  
	  /* lennnard jones */
	  for(j=0;j<3;j++) { f1[j] = p1->f[j]; f2[j] = p2->f[j]; }
	  add_lj_pair_force(p1,p2,ia_params,d,dist);
	  add_ljcos_pair_force(p1,p2,ia_params,d,dist);
#ifdef ROTATION  
	  add_gb_pair_force(p1,p2,ia_params,d,dist);
#endif
	  for(j=0;j<3;j++) { p1->f[j] -= (f1[j] = p1->f[j] - f1[j]); p2->f[j] -= (f2[j] = p2->f[j] - f2[j]); }
	  virials.node.e[type_num] += d[0]*f1[0] + d[1]*f1[1] + d[2]*f1[2];
	  
#ifdef ELECTROSTATICS
	  /* real space coulomb */
	  if(coulomb.method==COULOMB_P3M) 
	    virials.node.e[v_coulomb+1] += p3m_coulomb_pair_energy(p1,p2,d,dist2,dist);
	  else if(coulomb.method==COULOMB_DH)
	    virials.node.e[v_coulomb] += dh_coulomb_pair_energy(p1,p2,dist);
#endif
	} 
      }
    }
  }
#ifdef ELECTROSTATICS
  /* calculate k-space part of electrostatic interaction. */ 
  if(coulomb.method==COULOMB_P3M && (virials.ana_num == 0 || virials.ana_num >= v_coulomb) ) {
    virials.node.e[v_coulomb+2] = P3M_calc_kspace_forces(0,1);
    virials.node.e[v_coulomb] = virials.node.e[v_coulomb+1]+virials.node.e[v_coulomb+2];
  }
#endif
  /* rescale kinetic energy  &  sum virials over nodes */
  virials.node.e[1] /= (2.0*time_step*time_step);
  size=virials.n;

  MPI_Reduce(virials.node.e, virials.sum.e, size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(this_node==0) {
#ifdef ELECTROSTATICS
    if(coulomb.method==COULOMB_P3M) {
      for(i=virials.n_pre;i<virials.n-2;i++) 
	virials.sum.e[0] += virials.sum.e[i];
    } else 
#endif
      for(i=virials.n_pre;i<virials.n;i++)   
	virials.sum.e[0] += virials.sum.e[i];
    
  }
}
void calc_bins_sphere(int *new_bin,int *elements,double *volumes,double r_min,double r_max,int r_bins, double *center)
{
  int i,j,counter=0,*bins;
  double d[3],dist;
  double d_bin;
  d_bin = (r_max-r_min)/r_bins;
  bins = malloc(n_total_particles*sizeof(int));
  for(i=0; i<n_total_particles; i++) {
    get_mi_vector(d, center, partCfg[i].r.p);
    dist = sqrt(sqrlen(d))-r_min;
    bins[i] = floor(dist/d_bin);
    if(dist > r_max-r_min) bins[i] = -1;
  }
  for(i=0;i<r_bins;i++){
    new_bin[i] = 0;
    volumes[i] = 4./3.*3.1415926536*(pow(r_min+(i+1.)*d_bin,3.)-pow(r_min+i*d_bin,3.));
    for(j=0;j<n_total_particles;j++){
      if(bins[j]==i){
	elements[counter++] = partCfg[j].r.identity;
	new_bin[i] += 1;
      }
    }
    printf("vol %d %le\n",i,volumes[i]);
  }
}



int parse_bins(Tcl_Interp *interp, int argc, char **argv)
{
  /* bin particles for pressure calculation */
  /******************************************/
  /** Computes bins for pressure calculations, gives back lists
      with particles and bin volumes for each bin in spherical geometry**/
  char buffer[1000*TCL_INTEGER_SPACE];
  double r_min=0, r_max=-1.0, center[3];
  int r_bins=-1, i,j,k;
  int *new_bin;
  int *elements;
  double *volumes;

  if (ARG0_IS_S("sphere")) { 
    argc--; argv++;
    if(argc < 6) { Tcl_AppendResult(interp,"Too few arguments! Usage: 'analyze bins sphere <r_min> <r_max> <r_bins> <center1> <center2> <center3> '",(char *)NULL); return (TCL_ERROR); }
    if( argc>0 ) { if (!ARG0_IS_D(r_min)) return (TCL_ERROR); argc--; argv++; }
    if( argc>0 ) { if (!ARG0_IS_D(r_max)) return (TCL_ERROR); argc--; argv++; }
    if( argc>0 ) { if (!ARG0_IS_I(r_bins)) return (TCL_ERROR); argc--; argv++; }
    if( argc>0 ) { if (!ARG0_IS_D(center[0])) return (TCL_ERROR); argc--; argv++; }
    if( argc>0 ) { if (!ARG0_IS_D(center[1])) return (TCL_ERROR); argc--; argv++; }
    if( argc>0 ) { if (!ARG0_IS_D(center[2])) return (TCL_ERROR); argc--; argv++; }
    
    
    /* if not given use default */
    if(r_max == -1.0) r_max = min_box_l/2.0;
    if(r_bins == -1) r_bins = n_total_particles / 20;
    
    /* give back what you do */
    
    elements = malloc(n_total_particles*sizeof(int));
    new_bin = malloc(r_bins*sizeof(int));
    volumes = malloc(r_bins*sizeof(double));
    updatePartCfg();
    calc_bins_sphere(new_bin,elements,volumes, r_min, r_max, r_bins, center);
    /* append result */
    {
      Tcl_AppendResult(interp, " { ", (char *)NULL);
      sprintf(buffer,"%le",volumes[0]);
      Tcl_AppendResult(interp, buffer, (char *)NULL);
      Tcl_AppendResult(interp, " { ", (char *)NULL);
      /* i->particles, j->bin, k->particle/bin */ 
      for(i=0,j=0,k=0;i<n_total_particles;i++,k++){
	if(k==new_bin[j] || new_bin[j] == 0){
	  /* if all bins are full, rest of particles are outside r_min/r_max */
	  k=0,j++;
	  if(j==r_bins) break;
	  Tcl_AppendResult(interp, "} } { ", (char *)NULL);
	  sprintf(buffer,"%le",volumes[j]);
	  Tcl_AppendResult(interp, buffer, (char *)NULL);
	  Tcl_AppendResult(interp, " { ", (char *)NULL);
	}
	sprintf(buffer,"%d ",elements[i]);
	Tcl_AppendResult(interp, buffer, (char *)NULL);
      }
      Tcl_AppendResult(interp, "}\n", (char *)NULL);
    }
    free(new_bin);
    free(elements);
    return (TCL_OK);
  }
  else {
    sprintf(buffer, "The feature 'analyze bins %s' you requested is not yet implemented.",argv[0]);
    Tcl_AppendResult(interp,buffer,(char *)NULL);
    return (TCL_ERROR);
  }
  return (TCL_ERROR);
}

/* Calculate the pressure from the virials */
/*******************************************/
void calc_pressure() {
  /** Remarks:
      <ul><li> The ideal gas pressure P_ig is assumed to be the pressure which the system 
      would have if all interactions had been switched off.
      <li> This routine does not work for the forrest of rods.
      <li> Until now it only works for hydrophilic LJ.
      </ul>  */
  double p_total=0.0, volume;
  int    i,j,p;

  if (parameter_changed || interactions_changed || topology_changed || particle_changed) {
    mpi_integrate(0);
  }

  /* Derive Virials for all pairwise contribution + electrostatics */
  mpi_gather_stats(2, NULL); 
  virials.node.e[1] = virials.sum.e[1];

  /* Ideal Gas */
  volume = box_l[0]*box_l[1]*box_l[2];
  virials.sum.e[1] = n_total_particles*temperature/volume;
  p_total = virials.sum.e[1];

  /* Contribution of bonded interactions */
  for(i=0; i<n_bonded_ia; i++) {
    virials.sum.e[virials.n_pre+i] /= 3.0*volume;
    virials.node.e[virials.n_pre+i] = SQR(virials.sum.e[virials.n_pre+i]);
    p_total += virials.sum.e[virials.n_pre+i];
  }

  /* Contribution of non-bonded interactions */
  p = virials.n_pre+virials.n_bonded;
  for (i = 0; i < n_particle_types; i++) {
    for (j = i; j < n_particle_types; j++) {
      if (checkIfParticlesInteract(i, j)) {
	virials.sum.e[p]  /= 3.0*volume;
	virials.node.e[p]  = SQR(virials.sum.e[p]);
	p_total += virials.sum.e[p]; }
      else {
	virials.sum.e[p] = virials.node.e[p] = 0.0; 
      }
      p++;
    }
  }


#ifdef ELECTROSTATICS
  /* Contribution of electrostatics (if any) */
  if(coulomb.bjerrum > 0.0) {
    virials.sum.e[p]  /= 3.0*volume;
    virials.node.e[p]  = SQR(virials.sum.e[p]);
    p_total += virials.sum.e[p]; 
  }
#endif

  /* Check and return */
  virials.sum.e[0] = virials.sum.e[0]/(3.0*volume) + virials.sum.e[1];
  if( fabs(p_total - virials.sum.e[0]) > ROUND_ERROR_PREC*SQR(n_total_particles)*fabs(p_total)) {
    fprintf(stderr, "Sanity check failed: Derived pressure (%e) differs from given one (%e) by %e!\n",p_total,virials.sum.e[0],virials.sum.e[0]-p_total);
    errexit();
  }
  virials.node.e[0] = SQR(p_total);
}



/**********************/
/* Tensorial Pressure */
/**********************/

int parse_and_print_p_IK1(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze p_IK1 <bin_volume> { <ind_list> } <all>' */
  /*****************************************************/
  char buffer[9*TCL_DOUBLE_SPACE + 256];
  int i,j,p, flag=0;
  double volume;
  IntList p1;

  if (n_total_particles == 0) { Tcl_AppendResult(interp, "(no particles)",(char *)NULL); return (TCL_OK); }
  init_p_tensor();
  init_intlist(&p1);

  if(argc < 3) { Tcl_AppendResult(interp,"Too few arguments! Usage: 'analyze p_IK1 <bin_volume> { <ind_list> } <all>'",(char *)NULL); return (TCL_ERROR); }
  if ((!ARG0_IS_D(volume)) || (!ARG1_IS_INTLIST(p1)) || (!ARG_IS_I(2, flag))) { 
    Tcl_ResetResult(interp); Tcl_AppendResult(interp,"usage: 'analyze p_IK1 <bin_volume> { <ind_list> } <all>'",(char *)NULL); return (TCL_ERROR); 
  }

  p_tensor.ana_num=0;
  calc_p_tensor(volume,&p1,flag);

  sprintf(buffer,"%f %f { total ",p_tensor.node.e[0],p_tensor.node.e[1]); Tcl_AppendResult(interp, buffer, (char *)NULL);
  for(j=0; j<9; j++) { sprintf(buffer,"%f ",p_tensor.sum.e[j]); Tcl_AppendResult(interp, buffer, (char *)NULL); }
  Tcl_AppendResult(interp, "} ", (char *)NULL); 

  Tcl_AppendResult(interp, "{ ideal ", (char *)NULL);
  for(j=0; j<9; j++) { sprintf(buffer,"%f ",p_tensor.sum.e[9+j]); Tcl_AppendResult(interp, buffer, (char *)NULL); }
  Tcl_AppendResult(interp, "} ", (char *)NULL); 

  p = p_tensor.n_pre;
  Tcl_AppendResult(interp, "{ bonded ", (char *)NULL); 
  for(j=0; j<9; j++) { sprintf(buffer,"%f ",p_tensor.sum.e[p+j]); Tcl_AppendResult(interp, buffer, (char *)NULL); }
  for(i=0;i<n_bonded_ia;i++) {
    switch (bonded_ia_params[i].type) {
    case BONDED_IA_FENE:
      Tcl_AppendResult(interp, "{ FENE ", (char *)NULL);
      for(j=0; j<9; j++) { sprintf(buffer,"%f ",p_tensor.node.e[p+i*9+j]); Tcl_AppendResult(interp, buffer, (char *)NULL); }
      Tcl_AppendResult(interp, "} ", (char *)NULL);
      break;
    case BONDED_IA_ANGLE:
      Tcl_AppendResult(interp, "{ ANGLE ", (char *)NULL);
      for(j=0; j<9; j++) { sprintf(buffer,"%f ",p_tensor.node.e[p+i*9+j]); Tcl_AppendResult(interp, buffer, (char *)NULL); }
      Tcl_AppendResult(interp, "} ", (char *)NULL);
      break;
    case BONDED_IA_HARMONIC:
      Tcl_AppendResult(interp, "{ HARMONIC ", (char *)NULL);
      for(j=0; j<9; j++) { sprintf(buffer,"%f ",p_tensor.node.e[p+i*9+j]); Tcl_AppendResult(interp, buffer, (char *)NULL); }
      Tcl_AppendResult(interp, "} ", (char *)NULL);
      break;
    default: break; }
  }
  Tcl_AppendResult(interp, "} ", (char *)NULL); 

  p = p_tensor.n_pre+p_tensor.n_bonded;
  Tcl_AppendResult(interp, "{ nonbonded ", (char *)NULL); 
  for(j=0; j<9; j++) { sprintf(buffer,"%f ",p_tensor.sum.e[p+j]); Tcl_AppendResult(interp, buffer, (char *)NULL); }
  Tcl_AppendResult(interp, "} ", (char *)NULL); 

#ifdef ELECTROSTATICS
  p = p_tensor.n_pre+p_tensor.n_bonded+p_tensor.n_non_bonded;
  if(coulomb.bjerrum > 0.0) {
    Tcl_AppendResult(interp, "{ coulomb ", (char *)NULL); 
    for(j=0; j<9; j++) { sprintf(buffer,"%f ",p_tensor.sum.e[p+j]); Tcl_AppendResult(interp, buffer, (char *)NULL); }
    Tcl_AppendResult(interp, "} ", (char *)NULL); 
  }
#endif

  p_tensor.init_status=1;
  return (TCL_OK);
}


/* Initialize the p_tensor */
/***************************/
void init_p_tensor() {
  if (p_tensor.init_status != 0 && ! interactions_changed)
    return;
  p_tensor.n_pre        = 9+9;
  p_tensor.n_bonded     = 9*n_bonded_ia;
  p_tensor.n_non_bonded = 9*(n_particle_types*(n_particle_types+1))/2;

#ifdef ELECTROSTATICS
  switch (coulomb.method) {
  case COULOMB_NONE: p_tensor.n_coulomb = 0; break;
  case COULOMB_P3M: p_tensor.n_coulomb = 9; break;
  case COULOMB_DH:  p_tensor.n_coulomb = 9; break;
  default:
    fprintf(stderr, "%d: init_p_tensor: cannot calculate p_tensor for coulomb method %d\n",
	    this_node, coulomb.method);
    errexit();
  }
#endif
  p_tensor.n = p_tensor.n_pre+p_tensor.n_bonded+p_tensor.n_non_bonded+p_tensor.n_coulomb;
  realloc_doublelist(&(p_tensor.node),p_tensor.n);
  realloc_doublelist(&(p_tensor.sum),p_tensor.n);
  p_tensor.init_status = 0;
}


/* Derive the p_tensor */
/***********************/
void calc_p_tensor(double volume, IntList *p_list, int flag) {
  Particle *p1, *p2, *p3;
  int *p1_list, n_p1;
  int i,j,k,l, pp, indi,indj,startj,endj, type_num, type1,type2;
  double d[3],dist,dist2, f1[3],f2[3],f3[3];
  
  p1=malloc(1*sizeof(Particle)); p2=malloc(1*sizeof(Particle)); p3=malloc(1*sizeof(Particle)); 
  p1_list = p_list->e; n_p1 = p_list->n;
  for(i=0;i<p_tensor.n;i++) {
    p_tensor.node.e[i] = 0.0;
    p_tensor.sum.e[i]  = 0.0;
  }

  if (parameter_changed || interactions_changed || topology_changed || particle_changed) {
    mpi_integrate(0);
  }

  for(indi=0; indi<n_p1; indi++) {
    if (get_particle_data(p1_list[indi], p1) != TCL_OK) { fprintf(stderr,"The particle %d you requested does not exist! ",p1_list[indi]); errexit(); }

    /* ideal gas contribution (the rescaling of the velocities by '/=time_step' each will be done later) */
    for(k=0;k<3;k++) { for(l=0;l<3;l++) { 
      p_tensor.sum.e[9+ k*3 + l] += (p1->v[k])*(p1->v[l]);
    } }

    /* bonded interactions */
    i=0;
    while(i < p1->bl.n) {
      if((flag==1) || (intlist_contains(p_list,p1->bl.e[i+1])==1)) {
	get_particle_data(p1->bl.e[i+1], p2);
	f1[0] = p1->f[0]; f1[1] = p1->f[1]; f1[2] = p1->f[2];
	f2[0] = p2->f[0]; f2[1] = p2->f[1]; f2[2] = p2->f[2];
	get_mi_vector(d, p1->r.p, p2->r.p);
	type_num = p1->bl.e[i];
	switch(bonded_ia_params[type_num].type) {
	case BONDED_IA_FENE:
	  add_fene_pair_force(p1,p2,type_num);
	  for(k=0;k<3;k++) { p1->f[k] -= (f1[k] = p1->f[k] - f1[k]); p2->f[k] -= (f2[k] = p2->f[k] - f2[k]); }
	  for(k=0;k<3;k++) { for(l=0;l<3;l++) { 
	    p_tensor.node.e[p_tensor.n_pre+9*type_num + k*3 + l] += f1[k]*d[l];
	    p_tensor.sum.e[p_tensor.n_pre+ k*3 + l] += f1[k]*d[l];
	  } }
	  i+=2; break;
	case BONDED_IA_HARMONIC:
	  add_harmonic_pair_force(p1,p2,type_num);
	  for(k=0;k<3;k++) { p1->f[k] -= (f1[k] = p1->f[k] - f1[k]); p2->f[k] -= (f2[k] = p2->f[k] - f2[k]); }
	  for(k=0;k<3;k++) { for(l=0;l<3;l++) { 
	    p_tensor.node.e[p_tensor.n_pre+9*type_num + k*3 + l] += f1[k]*d[l];
	    p_tensor.sum.e[p_tensor.n_pre+ k*3 + l] += f1[k]*d[l];
	  } }
	  i+=2; break;
	case BONDED_IA_ANGLE:
	  get_particle_data(p1->bl.e[i+2], p3);
	  f3[0] = p3->f[0]; f3[1] = p3->f[1]; f3[2] = p3->f[2];
	  add_angle_force(p1,p2,p3,type_num);
	  for(k=0;k<3;k++) { p1->f[k] -= (f1[k] = p1->f[k] - f1[k]); p2->f[k] -= (f2[k] = p2->f[k] - f2[k]); }
	  for(k=0;k<3;k++) { for(l=0;l<3;l++) { 
	    p_tensor.node.e[p_tensor.n_pre+9*type_num + k*3 + l] += -f2[k]*d[l];
	    p_tensor.sum.e[p_tensor.n_pre+ k*3 + l] += -f2[k]*d[l];
	  } }
	  for(k=0;k<3;k++) { p3->f[k] -= (f3[k] = p3->f[k] - f3[k]); }
	  get_mi_vector(d, p1->r.p, p3->r.p);
	  for(k=0;k<3;k++) { for(l=0;l<3;l++) { 
	    p_tensor.node.e[p_tensor.n_pre+9*type_num + k*3 + l] += -f3[k]*d[l];
	    p_tensor.sum.e[p_tensor.n_pre+ k*3 + l] += -f3[k]*d[l];
	  } }
	  if (p3->bl.n > 0) { realloc_intlist(&(p3->bl),0); }
	  i+=3; break;
	default :
	  fprintf(stderr,"WARNING: Bonds of atom %d unknown\n",p1->r.identity);
	  i = p1->bl.n; break;
	}
      }
    }

    /* non-bonded interactions and electrostatics */
    if(flag==1) { startj=0; endj=n_total_particles; } else { startj=indi+1; endj=n_p1; }
    for(indj=startj; indj<endj; indj++) {
      if(flag==1) {
	if((indj == p1->r.identity) || (intlist_contains(p_list,indj)==1)) continue;
	get_particle_data(indj, p2); }
      else get_particle_data(p1_list[indj], p2);

      /* save current force information */
      for(j=0;j<3;j++) { f1[j] = p1->f[j]; f2[j] = p2->f[j]; }

      /* distance calculation */
      get_mi_vector(d, p1->r.p, p2->r.p);                 // for(j=0; j<3; j++) d[j] = p1->r.p[j] - p2->r.p[j];
      dist2 = SQR(d[0]) + SQR(d[1]) + SQR(d[2]);
      dist  = sqrt(dist2);

      /* non-bonded interactions */
      pp = p_tensor.n_pre+p_tensor.n_bonded;
      if (checkIfParticlesInteract(p1->r.type, p2->r.type)) { 
	ia_params = get_ia_param(p1->r.type,p2->r.type);

	/* derive index 'type_num' */
	if(p1->r.type > p2->r.type) { type1 = p2->r.type; type2 = p1->r.type; } else { type2 = p2->r.type; type1 = p1->r.type; }
	type_num = pp + 9*( ((2 * n_particle_types - 1 - type1) * type1) / 2  +  type2);
	
	/* lennnard jones */
	add_lj_pair_force(p1,p2,ia_params,d,dist);
	add_ljcos_pair_force(p1,p2,ia_params,d,dist);
#ifdef ROTATION  
	add_gb_pair_force(p1,p2,ia_params,d,dist);
#endif
	for(j=0;j<3;j++) { p1->f[j] -= (f1[j] = p1->f[j] - f1[j]); p2->f[j] -= (f2[j] = p2->f[j] - f2[j]); }
	for(k=0;k<3;k++) { for(l=0;l<3;l++) { 
	  p_tensor.node.e[type_num + k*3 + l] += f1[k]*d[l];
	  p_tensor.sum.e[pp+ k*3 + l] += f1[k]*d[l];
	} }
      }
	
#ifdef ELECTROSTATICS
      /* real space coulomb */
      pp = p_tensor.n_pre+p_tensor.n_bonded+p_tensor.n_non_bonded;
      if(coulomb.method==COULOMB_P3M) {
	add_p3m_coulomb_pair_force(p1,p2,d,dist2,dist);
	for(j=0;j<3;j++) { p1->f[j] -= (f1[j] = p1->f[j] - f1[j]); p2->f[j] -= (f2[j] = p2->f[j] - f2[j]); }
	for(k=0;k<3;k++) { for(l=0;l<3;l++) { 
	  p_tensor.sum.e[pp+ k*3 + l] += f1[k]*d[l];
	} }
      }
      else if(coulomb.method==COULOMB_DH) {
	add_dh_coulomb_pair_force(p1,p2,d,dist);
	for(j=0;j<3;j++) { p1->f[j] -= (f1[j] = p1->f[j] - f1[j]); p2->f[j] -= (f2[j] = p2->f[j] - f2[j]); }
	for(k=0;k<3;k++) { for(l=0;l<3;l++) { 
	  p_tensor.sum.e[pp+ k*3 + l] += f1[k]*d[l];
	} }
      }
#endif
    } 
  }

  /* Rescale entries and sum all contributions */
  for(i=9; i<2*9; i++) { 
    p_tensor.sum.e[i] /= 2.0*volume*SQR(time_step); 
    p_tensor.sum.e[i%9] += p_tensor.sum.e[i];
  }
  for(i=p_tensor.n_pre;i<p_tensor.n;i++) {
    p_tensor.node.e[i]  /= 3.0*volume;
  }
  pp=p_tensor.n_pre+p_tensor.n_bonded;
  for(i=p_tensor.n_pre; i<pp; i++) {
    p_tensor.sum.e[i]   /= 3.0*volume;
    p_tensor.sum.e[i%9] += p_tensor.sum.e[i];
  }
  pp+=p_tensor.n_non_bonded;
  for(i=pp-p_tensor.n_non_bonded; i<pp; i++) {
    p_tensor.sum.e[i]   /= 3.0*volume;
    p_tensor.sum.e[i%9] += p_tensor.sum.e[i];
  }
#ifdef ELECTROSTATICS
  pp+=p_tensor.n_coulomb;
  for(i=pp-p_tensor.n_coulomb; i<pp; i++) {
    p_tensor.sum.e[i]   /= 3.0*volume;
    p_tensor.sum.e[i%9] += p_tensor.sum.e[i];
  }
#endif

  /* Total Sum = Trace of 1st tensor */
  p_tensor.node.e[0] = p_tensor.sum.e[0] + p_tensor.sum.e[4] + p_tensor.sum.e[8];
  p_tensor.node.e[1] = SQR(p_tensor.node.e[0]);

  /* Clean up particles */
  if (p1->bl.n > 0) { realloc_intlist(&(p1->bl),0); }
  //  if (p2->bl.n > 0) { realloc_intlist(&(p2->bl),0); }
  //  free_particle(p1); free_particle(p2); free_particle(p3); 
  free(p1); free(p2); free(p3); 
}
