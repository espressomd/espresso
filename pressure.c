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
#include "cells.h"
#include "integrate.h"
#include "domain_decomposition.h"
#include "nsquare.h"

Observable_stat virials = {0, {NULL,0,0}, 0,0,0,0};
Observable_stat total_pressure = {0, {NULL,0,0}, 0,0,0,0};
Observable_stat p_tensor= {0, {NULL,0,0},0,0,0,0};

/************************************************************/
/* local prototypes                                         */
/************************************************************/

/** Calculate long range virials (P3M, MMM2d...). */
void calc_long_range_virials();

/** Initializes a virials Observable stat. */
void init_virials(Observable_stat *stat);

/** on the master node: calc energies only if necessary */
void master_pressure_calc();

/*******************/
/* Scalar Pressure */
/*******************/

void pressure_calc(double *result)
{
  int n;

  double volume = box_l[0]*box_l[1]*box_l[2];

  init_virials(&virials);

  if(resort_particles) {
    initialize_ghosts(DD_GLOBAL_EXCHANGE);
    resort_particles = 0;
  }

  switch (cell_structure.type) {
  case CELL_STRUCTURE_DOMDEC:
    if (rebuild_verletlist) build_verlet_lists();
    calculate_verlet_virials();
    break;
  case CELL_STRUCTURE_NSQUARE:
    nsq_calculate_virials();
  }
  /* rescale kinetic energy */
  virials.data.e[0] /= (2.0*time_step*time_step*volume);

  calc_long_range_virials();

  for (n = 1; n < virials.data.n; n++)
    virials.data.e[n] /= 3*volume;

  /* gather data */
  MPI_Reduce(virials.data.e, result, virials.data.n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
}

/************************************************************/

void calc_long_range_virials()
{
#ifdef ELECTROSTATICS  
  /* calculate k-space part of electrostatic interaction. */
  switch (coulomb.method) {
  case COULOMB_P3M:
    virials.coulomb[1] = P3M_calc_kspace_forces(0,1);
    break;
  }
#endif
}

/************************************************************/

void init_virials(Observable_stat *stat)
{
  int n_pre, n_non_bonded, n_coulomb;

  n_pre        = 1;
  n_non_bonded = (n_particle_types*(n_particle_types+1))/2;

  n_coulomb    = 0;
#ifdef ELECTROSTATICS
  if(coulomb.bjerrum != 0.0) {
    n_coulomb  = 1;
    if(coulomb.method==COULOMB_P3M)
      n_coulomb += 1;
  }
#endif

  obsstat_realloc_and_clear(stat, n_pre, n_bonded_ia, n_non_bonded, n_coulomb, 1);
  stat->init_status = 0;
}

/************************************************************/

void master_pressure_calc() {
  mpi_gather_stats(2, total_pressure.data.e);

  total_pressure.init_status=1;
}


/****************************************************************************************
 *                                 parser
 ****************************************************************************************/

static void print_detailed_pressure(Tcl_Interp *interp)
{
  char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE + 2];
  double value;
  int i, j;

  value = total_pressure.data.e[0];
  for (i = 1; i < total_pressure.data.n; i++)
    value += total_pressure.data.e[i];

  Tcl_PrintDouble(interp, value, buffer);
  Tcl_AppendResult(interp, "{ pressure ", buffer, " } ", (char *)NULL);

  Tcl_PrintDouble(interp, total_pressure.data.e[0], buffer);
  Tcl_AppendResult(interp, "{ ideal ", buffer, " } ", (char *)NULL);

  for(i=0;i<n_bonded_ia;i++) {
    sprintf(buffer, "%d ", i);
    Tcl_AppendResult(interp, "{ ", buffer, (char *)NULL);
    Tcl_PrintDouble(interp, *obsstat_bonded(&total_pressure, i), buffer);
    Tcl_AppendResult(interp,
		     get_name_of_bonded_ia(bonded_ia_params[i].type),
		     " ", buffer, " } ", (char *) NULL);
  }

  for (i = 0; i < n_particle_types; i++)
    for (j = i; j < n_particle_types; j++) {
      if (checkIfParticlesInteract(i, j)) {
	sprintf(buffer, "%d ", i);
	Tcl_AppendResult(interp, "{ ", buffer, (char *)NULL);
	sprintf(buffer, "%d ", j);
	Tcl_AppendResult(interp, " ", buffer, (char *)NULL);
	Tcl_PrintDouble(interp, *obsstat_nonbonded(&total_pressure, i, j), buffer);
	Tcl_AppendResult(interp, "nonbonded ", buffer, " } ", (char *)NULL);	    
      }
    }
  
#ifdef ELECTROSTATICS
  if(coulomb.bjerrum != 0.0) {
    /* total Coulomb pressure */
    value = total_pressure.coulomb[0];
    for (i = 1; i < total_pressure.n_coulomb; i++)
      value += total_pressure.coulomb[i];
    Tcl_PrintDouble(interp, value, buffer);
    Tcl_AppendResult(interp, "{ coulomb ", buffer, (char *)NULL);

    /* if it is split up, then print the split up parts */
    if (total_pressure.n_coulomb > 1) {
      for (i = 0; i < total_pressure.n_coulomb; i++) {
	Tcl_PrintDouble(interp, total_pressure.coulomb[i], buffer);
	Tcl_AppendResult(interp, " ", buffer, (char *)NULL);
      }
    }
    Tcl_AppendResult(interp, " }", (char *)NULL);
  }
#endif
}

/************************************************************/

int parse_and_print_pressure(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze pressure [{ fene <type_num> | harmonic <type_num> | lj <type1> <type2> | ljcos <type1> <type2> | gb <type1> <type2> | coulomb | ideal | total }]' */
  char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE + 2];
  int i, j;

  if (n_total_particles == 0) {
    Tcl_AppendResult(interp, "(no particles)",
		     (char *)NULL);
    return (TCL_OK);
  }

  if (total_pressure.init_status == 0) {
    init_virials(&total_pressure);
    master_pressure_calc();
  }

  if (argc == 0)
    print_detailed_pressure(interp);
  else {
    double value;
    if      (ARG0_IS_S("ideal"))
      value = total_pressure.data.e[0];
    else if (ARG0_IS_S("bonded") ||
	     ARG0_IS_S("fene") ||
	     ARG0_IS_S("subt_lj_harm") ||
	     ARG0_IS_S("subt_lj_fene") ||
	     ARG0_IS_S("harmonic")) {
      if(argc<2 || ! ARG1_IS_I(i)) {
	Tcl_ResetResult(interp);
	Tcl_AppendResult(interp, "wrong # or type of arguments for: analyze pressure bonded <type_num>",
			 (char *)NULL);
	return (TCL_ERROR);
      }
      if(i < 0 || i >= n_bonded_ia) {
	Tcl_AppendResult(interp, "bond type does not exist", (char *)NULL);
	return (TCL_ERROR);
      }
      value = *obsstat_bonded(&total_pressure, i);
    }
    else if (ARG0_IS_S("nonbonded") ||
	     ARG0_IS_S("lj") ||
	     ARG0_IS_S("lj-cos") ||
	     ARG0_IS_S("tabulated") ||
	     ARG0_IS_S("gb")) {
      if(argc<3 || ! ARG_IS_I(1, i) || ! ARG_IS_I(2, j)) {
	Tcl_ResetResult(interp);
	Tcl_AppendResult(interp, "wrong # or type of arguments for: analyze pressure nonbonded <type1> <type2>",
			 (char *)NULL);
	return (TCL_ERROR);
      }
      if(i < 0 || i >= n_particle_types || j < 0 || j >= n_particle_types) {
	Tcl_AppendResult(interp, "particle type does not exist", (char *)NULL);
	return (TCL_ERROR);
      }
      value = *obsstat_nonbonded(&total_pressure, i, j);
    }
    else if( ARG0_IS_S("coulomb")) {
#ifdef ELECTROSTATICS
      value = total_pressure.coulomb[0];
      for (i = 1; i < total_pressure.n_coulomb; i++)
	value += total_pressure.coulomb[i];
#else
      Tcl_AppendResult(interp, "ELECTROSTATICS not compiled (see config.h)\n", (char *)NULL);
#endif
    }
    else if (ARG0_IS_S("total")) {
      value = total_pressure.data.e[0];
      for (i = 1; i < total_pressure.data.n; i++)
	value += total_pressure.data.e[i];
    }
    else {
      Tcl_AppendResult(interp, "unknown feature of: analyze pressure",
		       (char *)NULL);
      return (TCL_ERROR);
    }
    Tcl_PrintDouble(interp, value, buffer);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
  }

  return (TCL_OK);
}

/**********************/
/* Tensorial Pressure */
/**********************/


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
	elements[counter++] = partCfg[j].p.identity;
	new_bin[i] += 1;
      }
    }
    printf("vol %d %le\n",i,volumes[i]);
  }
}

/************************************************************/

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
    updatePartCfg(WITHOUT_BONDS);
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

int parse_and_print_p_IK1(Tcl_Interp *interp, int argc, char **argv)
{

#if 0

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
    case BONDED_IA_SUBT_LJ_HARM:
      Tcl_AppendResult(interp, "{ SUBT_LJ_HARM ", (char *)NULL);
      for(j=0; j<9; j++) { sprintf(buffer,"%f ",p_tensor.node.e[p+i*9+j]); Tcl_AppendResult(interp, buffer, (char *)NULL); }
      Tcl_AppendResult(interp, "} ", (char *)NULL);
      break;
    case BONDED_IA_SUBT_LJ_FENE:
      Tcl_AppendResult(interp, "{ SUBT_LJ_FENE ", (char *)NULL);
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

#endif

  return (TCL_OK);

}


/* Initialize the p_tensor */
/***************************/
void init_p_tensor() {

#if 0

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
      p_tensor.sum.e[9+ k*3 + l] += (p1->m.v[k])*(p1->m.v[l]);
    } }

    /* bonded interactions */
    i=0;
    while(i < p1->bl.n) {
      if((flag==1) || (intlist_contains(p_list,p1->bl.e[i+1])==1)) {
	get_particle_data(p1->bl.e[i+1], p2);
	f1[0] = p1->f.f[0]; f1[1] = p1->f.f[1]; f1[2] = p1->f.f[2];
	f2[0] = p2->f.f[0]; f2[1] = p2->f.f[1]; f2[2] = p2->f.f[2];
	get_mi_vector(d, p1->r.p, p2->r.p);
	type_num = p1->bl.e[i];
	switch(bonded_ia_params[type_num].type) {
	case BONDED_IA_FENE:
	  add_fene_pair_force(p1,p2,type_num);
	  for(k=0;k<3;k++) { p1->f.f[k] -= (f1[k] = p1->f.f[k] - f1[k]); p2->f.f[k] -= (f2[k] = p2->f.f[k] - f2[k]); }
	  for(k=0;k<3;k++) { for(l=0;l<3;l++) { 
	    p_tensor.node.e[p_tensor.n_pre+9*type_num + k*3 + l] += f1[k]*d[l];
	    p_tensor.sum.e[p_tensor.n_pre+ k*3 + l] += f1[k]*d[l];
	  } }
	  i+=2; break;
	case BONDED_IA_HARMONIC:
	  add_harmonic_pair_force(p1,p2,type_num);
	  for(k=0;k<3;k++) { p1->f.f[k] -= (f1[k] = p1->f.f[k] - f1[k]); p2->f.f[k] -= (f2[k] = p2->f.f[k] - f2[k]); }
	  for(k=0;k<3;k++) { for(l=0;l<3;l++) { 
	    p_tensor.node.e[p_tensor.n_pre+9*type_num + k*3 + l] += f1[k]*d[l];
	    p_tensor.sum.e[p_tensor.n_pre+ k*3 + l] += f1[k]*d[l];
	  } }
	  i+=2; break;
	case BONDED_IA_SUBT_LJ_HARM:
	  add_subt_lj_harm_pair_force(p1,p2,type_num);
	  for(k=0;k<3;k++) { p1->f.f[k] -= (f1[k] = p1->f.f[k] - f1[k]); p2->f.f[k] -= (f2[k] = p2->f.f[k] - f2[k]); }
	  for(k=0;k<3;k++) { for(l=0;l<3;l++) { 
	    p_tensor.node.e[p_tensor.n_pre+9*type_num + k*3 + l] += f1[k]*d[l];
	    p_tensor.sum.e[p_tensor.n_pre+ k*3 + l] += f1[k]*d[l];
	  } }
	  i+=2; break;
	case BONDED_IA_SUBT_LJ_FENE:
	  add_subt_lj_fene_pair_force(p1,p2,type_num);
	  for(k=0;k<3;k++) { p1->f.f[k] -= (f1[k] = p1->f.f[k] - f1[k]); p2->f.f[k] -= (f2[k] = p2->f.f[k] - f2[k]); }
	  for(k=0;k<3;k++) { for(l=0;l<3;l++) { 
	    p_tensor.node.e[p_tensor.n_pre+9*type_num + k*3 + l] += f1[k]*d[l];
	    p_tensor.sum.e[p_tensor.n_pre+ k*3 + l] += f1[k]*d[l];
	  } }
	  i+=2; break;
	case BONDED_IA_ANGLE:
	  get_particle_data(p1->bl.e[i+2], p3);
	  f3[0] = p3->f.f[0]; f3[1] = p3->f.f[1]; f3[2] = p3->f.f[2];
	  add_angle_force(p1,p2,p3,type_num);
	  for(k=0;k<3;k++) { p1->f.f[k] -= (f1[k] = p1->f.f[k] - f1[k]); p2->f.f[k] -= (f2[k] = p2->f.f[k] - f2[k]); }
	  for(k=0;k<3;k++) { for(l=0;l<3;l++) { 
	    p_tensor.node.e[p_tensor.n_pre+9*type_num + k*3 + l] += -f2[k]*d[l];
	    p_tensor.sum.e[p_tensor.n_pre+ k*3 + l] += -f2[k]*d[l];
	  } }
	  for(k=0;k<3;k++) { p3->f.f[k] -= (f3[k] = p3->f.f[k] - f3[k]); }
	  get_mi_vector(d, p1->r.p, p3->r.p);
	  for(k=0;k<3;k++) { for(l=0;l<3;l++) { 
	    p_tensor.node.e[p_tensor.n_pre+9*type_num + k*3 + l] += -f3[k]*d[l];
	    p_tensor.sum.e[p_tensor.n_pre+ k*3 + l] += -f3[k]*d[l];
	  } }
	  if (p3->bl.n > 0) { realloc_intlist(&(p3->bl),0); }
	  i+=3; break;
	default :
	  fprintf(stderr,"WARNING: Bonds of atom %d unknown\n",p1->p.identity);
	  i = p1->bl.n; break;
	}
      }
    }

    /* non-bonded interactions and electrostatics */
    if(flag==1) { startj=0; endj=n_total_particles; } else { startj=indi+1; endj=n_p1; }
    for(indj=startj; indj<endj; indj++) {
      if(flag==1) {
	if((indj == p1->p.identity) || (intlist_contains(p_list,indj)==1)) continue;
	get_particle_data(indj, p2); }
      else get_particle_data(p1_list[indj], p2);

      /* save current force information */
      for(j=0;j<3;j++) { f1[j] = p1->f.f[j]; f2[j] = p2->f.f[j]; }

      /* distance calculation */
      get_mi_vector(d, p1->r.p, p2->r.p);                 // for(j=0; j<3; j++) d[j] = p1->r.p[j] - p2->r.p[j];
      dist2 = SQR(d[0]) + SQR(d[1]) + SQR(d[2]);
      dist  = sqrt(dist2);

      /* non-bonded interactions */
      pp = p_tensor.n_pre+p_tensor.n_bonded;
      if (checkIfParticlesInteract(p1->p.type, p2->p.type)) { 
	ia_params = get_ia_param(p1->p.type,p2->p.type);

	/* derive index 'type_num' */
	if(p1->p.type > p2->p.type) { type1 = p2->p.type; type2 = p1->p.type; } else { type2 = p2->p.type; type1 = p1->p.type; }
	type_num = pp + 9*( ((2 * n_particle_types - 1 - type1) * type1) / 2  +  type2);
	
	/* lennnard jones */
	add_lj_pair_force(p1,p2,ia_params,d,dist);
#ifdef LJCOS
	add_ljcos_pair_force(p1,p2,ia_params,d,dist);
#endif
#ifdef ROTATION  
	add_gb_pair_force(p1,p2,ia_params,d,dist);
#endif
	for(j=0;j<3;j++) { p1->f.f[j] -= (f1[j] = p1->f.f[j] - f1[j]); p2->f.f[j] -= (f2[j] = p2->f.f[j] - f2[j]); }
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
	for(j=0;j<3;j++) { p1->f.f[j] -= (f1[j] = p1->f.f[j] - f1[j]); p2->f.f[j] -= (f2[j] = p2->f.f[j] - f2[j]); }
	for(k=0;k<3;k++) { for(l=0;l<3;l++) { 
	  p_tensor.sum.e[pp+ k*3 + l] += f1[k]*d[l];
	} }
      }
      else if(coulomb.method==COULOMB_DH) {
	add_dh_coulomb_pair_force(p1,p2,d,dist);
	for(j=0;j<3;j++) { p1->f.f[j] -= (f1[j] = p1->f.f[j] - f1[j]); p2->f.f[j] -= (f2[j] = p2->f.f[j] - f2[j]); }
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

#endif

}

