/** \file statistics.c
    This is the place for analysation (so far...).
    Implementation of \ref statistics.h "statistics.h"
 */
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include "statistics.h"
#include "forces.h"
#include "communication.h"
#include "grid.h"
#include "integrate.h"
#include "particle_data.h"
#include "interaction_data.h"
#include "debug.h"
#include "cells.h"
#include "verlet.h"
#include "p3m.h"
/* include the force files */
#include "lj.h"
#include "fene.h"
#include "angle.h"
#include "debye_hueckel.h"

/** Particles' initial positions (needed for g1(t), g2(t), g3(t) in \ref #analyze) */
float *partCoord_g=NULL, *partCM_g=NULL;
int n_part_g = 0, n_chains_g = 0;

/** Previous particle configurations (needed for offline analysis and correlation analysis in \ref #analyze) */
double **configs = NULL; int n_configs = 0;

/** data for a system consisting of chains */
int chain_start = 0, chain_n_chains = 0, chain_length = 0;

Energy_stat energy = {0, {NULL,0,0}, {NULL,0,0}, 0,0,0,0,0,0};

static int get_reference_point(Tcl_Interp *interp, int *argc, char ***argv,
			       double *posx, double *posy, double *posz, int *pid)
{
  *pid = -1;

  if (*argc < 3) {
    Particle ref;
    if (Tcl_GetInt(interp, (*argv)[0], pid) == TCL_ERROR)
      return TCL_ERROR;

    if (get_particle_data(*pid, &ref) != TCL_OK) {
      Tcl_AppendResult(interp, "reference particle does not exist", (char *)NULL);
      return TCL_ERROR;
    }
    *posx = ref.r.p[0];
    *posy = ref.r.p[1];
    *posz = ref.r.p[2];
    realloc_intlist(&ref.bl, 0);
    
    (*argc)--;
    (*argv)++;

    return TCL_OK;
  }
  else {
    if (Tcl_GetDouble(interp, (*argv)[0], posx) == TCL_ERROR ||
	Tcl_GetDouble(interp, (*argv)[1], posy) == TCL_ERROR ||
	Tcl_GetDouble(interp, (*argv)[2], posz) == TCL_ERROR)
      return TCL_ERROR;

    (*argc) -= 3;
    (*argv) += 3;

    return TCL_OK;
  }
  return TCL_ERROR;
}

/** this function scans the arguments for a description of the chain structure,
    i.e. start of chains, number of chains and chain length. Since this structure
    requires the particles to be sorted, this is performed, too. */
static int prepare_chain_structure_info(Tcl_Interp *interp, int *argc, char ***argv)
{
  /* 'analyze set chains [<chain_start> <n_chains> <chain_length>]' */
  if (*argc > 0) {
    if (*argc < 3) {
      Tcl_AppendResult(interp, "chain structure info consists of <start> <n> <length>", (char *)NULL);    
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, (*argv)[0], &chain_start) == TCL_ERROR ||
	Tcl_GetInt(interp, (*argv)[1], &chain_n_chains) == TCL_ERROR ||
	Tcl_GetInt(interp, (*argv)[2], &chain_length) == TCL_ERROR)
      return TCL_ERROR;

    (*argc) -= 3;
    (*argv) += 3;
  }

  if (!sortPartCfg()) {
    Tcl_AppendResult(interp, "for analyze, store particles consecutively starting with 0.",
		     (char *) NULL);
    return (TCL_ERROR);      
  }

  return TCL_OK;
}
      
static int check_chain_structure_info(Tcl_Interp *interp)
{
  if (max_seen_particle <= chain_start + chain_n_chains*chain_length) {
    Tcl_AppendResult(interp, "not enough particles for chain structure", (char *)NULL);
    return TCL_ERROR;
  }
  return TCL_OK;
}

int analyze(ClientData data, Tcl_Interp *interp, int argc, char **argv)
{
  char buffer[50 + TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
  char *mode;
  int  i, j, p;
  double *buf = NULL;
  double result, posx, posy, posz, r_catch;

  if (argc < 2) {
    Tcl_AppendResult(interp, "Wrong # of args! Usage: analyze <what> ...", (char *)NULL);
    return (TCL_ERROR);
  }

  mode = argv[1];
  argc -= 2;
  argv += 2;

  if (!strncmp(mode, "set", strlen(mode))) {
    /* 'analyze set <structure info>' */
    /**********************************/
    if (argc == 0) {
      Tcl_AppendResult(interp, "which topology are you interested in?", (char *)NULL);
      return TCL_ERROR;
    }
    if (!strncmp(argv[0], "chains", strlen(argv[0]))) {
      /* 'analyze set chains [<chain_start> <n_chains> <chain_length>]' */
      /******************************************************************/
      argc--;
      argv++;
      if (argc == 0) {
 	sprintf(buffer, "%d ", chain_start);
	Tcl_AppendResult(interp, buffer, (char *)NULL);
 	sprintf(buffer, "%d ", chain_n_chains);
	Tcl_AppendResult(interp, buffer, (char *)NULL);
 	sprintf(buffer, "%d ", chain_length);
	Tcl_AppendResult(interp, buffer, (char *)NULL);
	return TCL_OK;
      }
      if (prepare_chain_structure_info(interp, &argc, &argv) == TCL_ERROR)
	return TCL_ERROR;
      return TCL_OK;
    }
    else {
      /* default */
      /***********/
      Tcl_AppendResult(interp, "The topology \"", argv[0],
		       "\" you requested is not implemented.", (char *)NULL);
      return (TCL_ERROR);
    }
  }
  else if (!strncmp(mode, "mindist", strlen(mode))) {
    /* 'analyze mindist' */
    /*********************/
    if (n_total_particles <= 1) {
      Tcl_AppendResult(interp, "(not enough particles)",
		       (char *)NULL);
      return (TCL_OK);
    }
    result = mindist();
    Tcl_PrintDouble(interp, result, buffer);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
    return TCL_OK;
  }
  else if (!strncmp(mode, "nbhood", strlen(mode))) {
    /* 'analyze nbhood { <partid> | <posx> <posy> <posz> } <r_catch>' */
    /******************************************************************/
    IntList il;

    if (n_total_particles == 0) {
      Tcl_AppendResult(interp, "(no particles)",
		       (char *)NULL);
      return (TCL_OK);
    }

    get_reference_point(interp, &argc, &argv, &posx, &posy, &posz, &p);

    if (argc != 1) {
      Tcl_AppendResult(interp, "usage: nbhood { <partid> | <posx> <posy> <posz> } <r_catch>",
		       (char *)NULL);
      return TCL_ERROR;
    }

    r_catch = atof(argv[0]);

    updatePartCfg();

    nbhood(posx, posy, posz, r_catch, &il);
    
    for (i = 0; i < il.n; i++) {
      sprintf(buffer, "%d ", il.e[i]);
      Tcl_AppendResult(interp, buffer, (char *)NULL);
    }
    realloc_intlist(&il, 0);
    return (TCL_OK);
  }
  else if (!strncmp(mode, "distto", strlen(mode))) {
    /* 'analyze distto { <part_id> | <posx> <posy> <posz> }' */
    /*********************************************************/

    if (n_total_particles == 0) {
      Tcl_AppendResult(interp, "(no particles)",
		       (char *)NULL);
      return (TCL_OK);
    }

    get_reference_point(interp, &argc, &argv, &posx, &posy, &posz, &p);
    if (argc != 0) {
      Tcl_AppendResult(interp, "usage: distto { <partid> | <posx> <posy> <posz> }",
		       (char *)NULL);
      return TCL_ERROR;
    }

    updatePartCfg();

    result = distto(posx, posy, posz, p);

    Tcl_PrintDouble(interp, result, buffer);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
    return (TCL_OK);
  }
  else if (!strncmp(mode, "energy", strlen(mode))) {
    /* 'analyze energy [{ fene <type_num> | lj <type1> <type2> | coulomb | kinetic }]' */
    /***********************************************************************************/
    /* checks */
    if (n_total_particles == 0) {
      Tcl_AppendResult(interp, "(no particles)",
		       (char *)NULL);
      return (TCL_OK);
    }

    /* check energy initialization */
    if(energy.init_status == 0 || interactions_changed) {
      init_energies();
      energy.init_status=0;
    }
 
    /* check integrator status */
    if (parameter_changed || interactions_changed || topology_changed || particle_changed) {
      mpi_integrate(0);
    }

    /* specify analysation */
    if(argc>0) {
      if(!strncmp(argv[0], "kinetic", strlen(argv[0]))) energy.ana_num=1;
      else if(!strncmp(argv[0], "fene", strlen(argv[0]))) {
	if(argc<2) {
	  Tcl_AppendResult(interp, "wrong # arguments for: analyze energy fene <type_num>",
			   (char *)NULL);
	  return (TCL_ERROR);
	}
	if(Tcl_GetInt(interp, argv[1], &i) == TCL_ERROR) return (TCL_ERROR);
	if(i >= energy.n_bonded) return (TCL_ERROR);
	energy.ana_num = energy.n_pre+i;
      }
      else if(!strncmp(argv[0], "lj", strlen(argv[0]))) {
	if(argc<3) {
	  Tcl_AppendResult(interp, "wrong # arguments for: analyze energy fene <type_num>",
			   (char *)NULL);
	  return (TCL_ERROR);
	}
	if(Tcl_GetInt(interp, argv[1], &i) == TCL_ERROR) return (TCL_ERROR);
	if(Tcl_GetInt(interp, argv[2], &j) == TCL_ERROR) return (TCL_ERROR);
	if(i >= n_particle_types || j >= n_particle_types) return (TCL_ERROR);
	energy.ana_num = energy.n_pre+energy.n_bonded + j -i;
	while(i>0) {
	  energy.ana_num += n_particle_types - (i-1);
	  i--;
	}
      }
      else if(!strncmp(argv[0], "coulomb", strlen(argv[0]))) {
	energy.ana_num = energy.n_pre+energy.n_bonded+energy.n_non_bonded;
      }
      else {
	Tcl_AppendResult(interp, "unknown feature of: analyze energy",
			 (char *)NULL);
	return (TCL_ERROR);
      }
    }
    else energy.ana_num=0;

    mpi_gather_stats(1, buf);

    if(energy.ana_num > 0) {
      Tcl_PrintDouble(interp, energy.sum.e[0], buffer);
      Tcl_AppendResult(interp, buffer, (char *)NULL);
    }
    else {
      Tcl_PrintDouble(interp, energy.sum.e[0], buffer);
      Tcl_AppendResult(interp, "{ energy ", buffer, " } ", (char *)NULL);
      Tcl_PrintDouble(interp, energy.sum.e[1], buffer);
      Tcl_AppendResult(interp, "{ kinetic ", buffer, " } ", (char *)NULL);
      for(i=0;i<n_bonded_ia;i++) {
	switch (bonded_ia_params[i].type) {
	case BONDED_IA_FENE:
	  sprintf(buffer, "%d ", i);
	  Tcl_AppendResult(interp, "{ ", buffer, (char *)NULL);
	  Tcl_PrintDouble(interp, energy.sum.e[energy.n_pre+i], buffer);
	  Tcl_AppendResult(interp, "FENE ", buffer, " } ", (char *) NULL);
	  break;
	case BONDED_IA_ANGLE:
	  sprintf(buffer, "%d ", i);
	  Tcl_AppendResult(interp, "{ ", buffer, (char *)NULL);
	  Tcl_PrintDouble(interp, energy.sum.e[energy.n_pre+i], buffer);
	  Tcl_AppendResult(interp, "angle ", buffer, " } ", (char *) NULL);
	  break;
	case BONDED_IA_DIHEDRAL:
	  sprintf(buffer, "%d ", i);
	  Tcl_AppendResult(interp, "{ ", buffer, (char *)NULL);
	  Tcl_PrintDouble(interp, energy.sum.e[energy.n_pre+i], buffer);
	  Tcl_AppendResult(interp, "dihedral", buffer, " } ", (char *) NULL);
	  break;
	default:
	  ;
	}
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
	    Tcl_AppendResult(interp, "lj ", buffer, " } ", (char *)NULL);
	  }
	  p++;
	}
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
    }

    energy.init_status=1;
    return (TCL_OK);
  }
  else if ( (!strncmp(mode, "re", strlen(mode))) || (!strncmp(mode, "<re>", strlen(mode))) ) {
    /* 'analyze { re | <re> } [<chain_start> <n_chains> <chain_length>]' */
    /*********************************************************************/
    if (prepare_chain_structure_info(interp, &argc, &argv) == TCL_ERROR)
      return TCL_ERROR;
    if (argc != 0) {
      Tcl_AppendResult(interp, "only chain structure info required", (char *)NULL);
      return TCL_ERROR;
    }
    if (!strncmp(mode, "re", strlen(mode))) result = calc_re(); 
    else if (n_configs == 0) {
      Tcl_AppendResult(interp, "no configurations found! ", (char *)NULL);
      Tcl_AppendResult(interp, "Use 'analyze append' to save some, or 'analyze re' to only look at current state!", (char *)NULL);
      return TCL_ERROR; }
    else result = calc_re_av();
    Tcl_PrintDouble(interp, result, buffer);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
    return (TCL_OK);
  }
  else if ( (!strncmp(mode, "rg", strlen(mode))) || (!strncmp(mode, "<rg>", strlen(mode))) ) {
    /* 'analyze { rg | <rg> } [<chain_start> <n_chains> <chain_length>]' */
    /*********************************************************************/
    if (prepare_chain_structure_info(interp, &argc, &argv) == TCL_ERROR)
      return TCL_ERROR;
    if (argc != 0) {
      Tcl_AppendResult(interp, "only chain structure info required", (char *)NULL);
      return TCL_ERROR;
    }
    if (!strncmp(mode, "rg", strlen(mode))) result = calc_rg(); 
    else if (n_configs == 0) {
      Tcl_AppendResult(interp, "no configurations found! ", (char *)NULL);
      Tcl_AppendResult(interp, "Use 'analyze append' to save some, or 'analyze rg' to only look at current state!", (char *)NULL);
      return TCL_ERROR; }
    else result = calc_rg_av();
    Tcl_PrintDouble(interp, result, buffer);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
    return (TCL_OK);
  }
  else if ( (!strncmp(mode, "rh", strlen(mode))) || (!strncmp(mode, "<rh>", strlen(mode))) ) {
    /* 'analyze { rh | <rh> } [<chain_start> <n_chains> <chain_length>]' */
    /**********************************************************/
    if (prepare_chain_structure_info(interp, &argc, &argv) == TCL_ERROR)
      return TCL_ERROR;
    if (argc != 0) {
      Tcl_AppendResult(interp, "only chain structure info required", (char *)NULL);
      return TCL_ERROR;
    }
    result = calc_rh();
    if (!strncmp(mode, "rh", strlen(mode))) result = calc_rh(); 
    else if (n_configs == 0) {
      Tcl_AppendResult(interp, "no configurations found! ", (char *)NULL);
      Tcl_AppendResult(interp, "Use 'analyze append' to save some, or 'analyze rh' to only look at current state!", (char *)NULL);
      return TCL_ERROR; }
    else result = calc_rh_av();
    Tcl_PrintDouble(interp, result, buffer);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
    return (TCL_OK);
  }
  else if (!strncmp(mode, "g123", strlen(mode))) {
    /* 'analyze g123 [-init] [<chain_start> <n_chains> <chain_length>]' */
    /********************************************************************/
    int init = 0;
    double g1, g2, g3;

    if (argc >= 1 && !strncmp(argv[0], "-init", strlen(argv[0]))) {
      init = 1;
      argc--;
      argv++;
    }

    if (prepare_chain_structure_info(interp, &argc, &argv) == TCL_ERROR)
      return TCL_ERROR;
    if (argc != 0) {
      Tcl_AppendResult(interp, "only chain structure info required", (char *)NULL);
      return TCL_ERROR;
    }

    if (init) {
      init_g123();
      return TCL_OK;
    }

    if (partCoord_g == NULL || partCM_g == NULL) {
      Tcl_AppendResult(interp, "please call with -init first", (char *)NULL);
      return TCL_ERROR;
    }
    if (chain_n_chains != n_chains_g ||
	n_total_particles != n_part_g) {
      fprintf(stderr, "%d %d %d %d\n", chain_n_chains, n_chains_g, n_total_particles, n_part_g);
      Tcl_AppendResult(interp, "initial config has different topology",
		       (char *)NULL);
      return TCL_ERROR;      
    }

    calc_g123(&g1, &g2, &g3);
    sprintf(buffer,"{ %f %f %f }",g1, g2, g3);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
    return (TCL_OK);
  }
  else if ( (!strncmp(mode, "<g1>", strlen(mode))) || (!strncmp(mode, "<g2>", strlen(mode))) || (!strncmp(mode, "<g3>", strlen(mode))) ){
    /* 'analyze { <g1> | <g2> | <g3> } [<chain_start> <n_chains> <chain_length>]' */
    /************************************************************/
    int i; double *gx; 

    if (prepare_chain_structure_info(interp, &argc, &argv) == TCL_ERROR) return TCL_ERROR;
    if (argc != 0) { Tcl_AppendResult(interp, "only chain structure info required", (char *)NULL); return TCL_ERROR; }
    if (!strncmp(mode, "<g1>", strlen(mode))) 
      calc_g1_av(&gx);
    else if (!strncmp(mode, "<g2>", strlen(mode))) 
      calc_g2_av(&gx);
    else if (!strncmp(mode, "<g3>", strlen(mode))) 
      calc_g3_av(&gx);
    for (i=0; i<n_configs; i++) { 
      sprintf(buffer,"%f ",gx[i]); Tcl_AppendResult(interp, buffer, (char *)NULL); 
    }
    free(gx);
    return (TCL_OK);
  }
  else if (!strncmp(mode, "append", strlen(mode))) {
    /* 'analyze append' */
    /********************/
    if (argc != 0) { 
      Tcl_AppendResult(interp, "Wrong # of args! Usage: analyze append", (char *)NULL); return TCL_ERROR; 
    }
    if (!sortPartCfg()) { Tcl_AppendResult(interp, "for analyze, store particles consecutively starting with 0.",(char *) NULL); return (TCL_ERROR); }
    analyze_append();
    sprintf(buffer,"%d",n_configs); Tcl_AppendResult(interp, buffer, (char *)NULL); return TCL_OK;
  }
  else if (!strncmp(mode, "push", strlen(mode))) {
    /* 'analyze push [<size>]' */
    /*****************************/
    if (!sortPartCfg()) { Tcl_AppendResult(interp, "for analyze, store particles consecutively starting with 0.",(char *) NULL); return (TCL_ERROR); }
    if (argc == 1) { 
      Tcl_GetInt(interp, argv[0], &i); argc--; argv++;
      if (n_configs < i) analyze_append(); else analyze_push();
      if (n_configs > i) for(j=0; j < n_configs-i; j++) analyze_remove(0);
    }
    else if (argc != 0) { Tcl_AppendResult(interp, "Wrong # of args! Usage: analyze push [<size>]", (char *)NULL); return TCL_ERROR; }
    else if (n_configs > 0) analyze_push();
    else analyze_append();
    sprintf(buffer,"%d",n_configs); Tcl_AppendResult(interp, buffer, (char *)NULL); return TCL_OK;
  }
  else if (!strncmp(mode, "replace", strlen(mode))) {
    /* 'analyze replace <index>' */
    /*****************************/
    if (argc != 1) {
      Tcl_AppendResult(interp, "Wrong # of args! Usage: analyze replace <index>", (char *)NULL); return TCL_ERROR; 
    }
    if (!sortPartCfg()) { Tcl_AppendResult(interp, "for analyze, store particles consecutively starting with 0.",(char *) NULL); return (TCL_ERROR); }
    Tcl_GetInt(interp, argv[0], &i); argc--; argv++;
    if((n_configs == 0) && (i==0)) analyze_append();
    else if ((n_configs == 0) && (i!=0)) {
      Tcl_AppendResult(interp, "Nice try, but there are no stored configurations that could be replaced!", (char *)NULL); return TCL_ERROR; }
    else if((i < 0) || (i > n_configs-1)) {
      sprintf(buffer,"Index %d out of range (must be in [0,%d])!",i,n_configs-1);
      Tcl_AppendResult(interp, buffer, (char *)NULL); return TCL_ERROR; }
    else analyze_replace(i);
    sprintf(buffer,"%d",n_configs); Tcl_AppendResult(interp, buffer, (char *)NULL); return TCL_OK;
  }
  else if (!strncmp(mode, "remove", strlen(mode))) {
    /* 'analyze remove [<index>]' */
    /******************************/
    if (!sortPartCfg()) { Tcl_AppendResult(interp, "for analyze, store particles consecutively starting with 0.",(char *) NULL); return (TCL_ERROR); }
    if (argc == 0) { for (i = n_configs-1; i >= 0; i--) analyze_remove(i); }
    else if (argc == 1) {
      Tcl_GetInt(interp, argv[0], &i); argc--; argv++;
      if(n_configs == 0) {
	Tcl_AppendResult(interp, "Nice try, but there are no stored configurations that could be removed!", (char *)NULL); return TCL_ERROR; }
      else if((i < 0) || (i > n_configs-1)) {
	sprintf(buffer,"Index %d out of range (must be in [0,%d])!",i,n_configs-1);
	Tcl_AppendResult(interp, buffer, (char *)NULL); return TCL_ERROR; }
      analyze_remove(i);
    }
    else {
      Tcl_AppendResult(interp, "Wrong # of args! Usage: analyze remove [<index>]", (char *)NULL); return TCL_ERROR; 
    }
    sprintf(buffer,"%d",n_configs); Tcl_AppendResult(interp, buffer, (char *)NULL); return TCL_OK;
  }
  else if (!strncmp(mode, "stored", strlen(mode))) {
    /* 'analyze stored' */
    /********************/
    if (argc != 0) {
      Tcl_AppendResult(interp, "Wrong # of args! Usage: analyze stored", (char *)NULL); return TCL_ERROR; 
    }
    sprintf(buffer,"%d",n_configs); Tcl_AppendResult(interp, buffer, (char *)NULL); return TCL_OK;
  }
  else {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "The operation \"", mode,
		     "\" you requested is not implemented.", (char *)NULL);
    return (TCL_ERROR);
  }
  return (TCL_ERROR);
}

double mindist()
{
  double *buf = NULL;
  double mindist, xt, yt, zt, dx, dy, dz;
  int i, j;

  /* minimal pair distance */
  /* if the integrator already ran, check his result */
  if (minimum_part_dist == -1)
    mindist = box_l[0] + box_l[1] + box_l[2];
  else {
    /* get data from integrator */ 
    buf = malloc(n_nodes*sizeof(double));
    mpi_gather_stats(0, buf);
    mindist = buf[0];
    for (i = 1; i < n_nodes; i++) {
      if (buf[i] < mindist)
	mindist = buf[i];
    }
    free(buf);
  }

  if (mindist >= box_l[0] + box_l[1] + box_l[2]) {
    /* ok, the integrator has not been started, or the distance
       is larger than the real space cutoffs, so calculate directly */
    mindist = box_l[0] + box_l[1] + box_l[2];
    mindist *= mindist;
    
    updatePartCfg();

    for (j=0; j<n_total_particles-1; j++) {
      xt = partCfg[j].r.p[0]; yt = partCfg[j].r.p[1]; zt = partCfg[j].r.p[2];
      for (i=j+1; i<n_total_particles; i++) {
	dx = xt - partCfg[i].r.p[0];   dx -= dround(dx/box_l[0])*box_l[0];
	dy = yt - partCfg[i].r.p[1];   dy -= dround(dy/box_l[1])*box_l[1];
	dz = zt - partCfg[i].r.p[2];   dz -= dround(dz/box_l[2])*box_l[2];
	mindist = dmin(mindist, SQR(dx)+SQR(dy)+SQR(dz));
      }
    }
    mindist = sqrt(mindist);
  }
  return mindist;
}

void nbhood(double xt, double yt, double zt, double r, IntList *il)
{
  double dx, dy, dz;
  int i;

  init_intlist(il);

  updatePartCfg();

  for (i = 0; i<n_total_particles; i++) {
    dx = xt - partCfg[i].r.p[0];   dx -= dround(dx/box_l[0])*box_l[0];
    dy = yt - partCfg[i].r.p[1];   dy -= dround(dy/box_l[1])*box_l[1];
    dz = zt - partCfg[i].r.p[2];   dz -= dround(dz/box_l[2])*box_l[2];
    if (sqrt(SQR(dx)+SQR(dy)+SQR(dz)) < r) {
      realloc_intlist(il, il->n + 1);
      il->e[il->n] = partCfg[i].r.identity;
      il->n++;
    }
  }
}

double distto(double xt, double yt, double zt, int pid)
{
  int i;
  double dx, dy, dz;
  double mindist;

  /* larger than possible */
  mindist=box_l[0] + box_l[1] + box_l[2];
  for (i=0; i<n_total_particles; i++) {
    if (pid != partCfg[i].r.identity) {
      dx = xt - partCfg[i].r.p[0];   dx -= dround(dx/box_l[0])*box_l[0];
      dy = yt - partCfg[i].r.p[1];   dy -= dround(dy/box_l[1])*box_l[1];
      dz = zt - partCfg[i].r.p[2];   dz -= dround(dz/box_l[2])*box_l[2];
      mindist = dmin(mindist, SQR(dx)+SQR(dy)+SQR(dz));
    }
  }
  return sqrt(mindist);
}

double calc_re()
{
  int i;
  double dx, dy, dz;
  double dist = 0;

  for (i=0; i<chain_n_chains; i++) {
    dx = partCfg[chain_start+i*chain_length + chain_length-1].r.p[0]
      - partCfg[chain_start+i*chain_length].r.p[0];
    dy = partCfg[chain_start+i*chain_length + chain_length-1].r.p[1]
      - partCfg[chain_start+i*chain_length].r.p[1];
    dz = partCfg[chain_start+i*chain_length + chain_length-1].r.p[2]
      - partCfg[chain_start+i*chain_length].r.p[2];
    dist += (SQR(dx) + SQR(dy) + SQR(dz));
  }
  return sqrt(dist/((double)chain_n_chains));
}

double calc_re_av()
{
  int i,j;
  double dx,dy,dz;
  double dist = 0.0;

  for (j=0; j<n_configs; j++) {
    for (i=0; i<chain_n_chains; i++) {
      dx = configs[j][3*(chain_start+i*chain_length + chain_length-1)]     
	- configs[j][3*(chain_start+i*chain_length)];
      dy = configs[j][3*(chain_start+i*chain_length + chain_length-1) + 1] 
	- configs[j][3*(chain_start+i*chain_length) + 1];
      dz = configs[j][3*(chain_start+i*chain_length + chain_length-1) + 2] 
	- configs[j][3*(chain_start+i*chain_length) + 2];
      dist += (SQR(dx) + SQR(dy) + SQR(dz));
    }
  }
  return sqrt(dist/((double)chain_n_chains*n_configs));
}

double calc_rg()
{
  int i, j;
  double dx, dy, dz;
  double r_CM_x,r_CM_y,r_CM_z, r_G=0.0, doubMPC;

  doubMPC = (double)chain_length;
  for (i=0; i<chain_n_chains; i++) {
    r_CM_x = r_CM_y = r_CM_z = 0.0;
    for (j=0; j<chain_length; j++) {
      r_CM_x += partCfg[chain_start+i*chain_length + j].r.p[0];
      r_CM_y += partCfg[chain_start+i*chain_length + j].r.p[1];
      r_CM_z += partCfg[chain_start+i*chain_length + j].r.p[2];
    }
    r_CM_x /= doubMPC;
    r_CM_y /= doubMPC;
    r_CM_z /= doubMPC;
    for (j=0; j<chain_length; ++j) {
      dx = partCfg[chain_start+i*chain_length + j].r.p[0] - r_CM_x;
      dy = partCfg[chain_start+i*chain_length + j].r.p[1] - r_CM_y;
      dz = partCfg[chain_start+i*chain_length + j].r.p[2] - r_CM_z;
      r_G += (SQR(dx) + SQR(dy) + SQR(dz));
    }
  }
  return sqrt(r_G / (double)(chain_length*chain_n_chains));
}

double calc_rg_av()
{
  int i, j, k;
  double dx,dy,dz;
  double r_CM_x,r_CM_y,r_CM_z, r_G=0.0, doubMPC;

  doubMPC = (double)chain_length;
  for (k=0; k<n_configs; k++) {
    for (i=0; i<chain_n_chains; i++) {
      r_CM_x = r_CM_y = r_CM_z = 0.0;
      for (j=0; j<chain_length; j++) {
	r_CM_x += configs[k][3*(chain_start+i*chain_length + j)];
	r_CM_y += configs[k][3*(chain_start+i*chain_length + j) + 1];
	r_CM_z += configs[k][3*(chain_start+i*chain_length + j) + 2];
      }
      r_CM_x /= doubMPC; r_CM_y /= doubMPC; r_CM_z /= doubMPC;
      for (j=0; j<chain_length; ++j) {
	dx = configs[k][3*(chain_start+i*chain_length + j)]     - r_CM_x;
	dy = configs[k][3*(chain_start+i*chain_length + j) + 1] - r_CM_y;
	dz = configs[k][3*(chain_start+i*chain_length + j) + 2] - r_CM_z;
	r_G += (SQR(dx) + SQR(dy) + SQR(dz));
      }
    }
  }
  return sqrt(r_G / (double)(chain_length*chain_n_chains*n_configs));
}

double calc_rh()
{
  int i, j, p;
  double dx, dy, dz;
  double rh=0.0, ri=0.0;

  for(p=0;p<chain_n_chains;p++) {
    ri=0.0;
    for(i=chain_start+chain_length*p;i<chain_start+chain_length*(p+1);i++) {
      for(j=i+1;j<chain_start+chain_length*(p+1);j++) {
	dx = partCfg[i].r.p[0]-partCfg[j].r.p[0]; dx*=dx;
	dy = partCfg[i].r.p[1]-partCfg[j].r.p[1]; dy*=dy;
	dz = partCfg[i].r.p[2]-partCfg[j].r.p[2]; dz*=dz;
	ri += 1.0/sqrt(dx+dy+dz);
      }
    }
    rh += 1.0/ri;
  }
  return rh * 0.5*(chain_length*(chain_length-1)) / (double)chain_n_chains;
}

double calc_rh_av()
{
  int i, j, p, k;
  double dx, dy, dz;
  double rh=0.0, ri=0.0;

  for(k=0; k<n_configs; k++) {
    for(p=0;p<chain_n_chains;p++) {
      ri=0.0;
      for(i=chain_start+chain_length*p;i<chain_start+chain_length*(p+1);i++) {
	for(j=i+1;j<chain_start+chain_length*(p+1);j++) {
	  dx = configs[k][3*i]  -configs[k][3*j];
	  dy = configs[k][3*i+1]-configs[k][3*j+1];
	  dz = configs[k][3*i+2]-configs[k][3*j+2];
	  ri += 1.0/sqrt(dx*dx + dy*dy + dz*dz);
	}
      }
      rh += 1.0/ri;
    }
  }
  return rh * 0.5*(chain_length*(chain_length-1)) / ((double)chain_n_chains*n_configs);
}

void init_g123()
{
  int i, j, p;
  double cm_tmp[3];

  /* Save particles' current positions 
     (which'll be used as initial position later on) */
  partCoord_g = (float *) realloc(partCoord_g, 3*n_total_particles*sizeof(float));
  partCM_g = (float *) realloc(partCM_g, 3*chain_n_chains*sizeof(float));
  n_part_g = n_total_particles;
  n_chains_g = chain_n_chains;
  for(j=0; j<chain_n_chains; j++) {
    cm_tmp[0] = cm_tmp[1] = cm_tmp[2] = 0.0;
    for(i=0; i<chain_length; i++) {
      p = chain_start+j*chain_length + i;
      partCoord_g[3*p]   = partCfg[p].r.p[0]; cm_tmp[0]+=partCfg[p].r.p[0];
      partCoord_g[3*p+1] = partCfg[p].r.p[1]; cm_tmp[1]+=partCfg[p].r.p[1];
      partCoord_g[3*p+2] = partCfg[p].r.p[2]; cm_tmp[2]+=partCfg[p].r.p[2];
    }
    partCM_g[3*j]   = cm_tmp[0]/(1.*chain_length);
    partCM_g[3*j+1] = cm_tmp[1]/(1.*chain_length);
    partCM_g[3*j+2] = cm_tmp[2]/(1.*chain_length);
  }
}

void calc_g123(double *_g1, double *_g2, double *_g3)
{
  /* - Mean square displacement of a monomer
     - Mean square displacement in the center of gravity of the chain itself
     - Motion of the center of mass */
  int i, j, p;
  double g1=0.0, g2=0.0, g3=0.0, cm_tmp[3];

  for(j=0; j<chain_n_chains; j++) {
    cm_tmp[0] = cm_tmp[1] = cm_tmp[2] = 0.0;
    for(i=0; i<chain_length; i++) {
      p = chain_start+j*chain_length + i;
      cm_tmp[0]+=partCfg[p].r.p[0];
      cm_tmp[1]+=partCfg[p].r.p[1];
      cm_tmp[2]+=partCfg[p].r.p[2];
    }
    cm_tmp[0] /= (1.*chain_length);
    cm_tmp[1] /= (1.*chain_length);
    cm_tmp[2] /= (1.*chain_length);
    for(i=0; i<chain_length; i++) {
      p = chain_start+j*chain_length + i;
      g1 += SQR(partCfg[p].r.p[0]-partCoord_g[3*p])
	+ SQR(partCfg[p].r.p[1]-partCoord_g[3*p+1])
	+ SQR(partCfg[p].r.p[2]-partCoord_g[3*p+2]);
      g2 += SQR( (partCfg[p].r.p[0]-partCoord_g[3*p])
		 - (cm_tmp[0]-partCM_g[3*j]  ) ) 
	  + SQR( (partCfg[p].r.p[1]-partCoord_g[3*p+1])
		 - (cm_tmp[1]-partCM_g[3*j+1]) ) 
	  + SQR( (partCfg[p].r.p[2]-partCoord_g[3*p+2])
		 - (cm_tmp[2]-partCM_g[3*j+2]) );
    }
    g3 += SQR(cm_tmp[0]-partCM_g[3*j])
      + SQR(cm_tmp[1]-partCM_g[3*j+1])
      + SQR(cm_tmp[2]-partCM_g[3*j+2]);
  }
  *_g1 = g1 / (1.*chain_n_chains*chain_length);
  *_g2 = g2 / (1.*chain_n_chains*chain_length);
  *_g3 = g3 / (1.*chain_n_chains);
}

void calc_g1_av(double **_g1) {
  int i, j, p, t,k;
  double *g1=NULL;
  *_g1 = g1 = realloc(g1,n_configs*sizeof(double));

  for(k=0; k < n_configs; k++) {
    g1[k] = 0.0;
    for(t=0; t < n_configs-k; t++) {
      for(j=0; j<chain_n_chains; j++) {
	for(i=0; i<chain_length; i++) {
	  p = chain_start+j*chain_length + i;
	  g1[k] += SQR(configs[t+k][3*p]-configs[t][3*p])
	    + SQR(configs[t+k][3*p+1]-configs[t][3*p+1])
	    + SQR(configs[t+k][3*p+2]-configs[t][3*p+2]);
	}
      }
    }
    g1[k] /= ((double)chain_n_chains*chain_length*(n_configs-k));
  }
}

void calc_g2_av(double **_g2) {
  int i, j, p, t,k;
  double *g2=NULL, cm_tmp[3];
  *_g2 = g2 = realloc(g2,n_configs*sizeof(double));

  for(k=0; k < n_configs; k++) {
    g2[k] = 0.0;
    for(t=0; t < n_configs-k; t++) {
      for(j=0; j<chain_n_chains; j++) {
	cm_tmp[0] = cm_tmp[1] = cm_tmp[2] = 0.0;
	for(i=0; i<chain_length; i++) {
	  p = chain_start+j*chain_length + i;
	  cm_tmp[0] += configs[t+k][3*p]   - configs[t][3*p];
	  cm_tmp[1] += configs[t+k][3*p+1] - configs[t][3*p+1];
	  cm_tmp[2] += configs[t+k][3*p+2] - configs[t][3*p+2];
	}
	cm_tmp[0] /= (1.*chain_length);	cm_tmp[1] /= (1.*chain_length);	cm_tmp[2] /= (1.*chain_length);
	for(i=0; i<chain_length; i++) {
	  p = chain_start+j*chain_length + i;
	  g2[k] += SQR( (configs[t+k][3*p]-configs[t][3*p]) - cm_tmp[0] )
	    + SQR( (configs[t+k][3*p+1]-configs[t][3*p+1])  - cm_tmp[1] ) 
	    + SQR( (configs[t+k][3*p+2]-configs[t][3*p+2])  - cm_tmp[2] );
	}
      }
    }
    g2[k] /= ((double)chain_n_chains*chain_length*(n_configs-k));
  }
}

void calc_g3_av(double **_g3) {
  int i, j, p, t,k;
  double *g3=NULL, cm_tmp[3];
  *_g3 = g3 = realloc(g3,n_configs*sizeof(double));

  for(k=0; k < n_configs; k++) {
    g3[k] = 0.0;
    for(t=0; t < n_configs-k; t++) {
      for(j=0; j<chain_n_chains; j++) {
	cm_tmp[0] = cm_tmp[1] = cm_tmp[2] = 0.0;
	for(i=0; i<chain_length; i++) {
	  p = chain_start+j*chain_length + i;
	  cm_tmp[0] += configs[t+k][3*p]   - configs[t][3*p];
	  cm_tmp[1] += configs[t+k][3*p+1] - configs[t][3*p+1];
	  cm_tmp[2] += configs[t+k][3*p+2] - configs[t][3*p+2];
	}
	g3[k] += SQR(cm_tmp[0] / (1.*chain_length)) 
	  + SQR(cm_tmp[1] / (1.*chain_length)) 
	  + SQR(cm_tmp[2] / (1.*chain_length));
      }
    }
    g3[k] /= ((double)chain_n_chains*chain_length*(n_configs-k));
  }
}

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
	/* bonded interaction energies */
	i=0;
	while(i<p1->bl.n) {
	  type_num = p1->bl.e[i];
	  switch(bonded_ia_params[type_num].type) {
	  case BONDED_IA_FENE:
	    energy.node.e[type_num + 2] +=
	      fene_pair_energy(p1, checked_particle_ptr(p1->bl.e[i+1]), type_num);
	    i+=2; break;
	  case BONDED_IA_ANGLE:
	    energy.node.e[type_num + 2] +=
	      angle_energy(p1, checked_particle_ptr(p1->bl.e[i+1]),
			   checked_particle_ptr(p1->bl.e[i+2]), type_num);
	    i+=3; break;
	  default :
	    fprintf(stderr,"WARNING: Bonds of atom %d unknown\n",p1->r.identity);
	    i = p1->bl.n; 
	    break;
	  }
	}
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
	  type_num = s_non_bonded + type2 - type1;
	  while(type1>0) {
	    type_num += n_particle_types - (type1-1) ;
	    type1--;
	  }
	  /* distance calculation */
	  for(j=0; j<3; j++) d[j] = p1->r.p[j] - p2->r.p[j];
	  dist2 = SQR(d[0]) + SQR(d[1]) + SQR(d[2]);
	  dist  = sqrt(dist2);
	  
	  /* lennnard jones */
	  energy.node.e[type_num] += lj_pair_energy(p1,p2,ia_params,d,dist);
	  
	  /* real space coulomb */
	  if(coulomb.method==COULOMB_P3M) 
	    energy.node.e[s_coulomb+1] += p3m_coulomb_pair_energy(p1,p2,d,dist2,dist);
	  else if(coulomb.method==COULOMB_DH)
	    energy.node.e[s_coulomb] += dh_coulomb_pair_energy(p1,p2,dist);

	  /* minimal particle distance calculation */
	  if (dist < minimum_part_dist)
	    minimum_part_dist = dist;
	} 
      }
    }
  }

  /* calculate k-space part of electrostatic interaction. */ 
  if(coulomb.method==COULOMB_P3M && (energy.ana_num == 0 || energy.ana_num >= s_coulomb) ) {
    energy.node.e[s_coulomb+2] = P3M_calc_kspace_forces(0,1);
    energy.node.e[s_coulomb] = energy.node.e[s_coulomb+1]+energy.node.e[s_coulomb+2];
  }

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
    if(coulomb.method==COULOMB_P3M) {
      for(i=1;i<energy.n-2;i++)
	energy.sum.e[0] += energy.sum.e[i];
    } else {
      for(i=1;i<energy.n;i++)
	energy.sum.e[0] += energy.sum.e[i];
    }
  }

}

void init_energies()
{
  energy.n_pre        = 2;
  energy.n_bonded     = n_bonded_ia;
  energy.n_non_bonded = (n_particle_types*(n_particle_types+1))/2;
  if(coulomb.bjerrum > 0.0)         energy.n_coulomb =  1;
  if(coulomb.method==COULOMB_P3M) energy.n_coulomb += 2;
  energy.n = energy.n_pre+energy.n_bonded+energy.n_non_bonded+energy.n_coulomb;
  realloc_doublelist(&(energy.node),energy.n);
  realloc_doublelist(&(energy.sum),energy.n);
}

void analyze_append() {
  int i;
  configs = realloc(configs,(n_configs+1)*sizeof(double *));
  configs[n_configs] = (double *) malloc(3*n_total_particles*sizeof(double));
  for(i=0; i<n_total_particles; i++) {
    configs[n_configs][3*i]   = partCfg[i].r.p[0];
    configs[n_configs][3*i+1] = partCfg[i].r.p[1];
    configs[n_configs][3*i+2] = partCfg[i].r.p[2];
  }
  n_configs++;
}

void analyze_push() {
  int i;
  free(configs[0]);
  for(i=0; i<n_configs-1; i++) {
    configs[i]=configs[i+1];
  }
  configs[n_configs-1] = (double *) malloc(3*n_total_particles*sizeof(double));
  for(i=0; i<n_total_particles; i++) {
    configs[n_configs-1][3*i]   = partCfg[i].r.p[0];
    configs[n_configs-1][3*i+1] = partCfg[i].r.p[1];
    configs[n_configs-1][3*i+2] = partCfg[i].r.p[2];
  }
}

void analyze_replace(int ind) {
  int i;
  for(i=0; i<n_total_particles; i++) {
    configs[ind][3*i]   = partCfg[i].r.p[0];
    configs[ind][3*i+1] = partCfg[i].r.p[1];
    configs[ind][3*i+2] = partCfg[i].r.p[2];
  }
}

void analyze_remove(int ind) {
  int i;
  free(configs[ind]);
  for(i=ind; i<n_configs-1; i++) {
    configs[i]=configs[i+1];
  }
  n_configs--;
  configs = realloc(configs,n_configs*sizeof(double *));
}
