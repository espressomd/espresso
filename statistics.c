/** \file statistics.c
    This is the place for analysation (so far...).
    Implementation of \ref statistics.h "statistics.h"
 */
#include <stdlib.h>
#include <string.h>
#include "statistics.h"
#include "forces.h"
#include "communication.h"
#include "grid.h"
#include "integrate.h"
#include "particle_data.h"
#include "debug.h"

/** Particles' initial positions (needed for g1(t), g2(t), g3(t) in \ref analyze) */
float *partCoord_g=NULL, *partCM_g=NULL;

/** data for a system consisting of chains */
int chain_start = 0, chain_n_chains = 0, chain_length = 0;

static int get_chain_structure_info(Tcl_Interp *interp, int *argc, char ***argv)
{
  if (*argc == 0)
    return TCL_OK;
  if (*argc < 3) {
    Tcl_AppendResult(interp, "chain structure info consists of <start> <n> <length>", (char *)NULL);    
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, (*argv)[0], &chain_start) == TCL_ERROR ||
      Tcl_GetInt(interp, (*argv)[1], &chain_n_chains) == TCL_ERROR ||
      Tcl_GetInt(interp, (*argv)[2], &chain_length) == TCL_ERROR)
    return TCL_ERROR;

  *argc -= 3;
  *argv += 3;

  return TCL_OK;
}
      
int analyze(ClientData data, Tcl_Interp *interp, int argc, char **argv)
{
  char buffer[50 + TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
  char *mode;
  int  i, j, p;
  double *buf;
  double mindist, xt, yt, zt, r_catch;
  double dx, dy, dz;

  if (argc < 2) {
    Tcl_AppendResult(interp, "Wrong # of args! Usage: analyze <what> ...", (char *)NULL);
    return (TCL_ERROR);
  }

  mode = argv[1];
  argc -= 2;
  argv += 2;

  if (!strncmp(mode, "set", strlen(mode))) {
    if (argc == 0) {
      Tcl_AppendResult(interp, "which topology are you interested in?", (char *)NULL);
      return TCL_ERROR;
    }
    if (!strncmp(argv[0], "chains", strlen(argv[0]))) {
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
      if (get_chain_structure_info(interp, &argc, &argv) == TCL_ERROR)
	return TCL_ERROR;
      return TCL_OK;
    }
    else {
      Tcl_AppendResult(interp, "The topology \"", argv[0],
		       "\" you requested is not implemented.", (char *)NULL);
      return (TCL_ERROR);
    }
  }
  else if (!strncmp(mode, "mindist", strlen(mode))) {
    if (n_total_particles <= 1) {
      Tcl_AppendResult(interp, "(not enough particles)",
		       (char *)NULL);
      return (TCL_OK);
    }
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
    Tcl_PrintDouble(interp, mindist, buffer);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
    return TCL_OK;
  }
  else if (!strncmp(mode, "nbhood", strlen(mode))) {
    Particle ref;

    if (n_total_particles == 0) {
      Tcl_AppendResult(interp, "(no particles)",
		       (char *)NULL);
      return (TCL_OK);
    }

    updatePartCfg();

    r_catch = atof(argv[1]);
    p = atoi(argv[0]);

    if (get_particle_data(p, &ref) != TCL_OK) {
      Tcl_AppendResult(interp, "particle does not exist", (char *)NULL);
      return TCL_ERROR;
    }

    xt = ref.r.p[0]; yt = ref.r.p[1]; zt = ref.r.p[2];

    for (i = 0; i<n_total_particles; i++) {
      if (partCfg[i].r.identity != p) {
	dx = xt - partCfg[i].r.p[0];   dx -= dround(dx/box_l[0])*box_l[0];
	dy = yt - partCfg[i].r.p[1];   dy -= dround(dy/box_l[1])*box_l[1];
	dz = zt - partCfg[i].r.p[2];   dz -= dround(dz/box_l[2])*box_l[2];
	if (sqrt(SQR(dx)+SQR(dy)+SQR(dz)) < r_catch) {
	  sprintf(buffer, "%d ", partCfg[i].r.identity);
	  Tcl_AppendResult(interp, buffer, (char *)NULL);
	}
      }
    }
    return (TCL_OK);
  }
  else if (!strncmp(mode, "distto", strlen(mode))) {
    xt = atof(argv[0]); yt = atof(argv[1]); zt = atof(argv[2]);

    if (n_total_particles == 0) {
      Tcl_AppendResult(interp, "(no particles)",
		       (char *)NULL);
      return (TCL_OK);
    }

    updatePartCfg();

    /* larger than possible */
    mindist=box_l[0] + box_l[1] + box_l[2];
    for (i=0; i<n_total_particles; i++) {
      dx = xt - partCfg[i].r.p[0];   dx -= dround(dx/box_l[0])*box_l[0];
      dy = yt - partCfg[i].r.p[1];   dy -= dround(dy/box_l[1])*box_l[1];
      dz = zt - partCfg[i].r.p[2];   dz -= dround(dz/box_l[2])*box_l[2];
      mindist = dmin(mindist, SQR(dx)+SQR(dy)+SQR(dz));
    }
    Tcl_PrintDouble(interp, sqrt(mindist), buffer);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
    return (TCL_OK);
  }
  else if (!strncmp(mode, "re", strlen(mode))) {
    double dist = 0;

    /* Averaged quadratic end-to-end-distance of the polymer chains */
    if (get_chain_structure_info(interp, &argc, &argv) == TCL_ERROR)
      return TCL_ERROR;
    if (argc != 0) {
      Tcl_AppendResult(interp, "only chain structure info required", (char *)NULL);
      return TCL_ERROR;
    }
    if (n_total_particles < chain_n_chains*chain_length) {
      Tcl_AppendResult(interp, "not enough particles for chain structure", (char *)NULL);
      return TCL_ERROR;
    }

    if (!sortPartCfg()) {
      Tcl_AppendResult(interp, "for analyze, store particles consecutively starting with 0.",
		       (char *) NULL);
      return (TCL_ERROR);      
    }
    
    for (i=0; i<chain_n_chains; i++) {
      dx = partCfg[chain_start+i*chain_length + chain_length-1].r.p[0] - partCfg[chain_start+i*chain_length].r.p[0];
      dy = partCfg[chain_start+i*chain_length + chain_length-1].r.p[1] - partCfg[chain_start+i*chain_length].r.p[1];
      dz = partCfg[chain_start+i*chain_length + chain_length-1].r.p[2] - partCfg[chain_start+i*chain_length].r.p[2];
      dist += (SQR(dx) + SQR(dy) + SQR(dz));
    }
    Tcl_PrintDouble(interp, sqrt(dist/((double)chain_n_chains)), buffer);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
    return (TCL_OK);
  }
  else if (!strncmp(mode, "rg", strlen(mode))) {
    double r_CM_x=0.0, r_CM_y=0.0, r_CM_z=0.0, r_G=0.0, doubMPC;

    /* Averaged radius of gyration */
    if (get_chain_structure_info(interp, &argc, &argv) == TCL_ERROR)
      return TCL_ERROR;
    if (argc != 0) {
      Tcl_AppendResult(interp, "only chain structure info required", (char *)NULL);
      return TCL_ERROR;
    }
    if (n_total_particles < chain_n_chains*chain_length) {
      Tcl_AppendResult(interp, "not enough particles for chain structure", (char *)NULL);
      return TCL_ERROR;
    }

    if (!sortPartCfg()) {
      Tcl_AppendResult(interp, "for analyze, store particles consecutively starting with 0.",
		       (char *) NULL);
      return (TCL_ERROR);      
    }

    doubMPC = (double)chain_length;
    for (i=0; i<chain_n_chains; i++) {
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
    r_G /= (double)(chain_length*chain_n_chains);
    Tcl_PrintDouble(interp, sqrt(r_G), buffer);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
    return (TCL_OK);
  }
  else if (!strncmp(mode, "rh", strlen(mode))) {
    double rh=0.0, ri=0.0;

    /* Averaged hydrodynamic radius */
    if (get_chain_structure_info(interp, &argc, &argv) == TCL_ERROR)
      return TCL_ERROR;
    if (argc != 0) {
      Tcl_AppendResult(interp, "only chain structure info required", (char *)NULL);
      return TCL_ERROR;
    }
    if (n_total_particles < chain_n_chains*chain_length) {
      Tcl_AppendResult(interp, "not enough particles for chain structure", (char *)NULL);
      return TCL_ERROR;
    }

    if (!sortPartCfg()) {
      Tcl_AppendResult(interp, "for analyze, store particles consecutively starting with 0.",
		       (char *) NULL);
      return (TCL_ERROR);      
    }

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
    rh *= 0.5*(chain_length*(chain_length-1));
    rh /= (double)chain_n_chains;
    Tcl_PrintDouble(interp, rh, buffer);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
    return (TCL_OK);
  }
  else if (!strncmp(mode, "g123", strlen(mode))) {
    double g1=0.0, g2=0.0, g3=0.0, cm_tmp[3];
    int init = 0;
    /* - Mean square displacement of a monomer
       - Mean square displacement in the center of gravity of the chain itself
       - Motion of the center of mass */

    if (!sortPartCfg()) {
      Tcl_AppendResult(interp, "for analyze, store particles consecutively starting with 0.",
		       (char *) NULL);
      return (TCL_ERROR);      
    }

    if (argc == 1 && !strncmp(argv[0], "-init", strlen(argv[0]))) {
      init = 1;
      argc--; argv++;
    }

    if (get_chain_structure_info(interp, &argc, &argv) == TCL_ERROR)
      return TCL_ERROR;
    if (argc != 0) {
      Tcl_AppendResult(interp, "only chain structure info required", (char *)NULL);
      return TCL_ERROR;
    }

    if (n_total_particles < chain_n_chains*chain_length) {
      Tcl_AppendResult(interp, "not enough particles for chain structure", (char *)NULL);
      return TCL_ERROR;
    }

    if (init) {
      /* Save particles' current positions 
	 (which'll be used as initial position later on) */
      partCoord_g = (float *) realloc(partCoord_g, 3*n_total_particles*sizeof(float));
      partCM_g = (float *) realloc(partCM_g, 3*chain_n_chains*sizeof(float));
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
      return TCL_OK;
    }

    if (partCoord_g == NULL || partCM_g == NULL) {
      Tcl_AppendResult(interp, "please call with -init first", (char *)NULL);
      return TCL_ERROR;
    }

    for(j=0; j<chain_n_chains; j++) {
      cm_tmp[0] = cm_tmp[1] = cm_tmp[2] = 0.0;
      for(i=0; i<chain_length; i++) {
	p = chain_start+j*chain_length + i;
	cm_tmp[0]+=partCfg[p].r.p[0]; cm_tmp[1]+=partCfg[p].r.p[1]; cm_tmp[2]+=partCfg[p].r.p[2];
      }
      cm_tmp[0] /= (1.*chain_length); cm_tmp[1] /= (1.*chain_length); cm_tmp[2] /= (1.*chain_length);
      for(i=0; i<chain_length; i++) {
	p = chain_start+j*chain_length + i;
	g1 += SQR(partCfg[p].r.p[0]-partCoord_g[3*p]) + SQR(partCfg[p].r.p[1]-partCoord_g[3*p+1])
	  + SQR(partCfg[p].r.p[2]-partCoord_g[3*p+2]);
	g2 += SQR( (partCfg[p].r.p[0]-partCoord_g[3*p]) - (cm_tmp[0]-partCM_g[3*j]  ) ) 
	  + SQR( (partCfg[p].r.p[1]-partCoord_g[3*p+1]) - (cm_tmp[1]-partCM_g[3*j+1]) ) 
	  + SQR( (partCfg[p].r.p[2]-partCoord_g[3*p+2]) - (cm_tmp[2]-partCM_g[3*j+2]) );
	g3 += SQR(cm_tmp[0]-partCM_g[3*j]) + SQR(cm_tmp[1]-partCM_g[3*j+1]) + SQR(cm_tmp[2]-partCM_g[3*j+2]);
      }
    }
    g1 /= (1.*chain_n_chains*chain_length); g2 /= (1.*chain_n_chains*chain_length); g3 /= (1.*chain_n_chains);
    sprintf(buffer,"{%f %f %f}",g1, g2, g3);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
    return (TCL_OK);
  }
  else {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "The operation \"", mode, "\" you requested is not implemented.", (char *)NULL);
    return (TCL_ERROR);
  }
}
