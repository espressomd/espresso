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

/** Particles' initial positions (needed for g1(t), g2(t), g3(t) in \ref analyze) */
float *partCoord_g=NULL, *partCM_g=NULL;


int mindist(ClientData data, Tcl_Interp *interp, int argc, char **argv)
{
  char buffer[TCL_DOUBLE_SPACE + 1];
  double *buf;
  double mindist;
  int i;

  if (argc != 1) {
    Tcl_AppendResult(interp, "mindist takes no arguments", (char *)NULL);
    return (TCL_ERROR);
  }

  if (minimum_part_dist == -1) {
    Tcl_AppendResult(interp, "(not yet set)", (char *)NULL);
    return (TCL_OK);
  }

  buf = malloc(n_nodes*sizeof(double));
  mpi_gather_stats(0, buf);
  mindist = buf[0];
  for (i = 1; i < n_nodes; i++) {
    if (buf[i] < mindist)
      mindist = buf[i];
  }
  free(buf);

  if (mindist >= box_l[0] + box_l[1] + box_l[2]) {
    Tcl_PrintDouble(interp, max_range, buffer);
    Tcl_AppendResult(interp, "> ", buffer, (char *)NULL);
    return (TCL_OK);
  }

  Tcl_PrintDouble(interp, mindist, buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return (TCL_OK);
}


int analyze(ClientData data, Tcl_Interp *interp, int argc, char **argv)
{
  char buffer[50 + TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
  Particle *partCfg;
  int N_P, MPC, N_CI, N_pS, N_nS;
  int mode, arg_i, i, j, p;
  double dx, dy, dz, dist = 0.0;
  double r_CM_x=0.0, r_CM_y=0.0, r_CM_z=0.0, r_G=0.0, doubMPC;
  double rh=0.0, ri=0.0;
  double g1=0.0, g2=0.0, g3=0.0, cm_tmp[3];

  /* If no further arguments are given, perform a self-consistency-test */
  if (argc == 1) {
    sprintf(buffer,"%d",n_total_particles);
    Tcl_AppendResult(interp, "Hello World! You have ", buffer, " particles!", (char *)NULL);
    mpi_gather_stats(2, &partCfg);
    for(i=0; i<n_total_particles; i++) {
      if(partCfg[i].r.identity!=i) { sprintf(buffer,"%d != %d \n",partCfg[i].r.identity,i); Tcl_AppendResult(interp, buffer, (char *)NULL); }
      sprintf(buffer,"%d: %d %d %f (%f, %f, %f)\n",i,partCfg[i].r.identity,partCfg[i].r.type,partCfg[i].r.q,partCfg[i].r.p[0],partCfg[i].r.p[1],partCfg[i].r.p[2]);
      Tcl_AppendResult(interp,buffer, (char *)NULL);
    }
    return (TCL_OK);
  }

  /* Now check for the operation the user wants to be done */
  mode = atol(argv[1]);
  if (mode<0) {
    Tcl_ResetResult(interp);
    sprintf(buffer,"The operation %d you requested does not exist.",mode);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
    return (TCL_ERROR);
  }

  /* The user has to provide some informations on the topology, such as
     number of polymer chains N_P and monomers per chain MPC;
     let's check if he complied! */
  if (argc < 4) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "Wrong # of args! Usage: analyze <what> <N_P> <MPC> [<N_CI> [<N_pS> <N_nS>]] [-g]", (char *)NULL);
    return (TCL_ERROR);
  }
  N_P = atol(argv[2]);
  MPC = atol(argv[3]);
  arg_i=4;
  if ((argc > 4) && (strncmp(argv[4], "-g", strlen(argv[4])))) {
    N_CI = atol(argv[4]); 
    arg_i=5; }
  else {
    N_CI = 0; }
  if ((argc > 6) && (strncmp(argv[5], "-g", strlen(argv[5]))) && (strncmp(argv[6], "-g", strlen(argv[6])))) {
    N_pS = atol(argv[5]); N_nS = atol(argv[6]); 
    arg_i=7; }
  else {
    N_pS = 0; N_nS = 0; }
  if (N_P*MPC + N_CI + N_pS + N_nS != n_total_particles) {
    Tcl_ResetResult(interp);
    sprintf(buffer,"The topology you specified does not match the current configuration (expected: %d particles, got: %d)!",n_total_particles,N_P*MPC+N_CI+N_pS+N_nS);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
    return (TCL_ERROR);
  }

  /* For mode==3 it is necessary to know the particles' initial positions;
     these are either saved automatically on the procedures' first run (since *partCoord_g==NULL at that time),
     or they are set to the current configuration if the [-g] option is specified. */
  if ((argc > arg_i) && !(strncmp(argv[arg_i], "-g", strlen(argv[arg_i])))) {
    if (mode == 3) {
      mode = -3; }
    else {
      Tcl_ResetResult(interp);
      sprintf(buffer,"The [-g] option you specified is only supported by mode 3, not by mode %d!",mode);
      Tcl_AppendResult(interp, buffer, (char *)NULL);
      return (TCL_ERROR);
    }
  } else if((mode==3) && ((partCoord_g==NULL) || (partCM_g==NULL))) {
    mode = -3; 
  }
      

  /* Get the complete informations on all particles 
     (this is supposed to be done by on_conf_change later on) */
  mpi_gather_stats(2, &partCfg);


  switch(mode) {
  case 0:
    /* Averaged quadratic end-to-end-distance of the polymer chains */
    for (i=0; i<N_P; i++) {
      dx = partCfg[(i+1)*MPC-1].r.p[0] - partCfg[i*MPC].r.p[0];
      dy = partCfg[(i+1)*MPC-1].r.p[1] - partCfg[i*MPC].r.p[0];
      dz = partCfg[(i+1)*MPC-1].r.p[2] - partCfg[i*MPC].r.p[0];
      dist += (SQR(dx) + SQR(dy) + SQR(dz));
    }
    Tcl_PrintDouble(interp, sqrt(dist/((double)N_P)), buffer);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
    return (TCL_OK);
    break;
  case 1:
    /* Averaged radius of gyration */
    doubMPC = (double)MPC;
    for (i=0; i<N_P; i++) {
      for (j=0; j<MPC; j++) {
	r_CM_x += partCfg[i*MPC+j].r.p[0];
	r_CM_y += partCfg[i*MPC+j].r.p[1];
	r_CM_z += partCfg[i*MPC+j].r.p[2];
      }
      r_CM_x /= doubMPC;
      r_CM_y /= doubMPC;
      r_CM_z /= doubMPC;
      for (j=0; j<MPC; ++j) {
	dx = partCfg[i*MPC+j].r.p[0]-r_CM_x;
	dy = partCfg[i*MPC+j].r.p[1]-r_CM_y;
	dz = partCfg[i*MPC+j].r.p[2]-r_CM_z;
	r_G += (SQR(dx) + SQR(dy) + SQR(dz));
      }
    }
    r_G /= (double)(MPC*N_P);
    Tcl_PrintDouble(interp, sqrt(r_G), buffer);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
    return (TCL_OK);
    break;
  case 2:
    /* Averaged hydrodynamic radius */
    for(p=0;p<N_P;p++) {
      ri=0.0;
      for(i=MPC*p;i<MPC*(p+1);i++) {
	for(j=i+1;j<MPC*(p+1);j++) {
	  dx = partCfg[i].r.p[0]-partCfg[j].r.p[0]; dx*=dx;
	  dy = partCfg[i].r.p[1]-partCfg[j].r.p[1]; dy*=dy;
	  dz = partCfg[i].r.p[2]-partCfg[j].r.p[2]; dz*=dz;
	  ri += 1.0/sqrt(dx+dy+dz);
	}
      }
      rh += 1.0/ri;
    }
    rh *= 0.5*(MPC*(MPC-1));
    rh /= (double)N_P;
    Tcl_PrintDouble(interp, rh, buffer);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
    return (TCL_OK);
    break;
  case -3:
    /* Save particles' current positions 
       (which'll be used as initial position later on) */
    partCoord_g = (float *) realloc(partCoord_g, 3*n_total_particles*sizeof(float));
    partCM_g = (float *) realloc(partCM_g, 3*N_P*sizeof(float));
    for(j=0; j<N_P; j++) {
      cm_tmp[0] = cm_tmp[1] = cm_tmp[2] = 0.0;
      for(i=0; i<MPC; i++) {
	p = j*MPC+i;
	partCoord_g[3*p]   = partCfg[p].r.p[0]; cm_tmp[0]+=partCfg[p].r.p[0];
	partCoord_g[3*p+1] = partCfg[p].r.p[1]; cm_tmp[1]+=partCfg[p].r.p[1];
	partCoord_g[3*p+2] = partCfg[p].r.p[2]; cm_tmp[2]+=partCfg[p].r.p[2];
      }
      partCM_g[3*j]   = cm_tmp[0]/(1.*MPC);
      partCM_g[3*j+1] = cm_tmp[1]/(1.*MPC);
      partCM_g[3*j+2] = cm_tmp[2]/(1.*MPC);
    }
  case 3:
    /* - Mean square displacement of a monomer
       - Mean square displacement in the center of gravity of the chain itself
       - Motion of the center of mass */
    for(j=0; j<N_P; j++) {
      cm_tmp[0] = cm_tmp[1] = cm_tmp[2] = 0.0;
      for(i=0; i<MPC; i++) {
	p = j*MPC+i;
	cm_tmp[0]+=partCfg[p].r.p[0]; cm_tmp[1]+=partCfg[p].r.p[1]; cm_tmp[2]+=partCfg[p].r.p[2];
      }
      cm_tmp[0] /= (1.*MPC); cm_tmp[1] /= (1.*MPC); cm_tmp[2] /= (1.*MPC);
      for(i=0; i<MPC; i++) {
	p = j*MPC+i;
	g1 += SQR(partCfg[p].r.p[0]-partCoord_g[3*p]) + SQR(partCfg[p].r.p[1]-partCoord_g[3*p+1]) + SQR(partCfg[p].r.p[2]-partCoord_g[3*p+2]);
	g2 += SQR( (partCfg[p].r.p[0]-partCoord_g[3*p]) - (cm_tmp[0]-partCM_g[3*j]  ) ) 
	  + SQR( (partCfg[p].r.p[1]-partCoord_g[3*p+1]) - (cm_tmp[1]-partCM_g[3*j+1]) ) 
	  + SQR( (partCfg[p].r.p[2]-partCoord_g[3*p+2]) - (cm_tmp[2]-partCM_g[3*j+2]) );
	g3 += SQR(cm_tmp[0]-partCM_g[3*j]) + SQR(cm_tmp[1]-partCM_g[3*j+1]) + SQR(cm_tmp[2]-partCM_g[3*j+2]);
      }
    }
    g1 /= (1.*N_P*MPC); g2 /= (1.*N_P*MPC); g3 /= (1.*N_P*MPC);
    sprintf(buffer,"{%f %f %f}",g1, g2, g3);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
    return (TCL_OK);
    break;
  default:
    Tcl_ResetResult(interp);
    sprintf(buffer,"The operation %d you requested is not implemented.",mode);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
    return (TCL_ERROR);
  }

  return (TCL_OK);
}
