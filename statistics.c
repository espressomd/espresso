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
  Particle *gdata;
  int N_P, MPC, N_CI, N_pS, N_nS, N_tot;
  int mode, i, j, p;
  double dx, dy, dz, dist = 0.0;
  double r_CM_x=0.0, r_CM_y=0.0, r_CM_z=0.0, r_G=0.0, doubMPC;
  double rh=0.0, ri=0.0;

  /* If no further arguments are given, perform a self-consistency-test */
  if (argc == 1) {
    sprintf(buffer,"%d",N_tot);
    Tcl_AppendResult(interp, "Hello World! You have ", buffer, " particles!", (char *)NULL);
    mpi_gather_stats(2, &gdata);
    for(i=0; i<n_total_particles; i++) {
      if(gdata[i].r.identity!=i) { sprintf(buffer,"%d != %d \n",gdata[i].r.identity,i); Tcl_AppendResult(interp, buffer, (char *)NULL); }
      sprintf(buffer,"%d: %d %d %f (%f, %f, %f)\n",i,gdata[i].r.identity,gdata[i].r.type,gdata[i].r.q,gdata[i].r.p[0],gdata[i].r.p[1],gdata[i].r.p[2]);
      Tcl_AppendResult(interp,buffer, (char *)NULL);
    }
    return (TCL_OK);
  }

  /* The user has to provide some informations on the topology, such as
  number of polymer chains N_P and monomers per chain MPC;
  let's check if he complied! */
  if (argc < 4) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "Wrong # of args! Usage: analyze <what> <N_P> <MPC> [<N_CI> [<N_pS> <N_nS>]]", (char *)NULL);
    return (TCL_ERROR);
  }
  N_P = atol(argv[2]);
  MPC = atol(argv[3]);
  if (argc > 4)
    N_CI = atol(argv[4]);
  else
    N_CI = 0;
  if (argc > 6) {
    N_pS = atol(argv[5]); N_nS = atol(argv[6]); }
  else {
    N_pS = 0; N_nS = 0; }
  if (N_P*MPC + N_CI + N_pS + N_nS != N_tot) {
    Tcl_ResetResult(interp);
    sprintf(buffer,"The topology you specified does not match the current configuration (expected: %d particles, got: %d)!",N_tot,N_P*MPC+N_CI+N_pS+N_nS);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
    return (TCL_ERROR);
  }

  /* Now check for the operation the user want to be done */
  mode = atol(argv[1]);
  if (mode<0) {
    Tcl_ResetResult(interp);
    sprintf(buffer,"The operation %d you requested does not exist.",mode);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
    return (TCL_ERROR);
  }

  /* Get the complete informations on all particles */
  mpi_gather_stats(2, &gdata);

  switch(mode) {
  case 0:
    /* Averaged quadratic end-to-end-distance of the polymer chains */
    for (i=0; i<N_P; i++) {
      dx = gdata[(i+1)*MPC-1].r.p[0] - gdata[i*MPC].r.p[0];
      dy = gdata[(i+1)*MPC-1].r.p[1] - gdata[i*MPC].r.p[0];
      dz = gdata[(i+1)*MPC-1].r.p[2] - gdata[i*MPC].r.p[0];
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
	r_CM_x += gdata[i*MPC+j].r.p[0];
	r_CM_y += gdata[i*MPC+j].r.p[1];
	r_CM_z += gdata[i*MPC+j].r.p[2];
      }
      r_CM_x /= doubMPC;
      r_CM_y /= doubMPC;
      r_CM_z /= doubMPC;
      for (j=0; j<MPC; ++j) {
	dx = gdata[i*MPC+j].r.p[0]-r_CM_x;
	dy = gdata[i*MPC+j].r.p[1]-r_CM_y;
	dz = gdata[i*MPC+j].r.p[2]-r_CM_z;
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
	  dx = gdata[i].r.p[0]-gdata[j].r.p[0]; dx*=dx;
	  dy = gdata[i].r.p[1]-gdata[j].r.p[1]; dy*=dy;
	  dz = gdata[i].r.p[2]-gdata[j].r.p[2]; dz*=dz;
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
  default:
    Tcl_ResetResult(interp);
    sprintf(buffer,"The operation %d you requested is not implemented.",mode);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
    return (TCL_ERROR);
  }

  return (TCL_OK);
}
