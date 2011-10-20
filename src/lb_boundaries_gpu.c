/* $Id$
 *
 * This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
 * It is therefore subject to the ESPResSo license agreement which you
 * accepted upon receiving the distribution and by which you are
 * legally bound while utilizing this file in any form or way.
 * There is NO WARRANTY, not even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
 * You should have received a copy of that license along with this
 * program; if not, refer to http://www.espresso.mpg.de/license.html
 * where its current version can be found, or write to
 * Max-Planck-Institute for Polymer Research, Theory Group, 
 * PO Box 3148, 55021 Mainz, Germany. 
 * Copyright (c) 2002-2007; all rights reserved unless otherwise stated.
 */

/** \file lb_boundaries_gpu.c
 *
 * Boundary conditions for Lattice Boltzmann fluid dynamics.
 * Header file for \ref lb_boundaries_gpu.h.
 *
 */
#include "utils.h"
#include "constraint.h"
#include "lb_boundaries_gpu.h"
#include "lbgpu.h"
#include "interaction_data.h"
#include "communication.h"

#ifdef LB_BOUNDARIES_GPU
int n_lb_boundaries_gpu       = 0;
LB_boundary_gpu *lb_boundaries_gpu = NULL;


int printLbBoundaryToResult(Tcl_Interp *interp, int i)
{
  LB_boundary_gpu *lbb = &lb_boundaries_gpu[i];
  char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
  sprintf(buffer, "%d ", i);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  
  switch (lbb->type) {
  case LB_BOUNDARY_WAL:
    Tcl_PrintDouble(interp, lbb->c.wal.n[0], buffer);
    Tcl_AppendResult(interp, "wall normal ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, lbb->c.wal.n[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, lbb->c.wal.n[2], buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    Tcl_PrintDouble(interp, lbb->c.wal.d, buffer);
    Tcl_AppendResult(interp, " dist ", buffer, (char *) NULL);
    break;
  case LB_BOUNDARY_SPH:
    Tcl_PrintDouble(interp, lbb->c.sph.pos[0], buffer);
    Tcl_AppendResult(interp, "sphere center ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, lbb->c.sph.pos[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, lbb->c.sph.pos[2], buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    Tcl_PrintDouble(interp, lbb->c.sph.rad, buffer);
    Tcl_AppendResult(interp, " radius ", buffer, (char *) NULL);
    Tcl_PrintDouble(interp, lbb->c.sph.direction, buffer);
    Tcl_AppendResult(interp, " direction ", buffer, (char *) NULL);
    break;
  case LB_BOUNDARY_CYL:
    Tcl_PrintDouble(interp, lbb->c.cyl.pos[0], buffer);
    Tcl_AppendResult(interp, "cylinder center ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, lbb->c.cyl.pos[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, lbb->c.cyl.pos[2], buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    Tcl_PrintDouble(interp, lbb->c.cyl.axis[0], buffer);
    Tcl_AppendResult(interp, " axis ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, lbb->c.cyl.axis[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, lbb->c.cyl.axis[2], buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    Tcl_PrintDouble(interp, lbb->c.cyl.rad, buffer);
    Tcl_AppendResult(interp, " radius ", buffer, (char *) NULL);
    Tcl_PrintDouble(interp, lbb->c.cyl.length, buffer);
    Tcl_AppendResult(interp, " length ", buffer, (char *) NULL);
    Tcl_PrintDouble(interp, lbb->c.cyl.direction, buffer);
    Tcl_AppendResult(interp, " direction ", buffer, (char *) NULL);
    break;

  default:
    sprintf(buffer, "%d", lbb->type);
    Tcl_AppendResult(interp, "unknown lb_boundary type ", buffer, ".", (char *) NULL);
    return (TCL_OK);
  }

  return (TCL_OK);
}

int lb_boundary_print_all(Tcl_Interp *interp)
{
  int i;
  if(n_lb_boundaries_gpu>0) Tcl_AppendResult(interp, "{", (char *)NULL);
  for (i = 0; i < n_lb_boundaries_gpu; i++) {
    if(i>0) Tcl_AppendResult(interp, " {", (char *)NULL);
    printLbBoundaryToResult(interp, i);
    Tcl_AppendResult(interp, "}", (char *)NULL);
  }
  return (TCL_OK);
}
/**(re-)allocte memory for boundary structs */
LB_boundary_gpu *generate_lb_boundary()
{
  n_lb_boundaries_gpu++;
  lb_boundaries_gpu = realloc(lb_boundaries_gpu, n_lb_boundaries_gpu*sizeof(LB_boundary_gpu));
  lb_boundaries_gpu[n_lb_boundaries_gpu-1].type = LB_BOUNDARY_NONE;
  
  return &lb_boundaries_gpu[n_lb_boundaries_gpu-1];
}

int lb_boundary_wall(LB_boundary_gpu *lbb, Tcl_Interp *interp,
		    int argc, char **argv)
{
  int i;
  double norm;

  lbb->type = LB_BOUNDARY_WAL;
  /* invalid entries to start of */
  lbb->c.wal.n[0] = 
  lbb->c.wal.n[1] = 
  lbb->c.wal.n[2] = 0;
  lbb->c.wal.d = 0;
			  
  while (argc > 0) {
    if(!strncmp(argv[0], "normal", strlen(argv[0]))) {
      if(argc < 4) {
	Tcl_AppendResult(interp, "lb_boundary wall normal <nx> <ny> <nz> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if(Tcl_GetDouble(interp, argv[1], &(lbb->c.wal.n[0])) == TCL_ERROR ||
	Tcl_GetDouble(interp, argv[2], &(lbb->c.wal.n[1])) == TCL_ERROR ||
        Tcl_GetDouble(interp, argv[3], &(lbb->c.wal.n[2])) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 4; argv += 4;
    }
    else if(!strncmp(argv[0], "dist", strlen(argv[0]))) {
      if(argc < 1) {
	Tcl_AppendResult(interp, "lb_boundary wall dist <d> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if(Tcl_GetDouble(interp, argv[1], &(lbb->c.wal.d)) == TCL_ERROR)
	return (TCL_ERROR);
        argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "type", strlen(argv[0]))) {
      if(argc < 1) {
	Tcl_AppendResult(interp, "lb_boundary wall type <t> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      argc -= 2; argv += 2;
    }
    else
      break;
  }
  /* length of the normal vector */
  norm = SQR(lbb->c.wal.n[0])+SQR(lbb->c.wal.n[1])+SQR(lbb->c.wal.n[2]);
  if (norm < 1e-10) {
    Tcl_AppendResult(interp, "usage: lb_boundary wall normal <nx> <ny> <nz> dist <d> type <t>",
		     (char *) NULL);
    return (TCL_ERROR);    
  }
  /* normalize the normal vector */
  for (i=0;i<3;i++) lbb->c.wal.n[i] /= sqrt(norm);

  return (TCL_OK);
}

int lb_boundary_sphere(LB_boundary_gpu *lbb, Tcl_Interp *interp,
		    int argc, char **argv)
{
  lbb->type = LB_BOUNDARY_SPH;

  /* invalid entries to start of */
  lbb->c.sph.pos[0] = 
  lbb->c.sph.pos[1] = 
  lbb->c.sph.pos[2] = 0;
  lbb->c.sph.rad = -1;
  lbb->c.sph.direction = -1;

  while (argc > 0) {
    if(!strncmp(argv[0], "center", strlen(argv[0]))) {
      if(argc < 4) {
	Tcl_AppendResult(interp, "lb_boundary sphere center <x> <y> <z> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if(Tcl_GetDouble(interp, argv[1], &(lbb->c.sph.pos[0])) == TCL_ERROR ||
	Tcl_GetDouble(interp, argv[2], &(lbb->c.sph.pos[1])) == TCL_ERROR ||
	Tcl_GetDouble(interp, argv[3], &(lbb->c.sph.pos[2])) == TCL_ERROR)
	return (TCL_ERROR);
        argc -= 4; argv += 4;
    }
    else if(!strncmp(argv[0], "radius", strlen(argv[0]))) {
      if(argc < 1) {
	Tcl_AppendResult(interp, "lb_boundary sphere radius <r> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if(Tcl_GetDouble(interp, argv[1], &(lbb->c.sph.rad)) == TCL_ERROR)
	return (TCL_ERROR);
        argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "direction", strlen(argv[0]))) {
      if(argc < 1) {
	Tcl_AppendResult(interp, "-1/1 or inside/outside is expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if(!strncmp(argv[1], "inside", strlen(argv[1])))
	lbb->c.sph.direction = -1;
      else if(!strncmp(argv[1], "outside", strlen(argv[1])))
	lbb->c.sph.direction = 1;
      else if(Tcl_GetDouble(interp, argv[1], &(lbb->c.sph.direction)) == TCL_ERROR)
	return (TCL_ERROR);
        argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "type", strlen(argv[0]))) {
      if(argc < 1) {
	Tcl_AppendResult(interp, "lb_boundary sphere type <t> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      argc -= 2; argv += 2;
    }
    else
    break;
  }

  if (lbb->c.sph.rad < 0.) {
    Tcl_AppendResult(interp, "usage: lb_boundary sphere center <x> <y> <z> radius <d> direction <direction> type <t>",
		     (char *) NULL);
    return (TCL_ERROR);    
  }

  return (TCL_OK);
}

int lb_boundary_cylinder(LB_boundary_gpu *lbb, Tcl_Interp *interp,
		    int argc, char **argv)
{
  double axis_len;
  int i;

  lbb->type = LB_BOUNDARY_CYL;
  /* invalid entries to start of */
  lbb->c.cyl.pos[0] = 
  lbb->c.cyl.pos[1] = 
  lbb->c.cyl.pos[2] = 0;
  lbb->c.cyl.axis[0] = 
  lbb->c.cyl.axis[1] = 
  lbb->c.cyl.axis[2] = 0;
  lbb->c.cyl.rad = 0;
  lbb->c.cyl.length = 0;
  lbb->c.cyl.direction = 0;
  while (argc > 0) {
    if(!strncmp(argv[0], "center", strlen(argv[0]))) {
      if(argc < 4) {
	Tcl_AppendResult(interp, "lb_boundary cylinder center <x> <y> <z> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if(Tcl_GetDouble(interp, argv[1], &(lbb->c.cyl.pos[0])) == TCL_ERROR ||
	Tcl_GetDouble(interp, argv[2], &(lbb->c.cyl.pos[1])) == TCL_ERROR ||
	Tcl_GetDouble(interp, argv[3], &(lbb->c.cyl.pos[2])) == TCL_ERROR)
	return (TCL_ERROR);
        argc -= 4; argv += 4;
    }
    else if(!strncmp(argv[0], "axis", strlen(argv[0]))) {
      if(argc < 4) {
	Tcl_AppendResult(interp, "lb_boundary cylinder axis <rx> <ry> <rz> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if(Tcl_GetDouble(interp, argv[1], &(lbb->c.cyl.axis[0])) == TCL_ERROR ||
	Tcl_GetDouble(interp, argv[2], &(lbb->c.cyl.axis[1])) == TCL_ERROR ||
	Tcl_GetDouble(interp, argv[3], &(lbb->c.cyl.axis[2])) == TCL_ERROR)
	return (TCL_ERROR);
        argc -= 4; argv += 4;    
    }
    else if(!strncmp(argv[0], "radius", strlen(argv[0]))) {
      if(argc < 1) {
	Tcl_AppendResult(interp, "lb_boundary cylinder radius <rad> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if(Tcl_GetDouble(interp, argv[1], &(lbb->c.cyl.rad)) == TCL_ERROR)
	return (TCL_ERROR);
        argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "length", strlen(argv[0]))) {
      if(argc < 1) {
	Tcl_AppendResult(interp, "lb_boundary cylinder length <len> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if(Tcl_GetDouble(interp, argv[1], &(lbb->c.cyl.length)) == TCL_ERROR)
	return (TCL_ERROR);
        argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "direction", strlen(argv[0]))) {
      if(argc < 1) {
	Tcl_AppendResult(interp, "lb_boundary cylinder direction <dir> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if(!strncmp(argv[1], "inside", strlen(argv[1])))
	lbb->c.cyl.direction = -1;
      else if (!strncmp(argv[1], "outside", strlen(argv[1])))
	lbb->c.cyl.direction = 1;
      else if (Tcl_GetDouble(interp, argv[1], &(lbb->c.cyl.direction)) == TCL_ERROR)
	return (TCL_ERROR);
        argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "type", strlen(argv[0]))) {
      if(argc < 1) {
	Tcl_AppendResult(interp, "lb_boundary cylinder type <t> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      argc -= 2; argv += 2;
    }
    else
    break;
  }

  axis_len=0.;
  for(i=0;i<3;i++)
    axis_len += SQR(lbb->c.cyl.axis[i]);

  if(lbb->c.cyl.rad < 0. || axis_len < 1e-30 ||
    lbb->c.cyl.direction == 0 || lbb->c.cyl.length <= 0) {
    Tcl_AppendResult(interp, "usage: lb_boundary cylinder center <x> <y> <z> axis <rx> <ry> <rz> radius <rad> length <length> direction <direction> type <t>",
		     (char *) NULL);
    return (TCL_ERROR);    
  }

  /*normalize the axis vector */
  axis_len = sqrt (axis_len);
  for(i=0;i<3;i++) {
    lbb->c.cyl.axis[i] /= axis_len;
  } 
  return (TCL_OK);
}

#endif /* LB_BOUNDARIES_GPU */
int tclcommand_lbboundary_gpu(Tcl_Interp *interp,
	       int argc, char **argv)
{
#ifdef LB_BOUNDARIES_GPU
  int status = TCL_ERROR;
  int c_num;

  if (argc < 2) return lb_boundary_print_all(interp);
  
  if(!strncmp(argv[1], "wall", strlen(argv[1]))) {
    status = lb_boundary_wall(generate_lb_boundary(),interp, argc - 2, argv + 2);
    mpi_bcast_lbboundary(-3);   
  }
  else if(!strncmp(argv[1], "sphere", strlen(argv[1]))) {
    status = lb_boundary_sphere(generate_lb_boundary(),interp, argc - 2, argv + 2); 
    mpi_bcast_lbboundary(-3);
  }
  else if(!strncmp(argv[1], "cylinder", strlen(argv[1]))) {
    status = lb_boundary_cylinder(generate_lb_boundary(),interp, argc - 2, argv + 2);
    mpi_bcast_lbboundary(-3);  
  }
  else if(!strncmp(argv[1], "delete", strlen(argv[1]))) {
    if(argc < 3) { 
      n_lb_boundaries_gpu = 0;     
      status = TCL_OK;
    }
    else {
      //TODO if(Tcl_GetInt(interp, argv[2], &(c_num)) == TCL_ERROR) return (TCL_ERROR);
      // if(c_num < 0 || c_num >= n_lb_boundaries_gpu) {
      Tcl_AppendResult(interp, "Cannot delete individual lb boundaries",(char *) NULL);
      status = TCL_ERROR;    
    }
    mpi_bcast_lbboundary(-3);
  }
  else if (argc == 2 && Tcl_GetInt(interp, argv[1], &c_num) == TCL_OK) {
    printLbBoundaryToResult(interp, c_num);
    status = TCL_OK;
  }
  else {
    Tcl_AppendResult(interp, "possible lb_boundaries: wall sphere cylinder lb_boundary delete {c} to delete lb_boundary or lb_boundaries",(char *) NULL);
    return status;
  }
 
  return status;
#else /* !defined(LB_BOUNDARIES) */
  Tcl_AppendResult(interp, "LB_BOUNDARIES_GPU not compiled in!" ,(char *) NULL);
  return (TCL_ERROR);
#endif /* LB_BOUNDARIES_GPU */
}

#ifdef LB_BOUNDARIES_GPU
/** Initialize boundary conditions for all constraints in the system. */
void lb_init_boundaries_gpu() {
  int n, x, y, z;
  char *errtxt;
  double pos[3], dist, dist_tmp, dist_vec[3];
  int number_of_boundnodes = 0;
  int *host_boundindex = (int*)malloc(sizeof(int));
  size_t size_of_index;
  dist_tmp = 0.;

  if (lbpar_gpu.agrid <= 0) return;

  for(z=0; z<lbpar_gpu.dim_z; z++) {
    for(y=0; y<lbpar_gpu.dim_y; y++) {
      for (x=0; x<lbpar_gpu.dim_x; x++) {	    
	       pos[0] = x*lbpar_gpu.agrid;
	       pos[1] = y*lbpar_gpu.agrid;
	       pos[2] = z*lbpar_gpu.agrid;
	     
	       dist = 0.;

        for (n=0;n<n_lb_boundaries_gpu;n++) {
          switch (lb_boundaries_gpu[n].type) {
            case LB_BOUNDARY_WAL:
              calculate_wall_dist((Particle*) NULL, pos, (Particle*) NULL, &lb_boundaries_gpu[n].c.wal, &dist_tmp, dist_vec);
              break;
            case LB_BOUNDARY_SPH:
              calculate_sphere_dist((Particle*) NULL, pos, (Particle*) NULL, &lb_boundaries_gpu[n].c.sph, &dist_tmp, dist_vec);
              break;
            case LB_BOUNDARY_CYL:
              calculate_cylinder_dist((Particle*) NULL, pos, (Particle*) NULL, &lb_boundaries_gpu[n].c.cyl, &dist_tmp, dist_vec);
              break;
            default:
              errtxt = runtime_error(128);
              ERROR_SPRINTF(errtxt, "{109 lb_boundary type %d not implemented in lb_init_boundaries()\n", lb_boundaries_gpu[n].type);
          }
          
          if (dist > dist_tmp || n == 0) {
            dist = dist_tmp;
          }
        }       
        
  	if (dist <= 0 && n_lb_boundaries_gpu > 0) {
   	  size_of_index = (number_of_boundnodes+1)*sizeof(int);
          host_boundindex = realloc(host_boundindex, size_of_index);
          host_boundindex[number_of_boundnodes] = x + lbpar_gpu.dim_x*y + lbpar_gpu.dim_x*lbpar_gpu.dim_y*z; 
          number_of_boundnodes++;   
        }
      }
    }
  }

  /**call of cuda fkt*/
  lb_init_boundaries_GPU(number_of_boundnodes, host_boundindex);

  LB_TRACE (fprintf(stderr,"lb_init_boundaries \n"));
  LB_TRACE (fprintf(stderr,"boundnumbers %i %i \n", n_lb_boundaries_gpu, number_of_boundnodes));

  free(host_boundindex);
}

#endif /* LB_BOUNDARIES_GPU */
