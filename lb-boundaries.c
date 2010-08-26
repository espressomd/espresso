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

/** \file lb-boundaries.c
 *
 * Boundary conditions for Lattice Boltzmann fluid dynamics.
 * Header file for \ref lb-boundaries.h.
 *
 */
#include "utils.h"
#include "constraint.h"
#include "lb-boundaries.h"
#include "lb.h"
#include "interaction_data.h"

#ifdef LB_BOUNDARIES

int n_lb_boundaries       = 0;
LB_Boundary *lb_boundaries = NULL;

#ifdef LB_BOUNDARIES

//    LB_Boundary lb_boundary_par = { LB_BOUNDARY_NONE, 0.0 }; //do we need this?

/** Initialize a planar boundary specified by a wall constraint.
 * @param plane The \ref Constraint_wall struct describing the boundary.
 */

int printLbBoundaryToResult(Tcl_Interp *interp, int i)
{
  LB_Boundary *lbb = &lb_boundaries[i];
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
  if(n_lb_boundaries>0) Tcl_AppendResult(interp, "{", (char *)NULL);
  for (i = 0; i < n_lb_boundaries; i++) {
    if(i>0) Tcl_AppendResult(interp, " {", (char *)NULL);
    printLbBoundaryToResult(interp, i);
    Tcl_AppendResult(interp, "}", (char *)NULL);
  }
  return (TCL_OK);
}

/*void printLbBoundaryForceToResult(Tcl_Interp *interp, int con)
{
  double f[3];
  char buffer[TCL_DOUBLE_SPACE];

  mpi_get_lb_boundary_force(lbb, f);

  Tcl_PrintDouble(interp, f[0], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, f[1], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, f[2], buffer);
  Tcl_AppendResult(interp, buffer, (char *) NULL);
}*/

LB_Boundary *generate_lb_boundary()
{
  n_lb_boundaries++;
  lb_boundaries = realloc(lb_boundaries,n_lb_boundaries*sizeof(LB_Boundary));
  lb_boundaries[n_lb_boundaries-1].type = LB_BOUNDARY_NONE;
  
  return &lb_boundaries[n_lb_boundaries-1];
}

int lb_boundary_wall(LB_Boundary *lbb, Tcl_Interp *interp,
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
      if (Tcl_GetDouble(interp, argv[1], &(lbb->c.wal.n[0])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[2], &(lbb->c.wal.n[1])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[3], &(lbb->c.wal.n[2])) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 4; argv += 4;
    }
    else if(!strncmp(argv[0], "dist", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "lb_boundary wall dist <d> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(lbb->c.wal.d)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "type", strlen(argv[0]))) {
      if (argc < 1) {
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

int lb_boundary_sphere(LB_Boundary *lbb, Tcl_Interp *interp,
		      int argc, char **argv)
{
  lbb->type = LB_BOUNDARY_SPH;

  /* invalid entries to start of */
  lbb->c.sph.pos[0] = 
    lbb->c.sph.pos[1] = 
    lbb->c.sph.pos[2] = 0;
  lbb->c.sph.rad = 0;
  lbb->c.sph.direction = -1;

  while (argc > 0) {
    if(!strncmp(argv[0], "center", strlen(argv[0]))) {
      if(argc < 4) {
	Tcl_AppendResult(interp, "lb_boundary sphere center <x> <y> <z> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(lbb->c.sph.pos[0])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[2], &(lbb->c.sph.pos[1])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[3], &(lbb->c.sph.pos[2])) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 4; argv += 4;
    }
    else if(!strncmp(argv[0], "radius", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "lb_boundary sphere radius <r> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(lbb->c.sph.rad)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "direction", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "-1/1 or inside/outside is expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (!strncmp(argv[1], "inside", strlen(argv[1])))
	lbb->c.sph.direction = -1;
      else if (!strncmp(argv[1], "outside", strlen(argv[1])))
	lbb->c.sph.direction = 1;
      else if (Tcl_GetDouble(interp, argv[1], &(lbb->c.sph.direction)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "type", strlen(argv[0]))) {
      if (argc < 1) {
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

int lb_boundary_cylinder(LB_Boundary *lbb, Tcl_Interp *interp,
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
      if (Tcl_GetDouble(interp, argv[1], &(lbb->c.cyl.pos[0])) == TCL_ERROR ||
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
      if (Tcl_GetDouble(interp, argv[1], &(lbb->c.cyl.axis[0])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[2], &(lbb->c.cyl.axis[1])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[3], &(lbb->c.cyl.axis[2])) == TCL_ERROR)
	return (TCL_ERROR);

      argc -= 4; argv += 4;    
    }
    else if(!strncmp(argv[0], "radius", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "lb_boundary cylinder radius <rad> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(lbb->c.cyl.rad)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "length", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "lb_boundary cylinder length <len> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(lbb->c.cyl.length)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "direction", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "lb_boundary cylinder direction <dir> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (!strncmp(argv[1], "inside", strlen(argv[1])))
	lbb->c.cyl.direction = -1;
      else if (!strncmp(argv[1], "outside", strlen(argv[1])))
	lbb->c.cyl.direction = 1;
      else if (Tcl_GetDouble(interp, argv[1], &(lbb->c.cyl.direction)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "type", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "lb_boundary cylinder type <t> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      argc -= 2; argv += 2;
    }
    else
      break;
  }

  axis_len=0.;
  for (i=0;i<3;i++)
    axis_len += SQR(lbb->c.cyl.axis[i]);

  if (lbb->c.cyl.rad < 0. || axis_len < 1e-30 ||
      lbb->c.cyl.direction == 0 || lbb->c.cyl.length <= 0) {
    Tcl_AppendResult(interp, "usage: lb_boundary cylinder center <x> <y> <z> axis <rx> <ry> <rz> radius <rad> length <length> direction <direction> type <t>",
		     (char *) NULL);
    return (TCL_ERROR);    
  }

  /*normalize the axis vector */
  axis_len = sqrt (axis_len);
  for (i=0;i<3;i++) {
    lbb->c.cyl.axis[i] /= axis_len;
  }
      
  return (TCL_OK);
}

#endif /* LB_BOUNDARIES */

int lb_boundary(ClientData _data, Tcl_Interp *interp,
	       int argc, char **argv)
{
#ifdef LB_BOUNDARIES
  int status, c_num;


  if (argc < 2) return lb_boundary_print_all(interp);
  
  if(!strncmp(argv[1], "wall", strlen(argv[1]))) {
    status = lb_boundary_wall(generate_lb_boundary(),interp, argc - 2, argv + 2);
    mpi_bcast_lb_boundary(-1);
  }
  else if(!strncmp(argv[1], "sphere", strlen(argv[1]))) {
    status = lb_boundary_sphere(generate_lb_boundary(),interp, argc - 2, argv + 2);
    mpi_bcast_lb_boundary(-1);
  }
  else if(!strncmp(argv[1], "cylinder", strlen(argv[1]))) {
    status = lb_boundary_cylinder(generate_lb_boundary(),interp, argc - 2, argv + 2);
    mpi_bcast_lb_boundary(-1);
  }
  else if(!strncmp(argv[1], "delete", strlen(argv[1]))) {
    if(argc < 3) {
      /* delete all */
      mpi_bcast_lb_boundary(-2);
      status = TCL_OK;
    }
    else {
      if(Tcl_GetInt(interp, argv[2], &(c_num)) == TCL_ERROR) return (TCL_ERROR);
      if(c_num < 0 || c_num >= n_lb_boundaries) {
	Tcl_AppendResult(interp, "Can not delete non existing lb_boundary",(char *) NULL);
	return (TCL_ERROR);
      }
      mpi_bcast_lb_boundary(c_num);
      status = TCL_OK;    
    }
  }
  else if (argc == 2 && Tcl_GetInt(interp, argv[1], &c_num) == TCL_OK) {
    printLbBoundaryToResult(interp, c_num);
    status = TCL_OK;
  }
  else {
    Tcl_AppendResult(interp, "possible lb_boundaries: wall sphere cylinder lb_boundary delete {c} to delete lb_boundary or lb_boundaries",(char *) NULL);
    return (TCL_ERROR);
  }

  return mpi_gather_runtime_errors(interp, status);

#else /* !defined(LB_BOUNDARIES) */
  Tcl_AppendResult(interp, "LB_BOUNDARIES not compiled in!" ,(char *) NULL);
  return (TCL_ERROR);
#endif /* LB_BOUNDARIES */
}

#ifdef LB_BOUNDARIES

static void lb_init_boundary_wall(Constraint_wall* wall) {
  int x, y, z, index, node_domain_position[3], offset[3];
  double pos[3], dist, dist_vec[3];
	
	map_node_array(this_node, node_domain_position); //constraint contains MD coordinates, therefore the conversion
	
	offset[0] = node_domain_position[0]*lblattice.grid[0];
	offset[1] = node_domain_position[1]*lblattice.grid[1];
	offset[2] = node_domain_position[2]*lblattice.grid[2];
  
  for (z=1; z<=lblattice.grid[2]; z++) {
    for (y=1; y<=lblattice.grid[1]; y++) {
	    for (x=1; x<=lblattice.grid[0]; x++) {
	    
	      pos[0] = (offset[0]+(x-1))*lblattice.agrid;
	      pos[1] = (offset[1]+(y-1))*lblattice.agrid;
	      pos[2] = (offset[2]+(z-1))*lblattice.agrid;
       
        calculate_wall_dist((Particle*) NULL, pos, (Particle*) NULL, wall, &dist, dist_vec);
                
  	    if (dist <= 0) {
   	      lbfields[get_linear_index(x,y,z,lblattice.halo_grid)].boundary = 1;     
 	      }
      }
    }
  }
}

static void lb_init_boundary_sphere(Constraint_sphere* sphere) {
  int x, y, z, index, node_domain_position[3], offset[3];
  double pos[3], dist, dist_vec[3];
	
	map_node_array(this_node, node_domain_position); //constraint contains MD coordinates, therefore the conversion
	
	offset[0] = node_domain_position[0]*lblattice.grid[0];
	offset[1] = node_domain_position[1]*lblattice.grid[1];
	offset[2] = node_domain_position[2]*lblattice.grid[2];
  
  for (z=1; z<=lblattice.grid[2]; z++) {
    for (y=1; y<=lblattice.grid[1]; y++) {
	    for (x=1; x<=lblattice.grid[0]; x++) {	    
	      pos[0] = (offset[0]+(x-1))*lblattice.agrid;
	      pos[1] = (offset[1]+(y-1))*lblattice.agrid;
	      pos[2] = (offset[2]+(z-1))*lblattice.agrid;
       
        calculate_sphere_dist((Particle*) NULL, pos, (Particle*) NULL, sphere, &dist, dist_vec);
        
  	    if (dist <= 0)
   	      lbfields[get_linear_index(x,y,z,lblattice.halo_grid)].boundary = 1;   
      }
    }
  }
}

static void lb_init_boundary_cylinder(Constraint_cylinder* cylinder) {
  int x, y, z, index, node_domain_position[3], offset[3];
  double pos[3], dist, dist_vec[3];
	
	map_node_array(this_node, node_domain_position);
	
	offset[0] = node_domain_position[0]*lblattice.grid[0];
	offset[1] = node_domain_position[1]*lblattice.grid[1];
	offset[2] = node_domain_position[2]*lblattice.grid[2];
  
  for (z=1; z<=lblattice.grid[2]; z++) {
    for (y=1; y<=lblattice.grid[1]; y++) {
	    for (x=1; x<=lblattice.grid[0]; x++) {	    
	      pos[0] = (offset[0]+(x-1))*lblattice.agrid;
	      pos[1] = (offset[1]+(y-1))*lblattice.agrid;
	      pos[2] = (offset[2]+(z-1))*lblattice.agrid;
       
        calculate_cylinder_dist((Particle*) NULL, pos, (Particle*) NULL, cylinder, &dist, dist_vec);
        
  	    if (dist <= 0) {
   	      lbfields[get_linear_index(x,y,z,lblattice.halo_grid)].boundary = 1;   
        }
      }
    }
  }
}

/** Initialize boundary conditions for all constraints in the system. */
void lb_init_boundaries() {
  int n, x, y, z, index, node_domain_position[3], offset[3];
  char *errtxt;
  double pos[3], dist, dist_tmp, dist_vec[3];
	
	map_node_array(this_node, node_domain_position);
	
	offset[0] = node_domain_position[0]*lblattice.grid[0];
	offset[1] = node_domain_position[1]*lblattice.grid[1];
	offset[2] = node_domain_position[2]*lblattice.grid[2];

  for (n=0;n<lblattice.halo_grid_volume;n++) {
    lbfields[n].boundary = 0;
  }
  
  for (z=1; z<=lblattice.grid[2]; z++) {
    for (y=1; y<=lblattice.grid[1]; y++) {
	    for (x=1; x<=lblattice.grid[0]; x++) {	    
	      pos[0] = (offset[0]+(x-1))*lblattice.agrid;
	      pos[1] = (offset[1]+(y-1))*lblattice.agrid;
	      pos[2] = (offset[2]+(z-1))*lblattice.agrid;
	      
	      dist = 0.;

        for (n=0;n<n_lb_boundaries;n++) {
          switch (lb_boundaries[n].type) {
            case LB_BOUNDARY_WAL:
              calculate_wall_dist((Particle*) NULL, pos, (Particle*) NULL, &lb_boundaries[n].c.wal, &dist_tmp, dist_vec);
              break;
            case LB_BOUNDARY_SPH:
              calculate_sphere_dist((Particle*) NULL, pos, (Particle*) NULL, &lb_boundaries[n].c.sph, &dist_tmp, dist_vec);
              break;
            case LB_BOUNDARY_CYL:
              calculate_cylinder_dist((Particle*) NULL, pos, (Particle*) NULL, &lb_boundaries[n].c.cyl, &dist_tmp, dist_vec);
              break;
            default:
              errtxt = runtime_error(128);
              ERROR_SPRINTF(errtxt, "{109 lb_boundary type %d not implemented in lb_init_boundaries()\n", lb_boundaries[n].type);
          }
          
          if (abs(dist) > abs(dist_tmp) || n == 0) {
            dist = dist_tmp;
          }
        }       
        
  	    if (dist <= 0 && n_lb_boundaries > 0) {
   	      lbfields[get_linear_index(x,y,z,lblattice.halo_grid)].boundary = 1;   
        }
      }
    }
  }
}

/** Initialize boundary conditions for all constraints in the system. */
/*void lb_init_boundaries() {
  int n;
  char *errtxt;
  
  //printf("executing lb_init_boundaries on node %d\n", this_node);

  for (n=0;n<lblattice.halo_grid_volume;n++) {
    lbfields[n].boundary = 0;
  }

  for (n=0;n<n_lb_boundaries;n++) {
    switch (lb_boundaries[n].type) {
    case LB_BOUNDARY_WAL:
      lb_init_boundary_wall(&lb_boundaries[n].c.wal);
      break;
    case LB_BOUNDARY_SPH:
      lb_init_boundary_sphere(&lb_boundaries[n].c.sph);
      break;
    case LB_BOUNDARY_CYL:
      lb_init_boundary_cylinder(&lb_boundaries[n].c.cyl);
      break;
    default:
      errtxt = runtime_error(128);
      ERROR_SPRINTF(errtxt, "{109 lb_boundary type %d not implemented in lb_init_boundaries()\n",lb_boundaries[n].type);
    }
  }
}*/


static int lbboundaries_parse_slip_reflection(Tcl_Interp *interp, int argc, char **argv) {
#if 0 //problems with slip_pref (georg, 03.08.10)
  if (argc <1) {
    Tcl_AppendResult(interp, "lbboundaries slip_reflection requires 1 argument", (char *)NULL);
    return TCL_ERROR;
  }
  if (!ARG0_IS_D(lb_boundary_par.slip_pref)) {
    Tcl_AppendResult(interp, "wrong argument for lbboundaries slip_reflection", (char *)NULL);
    return TCL_ERROR;
  }
  
  lb_boundary_par.type = LB_BOUNDARY_SLIP_REFLECTION;

//  mpi_bcast_lb_params(LBPAR_BOUNDARY);
  
  return TCL_OK;
#endif //if 0
}

static int lbboundaries_parse_partial_slip(Tcl_Interp *interp, int argc, char **argv) {
#if 0 //problems with slip_pref (georg, 03.08.10)
  if (argc <1) {
    Tcl_AppendResult(interp, "lbboundaries slip requires 1 argument", (char *)NULL);
    return TCL_ERROR;
  }
  if (!ARG0_IS_D(lb_boundary_par.slip_pref)) {
    Tcl_AppendResult(interp, "wrong argument for lbboundaries slip", (char *)NULL);
    return TCL_ERROR;
  }
  
  lb_boundary_par.type = LB_BOUNDARY_PARTIAL_SLIP;

//  mpi_bcast_lb_params(LBPAR_BOUNDARY);
  
  return TCL_OK;
#endif //if 0
}

#endif /* LB_BOUNDARIES */

/** Parser for the \ref lbfluid command. */
int lbboundaries_cmd(ClientData data, Tcl_Interp *interp, int argc, char **argv) {
#if 0 //problems with slip_pref (georg, 03.08.10)
#ifdef LB_BOUNDARIES
  int err = TCL_ERROR;

  if (argc < 2) {
    Tcl_AppendResult(interp, "too few arguments to \"lbboundaries\"", (char *)NULL);
    err = TCL_ERROR;
  }
  else if (ARG1_IS_S("off")) {
    err = TCL_ERROR;
  }
  else if (ARG1_IS_S("bounce_back")) {
    lb_boundary_par.type = LB_BOUNDARY_BOUNCE_BACK;
    err = TCL_OK;
  }
  else if (ARG1_IS_S("specular_reflections")) {
    lb_boundary_par.type = LB_BOUNDARY_SPECULAR_REFLECTION;
    err = TCL_OK;
  }
  else if (ARG1_IS_S("slip_reflection")) {
    err = lbboundaries_parse_slip_reflection(interp, argc-2, argv+2);
  }
  else if (ARG1_IS_S("partial_slip")) {
    err = lbboundaries_parse_partial_slip(interp, argc-2, argv+2);
  }
  else {
    Tcl_AppendResult(interp, "unkown boundary condition \"", argv[1], (char *)NULL);
    err = TCL_ERROR;
  }

  //mpi_bcast_lb_boundaries();

  return err;
#else /* !defined LB_BOUNDARIES */
  Tcl_AppendResult(interp, "LB_BOUNDARIES not compiled in!", NULL);
  return TCL_ERROR;
#endif /* LB_BOUNDARIES */
#endif //if 0
}
#endif /* LB_BOUNDARIES */

