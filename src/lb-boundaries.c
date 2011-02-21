/*
  Copyright (C) 2010,2011 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
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

// TCL Parser functions
int tclcommand_lbboundary(ClientData _data, Tcl_Interp *interp, int argc, char **argv);
int tclcommand_lbboundary_wall(LB_Boundary *lbb, Tcl_Interp *interp, int argc, char **argv);
int tclcommand_lbboundary_sphere(LB_Boundary *lbb, Tcl_Interp *interp, int argc, char **argv);
int tclcommand_lbboundary_cylinder(LB_Boundary *lbb, Tcl_Interp *interp, int argc, char **argv);
int tclcommand_printLbBoundaryToResult(Tcl_Interp *interp, int i);

/** Initialize a planar boundary specified by a wall constraint.
 * @param plane The \ref Constraint_wall struct describing the boundary.
 */
int tclcommand_printLbBoundaryToResult(Tcl_Interp *interp, int i)
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
    Tcl_AppendResult(interp, "unknown lbboundary type ", buffer, ".", (char *) NULL);
    return (TCL_OK);
  }

  return (TCL_OK);
}

int tclcommand_lbboundary_print_all(Tcl_Interp *interp)
{
  int i;
  if(n_lb_boundaries>0) Tcl_AppendResult(interp, "{", (char *)NULL);
  for (i = 0; i < n_lb_boundaries; i++) {
    if(i>0) Tcl_AppendResult(interp, " {", (char *)NULL);
    tclcommand_printLbBoundaryToResult(interp, i);
    Tcl_AppendResult(interp, "}", (char *)NULL);
  }
  return (TCL_OK);
}

LB_Boundary *generate_lbboundary()
{
  n_lb_boundaries++;
  lb_boundaries = realloc(lb_boundaries,n_lb_boundaries*sizeof(LB_Boundary));
  lb_boundaries[n_lb_boundaries-1].type = LB_BOUNDARY_BOUNCE_BACK;
  
  return &lb_boundaries[n_lb_boundaries-1];
}

int tclcommand_lbboundary_wall(LB_Boundary *lbb, Tcl_Interp *interp, int argc, char **argv)
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
	Tcl_AppendResult(interp, "lbboundary wall normal <nx> <ny> <nz> expected", (char *) NULL);
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
	Tcl_AppendResult(interp, "lbboundary wall dist <d> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(lbb->c.wal.d)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "type", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "lbboundary wall type <t> expected", (char *) NULL);
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
    Tcl_AppendResult(interp, "usage: lbboundary wall normal <nx> <ny> <nz> dist <d> type <t>",
		     (char *) NULL);
    return (TCL_ERROR);    
  }
  /* normalize the normal vector */
  for (i=0;i<3;i++) lbb->c.wal.n[i] /= sqrt(norm);

  return (TCL_OK);
}

int tclcommand_lbboundary_sphere(LB_Boundary *lbb, Tcl_Interp *interp, int argc, char **argv)
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
	Tcl_AppendResult(interp, "lbboundary sphere center <x> <y> <z> expected", (char *) NULL);
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
	Tcl_AppendResult(interp, "lbboundary sphere radius <r> expected", (char *) NULL);
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
	Tcl_AppendResult(interp, "lbboundary sphere type <t> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      argc -= 2; argv += 2;
    }
    else
      break;
  }

  if (lbb->c.sph.rad < 0.) {
    Tcl_AppendResult(interp, "usage: lbboundary sphere center <x> <y> <z> radius <d> direction <direction> type <t>",
		     (char *) NULL);
    return (TCL_ERROR);    
  }

  return (TCL_OK);
}

int tclcommand_lbboundary_cylinder(LB_Boundary *lbb, Tcl_Interp *interp, int argc, char **argv)
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
	Tcl_AppendResult(interp, "lbboundary cylinder center <x> <y> <z> expected", (char *) NULL);
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
	Tcl_AppendResult(interp, "lbboundary cylinder axis <rx> <ry> <rz> expected", (char *) NULL);
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
	Tcl_AppendResult(interp, "lbboundary cylinder radius <rad> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(lbb->c.cyl.rad)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "length", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "lbboundary cylinder length <len> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(lbb->c.cyl.length)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "direction", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "lbboundary cylinder direction <dir> expected", (char *) NULL);
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
	Tcl_AppendResult(interp, "lbboundary cylinder type <t> expected", (char *) NULL);
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
    Tcl_AppendResult(interp, "usage: lbboundary cylinder center <x> <y> <z> axis <rx> <ry> <rz> radius <rad> length <length> direction <direction> type <t>",
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

int tclcommand_lbboundary(ClientData _data, Tcl_Interp *interp, int argc, char **argv)
{
#ifdef LB_BOUNDARIES
  int status, c_num;


  if (argc < 2) return tclcommand_lbboundary_print_all(interp);
  
  if(!strncmp(argv[1], "wall", strlen(argv[1]))) {
    status = tclcommand_lbboundary_wall(generate_lbboundary(),interp, argc - 2, argv + 2);
    mpi_bcast_lbboundary(-1);
  }
  else if(!strncmp(argv[1], "sphere", strlen(argv[1]))) {
    status = tclcommand_lbboundary_sphere(generate_lbboundary(),interp, argc - 2, argv + 2);
    mpi_bcast_lbboundary(-1);
  }
  else if(!strncmp(argv[1], "cylinder", strlen(argv[1]))) {
    status = tclcommand_lbboundary_cylinder(generate_lbboundary(),interp, argc - 2, argv + 2);
    mpi_bcast_lbboundary(-1);
  }
  else if(!strncmp(argv[1], "delete", strlen(argv[1]))) {
    if(argc < 3) {
      /* delete all */
      mpi_bcast_lbboundary(-2);
      status = TCL_OK;
    }
    else {
      if(Tcl_GetInt(interp, argv[2], &(c_num)) == TCL_ERROR) return (TCL_ERROR);
      if(c_num < 0 || c_num >= n_lb_boundaries) {
	Tcl_AppendResult(interp, "Can not delete non existing lbboundary",(char *) NULL);
	return (TCL_ERROR);
      }
      mpi_bcast_lbboundary(c_num);
      status = TCL_OK;    
    }
  }
  else if (argc == 2 && Tcl_GetInt(interp, argv[1], &c_num) == TCL_OK) {
    tclcommand_printLbBoundaryToResult(interp, c_num);
    status = TCL_OK;
  }
  else {
    Tcl_AppendResult(interp, "possible lb_boundaries: wall sphere cylinder lbboundary delete {c} to delete lbboundary or lb_boundaries",(char *) NULL);
    return (TCL_ERROR);
  }

  return mpi_gather_runtime_errors(interp, status);

#else /* !defined(LB_BOUNDARIES) */
  Tcl_AppendResult(interp, "LB_BOUNDARIES not compiled in!" ,(char *) NULL);
  return (TCL_ERROR);
#endif /* LB_BOUNDARIES */
}
#ifdef LB_BOUNDARIES


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
              ERROR_SPRINTF(errtxt, "{109 lbboundary type %d not implemented in lb_init_boundaries()\n", lb_boundaries[n].type);
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

#endif /* LB_BOUNDARIES */

