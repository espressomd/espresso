/*
  Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
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
/** \file lb-boundaries_tcl.cpp
 *
 * Boundary conditions parser file for Lattice Boltzmann fluid dynamics.
 *
 */
#include "utils.hpp"
#include "parser.hpp"
#include "constraint.hpp"
#include "lb.hpp"
#include "interaction_data.hpp"
#include "lb-boundaries.hpp"
#include "communication.hpp"
#include <limits>

#if defined(LB_BOUNDARIES) || defined(LB_BOUNDARIES_GPU)

// TCL Parser functions
int tclcommand_lbboundary(ClientData _data, Tcl_Interp *interp, int argc, char **argv);
int tclcommand_lbboundary_wall(LB_Boundary *lbb, Tcl_Interp *interp, int argc, char **argv);
int tclcommand_lbboundary_sphere(LB_Boundary *lbb, Tcl_Interp *interp, int argc, char **argv);
int tclcommand_lbboundary_cylinder(LB_Boundary *lbb, Tcl_Interp *interp, int argc, char **argv);
int tclcommand_lbboundary_pore(LB_Boundary *lbb, Tcl_Interp *interp, int argc, char **argv);
int tclcommand_lbboundary_stomatocyte(LB_Boundary *lbb, Tcl_Interp *interp, int argc, char **argv);
int tclcommand_lbboundary_hollow_cone(LB_Boundary *lbb, Tcl_Interp *interp, int argc, char **argv);
int tclcommand_printLbBoundaryToResult(Tcl_Interp *interp, int i);

void LB_Boundary::Write_LB_Boundary_Tcl (Tcl_Interp *interp){
	char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
    sprintf(buffer, "%d", this->type);
    Tcl_AppendResult(interp, " type ", buffer, (char *) NULL);
    sprintf(buffer, "%E", this->slip_pref);
    Tcl_AppendResult(interp, " slip_pref ", buffer, (char *) NULL);
    for (int i=0; i<3; i++) {
			sprintf(buffer, "%E", this->force[i]);
    }
    Tcl_AppendResult(interp, " Force on Boundary ", buffer, (char *) NULL);
    for (int i=0; i<3; i++) {
			sprintf(buffer, "%E", this->velocity[i]);
    }
    Tcl_AppendResult(interp, " Surface Velocity ", buffer, (char *) NULL);
    this->_shape->Write_Shape_Tcl(interp);
}

int tclcommand_printLbBoundaryToResult(Tcl_Interp *interp, int i)
{
  LB_Boundary *lbb = &LB_Boundary::lb_boundaries[i];
  char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
  sprintf(buffer, "%d ", i);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  
  lbb->Write_LB_Boundary_Tcl(interp);

  return (TCL_OK);
}

int tclcommand_lbboundary_print_all(Tcl_Interp *interp)
{
  int i;
  
  if(LB_Boundary::n_lb_boundaries>0)
  	Tcl_AppendResult(interp, "{", (char *)NULL);
  	
  for (i = 0; i < LB_Boundary::n_lb_boundaries; i++) {
    if(i>0) Tcl_AppendResult(interp, " {", (char *)NULL);
    	tclcommand_printLbBoundaryToResult(interp, i);
    	
    Tcl_AppendResult(interp, "}", (char *)NULL);
  }
  
  return (TCL_OK);
}

LB_Boundary *generate_lbboundary()
{
  LB_Boundary::n_lb_boundaries++;

  LB_Boundary::lb_boundaries = (LB_Boundary*) realloc(LB_Boundary::lb_boundaries,LB_Boundary::n_lb_boundaries*sizeof(LB_Boundary));

  LB_Boundary::lb_boundaries[LB_Boundary::n_lb_boundaries-1].type = LB_BOUNDARY_BOUNCE_BACK;
  
  LB_Boundary::lb_boundaries[LB_Boundary::n_lb_boundaries-1].velocity[0]=
  LB_Boundary::lb_boundaries[LB_Boundary::n_lb_boundaries-1].velocity[1]=
  LB_Boundary::lb_boundaries[LB_Boundary::n_lb_boundaries-1].velocity[2]=0;
  
  LB_Boundary::lb_boundaries[LB_Boundary::n_lb_boundaries-1].force[0]=
  LB_Boundary::lb_boundaries[LB_Boundary::n_lb_boundaries-1].force[1]=
  LB_Boundary::lb_boundaries[LB_Boundary::n_lb_boundaries-1].force[2]=0;
  
#ifdef EK_BOUNDARIES
  if (ek_initialized)
  {
    LB_Boundary::lb_boundaries[LB_Boundary::n_lb_boundaries-1].charge_density = 0.0;
  }  
#endif
  
  return &LB_Boundary::lb_boundaries[LB_Boundary::n_lb_boundaries-1];
}

int tclcommand_lbboundary_wall(LB_Boundary *lbb, Tcl_Interp *interp, int argc, char **argv)
{
  int i;
  double norm;

  lbb->type = LB_BOUNDARY_WAL;
  /* create a new wall constraint */
  Constraint_wall* new_wall_constraint = new Constraint_wall;
  lbb->_shape = new_wall_constraint;

  
  /* invalid entries to start of */
  new_wall_constraint->n[0] =
  new_wall_constraint->n[1] =
  new_wall_constraint->n[2] = 0;
  
  new_wall_constraint->d = 0;
  
  while(argc > 0) {
    if(ARG_IS_S(0, "normal")) {
      if(argc < 4) {
				Tcl_AppendResult(interp, "lbboundary wall normal <nx> <ny> <nz> expected", (char *) NULL);
			  return (TCL_ERROR);
      }
    
      if(Tcl_GetDouble(interp, argv[1], &(new_wall_constraint->n[0])) == TCL_ERROR ||
	    	 Tcl_GetDouble(interp, argv[2], &(new_wall_constraint->n[1])) == TCL_ERROR ||
	    	 Tcl_GetDouble(interp, argv[3], &(new_wall_constraint->n[2])) == TCL_ERROR)
	      return (TCL_ERROR);
	      
      argc -= 4; argv += 4;
    }
    else if(ARG_IS_S(0, "distance")) {
      if (argc < 1) {
        Tcl_AppendResult(interp, "lbboundary wall dist <d> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      if (Tcl_GetDouble(interp, argv[1], &(new_wall_constraint->d)) == TCL_ERROR)
	      return (TCL_ERROR);
	      
      argc -= 2; argv += 2;
    }
    else if(ARG_IS_S(0, "type")) {
      if (argc < 1) {
	      Tcl_AppendResult(interp, "lbboundary wall type <t> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      argc -= 2; argv += 2;
    }
    else if(ARG_IS_S(0, "velocity")) {
      if(argc < 4) {
	      Tcl_AppendResult(interp, "lbboundary wall velocity <vx> <vy> <vz> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      if(Tcl_GetDouble(interp, argv[1], &(lbb->velocity[0])) == TCL_ERROR ||
      	 Tcl_GetDouble(interp, argv[2], &(lbb->velocity[1])) == TCL_ERROR ||
	       Tcl_GetDouble(interp, argv[3], &(lbb->velocity[2])) == TCL_ERROR)
	      return (TCL_ERROR);

      if (lattice_switch & LATTICE_LB_GPU) {	
#ifdef LB_GPU
        /* No velocity rescaling is required */
#endif
      } else {	
#ifdef LB
        lbb->velocity[0]*=lbpar.tau/lbpar.agrid;
        lbb->velocity[1]*=lbpar.tau/lbpar.agrid;
        lbb->velocity[2]*=lbpar.tau/lbpar.agrid;
#endif
			}
      
      argc -= 4; argv += 4;
    }
    else
      break;
  }
  
  /* length of the normal vector */
  norm = SQR(new_wall_constraint->n[0])+SQR(new_wall_constraint->n[1])+SQR(new_wall_constraint->n[2]);
  
  if (norm < 1e-10) {
    Tcl_AppendResult(interp, "usage: lbboundary wall normal <nx> <ny> <nz> dist <d> type <t>", (char *) NULL);
    return (TCL_ERROR);    
  }
  
  /* normalize the normal vector */
  for (i=0;i<3;i++) new_wall_constraint->n[i] /= sqrt(norm);

  return (TCL_OK);
}

int tclcommand_lbboundary_sphere(LB_Boundary *lbb, Tcl_Interp *interp, int argc, char **argv)
{
  lbb->type = LB_BOUNDARY_SPH;
  /* create a new sphere constraint */
  Constraint_sphere* new_sphere_constraint = new Constraint_sphere;
  lbb->_shape = new_sphere_constraint;

  /* invalid entries to start of */
  new_sphere_constraint->pos[0] =
  new_sphere_constraint->pos[1] =
  new_sphere_constraint->pos[2] = 0;
  new_sphere_constraint->rad = 0;
  new_sphere_constraint->direction = -1;

  while (argc > 0) {
    if(ARG_IS_S(0, "center")) {
      if(argc < 4) {
	      Tcl_AppendResult(interp, "lbboundary sphere center <x> <y> <z> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      if(Tcl_GetDouble(interp, argv[1], &(new_sphere_constraint->pos[0])) == TCL_ERROR ||
	       Tcl_GetDouble(interp, argv[2], &(new_sphere_constraint->pos[1])) == TCL_ERROR ||
	       Tcl_GetDouble(interp, argv[3], &(new_sphere_constraint->pos[2])) == TCL_ERROR)
	      return (TCL_ERROR);
	      
        argc -= 4; argv += 4;
    }
    else if(ARG_IS_S(0, "radius")) {
      if(argc < 1) {
	      Tcl_AppendResult(interp, "lbboundary sphere radius <r> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      if (Tcl_GetDouble(interp, argv[1], &(new_sphere_constraint->rad)) == TCL_ERROR)
	      return (TCL_ERROR);
	      
      argc -= 2; argv += 2;
    }
    else if(ARG_IS_S(0, "direction")) {
      if (argc < 1) {
	      Tcl_AppendResult(interp, "-1/1 or inside/outside is expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      if(ARG_IS_S(1, "inside"))
	      new_sphere_constraint->direction = -1;
      else if(ARG_IS_S(1, "outside"))
	      new_sphere_constraint->direction = 1;
      else if(Tcl_GetDouble(interp, argv[1], &(new_sphere_constraint->direction)) == TCL_ERROR)
	      return (TCL_ERROR);
	      
      argc -= 2; argv += 2;
    }
    else if(ARG_IS_S(0, "type")) {
      if (argc < 1) {
	      Tcl_AppendResult(interp, "lbboundary sphere type <t> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      argc -= 2; argv += 2;
    }
    else
      break;
  }

  if(new_sphere_constraint->rad < 0.) {
    Tcl_AppendResult(interp, "usage: lbboundary sphere center <x> <y> <z> radius <d> direction <direction> type <t>", (char *) NULL);
    return (TCL_ERROR);    
  }

  return (TCL_OK);
}

int tclcommand_lbboundary_cylinder(LB_Boundary *lbb, Tcl_Interp *interp, int argc, char **argv)
{
  double axis_len;
  int i;

  lbb->type = LB_BOUNDARY_CYL;
  /* create a new cylinder constraint */
  Constraint_cylinder* new_cylinder_constraint = new Constraint_cylinder;
  lbb->_shape = new_cylinder_constraint;

  /* invalid entries to start of */
  new_cylinder_constraint->pos[0] =
  new_cylinder_constraint->pos[1] =
  new_cylinder_constraint->pos[2] = 0;
  
  new_cylinder_constraint->axis[0] =
  new_cylinder_constraint->axis[1] =
  new_cylinder_constraint->axis[2] = 0;
  
  new_cylinder_constraint->rad = 0;
  
  new_cylinder_constraint->length = 0;
  
  new_cylinder_constraint->direction = 0;
  
  while (argc > 0) {
    if(ARG_IS_S(0, "center")) {
      if(argc < 4) {
	      Tcl_AppendResult(interp, "lbboundary cylinder center <x> <y> <z> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      if(Tcl_GetDouble(interp, argv[1], &(new_cylinder_constraint->pos[0])) == TCL_ERROR ||
	       Tcl_GetDouble(interp, argv[2], &(new_cylinder_constraint->pos[1])) == TCL_ERROR ||
	       Tcl_GetDouble(interp, argv[3], &(new_cylinder_constraint->pos[2])) == TCL_ERROR)
	      return (TCL_ERROR);
	      
      argc -= 4; argv += 4;
    }
    else if(ARG_IS_S(0, "axis")) {
      if(argc < 4) {
	      Tcl_AppendResult(interp, "lbboundary cylinder axis <rx> <ry> <rz> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      if(Tcl_GetDouble(interp, argv[1], &(new_cylinder_constraint->axis[0])) == TCL_ERROR ||
  	     Tcl_GetDouble(interp, argv[2], &(new_cylinder_constraint->axis[1])) == TCL_ERROR ||
    	   Tcl_GetDouble(interp, argv[3], &(new_cylinder_constraint->axis[2])) == TCL_ERROR)
	      return (TCL_ERROR);

      argc -= 4; argv += 4;    
    }
    else if(ARG_IS_S(0, "radius")) {
      if(argc < 1) {
	      Tcl_AppendResult(interp, "lbboundary cylinder radius <rad> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      if(Tcl_GetDouble(interp, argv[1], &(new_cylinder_constraint->rad)) == TCL_ERROR)
	      return (TCL_ERROR);
	      
      argc -= 2; argv += 2;
    }
    else if(ARG_IS_S(0, "length")) {
      if(argc < 1) {
	      Tcl_AppendResult(interp, "lbboundary cylinder length <len> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      if(Tcl_GetDouble(interp, argv[1], &(new_cylinder_constraint->length)) == TCL_ERROR)
	      return (TCL_ERROR);
	      
      argc -= 2; argv += 2;
    }
    else if(ARG_IS_S(0, "direction")) {
      if(argc < 1) {
	      Tcl_AppendResult(interp, "lbboundary cylinder direction <dir> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      if (ARG_IS_S(1, "inside"))
	      new_cylinder_constraint->direction = -1;
      else if (ARG_IS_S(1, "outside"))
	      new_cylinder_constraint->direction = 1;
      else if (Tcl_GetDouble(interp, argv[1], &(new_cylinder_constraint->direction)) == TCL_ERROR)
	      return (TCL_ERROR);
	      
      argc -= 2; argv += 2;
    }
    else if(ARG_IS_S(0, "type")) {
      if (argc < 1) {
	      Tcl_AppendResult(interp, "lbboundary cylinder type <t> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      argc -= 2; argv += 2;
    }
    else if(ARG_IS_S(0, "velocity")) {
      if(argc < 4) {
	      Tcl_AppendResult(interp, "lbboundary cylinder velocity <vx> <vy> <vz> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      if(Tcl_GetDouble(interp, argv[1], &(lbb->velocity[0])) == TCL_ERROR ||
      	 Tcl_GetDouble(interp, argv[2], &(lbb->velocity[1])) == TCL_ERROR ||
	       Tcl_GetDouble(interp, argv[3], &(lbb->velocity[2])) == TCL_ERROR)
	      return (TCL_ERROR);

      if (lattice_switch & LATTICE_LB_GPU) {	
#ifdef LB_GPU
        /* No velocity rescaling is required */
#endif
      } else {	
#ifdef LB
        lbb->velocity[0]*=lbpar.tau/lbpar.agrid;
        lbb->velocity[1]*=lbpar.tau/lbpar.agrid;
        lbb->velocity[2]*=lbpar.tau/lbpar.agrid;
#endif
			}
      
      argc -= 4; argv += 4;
    }
    else
      break;
  }

  axis_len=0.;
  
  for (i=0;i<3;i++)
    axis_len += SQR(new_cylinder_constraint->axis[i]);

  if(new_cylinder_constraint->rad < 0. || axis_len < 1e-30 ||
     new_cylinder_constraint->direction == 0 || new_cylinder_constraint->length <= 0) {
    Tcl_AppendResult(interp, "usage: lbboundary cylinder center <x> <y> <z> axis <rx> <ry> <rz> radius <rad> length <length> direction <direction> type <t>", (char *) NULL);
    return (TCL_ERROR);    
  }

  /*normalize the axis vector */
  axis_len = sqrt (axis_len);
  
  for (i=0;i<3;i++) {
    new_cylinder_constraint->axis[i] /= axis_len;
  }
      
  return (TCL_OK);
}

int tclcommand_lbboundary_spherocylinder(LB_Boundary *lbb, Tcl_Interp *interp, int argc, char **argv)
{
  double axis_len;
  int i;

  lbb->type = LB_BOUNDARY_SPHEROCYLINDER;
  /* create a new wall constraint */
  Constraint_spherocylinder* new_spherocylinder_constraint = new Constraint_spherocylinder;
  lbb->_shape = new_spherocylinder_constraint;

  
  /* invalid entries to start of */
  new_spherocylinder_constraint->pos[0] =
  new_spherocylinder_constraint->pos[1] =
  new_spherocylinder_constraint->pos[2] = 0;
  
  new_spherocylinder_constraint->axis[0] =
  new_spherocylinder_constraint->axis[1] =
  new_spherocylinder_constraint->axis[2] = 0;
  
  new_spherocylinder_constraint->rad = 0;
  
  new_spherocylinder_constraint->length = 0;
  
  new_spherocylinder_constraint->direction = 0;
  
  while (argc > 0) {
    if(ARG_IS_S(0, "center")) {
      if(argc < 4) {
	      Tcl_AppendResult(interp, "lbboundary spherocylinder center <x> <y> <z> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      if(Tcl_GetDouble(interp, argv[1], &(new_spherocylinder_constraint->pos[0])) == TCL_ERROR ||
	       Tcl_GetDouble(interp, argv[2], &(new_spherocylinder_constraint->pos[1])) == TCL_ERROR ||
	       Tcl_GetDouble(interp, argv[3], &(new_spherocylinder_constraint->pos[2])) == TCL_ERROR)
	      return (TCL_ERROR);
	      
      argc -= 4; argv += 4;
    }
    else if(ARG_IS_S(0, "axis")) {
      if(argc < 4) {
	      Tcl_AppendResult(interp, "lbboundary spherocylinder axis <rx> <ry> <rz> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      if(Tcl_GetDouble(interp, argv[1], &(new_spherocylinder_constraint->axis[0])) == TCL_ERROR ||
  	     Tcl_GetDouble(interp, argv[2], &(new_spherocylinder_constraint->axis[1])) == TCL_ERROR ||
    	   Tcl_GetDouble(interp, argv[3], &(new_spherocylinder_constraint->axis[2])) == TCL_ERROR)
	      return (TCL_ERROR);

      argc -= 4; argv += 4;    
    }
    else if(ARG_IS_S(0, "radius")) {
      if(argc < 1) {
	      Tcl_AppendResult(interp, "lbboundary spherocylinder radius <rad> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      if(Tcl_GetDouble(interp, argv[1], &(new_spherocylinder_constraint->rad)) == TCL_ERROR)
	      return (TCL_ERROR);
	      
      argc -= 2; argv += 2;
    }
    else if(ARG_IS_S(0, "length")) {
      if(argc < 1) {
	      Tcl_AppendResult(interp, "lbboundary spherocylinder length <len> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      if(Tcl_GetDouble(interp, argv[1], &(new_spherocylinder_constraint->length)) == TCL_ERROR)
	      return (TCL_ERROR);
	      
      argc -= 2; argv += 2;
    }
    else if(ARG_IS_S(0, "direction")) {
      if(argc < 1) {
	      Tcl_AppendResult(interp, "lbboundary spherocylinder direction <dir> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      if (ARG_IS_S(1, "inside"))
	      new_spherocylinder_constraint->direction = -1;
      else if (ARG_IS_S(1, "outside"))
	      new_spherocylinder_constraint->direction = 1;
      else if (Tcl_GetDouble(interp, argv[1], &(new_spherocylinder_constraint->direction)) == TCL_ERROR)
	      return (TCL_ERROR);
	      
      argc -= 2; argv += 2;
    }
    else if(ARG_IS_S(0, "type")) {
      if (argc < 1) {
	      Tcl_AppendResult(interp, "lbboundary spherocylinder type <t> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      argc -= 2; argv += 2;
    }
    else if(ARG_IS_S(0, "velocity")) {
      if(argc < 4) {
	      Tcl_AppendResult(interp, "lbboundary spherocylinder velocity <vx> <vy> <vz> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      if(Tcl_GetDouble(interp, argv[1], &(lbb->velocity[0])) == TCL_ERROR ||
      	 Tcl_GetDouble(interp, argv[2], &(lbb->velocity[1])) == TCL_ERROR ||
	       Tcl_GetDouble(interp, argv[3], &(lbb->velocity[2])) == TCL_ERROR)
	      return (TCL_ERROR);

      if (lattice_switch & LATTICE_LB_GPU) {	
#ifdef LB_GPU
        /* No velocity rescaling is required */
#endif
      } else {	
#ifdef LB
        lbb->velocity[0]*=lbpar.tau/lbpar.agrid;
        lbb->velocity[1]*=lbpar.tau/lbpar.agrid;
        lbb->velocity[2]*=lbpar.tau/lbpar.agrid;
#endif
			}
      
      argc -= 4; argv += 4;
    }
    else
      break;
  }

  axis_len=0.;
  
  for (i=0;i<3;i++)
    axis_len += SQR(new_spherocylinder_constraint->axis[i]);

  if(new_spherocylinder_constraint->rad < 0. || axis_len < 1e-30 ||
     new_spherocylinder_constraint->direction == 0 || new_spherocylinder_constraint->length <= 0) {
    Tcl_AppendResult(interp, "usage: lbboundary spherocylinder center <x> <y> <z> axis <rx> <ry> <rz> radius <rad> length <length> direction <direction> type <t>", (char *) NULL);
    return (TCL_ERROR);    
  }

  /*normalize the axis vector */
  axis_len = sqrt (axis_len);
  
  for (i=0;i<3;i++) {
    new_spherocylinder_constraint->axis[i] /= axis_len;
  }
      
  return (TCL_OK);
}

int tclcommand_lbboundary_rhomboid(LB_Boundary *lbb, Tcl_Interp *interp, int argc, char **argv)
{
  double triple_product;
  double tmp[3];
  
  lbb->type = LB_BOUNDARY_RHOMBOID;
  /* create a new rhomboid constraint */
  Constraint_rhomboid* new_rhomboid_constraint = new Constraint_rhomboid;
  lbb->_shape = new_rhomboid_constraint;

  
  new_rhomboid_constraint->pos[0] =
  new_rhomboid_constraint->pos[1] =
  new_rhomboid_constraint->pos[2] = 0;
  
  new_rhomboid_constraint->a[0] =
  new_rhomboid_constraint->a[1] =
  new_rhomboid_constraint->a[2] = 0;
  
  new_rhomboid_constraint->b[0] =
  new_rhomboid_constraint->b[1] =
  new_rhomboid_constraint->b[2] = 0;
  
  new_rhomboid_constraint->c[0] =
  new_rhomboid_constraint->c[1] =
  new_rhomboid_constraint->c[2] = 0;
  
  new_rhomboid_constraint->direction = 0;
  
  while (argc > 0) {
    if(ARG_IS_S(0, "a")) {
      if(argc < 4) {
				Tcl_AppendResult(interp, "lbboundary rhomboid a <ax> <ay> <az> expected", (char *) NULL);
				return TCL_ERROR;
      }
      
      if(Tcl_GetDouble(interp, argv[1], &(new_rhomboid_constraint->a[0])) == TCL_ERROR ||
	 			 Tcl_GetDouble(interp, argv[2], &(new_rhomboid_constraint->a[1])) == TCL_ERROR ||
	  		 Tcl_GetDouble(interp, argv[3], &(new_rhomboid_constraint->a[2])) == TCL_ERROR)
				return TCL_ERROR;
				
			argc -= 4; argv += 4;    
    }
    else if(ARG_IS_S(0, "b")) {
      if(argc < 4) {
				Tcl_AppendResult(interp, "lbboundary rhomboid b <bx> <by> <bz> expected", (char *) NULL);
				return TCL_ERROR;
      }
      
      if(Tcl_GetDouble(interp, argv[1], &(new_rhomboid_constraint->b[0])) == TCL_ERROR ||
	 			 Tcl_GetDouble(interp, argv[2], &(new_rhomboid_constraint->b[1])) == TCL_ERROR ||
	  		 Tcl_GetDouble(interp, argv[3], &(new_rhomboid_constraint->b[2])) == TCL_ERROR)
				return TCL_ERROR;
				
			argc -= 4; argv += 4;    
    }
    else if(ARG_IS_S(0, "c")) {
      if(argc < 4) {
				Tcl_AppendResult(interp, "lbboundary rhomboid c <cx> <cy> <cz> expected", (char *) NULL);
				return TCL_ERROR;
      }
      
      if(Tcl_GetDouble(interp, argv[1], &(new_rhomboid_constraint->c[0])) == TCL_ERROR ||
	 			 Tcl_GetDouble(interp, argv[2], &(new_rhomboid_constraint->c[1])) == TCL_ERROR ||
	  		 Tcl_GetDouble(interp, argv[3], &(new_rhomboid_constraint->c[2])) == TCL_ERROR)
				return TCL_ERROR;
				
			argc -= 4; argv += 4;    
    }
    else if(ARG_IS_S(0, "corner")) { //this has to come after c
      if(argc < 4) {
				Tcl_AppendResult(interp, "lbboundary rhomboid corner <x> <y> <z> expected", (char *) NULL);
				return TCL_ERROR;
      }
      
      if(Tcl_GetDouble(interp, argv[1], &(new_rhomboid_constraint->pos[0])) == TCL_ERROR ||
				 Tcl_GetDouble(interp, argv[2], &(new_rhomboid_constraint->pos[1])) == TCL_ERROR ||
	  		 Tcl_GetDouble(interp, argv[3], &(new_rhomboid_constraint->pos[2])) == TCL_ERROR)
				return TCL_ERROR;
				
      argc -= 4; argv += 4;
    }
    else if(ARG_IS_S(0, "velocity")) {
        if(argc < 4) {
            Tcl_AppendResult(interp, "lbboundary rhomboid velocity <vx> <vy> <vz> expected", (char *) NULL);
            return (TCL_ERROR);
        }
        
        if(Tcl_GetDouble(interp, argv[1], &(lbb->velocity[0])) == TCL_ERROR ||
           Tcl_GetDouble(interp, argv[2], &(lbb->velocity[1])) == TCL_ERROR ||
	       Tcl_GetDouble(interp, argv[3], &(lbb->velocity[2])) == TCL_ERROR)
            return (TCL_ERROR);
        
        if (lattice_switch & LATTICE_LB_GPU) {	
#ifdef LB_GPU
            /* No velocity rescaling is required */
#endif
        } else {	
#ifdef LB
            lbb->velocity[0]*=lbpar.tau/lbpar.agrid;
            lbb->velocity[1]*=lbpar.tau/lbpar.agrid;
            lbb->velocity[2]*=lbpar.tau/lbpar.agrid;
#endif
        }
        
        argc -= 4; argv += 4;
    }
    else if(ARG_IS_S(0, "direction")) {
      if (argc < 2) {
				Tcl_AppendResult(interp, "lbboundary rhomboid direction {inside|outside} expected", (char *) NULL);
				return (TCL_ERROR);
      }
      
      if(ARG_IS_S(1, "inside"))
				new_rhomboid_constraint->direction = -1;
      else if(ARG_IS_S(1, "outside"))
				new_rhomboid_constraint->direction = 1;
      else if(Tcl_GetDouble(interp, argv[1], &(new_rhomboid_constraint->direction)) == TCL_ERROR)
				return TCL_ERROR;
				
      argc -= 2; argv += 2;
    }
    else {
			Tcl_AppendResult(interp, "Error: Unknown parameter ", argv[0], " in lbboundary rhomboid", (char *) NULL);
			return TCL_ERROR;
    }
  }

  if( (new_rhomboid_constraint->a[0] == 0. && new_rhomboid_constraint->a[1] == 0. && new_rhomboid_constraint->a[2] == 0.) ||
  		(new_rhomboid_constraint->b[0] == 0. && new_rhomboid_constraint->b[1] == 0. && new_rhomboid_constraint->b[2] == 0.) ||
  		(new_rhomboid_constraint->c[0] == 0. && new_rhomboid_constraint->c[1] == 0. && new_rhomboid_constraint->c[2] == 0.) ||
  		new_rhomboid_constraint->direction == 0) {
    Tcl_AppendResult(interp, "usage: lbboundary rhomboid corner <x> <y> <z> a <ax> <ay> <az> b <bx> <by> <bz> c <cx> <cy> <cz> direction {inside|outside} type <t> [penetrable <0|1>] [reflecting <1|2>]", (char *) NULL);
    return TCL_ERROR;    
  }
                     
  //If the trihedron a, b, c is left handed, then inside and outside will be exchanged since all normals will be reversed. This compensates  for that, so that the user doesn't have to take care of the handedness.
  triple_product = new_rhomboid_constraint->a[0]*( new_rhomboid_constraint->b[1]*new_rhomboid_constraint->c[2] - new_rhomboid_constraint->b[2]*new_rhomboid_constraint->c[1] ) +
                new_rhomboid_constraint->a[1]*( new_rhomboid_constraint->b[2]*new_rhomboid_constraint->c[0] - new_rhomboid_constraint->b[0]*new_rhomboid_constraint->c[2] ) +
                new_rhomboid_constraint->a[2]*( new_rhomboid_constraint->b[0]*new_rhomboid_constraint->c[1] - new_rhomboid_constraint->b[1]*new_rhomboid_constraint->c[0] );
                
  if(triple_product < 0.)
  {    
    tmp[0] = new_rhomboid_constraint->a[0];
    tmp[1] = new_rhomboid_constraint->a[1];
    tmp[2] = new_rhomboid_constraint->a[2];
    
    new_rhomboid_constraint->a[0] = new_rhomboid_constraint->b[0];
    new_rhomboid_constraint->a[1] = new_rhomboid_constraint->b[1];
    new_rhomboid_constraint->a[2] = new_rhomboid_constraint->b[2];
    
    new_rhomboid_constraint->b[0] = tmp[0];
    new_rhomboid_constraint->b[1] = tmp[1];
    new_rhomboid_constraint->b[2] = tmp[2];
  }

  return TCL_OK;
}

int tclcommand_lbboundary_pore(LB_Boundary *lbb, Tcl_Interp *interp, int argc, char **argv)
{
  double axis_len;
  int i;

  lbb->type = LB_BOUNDARY_POR;
  /* create a new pore constraint */
  Constraint_pore* new_pore_constraint = new Constraint_pore;
  lbb->_shape = new_pore_constraint;

  /* invalid entries to start of */
  new_pore_constraint->pos[0] =
  new_pore_constraint->pos[1] =
  new_pore_constraint->pos[2] = 0;
  
  new_pore_constraint->axis[0] =
  new_pore_constraint->axis[1] =
  new_pore_constraint->axis[2] = 0;
  
  new_pore_constraint->rad_left = 0;
  
  new_pore_constraint->rad_right = 0;
  
  new_pore_constraint->length = 0;
  
  new_pore_constraint->reflecting = 0;
  
  new_pore_constraint->smoothing_radius = 1.;
  
  new_pore_constraint->outer_rad_left = std::numeric_limits<double>::max();
  new_pore_constraint->outer_rad_right = std::numeric_limits<double>::max();
  
  while (argc > 0) {
    if(ARG_IS_S(0, "center")) {
      if(argc < 4) {
	      Tcl_AppendResult(interp, "lbboundary pore center <x> <y> <z> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      if(Tcl_GetDouble(interp, argv[1], &(new_pore_constraint->pos[0])) == TCL_ERROR ||
	       Tcl_GetDouble(interp, argv[2], &(new_pore_constraint->pos[1])) == TCL_ERROR ||
	       Tcl_GetDouble(interp, argv[3], &(new_pore_constraint->pos[2])) == TCL_ERROR)
	      return (TCL_ERROR);
	      
      argc -= 4; argv += 4;
    }
    else if(ARG_IS_S(0, "axis")) {
      if(argc < 4) {
	      Tcl_AppendResult(interp, "lbboundary pore axis <rx> <ry> <rz> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      if(Tcl_GetDouble(interp, argv[1], &(new_pore_constraint->axis[0])) == TCL_ERROR ||
	       Tcl_GetDouble(interp, argv[2], &(new_pore_constraint->axis[1])) == TCL_ERROR ||
	       Tcl_GetDouble(interp, argv[3], &(new_pore_constraint->axis[2])) == TCL_ERROR)
	      return (TCL_ERROR);

      argc -= 4; argv += 4;    
    }
    else if(ARG_IS_S(0, "radius")) {
      if(argc < 1) {
    	  Tcl_AppendResult(interp, "lbboundary pore radius <rad> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      if(Tcl_GetDouble(interp, argv[1], &(new_pore_constraint->rad_left)) == TCL_ERROR)
	      return (TCL_ERROR);
	      
      new_pore_constraint->rad_right =  new_pore_constraint->rad_left;
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "outer_radius", strlen(argv[0]))) {
      if(argc < 1) {
    	  Tcl_AppendResult(interp, "lbboundary pore outer_radius <rad> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      if(Tcl_GetDouble(interp, argv[1], &(new_pore_constraint->outer_rad_left)) == TCL_ERROR)
	      return (TCL_ERROR);
	      
      new_pore_constraint->outer_rad_right =  new_pore_constraint->outer_rad_left;
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "smoothing_radius", strlen(argv[0]))) {
      if (argc < 1) {
	      Tcl_AppendResult(interp, "lbboundary pore smoothing_radius <smoothing_radius> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      if (Tcl_GetDouble(interp, argv[1], &(new_pore_constraint->smoothing_radius)) == TCL_ERROR)
	      return (TCL_ERROR);
	      
      argc -= 2; argv += 2;
    }
    else if(ARG_IS_S(0, "radii")) {
      if(argc < 1) {
	      Tcl_AppendResult(interp, "lbboundary pore radii <rad_left> <rad_right> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      if (Tcl_GetDouble(interp, argv[1], &(new_pore_constraint->rad_left)) == TCL_ERROR)
	      return (TCL_ERROR);
	      
      if (Tcl_GetDouble(interp, argv[2], &(new_pore_constraint->rad_right)) == TCL_ERROR)
	      return (TCL_ERROR);
	      
      argc -= 3; argv += 3;
    }
    else if(!strncmp(argv[0], "outer_radii", strlen(argv[0]))) {
      if(argc < 1) {
	      Tcl_AppendResult(interp, "lbboundary pore outer_radii <rad_left> <rad_right> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      if (Tcl_GetDouble(interp, argv[1], &(new_pore_constraint->outer_rad_left)) == TCL_ERROR)
	      return (TCL_ERROR);
	      
      if (Tcl_GetDouble(interp, argv[2], &(new_pore_constraint->outer_rad_right)) == TCL_ERROR)
	      return (TCL_ERROR);
	      
      argc -= 3; argv += 3;
    }
    else if(!strncmp(argv[0], "length", strlen(argv[0]))) {
      if (argc < 1) {
	      Tcl_AppendResult(interp, "lbboundary pore length <len/2> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      if (Tcl_GetDouble(interp, argv[1], &(new_pore_constraint->length)) == TCL_ERROR)
	      return (TCL_ERROR);

      argc -= 2; argv += 2;
    }
    else
      break;
  }

  axis_len=0.;
  
  for(i=0;i<3;i++)
    axis_len += SQR(new_pore_constraint->axis[i]);

  if(new_pore_constraint->rad_left < 0. || new_pore_constraint->rad_right < 0. || axis_len < 1e-30 || new_pore_constraint->length <= 0) {
    Tcl_AppendResult(interp, "usage: lbboundary pore center <x> <y> <z> axis <rx> <ry> <rz> radius <rad> length <length/2>", (char *) NULL);
    return (TCL_ERROR);
  }

  /*normalize the axis vector */
  axis_len = sqrt (axis_len);
  
  for (i=0;i<3;i++) {
    new_pore_constraint->axis[i] /= axis_len;
  }
  
  return (TCL_OK);
}

int tclcommand_lbboundary_stomatocyte(LB_Boundary *lbb, Tcl_Interp *interp, int argc, char **argv)
{
  /* DON'T PLAY WITH THIS CONSTRAINT UNLESS
     YOU KNOW WHAT IT IS THAT YOU ARE DOING */

  lbb->type = LB_BOUNDARY_STOMATOCYTE;
  /* create a new stomatocyte constraint */
  Constraint_stomatocyte* new_stomatocyte_constraint = new Constraint_stomatocyte;
  lbb->_shape = new_stomatocyte_constraint;

  /* invalid entries to start of */

  new_stomatocyte_constraint->position_x = -M_PI;
  new_stomatocyte_constraint->position_y = -M_PI;
  new_stomatocyte_constraint->position_z = -M_PI;
  new_stomatocyte_constraint->orientation_x = -M_PI;
  new_stomatocyte_constraint->orientation_y = -M_PI;
  new_stomatocyte_constraint->orientation_z = -M_PI;
  new_stomatocyte_constraint->outer_radius = -1.0;
  new_stomatocyte_constraint->inner_radius = -1.0;
  new_stomatocyte_constraint->layer_width = -1.0;
  new_stomatocyte_constraint->direction = 0;

  /* read the data */

  while ( argc > 0 )
  {
    if ( ARG_IS_S( 0, "center" ) ) 
    {
      if(argc < 4) 
      {
	      Tcl_AppendResult(interp, "lbboundary stomatocyte center <x> <y> <z> expected", (char *) NULL);
	      return (TCL_ERROR);
      }

      if ( !ARG_IS_D( 1, new_stomatocyte_constraint->position_x ) ||
	         !ARG_IS_D( 2, new_stomatocyte_constraint->position_y ) ||
	         !ARG_IS_D( 3, new_stomatocyte_constraint->position_z ) )
      {
	      return (TCL_ERROR);
      }

      argc -= 4; argv += 4;
    }
    else if ( ARG_IS_S( 0, "orientation" ) ) 
    {
      if(argc < 4) 
      {
	      Tcl_AppendResult(interp, "lbboundary stomatocyte orientation <ox> <oy> <oz> expected", (char *) NULL);
	      return (TCL_ERROR);
      }

      if ( !ARG_IS_D( 1, new_stomatocyte_constraint->orientation_x ) ||
	         !ARG_IS_D( 2, new_stomatocyte_constraint->orientation_y ) ||
	         !ARG_IS_D( 3, new_stomatocyte_constraint->orientation_z ) )
      {
	      return (TCL_ERROR);
      }

      argc -= 4; argv += 4;
    }
    else if ( ARG_IS_S( 0, "outer_radius" ) ) 
    {
      if(argc < 2) 
      {
	      Tcl_AppendResult(interp, "lbboundary stomatocyte outer_radius <Ro> expected", (char *) NULL);
	      return (TCL_ERROR);
      }

      if ( !ARG_IS_D(1, new_stomatocyte_constraint->outer_radius ) )
	      return (TCL_ERROR);

      argc -= 2; argv += 2;
    }
    else if ( ARG_IS_S( 0, "inner_radius" ) ) 
    {
      if(argc < 2) 
      {
	      Tcl_AppendResult(interp, "lbboundary stomatocyte inner_radius <Ri> expected", (char *) NULL);
	      return (TCL_ERROR);
      }

      if ( !ARG_IS_D( 1, new_stomatocyte_constraint->inner_radius ) )
	      return (TCL_ERROR);

      argc -= 2; argv += 2;
    }
    else if ( ARG_IS_S( 0, "layer_width" ) ) 
    {
      if(argc < 2) 
      {
	      Tcl_AppendResult(interp, "lbboundary stomatocyte layer_width <w> expected", (char *) NULL);
	      return (TCL_ERROR);
      }

      if ( !ARG_IS_D( 1, new_stomatocyte_constraint->layer_width ) )
	      return (TCL_ERROR);

      argc -= 2; argv += 2;
    }
    else if ( ARG_IS_S( 0, "direction" ) ) 
    {
      if ( argc < 2 ) 
      {
	      Tcl_AppendResult(interp, "lbboundary stomatocyte direction {-1|1} or {inside|outside} is expected", (char *) NULL);
	      return (TCL_ERROR);
      }

      if ( ARG_IS_S( 1, "inside" ) )
	      new_stomatocyte_constraint->direction = -1;
      else if ( ARG_IS_S( 1, "outside" ) )
	      new_stomatocyte_constraint->direction = 1;
      else if ( !ARG_IS_D( 1, new_stomatocyte_constraint->direction ) )
	      return (TCL_ERROR); 
      argc -= 2; argv += 2;
    }
    else
      break;
  }

  if ( new_stomatocyte_constraint->outer_radius < 0.0 ||
       new_stomatocyte_constraint->inner_radius < 0.0 ||
       new_stomatocyte_constraint->layer_width < 0.0 )
  {
    Tcl_AppendResult(interp, "stomatocyte radii and width have to be greater than zero",
		     (char *) NULL);
    return (TCL_ERROR);    
  }

  if ( new_stomatocyte_constraint->outer_radius < new_stomatocyte_constraint->inner_radius ||
       new_stomatocyte_constraint->inner_radius < new_stomatocyte_constraint->layer_width ||
       new_stomatocyte_constraint->outer_radius < new_stomatocyte_constraint->layer_width )
  {
    Tcl_AppendResult(interp, "stomatocyte requires layer_width < inner_radius < outer_radius",
		     (char *) NULL);
    return (TCL_ERROR);    
  }

  return (TCL_OK);
}


int tclcommand_lbboundary_hollow_cone(LB_Boundary *lbb, Tcl_Interp *interp, int argc, char **argv)
{
  /* DON'T PLAY WITH THIS CONSTRAINT UNLESS
     YOU KNOW WHAT IT IS THAT YOU ARE DOING */

  lbb->type = LB_BOUNDARY_HOLLOW_CONE;
  /* create a new hollow cone constraint */
  Constraint_hollow_cone* new_hollow_cone_constraint = new Constraint_hollow_cone;
  lbb->_shape = new_hollow_cone_constraint;

  /* invalid entries to start of */

  new_hollow_cone_constraint->position_x = -M_PI;
  new_hollow_cone_constraint->position_y = -M_PI;
  new_hollow_cone_constraint->position_z = -M_PI;
  new_hollow_cone_constraint->orientation_x = -M_PI;
  new_hollow_cone_constraint->orientation_y = -M_PI;
  new_hollow_cone_constraint->orientation_z = -M_PI;
  new_hollow_cone_constraint->outer_radius = -1.0;
  new_hollow_cone_constraint->inner_radius = -1.0;
  new_hollow_cone_constraint->width = -1.0;
  new_hollow_cone_constraint->opening_angle = -1.0;
  new_hollow_cone_constraint->direction = 0;

  /* read the data */

  while ( argc > 0 )
  {
    if ( ARG_IS_S( 0, "center" ) ) 
    {
      if(argc < 4) 
      {
	      Tcl_AppendResult(interp, "lbboundary hollow_cone center <x> <y> <z> expected", (char *) NULL);
	      return (TCL_ERROR);
      }

      if ( !ARG_IS_D( 1, new_hollow_cone_constraint->position_x ) ||
	         !ARG_IS_D( 2, new_hollow_cone_constraint->position_y ) ||
	         !ARG_IS_D( 3, new_hollow_cone_constraint->position_z ) )
      {
	      return (TCL_ERROR);
      }

      argc -= 4; argv += 4;
    }
    else if ( ARG_IS_S( 0, "orientation" ) ) 
    {
      if(argc < 4) 
      {
	      Tcl_AppendResult(interp, "lbboundary hollow_cone orientation <ox> <oy> <oz> expected", (char *) NULL);
	      return (TCL_ERROR);
      }

      if ( !ARG_IS_D( 1, new_hollow_cone_constraint->orientation_x ) ||
	         !ARG_IS_D( 2, new_hollow_cone_constraint->orientation_y ) ||
	         !ARG_IS_D( 3, new_hollow_cone_constraint->orientation_z ) )
      {
	      return (TCL_ERROR);
      }

      argc -= 4; argv += 4;
    }
    else if ( ARG_IS_S( 0, "outer_radius" ) ) 
    {
      if(argc < 2) 
      {
	      Tcl_AppendResult(interp, "lbboundary hollow_cone outer_radius <Ro> expected", (char *) NULL);
	      return (TCL_ERROR);
      }

      if ( !ARG_IS_D(1, new_hollow_cone_constraint->outer_radius ) )
	      return (TCL_ERROR);

      argc -= 2; argv += 2;
    }
    else if ( ARG_IS_S( 0, "inner_radius" ) ) 
    {
      if(argc < 2) 
      {
	      Tcl_AppendResult(interp, "lbboundary hollow_cone inner_radius <Ri> expected", (char *) NULL);
	      return (TCL_ERROR);
      }

      if ( !ARG_IS_D( 1, new_hollow_cone_constraint->inner_radius ) )
	      return (TCL_ERROR);

      argc -= 2; argv += 2;
    }
    else if ( ARG_IS_S( 0, "width" ) ) 
    {
      if(argc < 2) 
      {
	      Tcl_AppendResult(interp, "lbboundary hollow_cone width <w> expected", (char *) NULL);
	      return (TCL_ERROR);
      }

      if ( !ARG_IS_D( 1, new_hollow_cone_constraint->width ) )
	      return (TCL_ERROR);

      argc -= 2; argv += 2;
    }
    else if ( ARG_IS_S( 0, "opening_angle" ) ) 
    {
      if(argc < 2) 
      {
	      Tcl_AppendResult(interp, "lbboundary hollow_cone opening_angle <alpha> expected", (char *) NULL);
	      return (TCL_ERROR);
      }

      if ( !ARG_IS_D( 1, new_hollow_cone_constraint->opening_angle ) )
	      return (TCL_ERROR);

      argc -= 2; argv += 2;
    }
    else if ( ARG_IS_S( 0, "direction" ) ) 
    {
      if ( argc < 2 ) 
      {
	      Tcl_AppendResult(interp, "lbboundary hollow_cone direction {-1|1} or {inside|outside} is expected", (char *) NULL);
	      return (TCL_ERROR);
      }

      if ( ARG_IS_S( 1, "inside" ) )
	      new_hollow_cone_constraint->direction = -1;
      else if ( ARG_IS_S( 1, "outside" ) )
	      new_hollow_cone_constraint->direction = 1;
      else if ( !ARG_IS_D( 1, new_hollow_cone_constraint->direction ) )
	      return (TCL_ERROR); 
      argc -= 2; argv += 2;
    }
    else
      break;
  }

  if ( new_hollow_cone_constraint->outer_radius < 0.0 ||
       new_hollow_cone_constraint->inner_radius < 0.0 ||
       new_hollow_cone_constraint->width < 0.0 )
  {
    Tcl_AppendResult(interp, "hollow_cone radii and width have to be greater than zero",
		     (char *) NULL);
    return (TCL_ERROR);    
  }

  if ( new_hollow_cone_constraint->opening_angle < 0.0 ||
       new_hollow_cone_constraint->opening_angle > M_PI )
  {
    Tcl_AppendResult(interp, "hollow_cone requires 0.0 <= opening_angle <= Pi",
		     (char *) NULL);
    return (TCL_ERROR);    
  }

  if ( fabs( fmod( new_hollow_cone_constraint->outer_radius , 1.0 ) ) < 1.0e-05 ||
       fabs( fmod( new_hollow_cone_constraint->inner_radius , 1.0 ) ) < 1.0e-05 ||
       fabs( fmod( new_hollow_cone_constraint->width , 1.0 ) ) < 1.0e-05 )
  {
      fprintf( stderr, "Warning: Using (almost) exact integer values for the radii or width.\n");
      fprintf( stderr, "         can lead to numerical problems when the LB grid points coincide\n");
      fprintf( stderr, "         with the lattice, for specific values of the position and\n");
      fprintf( stderr, "         orientation. Consider adding or subtracting a small number\n");
      fprintf( stderr, "         to/from the specified sizes to overcome such problems.\n");
      fflush(stdout);
  }

  return (TCL_OK);
}


int tclcommand_lbboundary_box(LB_Boundary *lbb, Tcl_Interp *interp, int argc, char **argv)
{  
  lbb->type = LB_BOUNDARY_BOX;
  /* create a new boundary box constraint */
  Constraint_box* new_box_constraint = new Constraint_box;
  lbb->_shape = new_box_constraint;

  new_box_constraint->value = 0;

  return (TCL_OK);
}

#endif /* LB_BOUNDARIES or LB_BOUNDARIES_GPU */

int tclcommand_lbboundary(ClientData data, Tcl_Interp *interp, int argc, char **argv)
{
#if defined (LB_BOUNDARIES) || defined (LB_BOUNDARIES_GPU)
  int status = TCL_ERROR, c_num;
  
  if ( lattice_switch == LATTICE_OFF ) {
    fprintf (stderr ,"WARNING: Specifying boundaries before using lbfluid assumes a CPU implementation of the LB.\n");
    fprintf (stderr ,"WARNING: This will lead to unexpected behavior if a GPU LB fluid is later used since the boundaries wont exist.\n");
  }

  if (argc < 2)
    return tclcommand_lbboundary_print_all(interp);
  
  if(ARG_IS_S(1, "wall")) {
    status = tclcommand_lbboundary_wall(generate_lbboundary(),interp, argc - 2, argv + 2);
    if (lattice_switch & LATTICE_LB_GPU) {
        mpi_bcast_lbboundary(-3);
    } else 
        mpi_bcast_lbboundary(-1);
  }
  else if(ARG_IS_S(1, "sphere")) {
    status = tclcommand_lbboundary_sphere(generate_lbboundary(),interp, argc - 2, argv + 2);
    if (lattice_switch & LATTICE_LB_GPU) {
        mpi_bcast_lbboundary(-3);
    } else 
        mpi_bcast_lbboundary(-1);
  }
  else if(ARG_IS_S(1, "cylinder")) {
    status = tclcommand_lbboundary_cylinder(generate_lbboundary(),interp, argc - 2, argv + 2);
    if (lattice_switch & LATTICE_LB_GPU) {
        mpi_bcast_lbboundary(-3);
    } else 
        mpi_bcast_lbboundary(-1);
  }
  else if(ARG_IS_S(1, "spherocylinder")) {
    status = tclcommand_lbboundary_spherocylinder(generate_lbboundary(),interp, argc - 2, argv + 2);
    if (lattice_switch & LATTICE_LB_GPU) {
        mpi_bcast_lbboundary(-3);
    } else 
        mpi_bcast_lbboundary(-1);
  }
  else if(ARG_IS_S(1, "rhomboid")) {
    status = tclcommand_lbboundary_rhomboid(generate_lbboundary(),interp, argc - 2, argv + 2);
    if (lattice_switch & LATTICE_LB_GPU) {
        mpi_bcast_lbboundary(-3);
    } else 
        mpi_bcast_lbboundary(-1);
  }
  else if(ARG_IS_S(1, "pore")) {
    status = tclcommand_lbboundary_pore(generate_lbboundary(),interp, argc - 2, argv + 2);
    if (lattice_switch & LATTICE_LB_GPU) {
        mpi_bcast_lbboundary(-3);
    } else 
        mpi_bcast_lbboundary(-1);
  }
  else if(ARG_IS_S(1, "stomatocyte")) {
    status = tclcommand_lbboundary_stomatocyte(generate_lbboundary(),interp, argc - 2, argv + 2);
    if (lattice_switch & LATTICE_LB_GPU) {
        mpi_bcast_lbboundary(-3);
    } else 
        mpi_bcast_lbboundary(-1);
  }
  else if(ARG_IS_S(1, "hollow_cone")) {
    status = tclcommand_lbboundary_hollow_cone(generate_lbboundary(),interp, argc - 2, argv + 2);
    if (lattice_switch & LATTICE_LB_GPU) {
        mpi_bcast_lbboundary(-3);
    } else 
        mpi_bcast_lbboundary(-1);
  }
  else if(ARG_IS_S(1, "force")) {
    if(argc != 3 || Tcl_GetInt(interp, argv[2], &(c_num)) == TCL_ERROR) {
      Tcl_AppendResult(interp, "Usage: lbboundary force $n",(char *) NULL);
      return (TCL_ERROR);
    }
    
    if(c_num < 0 || c_num >= LB_Boundary::n_lb_boundaries) {
      Tcl_AppendResult(interp, "Error in lbboundary force: The selected boundary does not exist",(char *) NULL);
      return (TCL_ERROR);
    }

#if defined (LB) || defined (LB_GPU)
    char buffer[3*TCL_DOUBLE_SPACE+3];
    double force[3];

    status = lbboundary_get_force(c_num, force);
      
    for (int i = 0; i < 3; i++) {
      Tcl_PrintDouble(interp, force[i], buffer);
      Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
    }
#endif
  }
  else if(ARG_IS_S(1, "delete")) {
    if(argc < 3) {
      /* delete all */
        mpi_bcast_lbboundary(-2);
        status = TCL_OK;
    }
    else {
      if (lattice_switch & LATTICE_LB_GPU) {
        Tcl_AppendResult(interp, "Cannot delete individual lb boundaries",(char *) NULL);
        return(TCL_ERROR);
      } 

      if(Tcl_GetInt(interp, argv[2], &(c_num)) == TCL_ERROR) return (TCL_ERROR);
      if(c_num < 0 || c_num >= LB_Boundary::n_lb_boundaries) {
	Tcl_AppendResult(interp, "Can not delete non existing lbboundary",(char *) NULL);
	return (TCL_ERROR);
      }
      mpi_bcast_lbboundary(c_num);
      status = TCL_OK;    
    }
  }
  else if(argc == 2 && Tcl_GetInt(interp, argv[1], &c_num) == TCL_OK) {
    tclcommand_printLbBoundaryToResult(interp, c_num);
    status = TCL_OK;
  }
  else {
    Tcl_AppendResult(interp, "possible lbboundary parameters: wall, sphere, cylinder, rhomboid, pore, stomatocyte, hollow_cone, delete {c} to delete lbboundary",(char *) NULL);
    return (TCL_ERROR);
  }

  return gather_runtime_errors(interp, status);

#else /* !defined(LB_BOUNDARIES) || !defined(LB_BOUNDARIES_GPU */
  Tcl_AppendResult(interp, "LB_BOUNDARIES/LB_BOUNDARIES_GPU not compiled in!" ,(char *) NULL);
  return (TCL_ERROR);
#endif /* LB_BOUNDARIES or LB_BOUNDARIES_GPU */
}

