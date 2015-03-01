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
/** \file constraint.cpp
    Implementation of \ref constraint.hpp "constraint.h", here it's just the parsing stuff.
*/
#include <limits>
#include "constraint.hpp"
#include "communication.hpp"
#include "parser.hpp"
#include "TclOutputHelper.hpp"

#ifdef CONSTRAINTS
void Constraint::Write_Constraint_Tcl (Tcl_Interp *interp){
	char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
    sprintf(buffer, "%d", this->part_rep.p.type);
    Tcl_AppendResult(interp, " type ", buffer, (char *) NULL);
    sprintf(buffer, "%d", this->_shape->penetrable);
    Tcl_AppendResult(interp, " penetrable ", buffer, (char *) NULL);
    sprintf(buffer, "%d", this->_shape->reflecting);
    Tcl_AppendResult(interp, " reflecting flag ", buffer, (char *) NULL);
    this->_shape->Write_Shape_Tcl(interp);
}
void Constraint_wall::Write_Shape_Tcl (Tcl_Interp *interp){
    char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
    Tcl_PrintDouble(interp, this->n[0], buffer);
    Tcl_AppendResult(interp, "wall normal ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->n[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->n[2], buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    Tcl_PrintDouble(interp, this->d, buffer);
    Tcl_AppendResult(interp, " dist ", buffer, (char *) NULL);
    sprintf(buffer, "%d", this->only_positive);
    Tcl_AppendResult(interp, " only_positive ", buffer, (char *) NULL);
    sprintf(buffer, "%d", this->tunable_slip);
    Tcl_AppendResult(interp, " tunable_slip ", buffer, (char *) NULL);
}
void Constraint_sphere::Write_Shape_Tcl (Tcl_Interp *interp){
    char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
    Tcl_PrintDouble(interp, this->pos[0], buffer);
    Tcl_AppendResult(interp, "sphere center ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->pos[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->pos[2], buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    Tcl_PrintDouble(interp, this->rad, buffer);
    Tcl_AppendResult(interp, " radius ", buffer, (char *) NULL);
    Tcl_PrintDouble(interp, this->direction, buffer);
    Tcl_AppendResult(interp, " direction ", buffer, (char *) NULL);
}
void Constraint_cylinder::Write_Shape_Tcl (Tcl_Interp *interp){
    char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
    Tcl_PrintDouble(interp, this->pos[0], buffer);
    Tcl_AppendResult(interp, "cylinder center ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->pos[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->pos[2], buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    Tcl_PrintDouble(interp, this->axis[0], buffer);
    Tcl_AppendResult(interp, " axis ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->axis[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->axis[2], buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    Tcl_PrintDouble(interp, this->rad, buffer);
    Tcl_AppendResult(interp, " radius ", buffer, (char *) NULL);
    Tcl_PrintDouble(interp, this->length, buffer);
    Tcl_AppendResult(interp, " length ", buffer, (char *) NULL);
    Tcl_PrintDouble(interp, this->direction, buffer);
    Tcl_AppendResult(interp, " direction ", buffer, (char *) NULL);
}
void Constraint_rhomboid::Write_Shape_Tcl (Tcl_Interp *interp){
    char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
    Tcl_PrintDouble(interp, this->pos[0], buffer);
    Tcl_AppendResult(interp, "rhomboid corner ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->pos[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->pos[2], buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    Tcl_PrintDouble(interp, this->a[0], buffer);
    Tcl_AppendResult(interp, " a ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->a[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->a[2], buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    Tcl_PrintDouble(interp, this->b[0], buffer);
    Tcl_AppendResult(interp, " b ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->b[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->b[2], buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    Tcl_PrintDouble(interp, this->c[0], buffer);
    Tcl_AppendResult(interp, " c ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->c[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->c[2], buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    Tcl_PrintDouble(interp, this->direction, buffer);
    Tcl_AppendResult(interp, " direction ", buffer, (char *) NULL);
}
void Constraint_rod::Write_Shape_Tcl (Tcl_Interp *interp){
    char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
    Tcl_PrintDouble(interp, this->pos[0], buffer);
    Tcl_AppendResult(interp, "rod center ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->pos[1], buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    Tcl_PrintDouble(interp, this->lambda, buffer);
    Tcl_AppendResult(interp, " lambda ", buffer, (char *) NULL);
}
void Constraint_plate::Write_Shape_Tcl (Tcl_Interp *interp){
    char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
    Tcl_PrintDouble(interp, this->pos, buffer);
    Tcl_AppendResult(interp, "plate height ", buffer, (char *) NULL);
    Tcl_PrintDouble(interp, this->sigma, buffer);
    Tcl_AppendResult(interp, " sigma ", buffer, (char *) NULL);
}
void Constraint_maze::Write_Shape_Tcl (Tcl_Interp *interp){
    char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
    Tcl_PrintDouble(interp, this->nsphere, buffer);
    Tcl_AppendResult(interp, "maze nsphere ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->dim, buffer);
    Tcl_AppendResult(interp, " dim ", buffer, (char *) NULL);
    Tcl_PrintDouble(interp, this->sphrad, buffer);
    Tcl_AppendResult(interp, " sphrad ", buffer, (char *) NULL);
    Tcl_PrintDouble(interp, this->cylrad, buffer);
    Tcl_AppendResult(interp, " cylrad ", buffer, (char *) NULL);
}
void Constraint_pore::Write_Shape_Tcl (Tcl_Interp *interp){
    char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
    Tcl_PrintDouble(interp, this->pos[0], buffer);
    Tcl_AppendResult(interp, "pore center ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->pos[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->pos[2], buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    Tcl_PrintDouble(interp, this->axis[0], buffer);
    Tcl_AppendResult(interp, " axis ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->axis[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->axis[2], buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    Tcl_PrintDouble(interp, this->outer_rad_left, buffer);
    Tcl_AppendResult(interp, " outer radius left", buffer, (char *) NULL);
    Tcl_PrintDouble(interp, this->outer_rad_right, buffer);
    Tcl_AppendResult(interp, " outer radius right", buffer, (char *) NULL);
    Tcl_PrintDouble(interp, this->length, buffer);
    Tcl_AppendResult(interp, " length ", buffer, (char *) NULL);
}
void Constraint_slitpore::Write_Shape_Tcl (Tcl_Interp *interp){
    char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
    Tcl_PrintDouble(interp, this->channel_width, buffer);
    Tcl_AppendResult(interp, " channel width ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->lower_smoothing_radius, buffer);
    Tcl_AppendResult(interp, " lower smoothing radius ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->pore_length, buffer);
    Tcl_AppendResult(interp, " pore length ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->pore_mouth, buffer);
    Tcl_AppendResult(interp, " pore mouth ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->pore_width, buffer);
    Tcl_AppendResult(interp, " pore width ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->upper_smoothing_radius, buffer);
    Tcl_AppendResult(interp, " upper_smoothing_radius ", buffer, " ", (char *) NULL);
}
void Constraint_spherocylinder::Write_Shape_Tcl (Tcl_Interp *interp){
    char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
    Tcl_PrintDouble(interp, this->pos[0], buffer);
    Tcl_AppendResult(interp, "cylinder center ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->pos[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->pos[2], buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    Tcl_PrintDouble(interp, this->axis[0], buffer);
    Tcl_AppendResult(interp, " axis ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->axis[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->axis[2], buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    Tcl_PrintDouble(interp, this->rad, buffer);
    Tcl_AppendResult(interp, " radius ", buffer, (char *) NULL);
    Tcl_PrintDouble(interp, this->length, buffer);
    Tcl_AppendResult(interp, " length ", buffer, (char *) NULL);
    Tcl_PrintDouble(interp, this->direction, buffer);
    Tcl_AppendResult(interp, " direction ", buffer, (char *) NULL);
}
void Constraint_stomatocyte::Write_Shape_Tcl (Tcl_Interp *interp){
    char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
    Tcl_PrintDouble(interp, this->position_x, buffer);
    Tcl_AppendResult(interp, "stomatocyte center ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->position_y, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->position_z, buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    Tcl_PrintDouble(interp, this->orientation_x, buffer);
    Tcl_AppendResult(interp, " orientation ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->orientation_y, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->orientation_z, buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    Tcl_PrintDouble(interp, this->outer_radius, buffer);
    Tcl_AppendResult(interp, " outer radius ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->inner_radius, buffer);
    Tcl_AppendResult(interp, " inner radius ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->layer_width, buffer);
    Tcl_AppendResult(interp, " layer width ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->direction, buffer);
    Tcl_AppendResult(interp, " direction ", buffer, (char *) NULL);
}
void Constraint_hollow_cone::Write_Shape_Tcl (Tcl_Interp *interp){
    char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
    Tcl_PrintDouble(interp, this->position_x, buffer);
    Tcl_AppendResult(interp, "hollow_cone center ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->position_y, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->position_z, buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    Tcl_PrintDouble(interp, this->orientation_x, buffer);
    Tcl_AppendResult(interp, " orientation ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->orientation_y, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->orientation_z, buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    Tcl_PrintDouble(interp, this->outer_radius, buffer);
    Tcl_AppendResult(interp, " outer radius ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->inner_radius, buffer);
    Tcl_AppendResult(interp, " inner radius ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->width, buffer);
    Tcl_AppendResult(interp, " width ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->opening_angle, buffer);
    Tcl_AppendResult(interp, " opening angle ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->direction, buffer);
    Tcl_AppendResult(interp, " direction ", buffer, (char *) NULL);
}
void Constraint_ext_magn_field::Write_Shape_Tcl (Tcl_Interp *interp){
    char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
    Tcl_PrintDouble(interp, this->ext_magn_field[0], buffer);
    Tcl_AppendResult(interp, "ext_magn_field ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->ext_magn_field[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->ext_magn_field[2], buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
}
void Constraint_plane::Write_Shape_Tcl (Tcl_Interp *interp){
    char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
    Tcl_PrintDouble(interp, this->pos[0], buffer);
    Tcl_AppendResult(interp, "plane cell ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->pos[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, this->pos[2], buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
}
static int tclprint_to_result_Constraint(Tcl_Interp *interp, int i)
{
  Constraint *con = &Constraint::constraints[i];
  char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
  sprintf(buffer, "%d ", i);
  Tcl_AppendResult(interp, buffer, (char *)NULL);

  con->Write_Constraint_Tcl(interp);

  return (TCL_OK);
}

int tclcommand_constraint_print(Tcl_Interp *interp)
{
  int i;
  if(Constraint::n_constraints>0) Tcl_AppendResult(interp, "{", (char *)NULL);
  for (i = 0; i < Constraint::n_constraints; i++) {
    if(i>0) Tcl_AppendResult(interp, " {", (char *)NULL);
    tclprint_to_result_Constraint(interp, i);
    Tcl_AppendResult(interp, "}", (char *)NULL);
  }
  return (TCL_OK);
}

static void tclprint_to_result_ConstraintForce(Tcl_Interp *interp, int con)
{
  double f[3];
  char buffer[TCL_DOUBLE_SPACE];

  mpi_get_constraint_force(con, f);

  Tcl_PrintDouble(interp, f[0], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, f[1], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, f[2], buffer);
  Tcl_AppendResult(interp, buffer, (char *) NULL);
}


static void tclprint_to_result_n_constraints(Tcl_Interp *interp)
{
  char buffer[TCL_INTEGER_SPACE];
  sprintf(buffer, "%d", Constraint::n_constraints);
  Tcl_AppendResult(interp, buffer, (char *) NULL);
}

static int tclcommand_constraint_parse_wall(Constraint *con, Tcl_Interp *interp,
		    int argc, char **argv)
{
  int i;
  double norm;
  con->type = CONSTRAINT_WAL;
  /* create a new wall constraint */
  Constraint_wall* new_wall_constraint = new Constraint_wall;

  con->_shape = new_wall_constraint;

  /* invalid entries to start of */
  new_wall_constraint->n[0] =
    new_wall_constraint->n[1] =
    new_wall_constraint->n[2] = 0;
  new_wall_constraint->d = 0;
  con->_shape->penetrable = 0;
  con->_shape->only_positive = 0;
  con->part_rep.p.type = -1;

  while (argc > 0) {
    if(!strncmp(argv[0], "normal", strlen(argv[0]))) {
      if(argc < 4) {
	Tcl_AppendResult(interp, "constraint wall normal <nx> <ny> <nz> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(new_wall_constraint->n[0])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[2], &(new_wall_constraint->n[1])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[3], &(new_wall_constraint->n[2])) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 4; argv += 4;
    }
    else if(!strncmp(argv[0], "dist", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint wall dist <d> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(new_wall_constraint->d)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "type", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint wall type <t> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetInt(interp, argv[1], &(con->part_rep.p.type)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "penetrable", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint penetrable <0/1> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetInt(interp, argv[1], &(con->_shape->penetrable)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "reflecting", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint wall reflecting {0|1} expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetInt(interp, argv[1], &(con->_shape->reflecting)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "only_positive", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint wall only_positive {0|1} expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetInt(interp, argv[1], &(con->_shape->only_positive)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "tunable_slip", strlen(argv[0]))) {
      if (argc < 1) {
  	Tcl_AppendResult(interp, "constraint wall tunable_slip {0|1} expected", (char *) NULL);
   	return (TCL_ERROR);
      }
    if (Tcl_GetInt(interp, argv[1], &(con->_shape->tunable_slip)) == TCL_ERROR)
   	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else
      break;
  }
  /* length of the normal vector */
  norm = SQR(new_wall_constraint->n[0])+SQR(new_wall_constraint->n[1])+SQR(new_wall_constraint->n[2]);
  if (norm < 1e-10 || con->part_rep.p.type < 0) {
    Tcl_AppendResult(interp, "usage: constraint wall normal <nx> <ny> <nz> dist <d> type <t> penetrable <0/1> reflecting <1/2>",
		     (char *) NULL);
    return (TCL_ERROR);    
  }
  /* normalize the normal vector */
  for (i=0;i<3;i++) new_wall_constraint->n[i] /= sqrt(norm);

  make_particle_type_exist(con->part_rep.p.type);

  return (TCL_OK);
}

static int tclcommand_constraint_parse_sphere(Constraint *con, Tcl_Interp *interp,
		      int argc, char **argv)
{
  con->type = CONSTRAINT_SPH;

  /* create a new sphere constraint */
  Constraint_sphere* new_sphere_constraint = new Constraint_sphere;
  con->_shape = new_sphere_constraint;
  new_sphere_constraint->pos[0] =
    new_sphere_constraint->pos[1] =
    new_sphere_constraint->pos[2] = 0;
  new_sphere_constraint->rad = 0;
  new_sphere_constraint->direction = -1;
  con->_shape->penetrable = 0;
  con->_shape->reflecting = 0;
  con->part_rep.p.type = -1;

  while (argc > 0) {
    if(!strncmp(argv[0], "center", strlen(argv[0]))) {
      if(argc < 4) {
	Tcl_AppendResult(interp, "constraint sphere center <x> <y> <z> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(new_sphere_constraint->pos[0])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[2], &(new_sphere_constraint->pos[1])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[3], &(new_sphere_constraint->pos[2])) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 4; argv += 4;
    }
    else if(!strncmp(argv[0], "radius", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint sphere radius <r> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(new_sphere_constraint->rad)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "direction", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "-1/1 or inside/outside is expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (!strncmp(argv[1], "inside", strlen(argv[1])))
	new_sphere_constraint->direction = -1;
      else if (!strncmp(argv[1], "outside", strlen(argv[1])))
	new_sphere_constraint->direction = 1;
      else if (Tcl_GetDouble(interp, argv[1], &(new_sphere_constraint->direction)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "type", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint sphere type <t> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetInt(interp, argv[1], &(con->part_rep.p.type)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "penetrable", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint penetrable <0/1> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetInt(interp, argv[1], &(con->_shape->penetrable)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "reflecting", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint sphere reflecting {0|1} expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetInt(interp, argv[1], &(con->_shape->reflecting)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else
      break;
  }

  if (new_sphere_constraint->rad < 0. || con->part_rep.p.type < 0) {
    Tcl_AppendResult(interp, "usage: constraint sphere center <x> <y> <z> radius <d> direction <direction> type <t> penetrable <0/1> reflecting <1/2>",
		     (char *) NULL);
    return (TCL_ERROR);    
  }

  make_particle_type_exist(con->part_rep.p.type);

  return (TCL_OK);
}

static int tclcommand_constraint_parse_cylinder(Constraint *con, Tcl_Interp *interp,
			int argc, char **argv)
{
  double axis_len;
  int i;

  con->type = CONSTRAINT_CYL;
  /* create a new cylinder constraint */
  Constraint_cylinder* new_cylinder_constraint = new Constraint_cylinder;
  con->_shape = new_cylinder_constraint;
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
  con->_shape->penetrable = 0;
  con->part_rep.p.type = -1;
  con->_shape->reflecting = 0;

  while (argc > 0) {
    if(!strncmp(argv[0], "center", strlen(argv[0]))) {
      if(argc < 4) {
	Tcl_AppendResult(interp, "constraint cylinder center <x> <y> <z> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(new_cylinder_constraint->pos[0])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[2], &(new_cylinder_constraint->pos[1])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[3], &(new_cylinder_constraint->pos[2])) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 4; argv += 4;
    }
    else if(!strncmp(argv[0], "axis", strlen(argv[0]))) {
      if(argc < 4) {
	Tcl_AppendResult(interp, "constraint cylinder axis <rx> <ry> <rz> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(new_cylinder_constraint->axis[0])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[2], &(new_cylinder_constraint->axis[1])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[3], &(new_cylinder_constraint->axis[2])) == TCL_ERROR)
	return (TCL_ERROR);

      argc -= 4; argv += 4;    
    }
    else if(!strncmp(argv[0], "radius", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint cylinder radius <rad> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(new_cylinder_constraint->rad)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "length", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint cylinder length <len> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(new_cylinder_constraint->length)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "direction", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint cylinder direction <dir> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (!strncmp(argv[1], "inside", strlen(argv[1])))
	new_cylinder_constraint->direction = -1;
      else if (!strncmp(argv[1], "outside", strlen(argv[1])))
	new_cylinder_constraint->direction = 1;
      else if (Tcl_GetDouble(interp, argv[1], &(new_cylinder_constraint->direction)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "type", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint cylinder type <t> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetInt(interp, argv[1], &(con->part_rep.p.type)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "penetrable", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint cylinder penetrable <0/1> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetInt(interp, argv[1], &(con->_shape->penetrable)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "reflecting", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint cylinder reflecting {0|1} expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetInt(interp, argv[1], &(con->_shape->reflecting)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else
      break;
  }

  axis_len=0.;
  for (i=0;i<3;i++)
    axis_len += SQR(new_cylinder_constraint->axis[i]);

  if (new_cylinder_constraint->rad < 0. || con->part_rep.p.type < 0 || axis_len < 1e-30 ||
      new_cylinder_constraint->direction == 0 || new_cylinder_constraint->length <= 0) {
    Tcl_AppendResult(interp, "usage: constraint cylinder center <x> <y> <z> axis <rx> <ry> <rz> radius <rad> length <length> direction <direction> type <t> penetrable <0/1> reflecting <1/2>",
		     (char *) NULL);
    return (TCL_ERROR);    
  }

  /*normalize the axis vector */
  axis_len = sqrt (axis_len);
  for (i=0;i<3;i++) {
    new_cylinder_constraint->axis[i] /= axis_len;
  }

  make_particle_type_exist(con->part_rep.p.type);
      
  return (TCL_OK);
}

static int tclcommand_constraint_parse_rhomboid(Constraint *con, Tcl_Interp *interp,
			int argc, char **argv)
{
  double triple_product;
  double tmp[3];
  
  con->type = CONSTRAINT_RHOMBOID;

  /* create a new rhomboid constraint */
  Constraint_rhomboid* new_rhomboid_constraint = new Constraint_rhomboid;
  con->_shape = new_rhomboid_constraint;
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
  
  con->_shape->penetrable = 0;
  
  con->part_rep.p.type = -1;
  
  con->_shape->reflecting = 0;
  
  while (argc > 0) {
    if(!strncmp(argv[0], "a", strlen(argv[0]))) {
      if(argc < 4) {
				Tcl_AppendResult(interp, "constraint rhomboid a <ax> <ay> <az> expected", (char *) NULL);
				return TCL_ERROR;
      }
      
      if(Tcl_GetDouble(interp, argv[1], &(new_rhomboid_constraint->a[0])) == TCL_ERROR ||
	 			 Tcl_GetDouble(interp, argv[2], &(new_rhomboid_constraint->a[1])) == TCL_ERROR ||
	  		 Tcl_GetDouble(interp, argv[3], &(new_rhomboid_constraint->a[2])) == TCL_ERROR)
				return TCL_ERROR;
				
			argc -= 4; argv += 4;    
    }
    else if(!strncmp(argv[0], "b", strlen(argv[0]))) {
      if(argc < 4) {
				Tcl_AppendResult(interp, "constraint rhomboid b <bx> <by> <bz> expected", (char *) NULL);
				return TCL_ERROR;
      }
      
      if(Tcl_GetDouble(interp, argv[1], &(new_rhomboid_constraint->b[0])) == TCL_ERROR ||
	 			 Tcl_GetDouble(interp, argv[2], &(new_rhomboid_constraint->b[1])) == TCL_ERROR ||
	  		 Tcl_GetDouble(interp, argv[3], &(new_rhomboid_constraint->b[2])) == TCL_ERROR)
				return TCL_ERROR;
				
			argc -= 4; argv += 4;    
    }
    else if(!strncmp(argv[0], "c", strlen(argv[0]))) {
      if(argc < 4) {
				Tcl_AppendResult(interp, "constraint rhomboid c <cx> <cy> <cz> expected", (char *) NULL);
				return TCL_ERROR;
      }
      
      if(Tcl_GetDouble(interp, argv[1], &(new_rhomboid_constraint->c[0])) == TCL_ERROR ||
	 			 Tcl_GetDouble(interp, argv[2], &(new_rhomboid_constraint->c[1])) == TCL_ERROR ||
	  		 Tcl_GetDouble(interp, argv[3], &(new_rhomboid_constraint->c[2])) == TCL_ERROR)
				return TCL_ERROR;
				
			argc -= 4; argv += 4;    
    }
    else if(!strncmp(argv[0], "corner", strlen(argv[0]))) { //this has to come after c
      if(argc < 4) {
				Tcl_AppendResult(interp, "constraint rhomboid corner <x> <y> <z> expected", (char *) NULL);
				return TCL_ERROR;
      }
      
      if(Tcl_GetDouble(interp, argv[1], &(new_rhomboid_constraint->pos[0])) == TCL_ERROR ||
				 Tcl_GetDouble(interp, argv[2], &(new_rhomboid_constraint->pos[1])) == TCL_ERROR ||
	  		 Tcl_GetDouble(interp, argv[3], &(new_rhomboid_constraint->pos[2])) == TCL_ERROR)
				return TCL_ERROR;
				
      argc -= 4; argv += 4;
    }
    else if(!strncmp(argv[0], "direction", strlen(argv[0]))) {
      if (argc < 2) {
				Tcl_AppendResult(interp, "constraint rhomboid direction {inside|outside} expected", (char *) NULL);
				return (TCL_ERROR);
      }
      
      if(!strncmp(argv[1], "inside", strlen(argv[1])))
				new_rhomboid_constraint->direction = -1;
      else if(!strncmp(argv[1], "outside", strlen(argv[1])))
				new_rhomboid_constraint->direction = 1;
      else if(Tcl_GetDouble(interp, argv[1], &(new_rhomboid_constraint->direction)) == TCL_ERROR)
				return TCL_ERROR;
				
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "type", strlen(argv[0]))) {
      if(argc < 2) {
				Tcl_AppendResult(interp, "constraint rhomboid type <t> expected", (char *) NULL);
				return TCL_ERROR;
      }
      
      if(Tcl_GetInt(interp, argv[1], &(con->part_rep.p.type)) == TCL_ERROR)
				return TCL_ERROR;
				
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "penetrable", strlen(argv[0]))) {
      if (argc < 2) {
				Tcl_AppendResult(interp, "constraint rhomboid penetrable {0|1} expected", (char *) NULL);
				return TCL_ERROR;
      }
      
      if (Tcl_GetInt(interp, argv[1], &(con->_shape->penetrable)) == TCL_ERROR)
				return TCL_ERROR;
				
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "reflecting", strlen(argv[0]))) {
      if (argc < 2) {
				Tcl_AppendResult(interp, "constraint rhomboid reflecting {0|1} expected", (char *) NULL);
				return TCL_ERROR;
      }
      
      if (Tcl_GetInt(interp, argv[1], &(con->_shape->reflecting)) == TCL_ERROR)
				return TCL_ERROR;
				
      argc -= 2; argv += 2;
    }
    else {
			Tcl_AppendResult(interp, "Error: Unknown parameter ", argv[0], " in constraint rhomboid", (char *) NULL);
			return TCL_ERROR;
    }
  }

  if( (new_rhomboid_constraint->a[0] == 0. && new_rhomboid_constraint->a[1] == 0. && new_rhomboid_constraint->a[2] == 0.) ||
  		(new_rhomboid_constraint->b[0] == 0. && new_rhomboid_constraint->b[1] == 0. && new_rhomboid_constraint->b[2] == 0.) ||
  		(new_rhomboid_constraint->c[0] == 0. && new_rhomboid_constraint->c[1] == 0. && new_rhomboid_constraint->c[2] == 0.) ||
  		con->part_rep.p.type < 0 || new_rhomboid_constraint->direction == 0) {
    Tcl_AppendResult(interp, "usage: constraint rhomboid corner <x> <y> <z> a <ax> <ay> <az> b <bx> <by> <bz> c <cx> <cy> <cz> direction {inside|outside} type <t> [penetrable <0|1>] [reflecting <1|2>]",
		     (char *) NULL);
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

  make_particle_type_exist(con->part_rep.p.type);

  return TCL_OK;
}

static int tclcommand_constraint_parse_pore(Constraint *con, Tcl_Interp *interp,
		    int argc, char **argv)
{
  double axis_len;
  int i;

  con->type = CONSTRAINT_PORE;
  /* create a new pore constraint */
  Constraint_pore* new_pore_constraint = new Constraint_pore;
  con->_shape = new_pore_constraint;

  /* invalid entries to start of */
  new_pore_constraint->pos[0] =
    new_pore_constraint->pos[1] =
    new_pore_constraint->pos[2] = 0;
  new_pore_constraint->axis[0] =
    new_pore_constraint->axis[1] =
    new_pore_constraint->axis[2] = 0;
  new_pore_constraint->rad_left = 0;
  new_pore_constraint->rad_right = 0;
  new_pore_constraint->outer_rad_left = 1e99;
  new_pore_constraint->outer_rad_right = 1e99;
  new_pore_constraint->length = 0;
  con->_shape->reflecting = 0;
  con->part_rep.p.type = -1;
  new_pore_constraint->smoothing_radius = 1.;
  new_pore_constraint->outer_rad_left = std::numeric_limits<double>::max();
  new_pore_constraint->outer_rad_right = std::numeric_limits<double>::max();
  
  while (argc > 0) {
    if(!strncmp(argv[0], "center", strlen(argv[0]))) {
      if(argc < 4) {
	Tcl_AppendResult(interp, "constraint pore center <x> <y> <z> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(new_pore_constraint->pos[0])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[2], &(new_pore_constraint->pos[1])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[3], &(new_pore_constraint->pos[2])) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 4; argv += 4;
    }
    else if(!strncmp(argv[0], "axis", strlen(argv[0]))) {
      if(argc < 4) {
	Tcl_AppendResult(interp, "constraint pore axis <rx> <ry> <rz> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(new_pore_constraint->axis[0])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[2], &(new_pore_constraint->axis[1])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[3], &(new_pore_constraint->axis[2])) == TCL_ERROR)
	return (TCL_ERROR);

      argc -= 4; argv += 4;    
    }
    else if(!strncmp(argv[0], "radius", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint pore radius <rad> expected", (char *) NULL);
	return (TCL_ERROR);
      }  
      if (Tcl_GetDouble(interp, argv[1], &(new_pore_constraint->rad_left)) == TCL_ERROR)
	return (TCL_ERROR);
      new_pore_constraint->rad_right =  new_pore_constraint->rad_left;
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "outer_radius", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint pore outer_radius <rad> expected", (char *) NULL);
	return (TCL_ERROR);
      }  
      if (Tcl_GetDouble(interp, argv[1], &(new_pore_constraint->outer_rad_left)) == TCL_ERROR)
	return (TCL_ERROR);
      new_pore_constraint->outer_rad_right =  new_pore_constraint->outer_rad_left;
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "smoothing_radius", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint pore smoothing_radius <smoothing_radius> expected", (char *) NULL);
	return (TCL_ERROR);
      }  
      if (Tcl_GetDouble(interp, argv[1], &(new_pore_constraint->smoothing_radius)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "radii", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint pore radii <rad_left> <rad_right> expected", (char *) NULL);
	return (TCL_ERROR);
      }  
      if (Tcl_GetDouble(interp, argv[1], &(new_pore_constraint->rad_left)) == TCL_ERROR)
	return (TCL_ERROR);
      if (Tcl_GetDouble(interp, argv[2], &(new_pore_constraint->rad_right)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 3; argv += 3;
    }
    else if(!strncmp(argv[0], "outer_radii", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint pore outer_radii <rad_left> <rad_right> expected", (char *) NULL);
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
	Tcl_AppendResult(interp, "constraint pore length <len/2> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(new_pore_constraint->length)) == TCL_ERROR)
	return (TCL_ERROR);
//      new_pore_constraint->length *= 2;
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "type", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint pore type <t> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetInt(interp, argv[1], &(con->part_rep.p.type)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "reflecting", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint pore reflecting {0|1} expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetInt(interp, argv[1], &(con->_shape->reflecting)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else
      break;
  }

  axis_len=0.;
  for (i=0;i<3;i++)
    axis_len += SQR(new_pore_constraint->axis[i]);

  if (new_pore_constraint->rad_left < 0. || new_pore_constraint->rad_right < 0. || con->part_rep.p.type < 0 || axis_len < 1e-30 ||
      new_pore_constraint->length <= 0) {
    Tcl_AppendResult(interp, "usage: constraint pore center <x> <y> <z> axis <rx> <ry> <rz> radius <rad> length <length/2> type <t>",
		     (char *) NULL);
    return (TCL_ERROR);
  }

  /*normalize the axis vector */
  axis_len = sqrt (axis_len);
  for (i=0;i<3;i++) {
    new_pore_constraint->axis[i] /= axis_len;
  }
  
  make_particle_type_exist(con->part_rep.p.type);

  return (TCL_OK);
}


static int tclcommand_constraint_parse_slitpore(Constraint *con, Tcl_Interp *interp,
		    int argc, char **argv)
{

  con->type = CONSTRAINT_SLITPORE;
  /* create a new rhomboid constraint */
  Constraint_slitpore* new_slitpore_constraint = new Constraint_slitpore;
  con->_shape = new_slitpore_constraint;

  /* invalid entries to start of */
  new_slitpore_constraint->pore_mouth = 0;
  new_slitpore_constraint->channel_width = 0;
  new_slitpore_constraint->pore_width = 0;
  new_slitpore_constraint->pore_length = 0;
  new_slitpore_constraint->upper_smoothing_radius = 0;
  new_slitpore_constraint->lower_smoothing_radius = 0;
  con->_shape->reflecting = 0;
  con->part_rep.p.type = -1;
  while (argc > 0) {
    if(!strncmp(argv[0], "pore_mouth", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint slitpore mouth <mouth> expected", (char *) NULL);
	return (TCL_ERROR);
      }  
      if (Tcl_GetDouble(interp, argv[1], &(new_slitpore_constraint->pore_mouth)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "pore_width", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint slitpore pore_width <pore_width> expected", (char *) NULL);
	return (TCL_ERROR);
      }  
      if (Tcl_GetDouble(interp, argv[1], &(new_slitpore_constraint->pore_width)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "pore_length", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint slitpore pore_width <pore_length> expected", (char *) NULL);
	return (TCL_ERROR);
      }  
      if (Tcl_GetDouble(interp, argv[1], &(new_slitpore_constraint->pore_length)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "channel_width", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint slitpore channel_width <channel_width> expected", (char *) NULL);
	return (TCL_ERROR);
      }  
      if (Tcl_GetDouble(interp, argv[1], &(new_slitpore_constraint->channel_width)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "upper_smoothing_radius", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint slitpore upper_smoothing_radius <r> expected", (char *) NULL);
	return (TCL_ERROR);
      }  
      if (Tcl_GetDouble(interp, argv[1], &(new_slitpore_constraint->upper_smoothing_radius)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "lower_smoothing_radius", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint slitpore lower_smoothing_radius <r> expected", (char *) NULL);
	return (TCL_ERROR);
      }  
      if (Tcl_GetDouble(interp, argv[1], &(new_slitpore_constraint->lower_smoothing_radius)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "type", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint pore type <t> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetInt(interp, argv[1], &(con->part_rep.p.type)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "reflecting", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint pore reflecting {0|1} expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetInt(interp, argv[1], &(con->_shape->reflecting)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else
      break;
  }

  int error = 0;
  if (new_slitpore_constraint->channel_width <= 0.)  {
    Tcl_AppendResult(interp, "Error in contraint slitpore: Channel with must be > 0", (char *) NULL);
    error=1;
  }
  if ( new_slitpore_constraint->pore_width <= 0. ) {
    Tcl_AppendResult(interp, "Error in contraint slitpore: Pore width must be > 0", (char *) NULL);
    error=1;
  }
  if (  new_slitpore_constraint->pore_length < 0. ) {
    Tcl_AppendResult(interp, "Error in contraint slitpore: Pore length must be > 0", (char *) NULL);
    error=1;
  }
  if ( con->part_rep.p.type < 0 ) {
    Tcl_AppendResult(interp, "Error in contraint slitpore: Type not set", (char *) NULL);
    error=1;
  }
 
  if (error)
    return (TCL_ERROR);

  make_particle_type_exist(con->part_rep.p.type);

  return (TCL_OK);
}


static int tclcommand_constraint_parse_rod(Constraint *con, Tcl_Interp *interp,
		   int argc, char **argv)
{
  con->type = CONSTRAINT_ROD;
  con->part_rep.p.type = -1;
  /* create a new rod constraint */
  Constraint_rod* new_rod_constraint = new Constraint_rod;
  con->_shape = new_rod_constraint;

  new_rod_constraint->pos[0] = new_rod_constraint->pos[1] = 0;
  new_rod_constraint->lambda = 0;
  while (argc > 0) {
    if(!strncmp(argv[0], "center", strlen(argv[0]))) {
      if(argc < 3) {
	Tcl_AppendResult(interp, "constraint rod center <px> <py> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(new_rod_constraint->pos[0])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[2], &(new_rod_constraint->pos[1])) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 3; argv += 3;
    }
    else if(!strncmp(argv[0], "radius", strlen(argv[0]))) {
      Tcl_AppendResult(interp, "constraint rod radius <r> is deprecated, please use a cylinder for LJ component", (char *) NULL);
      return (TCL_ERROR);
    }
    else if(!strncmp(argv[0], "type", strlen(argv[0]))) {
      Tcl_AppendResult(interp, "constraint rod type <t> is deprecated, please use a cylinder for LJ component", (char *) NULL);
      return (TCL_ERROR);
    }
    else if(!strncmp(argv[0], "lambda", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint rod lambda <l> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(new_rod_constraint->lambda)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else
      break;
  }
  return (TCL_OK);
}

static int tclcommand_constraint_parse_plate(Constraint *con, Tcl_Interp *interp,
		     int argc, char **argv)
{
  con->type = CONSTRAINT_PLATE;
  con->part_rep.p.type = -1;
  /* create a new rhomboid constraint */
  Constraint_plate* new_plate_constraint = new Constraint_plate;
  con->_shape = new_plate_constraint;

  new_plate_constraint->pos = 0;
  new_plate_constraint->sigma = 0;
  while (argc > 0) {
    if(!strncmp(argv[0], "height", strlen(argv[0]))) {
      if(argc < 2) {
	Tcl_AppendResult(interp, "constraint plate height <pz> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(new_plate_constraint->pos)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "sigma", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint rod sigma <s> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(new_plate_constraint->sigma)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else
      break;
  }
  return (TCL_OK);
}

static int tclcommand_constraint_parse_stomatocyte(Constraint *con, Tcl_Interp *interp,
		      int argc, char **argv)
{
  /* DON'T PLAY WITH THIS CONSTRAINT UNLESS
     YOU KNOW WHAT IT IS THAT YOU ARE DOING */

  con->type = CONSTRAINT_STOMATOCYTE;

  /* create a new stomatocyte constraint */
  Constraint_stomatocyte* new_stomatocyte_constraint = new Constraint_stomatocyte;
  con->_shape = new_stomatocyte_constraint;

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
  con->_shape->penetrable = 0;
  con->_shape->reflecting = 0;
  con->part_rep.p.type = -1;

  /* read the data */

  while ( argc > 0 )
  {
    if ( ARG_IS_S( 0, "center" ) ) 
    {
      if(argc < 4) 
      {
	      Tcl_AppendResult(interp, "constraint stomatocyte center <x> <y> <z> expected", (char *) NULL);
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
	      Tcl_AppendResult(interp, "constraint stomatocyte orientation <ox> <oy> <oz> expected", (char *) NULL);
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
	      Tcl_AppendResult(interp, "constraint stomatocyte outer_radius <Ro> expected", (char *) NULL);
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
	      Tcl_AppendResult(interp, "constraint stomatocyte inner_radius <Ri> expected", (char *) NULL);
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
	      Tcl_AppendResult(interp, "constraint stomatocyte layer_width <w> expected", (char *) NULL);
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
	      Tcl_AppendResult(interp, "constraint stomatocyte direction {-1|1} or {inside|outside} is expected", (char *) NULL);
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
    else if ( ARG_IS_S( 0, "type" ) ) 
    {
      if ( argc < 2 )
      {
	      Tcl_AppendResult(interp, "constraint stomatocyte type <t> expected", (char *) NULL);
	      return (TCL_ERROR);
      }

      if ( !ARG_IS_I( 1, con->part_rep.p.type ) )
	      return (TCL_ERROR);

      argc -= 2; argv += 2;
    }
    else if ( ARG_IS_S( 0, "penetrable" ) ) 
    {
      if ( argc < 2 ) 
      {
	      Tcl_AppendResult(interp, "constraint stomatocyte penetrable {0|1} expected", (char *) NULL);
	      return (TCL_ERROR);
      }

      if ( !ARG_IS_I( 1, con->_shape->penetrable ) )
	      return (TCL_ERROR);

      argc -= 2; argv += 2;
    }
    else if ( ARG_IS_S( 0, "reflecting" ) ) 
    {
      if (argc < 1) 
      {
	      Tcl_AppendResult(interp, "constraint stomatocyte reflecting {0|1} expected", (char *) NULL);
	      return (TCL_ERROR);
      }

      if ( !ARG_IS_I( 1, con->_shape->reflecting ) )
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

  make_particle_type_exist(con->part_rep.p.type);

  return (TCL_OK);
}

static int tclcommand_constraint_parse_hollow_cone(Constraint *con, Tcl_Interp *interp,
		      int argc, char **argv)
{
  con->type = CONSTRAINT_HOLLOW_CONE;
  /* create a new hollow cone constraint */
  Constraint_hollow_cone* new_hollow_cone_constraint = new Constraint_hollow_cone;
  con->_shape = new_hollow_cone_constraint;

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
  con->_shape->penetrable = 0;
  con->_shape->reflecting = 0;
  con->part_rep.p.type = -1;

  /* read the data */

  while ( argc > 0 )
  {
    if ( ARG_IS_S( 0, "center" ) ) 
    {
      if(argc < 4) 
      {
	      Tcl_AppendResult(interp, "constraint hollow_cone center <x> <y> <z> expected", (char *) NULL);
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
	      Tcl_AppendResult(interp, "constraint hollow_cone orientation <ox> <oy> <oz> expected", (char *) NULL);
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
	      Tcl_AppendResult(interp, "constraint hollow_cone outer_radius <Ro> expected", (char *) NULL);
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
	      Tcl_AppendResult(interp, "constraint hollow_cone inner_radius <Ri> expected", (char *) NULL);
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
	      Tcl_AppendResult(interp, "constraint hollow_cone width <w> expected", (char *) NULL);
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
	      Tcl_AppendResult(interp, "constraint hollow_cone opening_angle <alpha> expected", (char *) NULL);
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
	      Tcl_AppendResult(interp, "constraint hollow_cone direction {-1|1} or {inside|outside} is expected", (char *) NULL);
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
    else if ( ARG_IS_S( 0, "type" ) ) 
    {
      if ( argc < 2 )
      {
	      Tcl_AppendResult(interp, "constraint hollow_cone type <t> expected", (char *) NULL);
	      return (TCL_ERROR);
      }

      if ( !ARG_IS_I( 1, con->part_rep.p.type ) )
	      return (TCL_ERROR);

      argc -= 2; argv += 2;
    }
    else if ( ARG_IS_S( 0, "penetrable" ) ) 
    {
      if ( argc < 2 ) 
      {
	      Tcl_AppendResult(interp, "constraint hollow_cone penetrable {0|1} expected", (char *) NULL);
	      return (TCL_ERROR);
      }

      if ( !ARG_IS_I( 1, con->_shape->penetrable ) )
	      return (TCL_ERROR);

      argc -= 2; argv += 2;
    }
    else if ( ARG_IS_S( 0, "reflecting" ) ) 
    {
      if (argc < 1) 
      {
	      Tcl_AppendResult(interp, "constraint hollow_cone reflecting {0|1} expected", (char *) NULL);
	      return (TCL_ERROR);
      }

      if ( !ARG_IS_I( 1, con->_shape->reflecting ) )
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

  make_particle_type_exist(con->part_rep.p.type);

  return (TCL_OK);
}

static int tclcommand_constraint_parse_maze(Constraint *con, Tcl_Interp *interp,
		      int argc, char **argv)
{
  con->type = CONSTRAINT_MAZE;
  /* create a new maze constraint */
  Constraint_maze* new_maze_constraint = new Constraint_maze;
  con->_shape = new_maze_constraint;


  /* invalid entries to start of */
  new_maze_constraint->nsphere = 0.;
  new_maze_constraint->dim = -1.;
  new_maze_constraint->sphrad = 0.;
  new_maze_constraint->cylrad = -1.;
  con->_shape->penetrable = 0;
  con->part_rep.p.type = -1;

  while (argc > 0) {
    if(!strncmp(argv[0], "nsphere", strlen(argv[0]))) {
      if(argc < 4) {
	Tcl_AppendResult(interp, "constraint maze nsphere <n> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(new_maze_constraint->nsphere)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "dim", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint maze dim <d> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(new_maze_constraint->dim)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "sphrad", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint maze sphrad <r> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(new_maze_constraint->sphrad)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "cylrad", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint maze cylrad <r> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(new_maze_constraint->cylrad)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "type", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint maze type <t> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetInt(interp, argv[1], &(con->part_rep.p.type)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "penetrable", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint penetrable <0/1> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetInt(interp, argv[1], &(con->_shape->penetrable)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else
      break;
  }

  if (new_maze_constraint->sphrad < 0. || new_maze_constraint->cylrad < 0. || con->part_rep.p.type < 0 || new_maze_constraint->dim < 0) {
    Tcl_AppendResult(interp, "usage: constraint maze nsphere <n> dim <d> sphrad <r> cylrad <r> type <t> penetrable <0/1>",
		     (char *) NULL);
    return (TCL_ERROR);    
  }

  make_particle_type_exist(con->part_rep.p.type);

  return (TCL_OK);
}

//ER
int tclcommand_constraint_parse_ext_magn_field(Constraint *con, Tcl_Interp *interp,
		      int argc, char **argv)
{
  int i;
  con->type = CONSTRAINT_EXT_MAGN_FIELD;
  con->part_rep.p.type=-1;
  /* create a new rhomboid constraint */
  Constraint_ext_magn_field* new_emfield_constraint = new Constraint_ext_magn_field;
  con->_shape = new_emfield_constraint;


  for(i=0; i<3; i++)
     new_emfield_constraint->ext_magn_field[i] = 0.;

  if(argc < 3) {
      Tcl_AppendResult(interp, "usage: constraint ext_magn_field <x> <y> <z>", (char *) NULL);
      return (TCL_ERROR);
  }
  for(i=0; i<3; i++){
     if (Tcl_GetDouble(interp, argv[i], &(new_emfield_constraint->ext_magn_field[i])) == TCL_ERROR)
	return (TCL_ERROR);
  }
  argc -= 3; argv += 3;

  return (TCL_OK);
}
//end ER

static int tclcommand_constraint_parse_plane_cell(Constraint *con, Tcl_Interp *interp,
                      int argc, char **argv)
{
  Tcl_AppendResult(interp, "constraint plane cell deprecated, use constraint wall instead!", (char *) NULL);
  return (TCL_ERROR);
  con->type = CONSTRAINT_PLANE;

  /* create a new plane constraint */
  Constraint_plane* new_plane_constraint = new Constraint_plane;
  con->_shape = new_plane_constraint;

  /* invalid entries to start of */
  new_plane_constraint->pos[0] =
    new_plane_constraint->pos[1] =
    new_plane_constraint->pos[2] = 0;
  con->part_rep.p.type = -1;

  while (argc > 0) {
    if(!strncmp(argv[0], "cell", strlen(argv[0]))) {
      if(argc < 4) {
        Tcl_AppendResult(interp, "constraint plane cell <x> <y> <z> expected", (char *) NULL);
        return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(new_plane_constraint->pos[0])) == TCL_ERROR ||
          Tcl_GetDouble(interp, argv[2], &(new_plane_constraint->pos[1])) == TCL_ERROR ||
          Tcl_GetDouble(interp, argv[3], &(new_plane_constraint->pos[2])) == TCL_ERROR)
        return (TCL_ERROR);
      argc -= 4; argv += 4;
    }
    else if(!strncmp(argv[0], "type", strlen(argv[0]))) {
      if (argc < 1) {
        Tcl_AppendResult(interp, "constraint plane cell type <t> expected", (char *) NULL);
        return (TCL_ERROR);
      }
      if (Tcl_GetInt(interp, argv[1], &(con->part_rep.p.type)) == TCL_ERROR)
        return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else
      break;
  }

  if (con->part_rep.p.type < 0) {
    Tcl_AppendResult(interp, "usage: constraint plane cell <x> <y> <z> type <t>",
                     (char *) NULL);
    return (TCL_ERROR);    
  }

  make_particle_type_exist(con->part_rep.p.type);

  return (TCL_OK);
}

static int tclcommand_constraint_mindist_position(Tcl_Interp *interp, int argc, char **argv) {
	double pos[3];
	double vec[3];
	double mindist = 1e100;
	int n;
	char buffer[TCL_DOUBLE_SPACE];
	if (Constraint::n_constraints==0) {
		Tcl_AppendResult(interp, "Error in constraint mindist_position: no constraints defined\n", (char*) NULL);
		return TCL_ERROR;
	}

	Particle* p1=0;
	if (ARG_IS_D(0, pos[0]) && ARG_IS_D(1, pos[1]) && ARG_IS_D(2, pos[2])) {
		for(n=0;n<Constraint::n_constraints;n++) {
			double dist=1e100;
			if ( !Constraint::constraints[n]._shape->penetrable )
				Constraint::constraints[n]._shape->calculate_dist(p1, pos, &Constraint::constraints[n].part_rep, &dist, vec);
			mindist = dist<mindist ? dist : mindist;
		}
		Tcl_PrintDouble(interp, mindist, buffer);
		Tcl_AppendResult(interp, " ", buffer, " ", (char *) NULL);
		return TCL_OK;
	} else {
		Tcl_AppendResult(interp, "\nError in constraint mindist_position: could not read position\n", (char*) NULL);
		return TCL_ERROR;
	}
}

int tclcommand_constraint_mindist_position_vec(Tcl_Interp *interp, int argc, char **argv) {
	double pos[3];
	double vec[3];
	double mindist = 1e100;
	double minvec[3] = { 1e100, 1e100, 1e100};
	int n;
	char buffer[TCL_DOUBLE_SPACE];
	if (Constraint::n_constraints==0) {
		Tcl_AppendResult(interp, "Error in constraint mindist_position: no constraints defined\n", (char*) NULL);
		return TCL_ERROR;
	}
	if (argc < 3) {
		Tcl_AppendResult(interp, "Usage: constraint mindist_position_vec <px> <py> <py>\n", (char*) NULL);
		return TCL_ERROR;
	}

	Particle* p1=0;
	if (ARG_IS_D(0, pos[0]) && ARG_IS_D(1, pos[1]) && ARG_IS_D(2, pos[2])) {
		for(n=0;n<Constraint::n_constraints;n++) {
			double dist=1e100;
			if ( !Constraint::constraints[n]._shape->penetrable )
				Constraint::constraints[n]._shape->calculate_dist(p1, pos, &Constraint::constraints[n].part_rep, &dist, vec);

			if (dist<mindist) {
				mindist = dist;
				minvec[0]=vec[0];
				minvec[1]=vec[1];
				minvec[2]=vec[2];
			}

		}
		Tcl_PrintDouble(interp, minvec[0], buffer);
		Tcl_AppendResult(interp, " ", buffer, " ", (char *) NULL);
		Tcl_PrintDouble(interp, minvec[1], buffer);
		Tcl_AppendResult(interp, " ", buffer, " ", (char *) NULL);
		Tcl_PrintDouble(interp, minvec[2], buffer);
		Tcl_AppendResult(interp, " ", buffer, (char *) NULL);
		return TCL_OK;
	} else {
		Tcl_AppendResult(interp, "\nError in constraint mindist_position: could not read position\n", (char*) NULL);
		return TCL_ERROR;
	}
}
#endif

int tclcommand_constraint(ClientData _data, Tcl_Interp *interp,
	       int argc, char **argv)
{
#ifdef CONSTRAINTS
  int status, c_num;

  if (argc < 2) return tclcommand_constraint_print(interp);

  if(!strncmp(argv[1], "mindist_position", strlen(argv[1]))) {
    return tclcommand_constraint_mindist_position(interp, argc-2, argv+2);
  }
  if(!strncmp(argv[1], "mindist_position_vec", strlen(argv[1]))) {
    return tclcommand_constraint_mindist_position_vec(interp, argc-2, argv+2);
  }
  else if(!strncmp(argv[1], "wall", strlen(argv[1]))) {
    status = tclcommand_constraint_parse_wall(generate_constraint(),interp, argc - 2, argv + 2);
    mpi_bcast_constraint(-1);
  }
  else if(!strncmp(argv[1], "sphere", strlen(argv[1]))) {
    status = tclcommand_constraint_parse_sphere(generate_constraint(),interp, argc - 2, argv + 2);
    mpi_bcast_constraint(-1);
  }
  else if(!strncmp(argv[1], "cylinder", strlen(argv[1]))) {
    status = tclcommand_constraint_parse_cylinder(generate_constraint(),interp, argc - 2, argv + 2);
    mpi_bcast_constraint(-1);
  }
  else if(!strncmp(argv[1], "rhomboid", strlen(argv[1]))) {
    status = tclcommand_constraint_parse_rhomboid(generate_constraint(),interp, argc - 2, argv + 2);
    mpi_bcast_constraint(-1);
  }
  else if(!strncmp(argv[1], "rod", strlen(argv[1]))) {
    status = tclcommand_constraint_parse_rod(generate_constraint(),interp, argc - 2, argv + 2);
    mpi_bcast_constraint(-1);
  }
  else if(!strncmp(argv[1], "plate", strlen(argv[1]))) {
    status = tclcommand_constraint_parse_plate(generate_constraint(),interp, argc - 2, argv + 2);
    mpi_bcast_constraint(-1);
  }
  else if(!strncmp(argv[1], "maze", strlen(argv[1]))) {
    status = tclcommand_constraint_parse_maze(generate_constraint(),interp, argc - 2, argv + 2);
    mpi_bcast_constraint(-1);
  }
  else if(!strncmp(argv[1], "pore", strlen(argv[1]))) {
    status = tclcommand_constraint_parse_pore(generate_constraint(),interp, argc - 2, argv + 2);
    mpi_bcast_constraint(-1);
  }
  else if(!strncmp(argv[1], "slitpore", strlen(argv[1]))) {
    status = tclcommand_constraint_parse_slitpore(generate_constraint(),interp, argc - 2, argv + 2);
    mpi_bcast_constraint(-1);
  }
  else if(!strncmp(argv[1], "stomatocyte", strlen(argv[1]))) {
    status = tclcommand_constraint_parse_stomatocyte(generate_constraint(),interp, argc - 2, argv + 2);
    mpi_bcast_constraint(-1);
  }
  else if(!strncmp(argv[1], "hollow_cone", strlen(argv[1]))) {
    status = tclcommand_constraint_parse_hollow_cone(generate_constraint(),interp, argc - 2, argv + 2);
    mpi_bcast_constraint(-1);
  }
  //ER
  else if(!strncmp(argv[1], "ext_magn_field", strlen(argv[1]))) {
    status = tclcommand_constraint_parse_ext_magn_field(generate_constraint(),interp, argc - 2, argv + 2);
    mpi_bcast_constraint(-1);
  }
  else if(!strncmp(argv[1], "plane cell", strlen(argv[1]))) {
    status = tclcommand_constraint_parse_plane_cell(generate_constraint(),interp, argc - 2, argv + 2);
    mpi_bcast_constraint(-1);
  }
  //end ER
  else if(!strncmp(argv[1], "force", strlen(argv[1]))) {
    if(argc < 3) {
      Tcl_AppendResult(interp, "which particles force?",(char *) NULL);
      return (TCL_ERROR);
    }
    if(Tcl_GetInt(interp, argv[2], &(c_num)) == TCL_ERROR) return (TCL_ERROR);
    if(c_num < 0 || c_num >= Constraint::n_constraints) {
      Tcl_AppendResult(interp, "constraint does not exist",(char *) NULL);
      return (TCL_ERROR);
    }
    tclprint_to_result_ConstraintForce(interp, c_num);
    status  = TCL_OK;
  }
  else if(!strncmp(argv[1], "n_constraints", strlen(argv[1]))) {
    tclprint_to_result_n_constraints(interp);
    status=TCL_OK;
  }
  else if(!strncmp(argv[1], "delete", strlen(argv[1]))) {
    if(argc < 3) {
      /* delete all */
      mpi_bcast_constraint(-2);
      status = TCL_OK;
    }
    else {
      if(Tcl_GetInt(interp, argv[2], &(c_num)) == TCL_ERROR) return (TCL_ERROR);
      if(c_num < 0 || c_num >= Constraint::n_constraints) {
	Tcl_AppendResult(interp, "Can not delete non existing constraint",(char *) NULL);
	return (TCL_ERROR);
      }
      mpi_bcast_constraint(c_num);
      status = TCL_OK;    
    }
  }
  else if (argc == 2 && Tcl_GetInt(interp, argv[1], &c_num) == TCL_OK) {
    tclprint_to_result_Constraint(interp, c_num);
    status = TCL_OK;
  }
  else {
  //ER "ext_magn_field" was put in the next line //end ER
    Tcl_AppendResult(interp, "possible constraints: wall sphere cylinder maze pore stomatocyte hollow_cone ext_magn_field or constraint delete {c} to delete constraint(s)",(char *) NULL);
    return (TCL_ERROR);
  }

  return gather_runtime_errors(interp, status);

#else /* !defined(CONSTRAINTS) */
  Tcl_AppendResult(interp, "Constraints not compiled in!" ,(char *) NULL);
  return (TCL_ERROR);
#endif
}

