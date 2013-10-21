/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
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
/** \file constraint.c
    Implementation of \ref constraint.h "constraint.h", here it's just the parsing stuff.
*/
#include "constraint.hpp"
#include "communication.hpp"
#include "parser.hpp"

#ifdef CONSTRAINTS
static int tclprint_to_result_Constraint(Tcl_Interp *interp, int i)
{
  Constraint *con = &constraints[i];
  char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
  sprintf(buffer, "%d ", i);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  
  switch (con->type) {
  case CONSTRAINT_WAL:
    Tcl_PrintDouble(interp, con->c.wal.n[0], buffer);
    Tcl_AppendResult(interp, "wall normal ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.wal.n[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.wal.n[2], buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    Tcl_PrintDouble(interp, con->c.wal.d, buffer);
    Tcl_AppendResult(interp, " dist ", buffer, (char *) NULL);
    sprintf(buffer, "%d", con->part_rep.p.type);
    Tcl_AppendResult(interp, " type ", buffer, (char *) NULL);
    sprintf(buffer, "%d", con->c.wal.penetrable);
    Tcl_AppendResult(interp, " penetrable ", buffer, (char *) NULL);
    break;
  case CONSTRAINT_SPH:
    Tcl_PrintDouble(interp, con->c.sph.pos[0], buffer);
    Tcl_AppendResult(interp, "sphere center ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.sph.pos[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.sph.pos[2], buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    Tcl_PrintDouble(interp, con->c.sph.rad, buffer);
    Tcl_AppendResult(interp, " radius ", buffer, (char *) NULL);
    Tcl_PrintDouble(interp, con->c.sph.direction, buffer);
    Tcl_AppendResult(interp, " direction ", buffer, (char *) NULL);
    sprintf(buffer, "%d", con->part_rep.p.type);
    Tcl_AppendResult(interp, " type ", buffer, (char *) NULL);
    sprintf(buffer, "%d", con->c.sph.penetrable);
    Tcl_AppendResult(interp, " penetrable ", buffer, (char *) NULL);
    break;
  case CONSTRAINT_CYL:
    Tcl_PrintDouble(interp, con->c.cyl.pos[0], buffer);
    Tcl_AppendResult(interp, "cylinder center ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.cyl.pos[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.cyl.pos[2], buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    Tcl_PrintDouble(interp, con->c.cyl.axis[0], buffer);
    Tcl_AppendResult(interp, " axis ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.cyl.axis[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.cyl.axis[2], buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    Tcl_PrintDouble(interp, con->c.cyl.rad, buffer);
    Tcl_AppendResult(interp, " radius ", buffer, (char *) NULL);
    Tcl_PrintDouble(interp, con->c.cyl.length, buffer);
    Tcl_AppendResult(interp, " length ", buffer, (char *) NULL);
    Tcl_PrintDouble(interp, con->c.cyl.direction, buffer);
    Tcl_AppendResult(interp, " direction ", buffer, (char *) NULL);
    sprintf(buffer, "%d", con->part_rep.p.type);
    Tcl_AppendResult(interp, " type ", buffer, (char *) NULL);
    sprintf(buffer, "%d", con->c.cyl.penetrable);
    Tcl_AppendResult(interp, " penetrable ", buffer, (char *) NULL);
    break;
  case CONSTRAINT_RHOMBOID:
    Tcl_PrintDouble(interp, con->c.rhomboid.pos[0], buffer);
    Tcl_AppendResult(interp, "rhomboid corner ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.rhomboid.pos[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.rhomboid.pos[2], buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    Tcl_PrintDouble(interp, con->c.rhomboid.a[0], buffer);
    Tcl_AppendResult(interp, " a ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.rhomboid.a[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.rhomboid.a[2], buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    Tcl_PrintDouble(interp, con->c.rhomboid.b[0], buffer);
    Tcl_AppendResult(interp, " b ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.rhomboid.b[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.rhomboid.b[2], buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    Tcl_PrintDouble(interp, con->c.rhomboid.c[0], buffer);
    Tcl_AppendResult(interp, " c ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.rhomboid.c[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.rhomboid.c[2], buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    Tcl_PrintDouble(interp, con->c.rhomboid.direction, buffer);
    Tcl_AppendResult(interp, " direction ", buffer, (char *) NULL);
    sprintf(buffer, "%d", con->part_rep.p.type);
    Tcl_AppendResult(interp, " type ", buffer, (char *) NULL);
    sprintf(buffer, "%d", con->c.rhomboid.penetrable);
    Tcl_AppendResult(interp, " penetrable ", buffer, (char *) NULL);    
    break;
  case CONSTRAINT_ROD:
    Tcl_PrintDouble(interp, con->c.rod.pos[0], buffer);
    Tcl_AppendResult(interp, "rod center ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.rod.pos[1], buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    Tcl_PrintDouble(interp, con->c.rod.lambda, buffer);
    Tcl_AppendResult(interp, " lambda ", buffer, (char *) NULL);
    break;
  case CONSTRAINT_PLATE:
    Tcl_PrintDouble(interp, con->c.plate.pos, buffer);
    Tcl_AppendResult(interp, "plate height ", buffer, (char *) NULL);
    Tcl_PrintDouble(interp, con->c.plate.sigma, buffer);
    Tcl_AppendResult(interp, " sigma ", buffer, (char *) NULL);
    break;
  case CONSTRAINT_MAZE:
    Tcl_PrintDouble(interp, con->c.maze.nsphere, buffer);
    Tcl_AppendResult(interp, "maze nsphere ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.maze.dim, buffer);
    Tcl_AppendResult(interp, " dim ", buffer, (char *) NULL);
    Tcl_PrintDouble(interp, con->c.maze.sphrad, buffer);
    Tcl_AppendResult(interp, " sphrad ", buffer, (char *) NULL);
    Tcl_PrintDouble(interp, con->c.maze.cylrad, buffer);
    Tcl_AppendResult(interp, " cylrad ", buffer, (char *) NULL);
    sprintf(buffer, "%d", con->part_rep.p.type);
    Tcl_AppendResult(interp, " type ", buffer, (char *) NULL);
    sprintf(buffer, "%d", con->c.maze.penetrable);
    Tcl_AppendResult(interp, " penetrable ", buffer, (char *) NULL);
    break;
  case CONSTRAINT_PORE:
    Tcl_PrintDouble(interp, con->c.cyl.pos[0], buffer);
    Tcl_AppendResult(interp, "pore center ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.cyl.pos[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.cyl.pos[2], buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    Tcl_PrintDouble(interp, con->c.cyl.axis[0], buffer);
    Tcl_AppendResult(interp, " axis ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.cyl.axis[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.cyl.axis[2], buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    Tcl_PrintDouble(interp, con->c.cyl.rad, buffer);
    Tcl_AppendResult(interp, " radius ", buffer, (char *) NULL);
    Tcl_PrintDouble(interp, con->c.cyl.length, buffer);
    Tcl_AppendResult(interp, " length ", buffer, (char *) NULL);
    sprintf(buffer, "%d", con->part_rep.p.type);
    Tcl_AppendResult(interp, " type ", buffer, (char *) NULL);
    break;
  case CONSTRAINT_STOMATOCYTE:

    Tcl_PrintDouble(interp, con->c.stomatocyte.position_x, buffer);
    Tcl_AppendResult(interp, "stomatocyte center ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.stomatocyte.position_y, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.stomatocyte.position_z, buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    Tcl_PrintDouble(interp, con->c.stomatocyte.orientation_x, buffer);
    Tcl_AppendResult(interp, " orientation ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.stomatocyte.orientation_y, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.stomatocyte.orientation_z, buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    Tcl_PrintDouble(interp, con->c.stomatocyte.outer_radius, buffer);
    Tcl_AppendResult(interp, " outer radius ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.stomatocyte.inner_radius, buffer);
    Tcl_AppendResult(interp, " inner radius ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.stomatocyte.layer_width, buffer);
    Tcl_AppendResult(interp, " layer width ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.stomatocyte.direction, buffer);
    Tcl_AppendResult(interp, " direction ", buffer, (char *) NULL);
    sprintf(buffer, "%d", con->part_rep.p.type);
    Tcl_AppendResult(interp, " type ", buffer, (char *) NULL);
    sprintf(buffer, "%d", con->c.stomatocyte.penetrable);
    Tcl_AppendResult(interp, " penetrable ", buffer, (char *) NULL);
    Tcl_PrintDouble(interp, con->c.stomatocyte.reflecting, buffer);
    Tcl_AppendResult(interp, " reflecting ", buffer, (char *) NULL);
    break;
//ER
  case CONSTRAINT_EXT_MAGN_FIELD:
    Tcl_PrintDouble(interp, con->c.emfield.ext_magn_field[0], buffer);
    Tcl_AppendResult(interp, "ext_magn_field ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.emfield.ext_magn_field[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.emfield.ext_magn_field[2], buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    break; 
//end ER
  case CONSTRAINT_PLANE:
    Tcl_PrintDouble(interp, con->c.plane.pos[0], buffer);
    Tcl_AppendResult(interp, "plane cell ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.plane.pos[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.plane.pos[2], buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    sprintf(buffer, "%d", con->part_rep.p.type);
    Tcl_AppendResult(interp, " type ", buffer, (char *) NULL);
    break;

  default:
    sprintf(buffer, "%d", con->type);
    Tcl_AppendResult(interp, "unknown constraint type ", buffer, ".", (char *) NULL);
    return (TCL_OK);
  }

  return (TCL_OK);
}

int tclcommand_constraint_print(Tcl_Interp *interp)
{
  int i;
  if(n_constraints>0) Tcl_AppendResult(interp, "{", (char *)NULL);
  for (i = 0; i < n_constraints; i++) {
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

static int tclcommand_constraint_parse_wall(Constraint *con, Tcl_Interp *interp,
		    int argc, char **argv)
{
  int i;
  double norm;
  con->type = CONSTRAINT_WAL;
  /* invalid entries to start of */
  con->c.wal.n[0] = 
    con->c.wal.n[1] = 
    con->c.wal.n[2] = 0;
  con->c.wal.d = 0;
  con->c.wal.penetrable = 0;
  con->part_rep.p.type = -1;
  while (argc > 0) {
    if(!strncmp(argv[0], "normal", strlen(argv[0]))) {
      if(argc < 4) {
	Tcl_AppendResult(interp, "constraint wall normal <nx> <ny> <nz> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(con->c.wal.n[0])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[2], &(con->c.wal.n[1])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[3], &(con->c.wal.n[2])) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 4; argv += 4;
    }
    else if(!strncmp(argv[0], "dist", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint wall dist <d> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(con->c.wal.d)) == TCL_ERROR)
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
      if (Tcl_GetInt(interp, argv[1], &(con->c.wal.penetrable)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "reflecting", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint wall reflecting {0|1} expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetInt(interp, argv[1], &(con->c.wal.reflecting)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else
      break;
  }
  /* length of the normal vector */
  norm = SQR(con->c.wal.n[0])+SQR(con->c.wal.n[1])+SQR(con->c.wal.n[2]);
  if (norm < 1e-10 || con->part_rep.p.type < 0) {
    Tcl_AppendResult(interp, "usage: constraint wall normal <nx> <ny> <nz> dist <d> type <t> penetrable <0/1> reflecting <1/2>",
		     (char *) NULL);
    return (TCL_ERROR);    
  }
  /* normalize the normal vector */
  for (i=0;i<3;i++) con->c.wal.n[i] /= sqrt(norm);

  make_particle_type_exist(con->part_rep.p.type);

  return (TCL_OK);
}

static int tclcommand_constraint_parse_sphere(Constraint *con, Tcl_Interp *interp,
		      int argc, char **argv)
{
  con->type = CONSTRAINT_SPH;

  /* invalid entries to start of */
  con->c.sph.pos[0] = 
    con->c.sph.pos[1] = 
    con->c.sph.pos[2] = 0;
  con->c.sph.rad = 0;
  con->c.sph.direction = -1;
  con->c.sph.penetrable = 0;
  con->c.sph.reflecting = 0;
  con->part_rep.p.type = -1;

  while (argc > 0) {
    if(!strncmp(argv[0], "center", strlen(argv[0]))) {
      if(argc < 4) {
	Tcl_AppendResult(interp, "constraint sphere center <x> <y> <z> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(con->c.sph.pos[0])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[2], &(con->c.sph.pos[1])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[3], &(con->c.sph.pos[2])) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 4; argv += 4;
    }
    else if(!strncmp(argv[0], "radius", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint sphere radius <r> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(con->c.sph.rad)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "direction", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "-1/1 or inside/outside is expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (!strncmp(argv[1], "inside", strlen(argv[1])))
	con->c.sph.direction = -1;
      else if (!strncmp(argv[1], "outside", strlen(argv[1])))
	con->c.sph.direction = 1;
      else if (Tcl_GetDouble(interp, argv[1], &(con->c.sph.direction)) == TCL_ERROR)
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
      if (Tcl_GetInt(interp, argv[1], &(con->c.sph.penetrable)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "reflecting", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint sphere reflecting {0|1} expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetInt(interp, argv[1], &(con->c.sph.reflecting)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else
      break;
  }

  if (con->c.sph.rad < 0. || con->part_rep.p.type < 0) {
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
  /* invalid entries to start of */
  con->c.cyl.pos[0] = 
    con->c.cyl.pos[1] = 
    con->c.cyl.pos[2] = 0;
  con->c.cyl.axis[0] = 
    con->c.cyl.axis[1] = 
    con->c.cyl.axis[2] = 0;
  con->c.cyl.rad = 0;
  con->c.cyl.length = 0;
  con->c.cyl.direction = 0;
  con->c.cyl.penetrable = 0;
  con->part_rep.p.type = -1;
  con->c.cyl.reflecting = 0;
  while (argc > 0) {
    if(!strncmp(argv[0], "center", strlen(argv[0]))) {
      if(argc < 4) {
	Tcl_AppendResult(interp, "constraint cylinder center <x> <y> <z> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(con->c.cyl.pos[0])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[2], &(con->c.cyl.pos[1])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[3], &(con->c.cyl.pos[2])) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 4; argv += 4;
    }
    else if(!strncmp(argv[0], "axis", strlen(argv[0]))) {
      if(argc < 4) {
	Tcl_AppendResult(interp, "constraint cylinder axis <rx> <ry> <rz> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(con->c.cyl.axis[0])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[2], &(con->c.cyl.axis[1])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[3], &(con->c.cyl.axis[2])) == TCL_ERROR)
	return (TCL_ERROR);

      argc -= 4; argv += 4;    
    }
    else if(!strncmp(argv[0], "radius", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint cylinder radius <rad> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(con->c.cyl.rad)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "length", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint cylinder length <len> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(con->c.cyl.length)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "direction", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint cylinder direction <dir> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (!strncmp(argv[1], "inside", strlen(argv[1])))
	con->c.cyl.direction = -1;
      else if (!strncmp(argv[1], "outside", strlen(argv[1])))
	con->c.cyl.direction = 1;
      else if (Tcl_GetDouble(interp, argv[1], &(con->c.cyl.direction)) == TCL_ERROR)
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
      if (Tcl_GetInt(interp, argv[1], &(con->c.cyl.penetrable)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "reflecting", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint cylinder reflecting {0|1} expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetInt(interp, argv[1], &(con->c.cyl.reflecting)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else
      break;
  }

  axis_len=0.;
  for (i=0;i<3;i++)
    axis_len += SQR(con->c.cyl.axis[i]);

  if (con->c.cyl.rad < 0. || con->part_rep.p.type < 0 || axis_len < 1e-30 ||
      con->c.cyl.direction == 0 || con->c.cyl.length <= 0) {
    Tcl_AppendResult(interp, "usage: constraint cylinder center <x> <y> <z> axis <rx> <ry> <rz> radius <rad> length <length> direction <direction> type <t> penetrable <0/1> reflecting <1/2>",
		     (char *) NULL);
    return (TCL_ERROR);    
  }

  /*normalize the axis vector */
  axis_len = sqrt (axis_len);
  for (i=0;i<3;i++) {
    con->c.cyl.axis[i] /= axis_len;
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
  
  con->c.rhomboid.pos[0] = 
  con->c.rhomboid.pos[1] = 
  con->c.rhomboid.pos[2] = 0;
  
  con->c.rhomboid.a[0] = 
  con->c.rhomboid.a[1] = 
  con->c.rhomboid.a[2] = 0;
  
  con->c.rhomboid.b[0] = 
  con->c.rhomboid.b[1] = 
  con->c.rhomboid.b[2] = 0;
  
  con->c.rhomboid.c[0] = 
  con->c.rhomboid.c[1] = 
  con->c.rhomboid.c[2] = 0;
  
  con->c.rhomboid.direction = 0;
  
  con->c.rhomboid.penetrable = 0;
  
  con->part_rep.p.type = -1;
  
  con->c.rhomboid.reflecting = 0;
  
  while (argc > 0) {
    if(!strncmp(argv[0], "a", strlen(argv[0]))) {
      if(argc < 4) {
				Tcl_AppendResult(interp, "constraint rhomboid a <ax> <ay> <az> expected", (char *) NULL);
				return TCL_ERROR;
      }
      
      if(Tcl_GetDouble(interp, argv[1], &(con->c.rhomboid.a[0])) == TCL_ERROR ||
	 			 Tcl_GetDouble(interp, argv[2], &(con->c.rhomboid.a[1])) == TCL_ERROR ||
	  		 Tcl_GetDouble(interp, argv[3], &(con->c.rhomboid.a[2])) == TCL_ERROR)
				return TCL_ERROR;
				
			argc -= 4; argv += 4;    
    }
    else if(!strncmp(argv[0], "b", strlen(argv[0]))) {
      if(argc < 4) {
				Tcl_AppendResult(interp, "constraint rhomboid b <bx> <by> <bz> expected", (char *) NULL);
				return TCL_ERROR;
      }
      
      if(Tcl_GetDouble(interp, argv[1], &(con->c.rhomboid.b[0])) == TCL_ERROR ||
	 			 Tcl_GetDouble(interp, argv[2], &(con->c.rhomboid.b[1])) == TCL_ERROR ||
	  		 Tcl_GetDouble(interp, argv[3], &(con->c.rhomboid.b[2])) == TCL_ERROR)
				return TCL_ERROR;
				
			argc -= 4; argv += 4;    
    }
    else if(!strncmp(argv[0], "c", strlen(argv[0]))) {
      if(argc < 4) {
				Tcl_AppendResult(interp, "constraint rhomboid c <cx> <cy> <cz> expected", (char *) NULL);
				return TCL_ERROR;
      }
      
      if(Tcl_GetDouble(interp, argv[1], &(con->c.rhomboid.c[0])) == TCL_ERROR ||
	 			 Tcl_GetDouble(interp, argv[2], &(con->c.rhomboid.c[1])) == TCL_ERROR ||
	  		 Tcl_GetDouble(interp, argv[3], &(con->c.rhomboid.c[2])) == TCL_ERROR)
				return TCL_ERROR;
				
			argc -= 4; argv += 4;    
    }
    else if(!strncmp(argv[0], "corner", strlen(argv[0]))) { //this has to come after c
      if(argc < 4) {
				Tcl_AppendResult(interp, "constraint rhomboid corner <x> <y> <z> expected", (char *) NULL);
				return TCL_ERROR;
      }
      
      if(Tcl_GetDouble(interp, argv[1], &(con->c.rhomboid.pos[0])) == TCL_ERROR ||
				 Tcl_GetDouble(interp, argv[2], &(con->c.rhomboid.pos[1])) == TCL_ERROR ||
	  		 Tcl_GetDouble(interp, argv[3], &(con->c.rhomboid.pos[2])) == TCL_ERROR)
				return TCL_ERROR;
				
      argc -= 4; argv += 4;
    }
    else if(!strncmp(argv[0], "direction", strlen(argv[0]))) {
      if (argc < 2) {
				Tcl_AppendResult(interp, "constraint rhomboid direction {inside|outside} expected", (char *) NULL);
				return (TCL_ERROR);
      }
      
      if(!strncmp(argv[1], "inside", strlen(argv[1])))
				con->c.rhomboid.direction = -1;
      else if(!strncmp(argv[1], "outside", strlen(argv[1])))
				con->c.rhomboid.direction = 1;
      else if(Tcl_GetDouble(interp, argv[1], &(con->c.rhomboid.direction)) == TCL_ERROR)
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
      
      if (Tcl_GetInt(interp, argv[1], &(con->c.rhomboid.penetrable)) == TCL_ERROR)
				return TCL_ERROR;
				
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "reflecting", strlen(argv[0]))) {
      if (argc < 2) {
				Tcl_AppendResult(interp, "constraint rhomboid reflecting {0|1} expected", (char *) NULL);
				return TCL_ERROR;
      }
      
      if (Tcl_GetInt(interp, argv[1], &(con->c.rhomboid.reflecting)) == TCL_ERROR)
				return TCL_ERROR;
				
      argc -= 2; argv += 2;
    }
    else {
			Tcl_AppendResult(interp, "Error: Unknown parameter ", argv[0], " in constraint rhomboid", (char *) NULL);
			return TCL_ERROR;
    }
  }

  if( (con->c.rhomboid.a[0] == 0. && con->c.rhomboid.a[1] == 0. && con->c.rhomboid.a[2] == 0.) ||
  		(con->c.rhomboid.b[0] == 0. && con->c.rhomboid.b[1] == 0. && con->c.rhomboid.b[2] == 0.) ||
  		(con->c.rhomboid.c[0] == 0. && con->c.rhomboid.c[1] == 0. && con->c.rhomboid.c[2] == 0.) ||
  		con->part_rep.p.type < 0 || con->c.rhomboid.direction == 0) {
    Tcl_AppendResult(interp, "usage: constraint rhomboid corner <x> <y> <z> a <ax> <ay> <az> b <bx> <by> <bz> c <cx> <cy> <cz> direction {inside|outside} type <t> [penetrable <0|1>] [reflecting <1|2>]",
		     (char *) NULL);
    return TCL_ERROR;    
  }
                     
  //If the trihedron a, b, c is left handed, then inside and outside will be exchanged since all normals will be reversed. This compensates  for that, so that the user doesn't have to take care of the handedness.
  triple_product = con->c.rhomboid.a[0]*( con->c.rhomboid.b[1]*con->c.rhomboid.c[2] - con->c.rhomboid.b[2]*con->c.rhomboid.c[1] ) +
                   con->c.rhomboid.a[1]*( con->c.rhomboid.b[2]*con->c.rhomboid.c[0] - con->c.rhomboid.b[0]*con->c.rhomboid.c[2] ) + 
                   con->c.rhomboid.a[2]*( con->c.rhomboid.b[0]*con->c.rhomboid.c[1] - con->c.rhomboid.b[1]*con->c.rhomboid.c[0] );
                
  if(triple_product < 0.)
  {    
    tmp[0] = con->c.rhomboid.a[0];
    tmp[1] = con->c.rhomboid.a[1];
    tmp[2] = con->c.rhomboid.a[2];
    
    con->c.rhomboid.a[0] = con->c.rhomboid.b[0];
    con->c.rhomboid.a[1] = con->c.rhomboid.b[1];
    con->c.rhomboid.a[2] = con->c.rhomboid.b[2];
    
    con->c.rhomboid.b[0] = tmp[0];
    con->c.rhomboid.b[1] = tmp[1];
    con->c.rhomboid.b[2] = tmp[2];
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
  /* invalid entries to start of */
  con->c.pore.pos[0] = 
    con->c.pore.pos[1] = 
    con->c.pore.pos[2] = 0;
  con->c.pore.axis[0] = 
    con->c.pore.axis[1] = 
    con->c.pore.axis[2] = 0;
  con->c.pore.rad_left = 0;
  con->c.pore.rad_right = 0;
  con->c.pore.length = 0;
  con->c.pore.reflecting = 0;
  con->part_rep.p.type = -1;
  con->c.pore.smoothing_radius = 1.;
  while (argc > 0) {
    if(!strncmp(argv[0], "center", strlen(argv[0]))) {
      if(argc < 4) {
	Tcl_AppendResult(interp, "constraint pore center <x> <y> <z> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(con->c.pore.pos[0])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[2], &(con->c.pore.pos[1])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[3], &(con->c.pore.pos[2])) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 4; argv += 4;
    }
    else if(!strncmp(argv[0], "axis", strlen(argv[0]))) {
      if(argc < 4) {
	Tcl_AppendResult(interp, "constraint pore axis <rx> <ry> <rz> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(con->c.pore.axis[0])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[2], &(con->c.pore.axis[1])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[3], &(con->c.pore.axis[2])) == TCL_ERROR)
	return (TCL_ERROR);

      argc -= 4; argv += 4;    
    }
    else if(!strncmp(argv[0], "radius", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint pore radius <rad> expected", (char *) NULL);
	return (TCL_ERROR);
      }  
      if (Tcl_GetDouble(interp, argv[1], &(con->c.pore.rad_left)) == TCL_ERROR)
	return (TCL_ERROR);
      con->c.pore.rad_right =  con->c.pore.rad_left; 
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "smoothing_radius", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint pore smoothing_radius <smoothing_radius> expected", (char *) NULL);
	return (TCL_ERROR);
      }  
      if (Tcl_GetDouble(interp, argv[1], &(con->c.pore.smoothing_radius)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "radii", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint pore radii <rad_left> <rad_right> expected", (char *) NULL);
	return (TCL_ERROR);
      }  
      if (Tcl_GetDouble(interp, argv[1], &(con->c.pore.rad_left)) == TCL_ERROR)
	return (TCL_ERROR);
      if (Tcl_GetDouble(interp, argv[2], &(con->c.pore.rad_right)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 3; argv += 3;
    }
    else if(!strncmp(argv[0], "length", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint pore length <len/2> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(con->c.pore.length)) == TCL_ERROR)
	return (TCL_ERROR);
//      con->c.pore.length *= 2;
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
      if (Tcl_GetInt(interp, argv[1], &(con->c.pore.reflecting)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else
      break;
  }

  axis_len=0.;
  for (i=0;i<3;i++)
    axis_len += SQR(con->c.pore.axis[i]);

  if (con->c.pore.rad_left < 0. || con->c.pore.rad_right < 0. || con->part_rep.p.type < 0 || axis_len < 1e-30 ||
      con->c.pore.length <= 0) {
    Tcl_AppendResult(interp, "usage: constraint pore center <x> <y> <z> axis <rx> <ry> <rz> radius <rad> length <length/2> type <t>",
		     (char *) NULL);
    return (TCL_ERROR);
  }

  /*normalize the axis vector */
  axis_len = sqrt (axis_len);
  for (i=0;i<3;i++) {
    con->c.pore.axis[i] /= axis_len;
  }
  
  make_particle_type_exist(con->part_rep.p.type);

  return (TCL_OK);
}

static int tclcommand_constraint_parse_rod(Constraint *con, Tcl_Interp *interp,
		   int argc, char **argv)
{
  con->type = CONSTRAINT_ROD;
  con->part_rep.p.type = -1;
  con->c.rod.pos[0] = con->c.rod.pos[1] = 0;
  con->c.rod.lambda = 0;
  while (argc > 0) {
    if(!strncmp(argv[0], "center", strlen(argv[0]))) {
      if(argc < 3) {
	Tcl_AppendResult(interp, "constraint rod center <px> <py> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(con->c.rod.pos[0])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[2], &(con->c.rod.pos[1])) == TCL_ERROR)
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
      if (Tcl_GetDouble(interp, argv[1], &(con->c.rod.lambda)) == TCL_ERROR)
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
  con->c.plate.pos = 0;
  con->c.plate.sigma = 0;
  while (argc > 0) {
    if(!strncmp(argv[0], "height", strlen(argv[0]))) {
      if(argc < 2) {
	Tcl_AppendResult(interp, "constraint plate height <pz> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(con->c.plate.pos)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "sigma", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint rod sigma <s> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(con->c.plate.sigma)) == TCL_ERROR)
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

  /* invalid entries to start of */

  con->c.stomatocyte.position_x = -M_PI;
  con->c.stomatocyte.position_y = -M_PI;
  con->c.stomatocyte.position_z = -M_PI;
  con->c.stomatocyte.orientation_x = -M_PI;
  con->c.stomatocyte.orientation_y = -M_PI;
  con->c.stomatocyte.orientation_z = -M_PI;
  con->c.stomatocyte.outer_radius = -1.0;
  con->c.stomatocyte.inner_radius = -1.0;
  con->c.stomatocyte.layer_width = -1.0;
  con->c.stomatocyte.direction = 0;
  con->c.stomatocyte.penetrable = 0;
  con->c.stomatocyte.reflecting = 0;
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

      if ( !ARG_IS_D( 1, con->c.stomatocyte.position_x ) ||
	         !ARG_IS_D( 2, con->c.stomatocyte.position_y ) ||
	         !ARG_IS_D( 3, con->c.stomatocyte.position_z ) )
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

      if ( !ARG_IS_D( 1, con->c.stomatocyte.orientation_x ) ||
	         !ARG_IS_D( 2, con->c.stomatocyte.orientation_y ) ||
	         !ARG_IS_D( 3, con->c.stomatocyte.orientation_z ) )
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

      if ( !ARG_IS_D(1, con->c.stomatocyte.outer_radius ) )
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

      if ( !ARG_IS_D( 1, con->c.stomatocyte.inner_radius ) )
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

      if ( !ARG_IS_D( 1, con->c.stomatocyte.layer_width ) )
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
	      con->c.stomatocyte.direction = -1;
      else if ( ARG_IS_S( 1, "outside" ) )
	      con->c.stomatocyte.direction = 1;
      else if ( !ARG_IS_D( 1, con->c.stomatocyte.direction ) )
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

      if ( !ARG_IS_I( 1, con->c.stomatocyte.penetrable ) ) 
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

      if ( !ARG_IS_I( 1, con->c.stomatocyte.reflecting ) )
	      return (TCL_ERROR);

      argc -= 2; argv += 2;
    }
    else
      break;
  }

  if ( con->c.stomatocyte.outer_radius < 0.0 || 
       con->c.stomatocyte.inner_radius < 0.0 || 
       con->c.stomatocyte.layer_width < 0.0 ) 
  {
    Tcl_AppendResult(interp, "stomatocyte radii and width have to be greater than zero",
		     (char *) NULL);
    return (TCL_ERROR);    
  }

  if ( con->c.stomatocyte.outer_radius < con->c.stomatocyte.inner_radius || 
       con->c.stomatocyte.inner_radius < con->c.stomatocyte.layer_width ||
       con->c.stomatocyte.outer_radius < con->c.stomatocyte.layer_width ) 
  {
    Tcl_AppendResult(interp, "stomatocyte requires layer_width < inner_radius < outer_radius",
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

  /* invalid entries to start of */
  con->c.maze.nsphere = 0.;
  con->c.maze.dim = -1.;
  con->c.maze.sphrad = 0.;
  con->c.maze.cylrad = -1.;
  con->c.maze.penetrable = 0;
  con->part_rep.p.type = -1;

  while (argc > 0) {
    if(!strncmp(argv[0], "nsphere", strlen(argv[0]))) {
      if(argc < 4) {
	Tcl_AppendResult(interp, "constraint maze nsphere <n> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(con->c.maze.nsphere)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "dim", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint maze dim <d> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(con->c.maze.dim)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "sphrad", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint maze sphrad <r> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(con->c.maze.sphrad)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "cylrad", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint maze cylrad <r> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(con->c.maze.cylrad)) == TCL_ERROR)
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
      if (Tcl_GetInt(interp, argv[1], &(con->c.maze.penetrable)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else
      break;
  }

  if (con->c.maze.sphrad < 0. || con->c.maze.cylrad < 0. || con->part_rep.p.type < 0 || con->c.maze.dim < 0) {
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

  for(i=0; i<3; i++)
     con->c.emfield.ext_magn_field[i] = 0.;

  if(argc < 3) {
      Tcl_AppendResult(interp, "usage: constraint ext_magn_field <x> <y> <z>", (char *) NULL);
      return (TCL_ERROR);
  }
  for(i=0; i<3; i++){
     if (Tcl_GetDouble(interp, argv[i], &(con->c.emfield.ext_magn_field[i])) == TCL_ERROR)
	return (TCL_ERROR);
  }
  argc -= 3; argv += 3;

  return (TCL_OK);
}
//end ER

static int tclcommand_constraint_parse_plane_cell(Constraint *con, Tcl_Interp *interp,
                      int argc, char **argv)
{
  con->type = CONSTRAINT_PLANE;

  /* invalid entries to start of */
  con->c.plane.pos[0] = 
    con->c.plane.pos[1] = 
    con->c.plane.pos[2] = 0;
  con->part_rep.p.type = -1;

  while (argc > 0) {
    if(!strncmp(argv[0], "cell", strlen(argv[0]))) {
      if(argc < 4) {
        Tcl_AppendResult(interp, "constraint plane cell <x> <y> <z> expected", (char *) NULL);
        return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(con->c.plane.pos[0])) == TCL_ERROR ||
          Tcl_GetDouble(interp, argv[2], &(con->c.plane.pos[1])) == TCL_ERROR ||
          Tcl_GetDouble(interp, argv[3], &(con->c.plane.pos[2])) == TCL_ERROR)
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
  double dist=1e100;
  double mindist = 1e100;
  int n;
  char buffer[TCL_DOUBLE_SPACE];
  if (n_constraints==0) {
    Tcl_AppendResult(interp, "Error in constraint mindist_position: no constraints defined\n", (char*) NULL);
    return TCL_ERROR;
  }

  Particle* p1=0;
  if (ARG_IS_D(0, pos[0]) && ARG_IS_D(1, pos[1]) && ARG_IS_D(2, pos[2])) {
    for(n=0;n<n_constraints;n++) {
      switch(constraints[n].type) {
        case CONSTRAINT_WAL: 
          if ( !constraints[n].c.wal.penetrable )
	          calculate_wall_dist(p1, pos, &constraints[n].part_rep, &constraints[n].c.wal, &dist, vec); 
          else
            dist=1e100;
          break;
        case CONSTRAINT_SPH:
          if ( !constraints[n].c.sph.penetrable )
	          calculate_sphere_dist(p1, pos, &constraints[n].part_rep, &constraints[n].c.sph, &dist, vec); 
          else
            dist=1e100;
          break;
        case CONSTRAINT_CYL: 
          if ( !constraints[n].c.cyl.penetrable )
	          calculate_cylinder_dist(p1, pos, &constraints[n].part_rep, &constraints[n].c.cyl, &dist, vec); 
          else
            dist=1e100;
          break;
        case CONSTRAINT_RHOMBOID: 
          if ( !constraints[n].c.rhomboid.penetrable )
	          calculate_rhomboid_dist(p1, pos, &constraints[n].part_rep, &constraints[n].c.rhomboid, &dist, vec); 
          else
            dist=1e100;
          break;
        case CONSTRAINT_MAZE: 
          if ( !constraints[n].c.maze.penetrable )
	          calculate_maze_dist(p1, pos, &constraints[n].part_rep, &constraints[n].c.maze, &dist, vec); 
          else
            dist=1e100;
          break;
        case CONSTRAINT_STOMATOCYTE:
          if ( !constraints[n].c.stomatocyte.penetrable )
	          calculate_stomatocyte_dist(p1, pos, &constraints[n].part_rep, &constraints[n].c.stomatocyte, &dist, vec); 
          else
            dist=1e100;
          break;
        case CONSTRAINT_PORE: 
	        calculate_pore_dist(p1, pos, &constraints[n].part_rep, &constraints[n].c.pore, &dist, vec); 
          break;
        case CONSTRAINT_PLANE:
	        calculate_plane_dist(p1, pos, &constraints[n].part_rep, &constraints[n].c.plane, &dist, vec); 
          break;
      }
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
  double dist=1e100;
  double mindist = 1e100;
  double minvec[3] = { 1e100, 1e100, 1e100};
  int n;
  char buffer[TCL_DOUBLE_SPACE];
  if (n_constraints==0) {
    Tcl_AppendResult(interp, "Error in constraint mindist_position: no constraints defined\n", (char*) NULL);
    return TCL_ERROR;
  }
  if (argc < 3) {
    Tcl_AppendResult(interp, "Usage: constraint mindist_position_vec <px> <py> <py>\n", (char*) NULL);
    return TCL_ERROR;
  }

  Particle* p1=0;
  if (ARG_IS_D(0, pos[0]) && ARG_IS_D(1, pos[1]) && ARG_IS_D(2, pos[2])) {
    for(n=0;n<n_constraints;n++) {
      switch(constraints[n].type) {
        case CONSTRAINT_WAL: 
          if ( !constraints[n].c.wal.penetrable )
	          calculate_wall_dist(p1, pos, &constraints[n].part_rep, &constraints[n].c.wal, &dist, vec); 
          else
            dist=1e100;
          break;
        case CONSTRAINT_SPH:
          if ( !constraints[n].c.sph.penetrable )
	          calculate_sphere_dist(p1, pos, &constraints[n].part_rep, &constraints[n].c.sph, &dist, vec); 
          else
            dist=1e100;
          break;
        case CONSTRAINT_CYL: 
          if ( !constraints[n].c.cyl.penetrable )
	          calculate_cylinder_dist(p1, pos, &constraints[n].part_rep, &constraints[n].c.cyl, &dist, vec); 
          else
            dist=1e100;
          break;
        case CONSTRAINT_RHOMBOID: 
          if ( !constraints[n].c.rhomboid.penetrable )
	          calculate_rhomboid_dist(p1, pos, &constraints[n].part_rep, &constraints[n].c.rhomboid, &dist, vec); 
          else
            dist=1e100;
          break;
        case CONSTRAINT_MAZE: 
          if ( !constraints[n].c.maze.penetrable )
	          calculate_maze_dist(p1, pos, &constraints[n].part_rep, &constraints[n].c.maze, &dist, vec); 
          else
            dist=1e100;
          break;
        case CONSTRAINT_STOMATOCYTE:
          if ( !constraints[n].c.stomatocyte.penetrable )
	          calculate_stomatocyte_dist(p1, pos, &constraints[n].part_rep, &constraints[n].c.stomatocyte, &dist, vec); 
          else
            dist=1e100;
          break;
        case CONSTRAINT_PORE: 
	        calculate_pore_dist(p1, pos, &constraints[n].part_rep, &constraints[n].c.pore, &dist, vec); 
          break;
        case CONSTRAINT_PLANE:
	        calculate_plane_dist(p1, pos, &constraints[n].part_rep, &constraints[n].c.plane, &dist, vec); 
          break;
      }
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
  else if(!strncmp(argv[1], "stomatocyte", strlen(argv[1]))) {
    status = tclcommand_constraint_parse_stomatocyte(generate_constraint(),interp, argc - 2, argv + 2);
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
    if(c_num < 0 || c_num >= n_constraints) {
      Tcl_AppendResult(interp, "constraint does not exist",(char *) NULL);
      return (TCL_ERROR);
    }
    tclprint_to_result_ConstraintForce(interp, c_num);
    status  = TCL_OK;
  }
  else if(!strncmp(argv[1], "delete", strlen(argv[1]))) {
    if(argc < 3) {
      /* delete all */
      mpi_bcast_constraint(-2);
      status = TCL_OK;
    }
    else {
      if(Tcl_GetInt(interp, argv[2], &(c_num)) == TCL_ERROR) return (TCL_ERROR);
      if(c_num < 0 || c_num >= n_constraints) {
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
    Tcl_AppendResult(interp, "possible constraints: wall sphere cylinder maze pore stomatocyte ext_magn_field or constraint delete {c} to delete constraint(s)",(char *) NULL);
    return (TCL_ERROR);
  }

  return gather_runtime_errors(interp, status);

#else /* !defined(CONSTRAINTS) */
  Tcl_AppendResult(interp, "Constraints not compiled in!" ,(char *) NULL);
  return (TCL_ERROR);
#endif
}

