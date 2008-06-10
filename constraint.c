/** \file constraint.c
    Implementation of \ref constraint.h "constraint.h", here it's just the parsing stuff.
*/
#include "constraint.h"

#ifdef CONSTRAINTS
int printConstraintToResult(Tcl_Interp *interp, int i)
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
  default:
    sprintf(buffer, "%d", con->type);
    Tcl_AppendResult(interp, "unknown constraint type ", buffer, ".", (char *) NULL);
    return (TCL_OK);
  }

  return (TCL_OK);
}

int constraint_print_all(Tcl_Interp *interp)
{
  int i;
  if(n_constraints>0) Tcl_AppendResult(interp, "{", (char *)NULL);
  for (i = 0; i < n_constraints; i++) {
    if(i>0) Tcl_AppendResult(interp, " {", (char *)NULL);
    printConstraintToResult(interp, i);
    Tcl_AppendResult(interp, "}", (char *)NULL);
  }
  return (TCL_OK);
}

void printConstraintForceToResult(Tcl_Interp *interp, int con)
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

Constraint *generate_constraint()
{
  n_constraints++;
  constraints = realloc(constraints,n_constraints*sizeof(Constraint));
  constraints[n_constraints-1].type = CONSTRAINT_NONE;
  constraints[n_constraints-1].part_rep.p.identity = -n_constraints;
  
  return &constraints[n_constraints-1];
}

int constraint_wall(Constraint *con, Tcl_Interp *interp,
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
    else
      break;
  }
  /* length of the normal vector */
  norm = SQR(con->c.wal.n[0])+SQR(con->c.wal.n[1])+SQR(con->c.wal.n[2]);
  if (norm < 1e-10 || con->part_rep.p.type < 0) {
    Tcl_AppendResult(interp, "usage: constraint wall normal <nx> <ny> <nz> dist <d> type <t>",
		     (char *) NULL);
    return (TCL_ERROR);    
  }
  /* normalize the normal vector */
  for (i=0;i<3;i++) con->c.wal.n[i] /= sqrt(norm);

  make_particle_type_exist(con->part_rep.p.type);

  return (TCL_OK);
}

int constraint_sphere(Constraint *con, Tcl_Interp *interp,
		      int argc, char **argv)
{
  con->type = CONSTRAINT_SPH;

  /* invalid entries to start of */
  con->c.sph.pos[0] = 
    con->c.sph.pos[1] = 
    con->c.sph.pos[2] = 0;
  con->c.sph.rad = 0;
  con->c.sph.direction = -1;
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
    else
      break;
  }

  if (con->c.sph.rad < 0. || con->part_rep.p.type < 0) {
    Tcl_AppendResult(interp, "usage: constraint sphere center <x> <y> <z> radius <d> direction <direction> type <t>",
		     (char *) NULL);
    return (TCL_ERROR);    
  }

  make_particle_type_exist(con->part_rep.p.type);

  return (TCL_OK);
}

int constraint_cylinder(Constraint *con, Tcl_Interp *interp,
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
  con->part_rep.p.type = -1;
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
    else
      break;
  }

  axis_len=0.;
  for (i=0;i<3;i++)
    axis_len += SQR(con->c.cyl.axis[i]);

  if (con->c.cyl.rad < 0. || con->part_rep.p.type < 0 || axis_len < 1e-30 ||
      con->c.cyl.direction == 0 || con->c.cyl.length <= 0) {
    Tcl_AppendResult(interp, "usage: constraint cylinder center <x> <y> <z> axis <rx> <ry> <rz> radius <rad> length <length> direction <direction> type <t>",
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

int constraint_pore(Constraint *con, Tcl_Interp *interp,
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
  con->c.pore.rad = 0;
  con->c.pore.length = 0;
  con->part_rep.p.type = -1;
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
      if (Tcl_GetDouble(interp, argv[1], &(con->c.pore.rad)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "length", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint pore length <len/2> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(con->c.pore.length)) == TCL_ERROR)
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
    else
      break;
  }

  axis_len=0.;
  for (i=0;i<3;i++)
    axis_len += SQR(con->c.pore.axis[i]);

  if (con->c.pore.rad < 0. || con->part_rep.p.type < 0 || axis_len < 1e-30 ||
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

int constraint_rod(Constraint *con, Tcl_Interp *interp,
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

int constraint_plate(Constraint *con, Tcl_Interp *interp,
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


int constraint_maze(Constraint *con, Tcl_Interp *interp,
		      int argc, char **argv)
{
  con->type = CONSTRAINT_MAZE;

  /* invalid entries to start of */
  con->c.maze.nsphere = 0.;
  con->c.maze.dim = -1.;
  con->c.maze.sphrad = 0.;
  con->c.maze.cylrad = -1.;
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
    else
      break;
  }

  if (con->c.maze.sphrad < 0. || con->c.maze.cylrad < 0. || con->part_rep.p.type < 0 || con->c.maze.dim < 0) {
    Tcl_AppendResult(interp, "usage: constraint maze nsphere <n> dim <d> sphrad <r> cylrad <r> type <t>",
		     (char *) NULL);
    return (TCL_ERROR);    
  }

  make_particle_type_exist(con->part_rep.p.type);

  return (TCL_OK);
}

//ER
int constraint_ext_magn_field(Constraint *con, Tcl_Interp *interp,
		      int argc, char **argv)
{
  int i;
  con->type = CONSTRAINT_EXT_MAGN_FIELD;

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

#endif


int constraint(ClientData _data, Tcl_Interp *interp,
	       int argc, char **argv)
{
#ifdef CONSTRAINTS
  int status, c_num;

  if (argc < 2) return constraint_print_all(interp);
  
  if(!strncmp(argv[1], "wall", strlen(argv[1]))) {
    status = constraint_wall(generate_constraint(),interp, argc - 2, argv + 2);
    mpi_bcast_constraint(-1);
  }
  else if(!strncmp(argv[1], "sphere", strlen(argv[1]))) {
    status = constraint_sphere(generate_constraint(),interp, argc - 2, argv + 2);
    mpi_bcast_constraint(-1);
  }
  else if(!strncmp(argv[1], "cylinder", strlen(argv[1]))) {
    status = constraint_cylinder(generate_constraint(),interp, argc - 2, argv + 2);
    mpi_bcast_constraint(-1);
  }
  else if(!strncmp(argv[1], "rod", strlen(argv[1]))) {
    status = constraint_rod(generate_constraint(),interp, argc - 2, argv + 2);
    mpi_bcast_constraint(-1);
  }
  else if(!strncmp(argv[1], "plate", strlen(argv[1]))) {
    status = constraint_plate(generate_constraint(),interp, argc - 2, argv + 2);
    mpi_bcast_constraint(-1);
  }
  else if(!strncmp(argv[1], "maze", strlen(argv[1]))) {
    status = constraint_maze(generate_constraint(),interp, argc - 2, argv + 2);
    mpi_bcast_constraint(-1);
  }
  else if(!strncmp(argv[1], "pore", strlen(argv[1]))) {
    status = constraint_pore(generate_constraint(),interp, argc - 2, argv + 2);
    mpi_bcast_constraint(-1);
  }
  //ER
  else if(!strncmp(argv[1], "ext_magn_field", strlen(argv[1]))) {
    status = constraint_ext_magn_field(generate_constraint(),interp, argc - 2, argv + 2);
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
    printConstraintForceToResult(interp, c_num);
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
    printConstraintToResult(interp, c_num);
    status = TCL_OK;
  }
  else {
  //ER "ext_magn_field" was put in the next line //end ER
    Tcl_AppendResult(interp, "possible constraints: wall sphere cylinder maze pore ext_magn_field or constraint delete {c} to delete constraint(s)",(char *) NULL);
    return (TCL_ERROR);
  }

  return mpi_gather_runtime_errors(interp, status);

#else /* !defined(CONSTRAINTS) */
  Tcl_AppendResult(interp, "Constraints not compiled in!" ,(char *) NULL);
  return (TCL_ERROR);
#endif
}
