#include "tcl_datafield.h"
#include "global.h"
#include "communication.h"
#include "grid.h"
#include "binaryfile.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* cwz-build-comman: ssh chakotay "builtin cd /nhomes/janeway/axel/progs/tcl_md; make" 
   cwz-build-command: make
*/

/**************************************************************
 * function prototypes
 **************************************************************/

/** tcl procedure for datafield access */
int setmd(ClientData data, Tcl_Interp *interp,
	  int argc, char **argv);
/** tcl procedure for interaction access */
int inter(ClientData data, Tcl_Interp *interp,
	  int argc, char **argv);
/** tcl procedure for particle access */
int part(ClientData data, Tcl_Interp *interp,
	 int argc, char **argv);
/** tcl procedure for writing particle data */
int writemd(ClientData data, Tcl_Interp *interp,
	    int argc, char **argv);
/** tcl procedure for writing particle data */
int readmd(ClientData data, Tcl_Interp *interp,
	   int argc, char **argv);

void tcl_datafield_init(Tcl_Interp *interp)
{
  Tcl_CreateCommand(interp, "setmd", setmd, 0, NULL);
  Tcl_CreateCommand(interp, "inter", inter, 0, NULL);
  Tcl_CreateCommand(interp, "part", part, 0, NULL);
  Tcl_CreateCommand(interp, "writemd", writemd, 0, NULL);
  Tcl_CreateCommand(interp, "readmd", readmd, 0, NULL);
}

int setmd(ClientData data, Tcl_Interp *interp,
	  int argc, char **argv)
{
  double dbuffer[MAX_DIMENSION];
  int    ibuffer[MAX_DIMENSION];
  char   buffer[TCL_DOUBLE_SPACE + 5];
  int i, j;
  int status;

  if (argc < 2) {
    Tcl_AppendResult(interp, "wrong # args:  should be \"",
		     argv[0], " <variable> ?value? ?value? ...\"",
		     (char *) NULL);
    return (TCL_ERROR);
  }

  for (i = 0; fields[i].data != NULL; i++) {
    if (!strncmp(argv[1], fields[i].name, strlen(argv[1]))) {
      if (argc >= 3) {
	/* set */
	/* parse in data */
	if (argc != 2 + fields[i].dimension) {
	  sprintf(buffer, "%d", fields[i].dimension);	  
	  Tcl_AppendResult(interp, "\"", argv[1],
			   "\" has dimension ",
			   buffer, (char *) NULL);
	  sprintf(buffer, " not %d", argc - 2);	  
	  Tcl_AppendResult(interp, buffer, (char *) NULL);
	  return (TCL_ERROR);
	}

	/* get new value */
	for (j = 0; j < fields[i].dimension; j++) {
	  switch (fields[i].type) {
	  case TYPE_INT:
	    if (Tcl_GetInt(interp, argv[2 + j], &ibuffer[j]) == TCL_ERROR)
	      return (TCL_ERROR);
	    break;
	  case TYPE_DOUBLE:
	    if (Tcl_GetDouble(interp, argv[2 + j], &dbuffer[j]))
	      return (TCL_ERROR);
	    break;
	  default: ;
	  }
	}

	/* call changeproc */
	switch(fields[i].type) {
	case TYPE_INT:
	  status = fields[i].changeproc(interp, ibuffer);
	  break;
	case TYPE_DOUBLE:
	  status = fields[i].changeproc(interp, dbuffer);
	  break;
	default:
	  status = TCL_ERROR;
	}
	if (status != TCL_OK)
	  return TCL_ERROR;

	mpi_bcast_parameter(i);
      }

      /* get */
      for (j = 0; j < fields[i].dimension; j++) {
	switch (fields[i].type) {
	case TYPE_INT:
	  sprintf(buffer, "%d", ((int *)fields[i].data)[j]);
	  break;
	case TYPE_DOUBLE:
	  Tcl_PrintDouble(interp, ((double *)fields[i].data)[j], buffer);
	  break;
	default: ;
	}
	Tcl_AppendResult(interp, buffer, (char *) NULL);
	if (j < fields[i].dimension - 1)
	  Tcl_AppendResult(interp, " ", (char *) NULL);
      }

      return (TCL_OK);
    }
  }
  Tcl_AppendResult(interp, "unkown variable \"",
		   argv[1], "\"", (char *) NULL);
  return (TCL_ERROR);
}

int part(ClientData data, Tcl_Interp *interp,
	 int argc, char **argv)
{
  int part_num = -1;
  int node, j;
  char buffer[50 + TCL_DOUBLE_SPACE];

  if (argc < 2) {
    Tcl_AppendResult(interp, "wrong # args:  should be \"",
		     argv[0], " <part num> ?what? ?value?\"", (char *) NULL);
    return (TCL_ERROR);
  }
  if (!processor_grid_is_set())
    setup_processor_grid();

  part_num = atol(argv[1]);
  if (part_num < 0) {
    Tcl_AppendResult(interp, "illegal particle", (char *) NULL);
    return (TCL_ERROR);
  }

  if (!particle_node)
    build_particle_node();

  node = (part_num < n_total_particles) ? particle_node[part_num] : -1;
 
  /* print out particle information */
  if (argc == 2) {
    Particle part;
    /* retrieve particle data */
    if (node == -1) {
      sprintf(buffer, "particle %d does not exist", part_num);
      Tcl_AppendResult(interp, buffer, (char *) NULL);
      return (TCL_ERROR);
    }
    mpi_recv_part(node, part_num, &part);
    Tcl_PrintDouble(interp, part.p[0], buffer);
    Tcl_AppendResult(interp, "p ", buffer, " ", (char *)NULL);
    Tcl_PrintDouble(interp, part.p[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
    Tcl_PrintDouble(interp, part.p[2], buffer);
    Tcl_AppendResult(interp, buffer, " type ", (char *)NULL);
    sprintf(buffer, "%d", part.type);
    Tcl_AppendResult(interp, buffer, " q ", (char *)NULL);
    Tcl_PrintDouble(interp, part.q, buffer);
    Tcl_AppendResult(interp, buffer, " v ", (char *)NULL);
    Tcl_PrintDouble(interp, part.v[0], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
    Tcl_PrintDouble(interp, part.v[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
    Tcl_PrintDouble(interp, part.v[2], buffer);
    Tcl_AppendResult(interp, buffer, " f ", (char *)NULL);
    Tcl_PrintDouble(interp, part.f[0], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
    Tcl_PrintDouble(interp, part.f[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
    Tcl_PrintDouble(interp, part.f[2], buffer);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
    /* FIXME: print bonding structure here */
    return (TCL_OK);
  }
  
  /* set particle data */
  argc -= 2;
  argv += 2;
  while (argc > 0) {
    if (!strncmp(argv[0], "pos", strlen(argv[0]))) {
      double pos[3];
      if (argc < 4) {
	Tcl_AppendResult(interp, "pos requires 3 arguments", (char *) NULL);
	return (TCL_ERROR);
      }
      /* set position */
      for (j = 0; j < 3; j++) {
	if (Tcl_GetDouble(interp, argv[1 + j], &pos[j]) == TCL_ERROR)
	  return (TCL_ERROR);
      }
    
      if (node == -1) {
	/* spatial position */
	node = find_node(pos);
	mpi_attach_particle(part_num, node);
	map_particle_node(part_num, node);
      }

      mpi_send_pos(node, part_num, pos);

      argc -= 4;
      argv += 4;
    }
    else {
      if (node == -1) {
	Tcl_AppendResult(interp, "set particle position first", (char *)NULL);
	return (TCL_ERROR);
      }
      
      if (!strncmp(argv[0], "q", strlen(argv[0]))) {
	double q;
	if (argc < 2) {
	  Tcl_AppendResult(interp, "q requires 1 argument", (char *) NULL);
	  return (TCL_ERROR);
	}
	/* set charge */
	if (Tcl_GetDouble(interp, argv[1], &q) == TCL_ERROR)
	  return (TCL_ERROR);
	
	mpi_send_q(node, part_num, q);

	argc -= 2;
	argv += 2;
      }
      else if (!strncmp(argv[0], "v", strlen(argv[0]))) {
	double v[3];
	if (argc < 4) {
	  Tcl_AppendResult(interp, "v requires 3 arguments", (char *) NULL);
	  return (TCL_ERROR);
	}
	/* set v */
	if (Tcl_GetDouble(interp, argv[1], &v[0]) == TCL_ERROR)
	  return (TCL_ERROR);
	if (Tcl_GetDouble(interp, argv[2], &v[1]) == TCL_ERROR)
	  return (TCL_ERROR);
	if (Tcl_GetDouble(interp, argv[3], &v[2]) == TCL_ERROR)
	  return (TCL_ERROR);

	mpi_send_v(node, part_num, v);

	argc -= 4;
	argv += 4;
      }
      else if (!strncmp(argv[0], "f", strlen(argv[0]))) {
	double f[3];
	if (argc < 4) {
	  Tcl_AppendResult(interp, "f requires 3 arguments", (char *) NULL);
	  return (TCL_ERROR);
	}
	/* set v */
	if (Tcl_GetDouble(interp, argv[1], &f[0]) == TCL_ERROR)
	  return (TCL_ERROR);
	if (Tcl_GetDouble(interp, argv[2], &f[1]) == TCL_ERROR)
	  return (TCL_ERROR);
	if (Tcl_GetDouble(interp, argv[3], &f[2]) == TCL_ERROR)
	  return (TCL_ERROR);

	mpi_send_v(node, part_num, f);

	argc -= 4;
	argv += 4;
      }
      else if (!strncmp(argv[0], "type", strlen(argv[0]))) {
	int type;
	if (argc < 2) {
	  Tcl_AppendResult(interp, "type requires 1 argument", (char *) NULL);
	  return (TCL_ERROR);
	}
	/* set type */
	if (Tcl_GetInt(interp, argv[1], &type) == TCL_ERROR)
	  return (TCL_ERROR);

	if (type < 0) {
	  Tcl_AppendResult(interp, "invalid particle type", (char *) NULL);
	  return (TCL_ERROR);	  
	} 

	/* make sure type exists */
	realloc_ia_params(type);

	mpi_send_type(node, part_num, type);

	argc -= 2;
	argv += 2;
      }
      else {
	Tcl_AppendResult(interp, "unknown particle parameter \"", argv[0],"\"", (char *)NULL);
	return (TCL_ERROR);
      }
    }
  }

  return (TCL_OK);
}

int inter(ClientData _data, Tcl_Interp *interp,
	  int argc, char **argv)
{
  int i, j;
  IA_parameters *data, *data_sym;

  if (argc < 3) {
    Tcl_AppendResult(interp, "wrong # args:  should be \"",
		     argv[0], " <type 1> <type 2> ?interaction? ?values?\"",
		     (char *) NULL);
    return (TCL_ERROR);
  }

  if ((Tcl_GetInt(interp, argv[1], &i) == TCL_ERROR) ||
      (Tcl_GetInt(interp, argv[2], &j) == TCL_ERROR))
    return (TCL_ERROR);

  data     = safe_get_ia_param(i, j);
  data_sym = safe_get_ia_param(j, i);

  if (!data || !data_sym) {
    Tcl_AppendResult(interp, "particle types must be nonnegative",
		     (char *) NULL);
    return (TCL_ERROR);
  }

  if (argc == 3) {
    /* print interaction information */
    char buffer[TCL_DOUBLE_SPACE];
    Tcl_PrintDouble(interp, data->LJ_eps, buffer);
    Tcl_AppendResult(interp, "{lennard-jones ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->LJ_sig, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->LJ_cut, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->LJ_shift, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->LJ_offset, buffer);
    Tcl_AppendResult(interp, buffer, "}", (char *) NULL);
    return (TCL_OK);
  }

  /* set interaction parameters */
  argc -= 3;
  argv += 3;

  while (argc > 0) {
    if (!strncmp(argv[0], "lennard-jones", strlen(argv[0]))) {
      if (argc < 6) {
	Tcl_AppendResult(interp, "lennard-jones needs 4 parameters: "
			 "<lj_eps> <lj_sig> <lj_cut> <lj_shift> <lj_offset>",
			 (char *) NULL);
	return (TCL_ERROR);
      }
      
      if ((Tcl_GetDouble(interp, argv[1], &data->LJ_eps) == TCL_ERROR) ||
	  (Tcl_GetDouble(interp, argv[2], &data->LJ_sig)  == TCL_ERROR) ||
	  (Tcl_GetDouble(interp, argv[3], &data->LJ_cut)  == TCL_ERROR) ||
	  (Tcl_GetDouble(interp, argv[4], &data->LJ_shift)   == TCL_ERROR) ||
	  (Tcl_GetDouble(interp, argv[5], &data->LJ_offset)  == TCL_ERROR))
	return (TCL_ERROR);

      /* LJ should be symmetrically */
      data_sym->LJ_eps = data->LJ_eps;
      data_sym->LJ_sig = data->LJ_sig;
      data_sym->LJ_cut = data->LJ_cut;
      data_sym->LJ_shift   = data->LJ_shift;
      data_sym->LJ_offset  = data->LJ_offset;
      argc -= 6;
      argv += 6;

      mpi_bcast_ia_params(i, j);
    }
    else {
      Tcl_AppendResult(interp, "unknown interaction type \"", argv[3],
		       "\"", (char *)NULL);
      return (TCL_ERROR);
    }
  }

  return (TCL_OK);
}

int writemd(ClientData data, Tcl_Interp *interp,
	    int argc, char **argv)
{
  static int end_num = -1;
  char *row;
  int p, node, i;
  struct MDHeader header;
  int tcl_file_mode;
  Tcl_Channel channel;

  if (argc < 3) {
    Tcl_AppendResult(interp, "wrong # args:  should be \"",
		     argv[0], " <file> ?posx|posy|posz|q|vx|vy|vz|fx|fy|fz|type?* ...\"",
		     (char *) NULL);
    return (TCL_ERROR);
  }

  if ((channel = Tcl_GetChannel(interp, argv[1], &tcl_file_mode)) == NULL)
    return (TCL_ERROR);
  if (!(tcl_file_mode & TCL_WRITABLE)) {
    Tcl_AppendResult(interp, "\"", argv[1], "\" not writeable", (char *) NULL);
    return (TCL_ERROR);
  }

  /* tune channel to binary translation, e.g. none */
  Tcl_SetChannelOption(interp, channel, "-translation", "binary");

  /* assemble rows */
  argc -= 2;
  argv += 2;
  row = malloc(sizeof(char)*argc);
  for (i = 0; i < argc; i++) {
    if (!strncmp(*argv, "posx", strlen(*argv))) {
      row[i] = POSX;
    }
    else if (!strncmp(*argv, "posy", strlen(*argv))) {
      row[i] = POSY;
    }
    else if (!strncmp(*argv, "posz", strlen(*argv))) {
      row[i] = POSZ;
    }
    else if (!strncmp(*argv, "q", strlen(*argv))) {
      row[i] = Q;
    }
    else if (!strncmp(*argv, "vx", strlen(*argv))) {
      row[i] = VX;
    }
    else if (!strncmp(*argv, "vy", strlen(*argv))) {
      row[i] = VY;
    }
    else if (!strncmp(*argv, "vz", strlen(*argv))) {
      row[i] = VZ;
    }
    else if (!strncmp(*argv, "fx", strlen(*argv))) {
      row[i] = FX;
    }
    else if (!strncmp(*argv, "fy", strlen(*argv))) {
      row[i] = FY;
    }
    else if (!strncmp(*argv, "fz", strlen(*argv))) {
      row[i] = FZ;
    }
    else if (!strncmp(*argv, "type", strlen(*argv))) {
      row[i] = TYPE;
    }
    else {
      Tcl_AppendResult(interp, "no particle data field \"", *argv, "\"?",
		       (char *) NULL);
      return (TCL_ERROR);
    }
    argv++;
  }

  if (!particle_node)
    build_particle_node();

  /* write header and row data */
  memcpy(header.magic, MDMAGIC, 4*sizeof(char));
  header.n_rows = argc;
  Tcl_Write(channel, (char *)&header, sizeof(header));
  Tcl_Write(channel, row, header.n_rows*sizeof(char));

  for (p = 0; p < n_total_particles; p++) {
    node = particle_node[p];
    if (node != -1) {
      Particle data;
      /* fetch particle data */
      mpi_recv_part(node, p, &data);
      for (i = 0; i < 3; i++)
	data.p[i] += data.i[i]*box_l[i];

      /* write particle index */
      Tcl_Write(channel, (char *)&p, sizeof(int));

      for (i = 0; i < header.n_rows; i++) {
	switch (row[i]) {
	case POSX: Tcl_Write(channel, (char *)&data.p[0], sizeof(double)); break;
	case POSY: Tcl_Write(channel, (char *)&data.p[1], sizeof(double)); break;
	case POSZ: Tcl_Write(channel, (char *)&data.p[2], sizeof(double)); break;
	case VX:   Tcl_Write(channel, (char *)&data.v[0], sizeof(double)); break;
	case VY:   Tcl_Write(channel, (char *)&data.v[1], sizeof(double)); break;
	case VZ:   Tcl_Write(channel, (char *)&data.v[2], sizeof(double)); break;
	case FX:   Tcl_Write(channel, (char *)&data.f[0], sizeof(double)); break;
	case FY:   Tcl_Write(channel, (char *)&data.f[1], sizeof(double)); break;
	case FZ:   Tcl_Write(channel, (char *)&data.f[2], sizeof(double)); break;
	case Q:    Tcl_Write(channel, (char *)&data.q, sizeof(double)); break;
	case TYPE: Tcl_Write(channel, (char *)&data.type, sizeof(int)); break;
	}
      }
    }
  }
  /* end marker */
  Tcl_Write(channel, (char *)&end_num, sizeof(int));
  return TCL_OK;
}

int readmd(ClientData dummy, Tcl_Interp *interp,
	   int argc, char **argv)
{
  char *row;
  int pos_row[3] = { -1 }, v_row[3] = { -1 }, f_row[3] = { -1 };
  int av_pos = 0, av_v = 0, av_f = 0, av_q = 0, av_type = 0;
  int node, i;
  struct MDHeader header;
  Particle data;
  int tcl_file_mode;
  Tcl_Channel channel;

  if (argc != 2) {
    Tcl_AppendResult(interp, "wrong # args:  should be \"",
		     argv[0], " <file>\"",
		     (char *) NULL);
    return (TCL_ERROR);
  }

  if ((channel = Tcl_GetChannel(interp, argv[1], &tcl_file_mode)) == NULL)
    return (TCL_ERROR);

  /* tune channel to binary translation, e.g. none */
  Tcl_SetChannelOption(interp, channel, "-translation", "binary");

  Tcl_Read(channel, (char *)&header, sizeof(header));
  /* check token */
  if (strncmp(header.magic, MDMAGIC, 4) || header.n_rows < 0) {
    Tcl_AppendResult(interp, "data file \"", argv[1],
		     "\" does not contain tcl MD data",
		     (char *) NULL);
    return (TCL_ERROR);
  }

  if (!particle_node)
    build_particle_node();

  if (!processor_grid_is_set())
    setup_processor_grid();

  /* parse rows */
  row = malloc(header.n_rows*sizeof(char));
  for (i = 0; i < header.n_rows; i++) {
    Tcl_Read(channel, (char *)&row[i], sizeof(char));
    switch (row[i]) {
    case POSX: pos_row[0] = i; break;
    case POSY: pos_row[1] = i; break;
    case POSZ: pos_row[2] = i; break;
    case   VX:   v_row[0] = i; break;
    case   VY:   v_row[1] = i; break;
    case   VZ:   v_row[2] = i; break;
    case   FX:   f_row[0] = i; break;
    case   FY:   f_row[1] = i; break;
    case   FZ:   f_row[2] = i; break;
    case    Q:   av_q = 1; break;
    case TYPE:   av_type = 1; break;
    }
  }

  /* *_row[0] tells if * data is completely available -
   * otherwise we ignore it */
  if (pos_row[0] != -1 && pos_row[1] != -1 && pos_row[2] != -1) {
    av_pos = 1;
  }
  if (v_row[0] != -1 && v_row[1] != -1 && v_row[2] != -1) {
    av_v = 1;
  }
  if (f_row[0] != -1 && f_row[1] != -1 && f_row[2] != -1) {
    av_f = 1;
  }

  while (!Tcl_Eof(channel)) {
    Tcl_Read(channel, (char *)&data.identity, sizeof(int));
    if (data.identity == -1)
      break;

    /* printf("id=%d\n", data.identity); */

    if (data.identity < 0) {
      Tcl_AppendResult(interp, "illegal data format in data file \"", argv[1],
		       "\", perhaps wrong file?",
		       (char *) NULL);
      return (TCL_ERROR);
    }

    for (i = 0; i < header.n_rows; i++) {
      switch (row[i]) {
      case POSX: Tcl_Read(channel, (char *)&data.p[0], sizeof(double)); break;
      case POSY: Tcl_Read(channel, (char *)&data.p[1], sizeof(double)); break;
      case POSZ: Tcl_Read(channel, (char *)&data.p[2], sizeof(double)); break;
      case   VX: Tcl_Read(channel, (char *)&data.v[0], sizeof(double)); break;
      case   VY: Tcl_Read(channel, (char *)&data.v[1], sizeof(double)); break;
      case   VZ: Tcl_Read(channel, (char *)&data.v[2], sizeof(double)); break;
      case   FX: Tcl_Read(channel, (char *)&data.f[0], sizeof(double)); break;
      case   FY: Tcl_Read(channel, (char *)&data.f[1], sizeof(double)); break;
      case   FZ: Tcl_Read(channel, (char *)&data.f[2], sizeof(double)); break;
      case    Q: Tcl_Read(channel, (char *)&data.q, sizeof(double)); break;
      case TYPE: Tcl_Read(channel, (char *)&data.type, sizeof(int)); break;
      }
    }

    node = (data.identity < n_total_particles) ? particle_node[data.identity] : -1; 
    if (node == -1) {
      if (!av_pos) {
	Tcl_AppendResult(interp, "new particle without position data",
			 (char *) NULL);
	return (TCL_ERROR);
      }

      node = find_node(data.p);
      mpi_attach_particle(data.identity, node);
      map_particle_node(data.identity, node);
    }

    if (av_pos)
      mpi_send_pos(node, data.identity, data.p);
    if (av_q)
      mpi_send_q(node, data.identity, data.q);
    if (av_v)
      mpi_send_v(node, data.identity, data.v);
    if (av_f)
      mpi_send_f(node, data.identity, data.f);
    if (av_type)
      mpi_send_type(node, data.identity, data.type);
  }
  return TCL_OK;
}

