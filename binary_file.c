// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
/** \file binary_file.c
    Implementation of \ref binary_file.h "binary_file.h".
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "binary_file.h"
#include "global.h"
#include "communication.h"
#include "grid.h"
#include "interaction_data.h"
#include "debug.h"

/* cwz-build-comman: ssh chakotay "builtin cd /nhomes/janeway/axel/progs/Espresso; make" 
   cwz-build-command: make
*/

int writemd(ClientData data, Tcl_Interp *interp,
	    int argc, char **argv)
{
  static int end_num = -1;
  char *row;
  int p, i;
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
      free(row);
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

  for (p = 0; p <= max_seen_particle; p++) {
    Particle data;
    if (get_particle_data(p, &data)) {
      unfold_position(data.r.p, data.l.i);

      /* write particle index */
      Tcl_Write(channel, (char *)&p, sizeof(int));

      for (i = 0; i < header.n_rows; i++) {
	switch (row[i]) {
	case POSX: Tcl_Write(channel, (char *)&data.r.p[0], sizeof(double)); break;
	case POSY: Tcl_Write(channel, (char *)&data.r.p[1], sizeof(double)); break;
	case POSZ: Tcl_Write(channel, (char *)&data.r.p[2], sizeof(double)); break;
	case VX:   Tcl_Write(channel, (char *)&data.m.v[0], sizeof(double)); break;
	case VY:   Tcl_Write(channel, (char *)&data.m.v[1], sizeof(double)); break;
	case VZ:   Tcl_Write(channel, (char *)&data.m.v[2], sizeof(double)); break;
	case FX:   Tcl_Write(channel, (char *)&data.f.f[0], sizeof(double)); break;
	case FY:   Tcl_Write(channel, (char *)&data.f.f[1], sizeof(double)); break;
	case FZ:   Tcl_Write(channel, (char *)&data.f.f[2], sizeof(double)); break;
#ifdef ELECTROSTATICS
	case Q:    Tcl_Write(channel, (char *)&data.p.q, sizeof(double)); break;
#endif
	case TYPE: Tcl_Write(channel, (char *)&data.p.type, sizeof(int)); break;
	}
      }
      free_particle(&data);
    }
  }
  /* end marker */
  Tcl_Write(channel, (char *)&end_num, sizeof(int));
  free(row);
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
    Tcl_Read(channel, (char *)&data.p.identity, sizeof(int));
    if (data.p.identity == -1)
      break;

    /* printf("id=%d\n", data.identity); */

    if (data.p.identity < 0) {
      Tcl_AppendResult(interp, "illegal data format in data file \"", argv[1],
		       "\", perhaps wrong file?",
		       (char *) NULL);
      free(row);
      return (TCL_ERROR);
    }

    for (i = 0; i < header.n_rows; i++) {
      switch (row[i]) {
      case POSX: Tcl_Read(channel, (char *)&data.r.p[0], sizeof(double)); break;
      case POSY: Tcl_Read(channel, (char *)&data.r.p[1], sizeof(double)); break;
      case POSZ: Tcl_Read(channel, (char *)&data.r.p[2], sizeof(double)); break;
      case   VX: Tcl_Read(channel, (char *)&data.m.v[0], sizeof(double)); break;
      case   VY: Tcl_Read(channel, (char *)&data.m.v[1], sizeof(double)); break;
      case   VZ: Tcl_Read(channel, (char *)&data.m.v[2], sizeof(double)); break;
      case   FX: Tcl_Read(channel, (char *)&data.f.f[0], sizeof(double)); break;
      case   FY: Tcl_Read(channel, (char *)&data.f.f[1], sizeof(double)); break;
      case   FZ: Tcl_Read(channel, (char *)&data.f.f[2], sizeof(double)); break;
#ifdef ELECTROSTATICS
      case    Q: Tcl_Read(channel, (char *)&data.p.q, sizeof(double)); break;
#endif
      case TYPE: Tcl_Read(channel, (char *)&data.p.type, sizeof(int)); break;
      }
    }

    node = (data.p.identity <= max_seen_particle) ? particle_node[data.p.identity] : -1; 
    if (node == -1) {
      if (!av_pos) {
	Tcl_AppendResult(interp, "new particle without position data",
			 (char *) NULL);
	free(row);
	return (TCL_ERROR);
      }
    }

    if (av_pos)
      place_particle(data.p.identity, data.r.p);
#ifdef ELECTROSTATICS
    if (av_q)
      set_particle_q(data.p.identity, data.p.q);
#endif
    if (av_v)
      set_particle_v(data.p.identity, data.m.v);
    if (av_f)
      set_particle_f(data.p.identity, data.f.f);
    if (av_type)
      set_particle_type(data.p.identity, data.p.type);
  }

  free(row);
  return TCL_OK;
}

