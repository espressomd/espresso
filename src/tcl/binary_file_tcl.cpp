/*
  Copyright (C) 2010,2012,2013 The ESPResSo project
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
/** \file binary_file_tcl.cpp
    Implementation of \ref binary_file_tcl.hpp "binary_file_tcl.h".
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "utils.hpp"
#include "binary_file_tcl.hpp"
#include "global.hpp"
#include "communication.hpp"
#include "grid.hpp"
#include "interaction_data.hpp"

int tclcommand_writemd(ClientData data, Tcl_Interp *interp,
	    int argc, char **argv)
{
  static int end_num = -1;
  char *row;
  int p, i;
  struct MDHeader header;
  int tcl_file_mode;
  Tcl_Channel channel;

  if (argc < 3) {
    #if defined(ELECTROSTATICS) && defined(DIPOLES)
      Tcl_AppendResult(interp, "wrong # args:  should be \"",
	  	       argv[0], " <file> ?posx|posy|posz|q|mx|my|mz|vx|vy|vz|fx|fy|fz|type?* ...\"",
		       (char *) NULL);
    #else
      #ifdef ELECTROSTATICS
      Tcl_AppendResult(interp, "wrong # args:  should be \"",
	  	       argv[0], " <file> ?posx|posy|posz|q|vx|vy|vz|fx|fy|fz|type?* ...\"",
		       (char *) NULL);
      #endif
      
      #ifdef DIPOLES
      Tcl_AppendResult(interp, "wrong # args:  should be \"",
	  	       argv[0], " <file> ?posx|posy|posz|mx|my|mz|vx|vy|vz|fx|fy|fz|type?* ...\"",
		       (char *) NULL);
      #endif
    
    #endif		       
   		     
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
  row = (char*)malloc(sizeof(char)*argc);
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
#ifdef MASS
    else if (!strncmp(*argv, "mass", strlen(*argv))) {
      row[i] = MASSES;
    }
#endif
    else if (!strncmp(*argv, "q", strlen(*argv))) {
      row[i] = Q;
    }
#ifdef DIPOLES    
    else if (!strncmp(*argv, "mx", strlen(*argv))) {
      row[i] = MX;
    }
    else if (!strncmp(*argv, "my", strlen(*argv))) {
      row[i] = MY;
    }
    else if (!strncmp(*argv, "mz", strlen(*argv))) {
      row[i] = MZ;
    }    
#endif    
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
    if (get_particle_data(p, &data) == ES_OK) {
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
#ifdef MASS
	case MASSES: Tcl_Write(channel, (char *)&data.p.mass, sizeof(double)); break;
#endif
#ifdef ELECTROSTATICS
	case Q:    Tcl_Write(channel, (char *)&data.p.q, sizeof(double)); break;
#endif
#ifdef DIPOLES
	case MX:   Tcl_Write(channel, (char *)&data.r.dip[0], sizeof(double)); break;
	case MY:   Tcl_Write(channel, (char *)&data.r.dip[1], sizeof(double)); break;
	case MZ:   Tcl_Write(channel, (char *)&data.r.dip[2], sizeof(double)); break;
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

int tclcommand_readmd(ClientData dummy, Tcl_Interp *interp,
	   int argc, char **argv)
{
  char *row;
  int pos_row[3] = { -1 }, v_row[3] = { -1 }, 
  #ifdef DIPOLES 
  dip_row[3] = { -1 }, 
  #endif
  f_row[3] = { -1 };
  
  int av_pos = 0, av_v = 0, 
#ifdef DIPOLES 
    av_dip=0, 
#endif
#ifdef MASS
    av_mass=0,
#endif
#ifdef SHANCHEN
    av_solvation=0,
#endif
    av_f = 0,
#ifdef ELECTROSTATICS
    av_q = 0,
#endif
    av_type = 0;
  
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
  row = (char*)malloc(header.n_rows*sizeof(char));
  for (i = 0; i < header.n_rows; i++) {
    Tcl_Read(channel, (char *)&row[i], sizeof(char));
    switch (row[i]) {
    case POSX: pos_row[0] = i; break;
    case POSY: pos_row[1] = i; break;
    case POSZ: pos_row[2] = i; break;
    case   VX:   v_row[0] = i; break;
    case   VY:   v_row[1] = i; break;
    case   VZ:   v_row[2] = i; break;
#ifdef DIPOLES
    case   MX:   dip_row[0] = i; break;
    case   MY:   dip_row[1] = i; break;
    case   MZ:   dip_row[2] = i; break;
#endif
    case   FX:   f_row[0] = i; break;
    case   FY:   f_row[1] = i; break;
    case   FZ:   f_row[2] = i; break;
#ifdef MASS
    case MASSES: av_mass  = 1; break;
#endif
#ifdef SHANCHEN
    case SOLVATION: av_solvation = 1; break;
#endif
#ifdef ELECTROSTATICS
    case    Q:   av_q     = 1; break;
#endif
    case TYPE:   av_type  = 1; break;
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
  
  #ifdef DIPOLES
  if (dip_row[0] != -1 && dip_row[1] != -1 && dip_row[2] != -1) {
    av_dip = 1;
  }
  #endif


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
      case MASSES:
#ifdef MASS
          Tcl_Read(channel, (char *)&data.p.mass, sizeof(double)); break;
#else
          {
              double dummy_mass;
              Tcl_Read(channel, (char *)&dummy_mass, sizeof(double)); break;
          }
#endif
#ifdef ELECTROSTATICS
      case    Q: Tcl_Read(channel, (char *)&data.p.q, sizeof(double)); break;
#endif
#ifdef DIPOLES
      case   MX: Tcl_Read(channel, (char *)&data.r.dip[0], sizeof(double)); break;
      case   MY: Tcl_Read(channel, (char *)&data.r.dip[1], sizeof(double)); break;
      case   MZ: Tcl_Read(channel, (char *)&data.r.dip[2], sizeof(double)); break;
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
#ifdef MASS
    if (av_mass)
      set_particle_mass(data.p.identity, data.p.mass);
#endif
#ifdef SHANCHEN
    if (av_solvation)
      set_particle_solvation(data.p.identity, data.p.solvation);
#endif
#ifdef ELECTROSTATICS
    if (av_q)
      set_particle_q(data.p.identity, data.p.q);
#endif
#ifdef DIPOLES
    if (av_dip)
      set_particle_dip(data.p.identity, data.r.dip);
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

