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
/** \file imd.c
    Implementation of \ref imd.h "imd.h".
 */
#include <unistd.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cerrno>
#include <tcl.h>
#include "vmdsock.hpp"
#include "imd.hpp"
#include "particle_data.hpp"
#include "parser.hpp"
#include "grid.hpp"
#include "statistics_molecule.hpp"
#include <inttypes.h>

static void *initsock = 0;
static void *sock = 0;


/***************************** Espresso stuff ************************/

int tclcommand_imd_print_drain_socket(Tcl_Interp *interp)
{
  while (vmdsock_selread(sock,0) > 0)  {
    int32_t length;
    IMDType type = imd_recv_header(sock, &length);
    switch (type) {
    case IMD_MDCOMM:
      /* ignore forces for now */
      Tcl_AppendResult(interp, "IMD force feedback not yet implemented",
		       (char *) NULL);
      return (TCL_ERROR);
      /* Expect the msglength to give number of indicies, and the data
	 message to consist of first the indicies, then the coordinates
	 in xyz1 xyz2... format.
	 int32_t *addindices = new int32_t[length];
	 float *addforces = new float[3*length];

	 if (imd_recv_mdcomm(sock, length, addindices, addforces))
	 throw Error(IOERROR);
	 oversized tmp buffer
	 int32_t *tmpindices = new int32_t[n_atoms + length];
	 float *tmpforces = new float[3*(n_atoms + length)];
	 int32_t tmpatoms = n_atoms;
	 for (int i = 0; i < tmpatoms; i++) {
	 tmpindices[i]      = indices[i];
	 tmpforces[3*i    ] = forces[3*i];
	 tmpforces[3*i + 1] = forces[3*i + 1];
	 tmpforces[3*i + 2] = forces[3*i + 2];
	 }

	 for (int i = 0; i < length; i++) {
	 int32_t index = addindices[i];
	 check if there is a force for this atom already
	 int j;
	 for (j = 0; j < tmpatoms; j++)
	 if (tmpindices[j] == index)
	 break;
	 if (j == tmpatoms) {
	 tmpindices[j] = index;
	 tmpatoms++;
	 }
	 tmpforces[3*j    ] = addforces[3*i    ];
	 tmpforces[3*j + 1] = addforces[3*i + 1];
	 tmpforces[3*j + 2] = addforces[3*i + 2];
	 }
	 if (n_atoms > 0) {
	 delete[] indices;
	 delete[] forces;
	 }
	 n_atoms = tmpatoms;
	 indices = new int32_t[n_atoms];
	 forces  = new float[3*n_atoms];
	 cout << "now forces are" << endl;
	 for (int i = 0; i < n_atoms; i++) {
	 indices[i] = tmpindices[i];
	 forces[3*i    ] = tmpforces[3*i];
	 forces[3*i + 1] = tmpforces[3*i + 1];
	 forces[3*i + 2] = tmpforces[3*i + 2];
	 
	 cout << indices[i] << " force "
	 << forces[3*i    ] << " "
	 << forces[3*i + 1] << " "
	 << forces[3*i + 2] << endl;
	 
	 }
	 cout << "end" << endl;
	 */
      break;
    case IMD_TRATE:
      transfer_rate = length;
      break;
    case IMD_IOERROR:
      vmdsock_destroy(sock);
      sock = 0;
      Tcl_AppendResult(interp, "IMD reports IO error",
		       (char *) NULL);
      return (TCL_ERROR);
    case IMD_DISCONNECT:
    case IMD_KILL:
      vmdsock_destroy(sock);
      sock = 0;
      Tcl_AppendResult(interp, "no connection",
		       (char *) NULL);
      return (TCL_OK);
    case IMD_ENERGIES:
    case IMD_FCOORDS:
      Tcl_AppendResult(interp, "IMD protocol failure: unexpected token.",
		       (char *) NULL);
      return (TCL_ERROR);
      break;
    default: ;
    }
  }
  return (TCL_OK);
}

int tclcommand_imd_print_check_connect(Tcl_Interp *interp)
{
  /* handshaking */
  if (vmdsock_selread(initsock, 0) > 0) {
    int32_t length;

    sock = vmdsock_accept(initsock);
    if (imd_handshake(sock)) {
      Tcl_AppendResult(interp, "IMD handshake failed. Wrong VMD version ?",
		       (char *) NULL);
      vmdsock_destroy(sock);
      sock = 0;
      return (TCL_ERROR);
    }

    sleep(1);
    if ((vmdsock_selread(sock, 0) != 1) ||
	(imd_recv_header(sock, &length) != IMD_GO)) {
      Tcl_AppendResult(interp, "No go from VMD. Wrong VMD version ?",
		       (char *) NULL);
      vmdsock_destroy(sock);
      sock = 0;
      return (TCL_ERROR);
    }

    sleep(1);
  }

  return (TCL_OK);
}

int tclcommand_imd_parse_pos(Tcl_Interp *interp, int argc, char **argv)
{
  enum flag {NONE, UNFOLDED, FOLD_CHAINS};
  double shift[3] = {0.0,0.0,0.0};
  //double part_selected=n_total_particles;

  float *coord;
  int flag = NONE;
  int i, j;

  // Determine how many arguments we have and set the value of flag
  switch (argc)
    {
    case 2:
      flag = NONE; 
      break;
    case 3:
      {
	if (ARG_IS_S(2,"-unfolded"))
	  {flag = UNFOLDED;}
	else if (ARG_IS_S(2,"-fold_chains"))
	  {flag = FOLD_CHAINS;}
	else{
	  Tcl_AppendResult(interp, "wrong flag to",argv[0],
			   " positions: should be \" -fold_chains or -unfolded \"",
			   (char *) NULL);
	  return (TCL_ERROR);
	}
      }
      break;
    default:
      Tcl_AppendResult(interp, "wrong # args:  should be \"",
		       argv[0], " positions [-flag]\"",
		       (char *) NULL);
      return (TCL_ERROR);
    }

  if (!initsock) {
    Tcl_AppendResult(interp, "no connection",
		     (char *) NULL);
    return (TCL_OK);
  }
  if (!sock) {
    if (tclcommand_imd_print_check_connect(interp) == TCL_ERROR)
      return (TCL_ERROR);
    
    /* no VMD is ok, but tell the user */
    if (!sock) {
      Tcl_AppendResult(interp, "no connection",
		       (char *) NULL);
      return (TCL_OK);
    }
  }

  if (tclcommand_imd_print_drain_socket(interp) == TCL_ERROR)
    return (TCL_ERROR);
  
  /* we do not consider a non connected VMD as error, but tell the user */
  if (!sock) {
    Tcl_AppendResult(interp, "no connection",
		     (char *) NULL);
    return (TCL_OK);
  }

  if (!(vmdsock_selwrite(sock, 60) > 0)) {
    Tcl_AppendResult(interp, "could not write to IMD socket.",
		     (char *) NULL);
    return (TCL_ERROR);
  }

  if (n_total_particles != max_seen_particle + 1) {
    Tcl_AppendResult(interp, "for IMD, store particles consecutively starting with 0.",
		     (char *) NULL);
    return (TCL_ERROR);      
  }

  updatePartCfg(WITH_BONDS);
  coord = (float*)malloc(n_total_particles*3*sizeof(float));
  /* sort partcles according to identities */
  for (i = 0; i < n_total_particles; i++) {
    int dummy[3] = {0,0,0};
    double tmpCoord[3];
    tmpCoord[0] = partCfg[i].r.p[0];
    tmpCoord[1] = partCfg[i].r.p[1];
    tmpCoord[2] = partCfg[i].r.p[2];
    if (flag == NONE)  {   // perform folding by particle
      fold_position(tmpCoord, dummy);
    }
    j = 3*partCfg[i].p.identity;
    coord[j    ] = tmpCoord[0];
    coord[j + 1] = tmpCoord[1];
    coord[j + 2] = tmpCoord[2];
  }


  // Use information from the analyse set command to fold chain molecules
  if ( flag == FOLD_CHAINS ){
    if(analyze_fold_molecules(coord, shift ) != TCL_OK){
      Tcl_AppendResult(interp, "could not fold chains: \"analyze set chains <chain_start> <n_chains> <chain_length>\" must be used first",
		       (char *) NULL);
      return (TCL_ERROR);   
    }
  }

 
  if (imd_send_fcoords(sock, n_total_particles, coord)) {
    Tcl_AppendResult(interp, "could not write to IMD socket.",
		     (char *) NULL);
    return (TCL_ERROR);      
  }
  free(coord);
  
  Tcl_AppendResult(interp, "connected",
		   (char *) NULL);
  return (TCL_OK);
}

int tclcommand_imd(ClientData data, Tcl_Interp *interp,
	int argc, char **argv)
{
  if (argc < 2) {
    Tcl_AppendResult(interp, "wrong # args:  should be \"",
		     argv[0], " connect|disconnect|listen|positions|energies ?values?\"",
		     (char *) NULL);
    return (TCL_ERROR);
  }

  if (ARG1_IS_S("connect")) {
    /* connect to vmd */
    int port = 12346;

    if (argc > 3) {
      Tcl_AppendResult(interp, "wrong # args:  should be \"",
		       argv[0], " connect ?port?\"",
		       (char *) NULL);
      return (TCL_ERROR);
    }
    if (argc == 3)
      if (!ARG_IS_I(2, port))
	return (TCL_ERROR);

    if (sock)
      vmdsock_destroy(sock);
    if (initsock)
      vmdsock_destroy(initsock);
    sock = 0;
    initsock = 0;

    vmdsock_init();
    initsock = vmdsock_create();
    if (vmdsock_bind(initsock, port) != 0) {
      Tcl_AppendResult(interp, "IMD bind failed. Port already in use ?",
		       (char *) NULL);
      vmdsock_destroy(initsock);
      initsock = 0;
      return (TCL_ERROR);
    }

    if (vmdsock_listen(initsock)) {
      Tcl_AppendResult(interp, "IMD listen failed. Port already in use ?",
		       (char *) NULL);
      vmdsock_destroy(initsock);
      initsock = 0;
      return (TCL_ERROR);
    }

    return (TCL_OK);
  }
  if (ARG1_IS_S("disconnect")) {
    if (argc > 2) {
      Tcl_AppendResult(interp, "wrong # args:  should be \"",
		       argv[0], " disconnect\"",
		       (char *) NULL);
      return (TCL_ERROR);
    }

    if (sock)
      vmdsock_destroy(sock);
    if (initsock)
      vmdsock_destroy(initsock);
    sock = 0;
    initsock = 0;

    Tcl_AppendResult(interp, "no connection",
		     (char *) NULL);
    return (TCL_OK);
  }

  if (ARG1_IS_S("listen")) {
    /* wait until vmd connects */
    int cnt = 3600;
    
    if (argc != 3) {
      Tcl_AppendResult(interp, "wrong # args:  should be \"",
		       argv[0], " listen <secs>\"",
		       (char *) NULL);
      return (TCL_ERROR);
    } 
   
    if (!ARG_IS_I(2, cnt))
      return (TCL_ERROR);

    while (initsock && !sock && cnt--) {
      if (tclcommand_imd_print_check_connect(interp) == TCL_ERROR)
	return (TCL_ERROR);
      sleep(1);
    }

    if (!sock)
      Tcl_AppendResult(interp, "no connection",
		       (char *) NULL);
    else {
      if (tclcommand_imd_print_drain_socket(interp) == TCL_ERROR)
	return (TCL_ERROR);
      Tcl_AppendResult(interp, "connected",
		       (char *) NULL);
    }
    return (TCL_OK);
  }

  if (ARG1_IS_S("positions")) 
    return tclcommand_imd_parse_pos(interp, argc, argv);
  
  if (ARG1_IS_S("energies")) {
    Tcl_AppendResult(interp, "Sorry. imd energies not yet implemented",
		     (char *) NULL);
    return (TCL_ERROR);      
  }

  Tcl_AppendResult(interp, "imd: unkown job.",
		   (char *) NULL);
  return (TCL_ERROR);      
}
