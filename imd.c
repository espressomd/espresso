/*
  Copyright (C) 2010 The ESPResSo project
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
/** \file imd.c
    Implementation of \ref imd.h "imd.h".
 */
#include <limits.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "utils.h"
#include "vmdsock.h"
#include "imd.h"
#include "communication.h"
#include "particle_data.h"
#include "parser.h"
#include "statistics_chain.h"
#include "statistics_molecule.h"

int transfer_rate = 0;

#include <inttypes.h>

typedef enum {
  IMD_DISCONNECT,
  IMD_ENERGIES, 
  IMD_FCOORDS,   
  IMD_GO,
  IMD_HANDSHAKE, 
  IMD_KILL,      
  IMD_MDCOMM,    
  IMD_PAUSE,
  IMD_TRATE,
  IMD_IOERROR
} IMDType;

typedef struct {
  int32_t tstep;
  float T;
  float Etot;
  float Epot;
  float Evdw;
  float Eelec;
  float Ebond;
  float Eangle;
  float Edihe;
  float Eimpr;
} IMDEnergies;

/* Send simple messages - these consist of a header with no subsequent data */
int   imd_disconnect(void *);
int   imd_pause(void *);
int   imd_kill(void *);
int   imd_handshake(void *);
int   imd_trate(void *, int32_t);

/* Send data */
int   imd_send_mdcomm(void *, int32_t, const int32_t *, const float *);
int   imd_send_energies(void *, const IMDEnergies *);
int   imd_send_fcoords(void *, int32_t, const float *);

/* Receive header and data */
int imd_recv_handshake(void *);

IMDType imd_recv_header(void *, int32_t *);
int imd_recv_mdcomm(void *, int32_t, int32_t *, float *);
int imd_recv_energies(void *, IMDEnergies *);
int imd_recv_fcoords(void *, int32_t, float *);


typedef struct {
  int32_t type;
  int32_t length;
} IMDheader;

#define HEADERSIZE 8
#define IMDVERSION 2

static void *initsock = 0;
static void *sock = 0;

static void swap4(char *data, int ndata) {
  int i;
  char *dataptr;
  char b0, b1;

  dataptr = data;
  for (i=0; i<ndata; i+=4) {
    b0 = dataptr[0];
    b1 = dataptr[1];
    dataptr[0] = dataptr[3];
    dataptr[1] = dataptr[2];
    dataptr[2] = b1;
    dataptr[3] = b0;
    dataptr += 4;
  }
}

static int32_t imd_htonl(int32_t h) {
  int32_t n;
  ((char *)&n)[0] = (h >> 24) & 0x0FF;
  ((char *)&n)[1] = (h >> 16) & 0x0FF;
  ((char *)&n)[2] = (h >> 8) & 0x0FF;
  ((char *)&n)[3] = h & 0x0FF;
  return n;
}

static int32_t imd_ntohl(int32_t n) {
  int32_t h = 0;

  h  = ((char *)&n)[0] << 24;
  h |= ((char *)&n)[1] << 16;
  h |= ((char *)&n)[2] << 8;
  h |= ((char *)&n)[3];

  return h;
}

static void fill_header(IMDheader *header, IMDType type, int32_t length) {
  header->type = imd_htonl((int32_t)type);
  header->length = imd_htonl(length);
}

static void swap_header(IMDheader *header) {
  header->type = imd_ntohl(header->type);
  header->length= imd_ntohl(header->length);
}

static int32_t imd_readn(void *s, char *ptr, int32_t n) {
  int32_t nleft;
  int32_t nread;
 
  nleft = n;
  while (nleft > 0) {
    if ((nread = vmdsock_read(s, ptr, nleft)) < 0) {
      if (errno == EINTR)
        nread = 0;         /* and call read() again */
      else
        return -1;
    } else if (nread == 0)
      break;               /* EOF */
    nleft -= nread;
    ptr += nread;
  }
  return n-nleft;
}

static int32_t imd_writen(void *s, const char *ptr, int32_t n) {
  int32_t nleft;
  int32_t nwritten;

  nleft = n;
  while (nleft > 0) {
    if ((nwritten = vmdsock_write(s, ptr, nleft)) <= 0) {
      if (errno == EINTR)
        nwritten = 0;
      else
        return -1;
    }
    nleft -= nwritten;
    ptr += nwritten;
  }
  return n;
}
 

int imd_disconnect(void *s) {
  IMDheader header;
  fill_header(&header, IMD_DISCONNECT, 0);
  return (imd_writen(s, (char *)&header, HEADERSIZE) != HEADERSIZE);
}

int imd_pause(void *s) {
  IMDheader header;
  fill_header(&header, IMD_PAUSE, 0);
  return (imd_writen(s, (char *)&header, HEADERSIZE) != HEADERSIZE);
}

int imd_kill(void *s) {
  IMDheader header;
  fill_header(&header, IMD_KILL, 0);
  return (imd_writen(s, (char *)&header, HEADERSIZE) != HEADERSIZE);
}

static int imd_go(void *s) {
  IMDheader header;
  fill_header(&header, IMD_GO, 0);
  return (imd_writen(s, (char *)&header, HEADERSIZE) != HEADERSIZE);
}


int imd_handshake(void *s) {
  IMDheader header;
  fill_header(&header, IMD_HANDSHAKE, 1);
  header.length = IMDVERSION;   /* Not byteswapped! */
  return (imd_writen(s, (char *)&header, HEADERSIZE) != HEADERSIZE);
}

int imd_trate(void *s, int32_t rate) {
  IMDheader header;
  fill_header(&header, IMD_TRATE, rate);
  return (imd_writen(s, (char *)&header, HEADERSIZE) != HEADERSIZE);
}

/* Data methods */

int imd_send_mdcomm(void *s,int32_t n,const int32_t *indices,const float *forces) {
  int32_t size = HEADERSIZE+16*n;
  char *buf = malloc(sizeof(char)*size);
  int rc;

  fill_header((IMDheader *)buf, IMD_MDCOMM, n);
  memcpy((void *)(buf+HEADERSIZE), (const void *)indices, 4*n);
  memcpy((void *)(buf+HEADERSIZE+4*n), (const void *)forces, 12*n);
  rc = (imd_writen(s, buf, size) != size);
  free(buf);
  return rc;
}

int imd_send_energies(void *s, const IMDEnergies *energies) {
  int32_t size = HEADERSIZE+sizeof(IMDEnergies);
  char *buf = malloc(sizeof(char)*size);
  int rc;

  fill_header((IMDheader *)buf, IMD_ENERGIES, 1);
  memcpy((void *)(buf+HEADERSIZE), (const void *)energies, sizeof(IMDEnergies));
  rc = (imd_writen(s, buf, size) != size);
  free(buf);
  return rc;
}

int imd_send_fcoords(void *s, int32_t n, const float *coords) {
  int32_t size = HEADERSIZE+12*n;
  char *buf = malloc(sizeof(char)*size);
  int rc;

  fill_header((IMDheader *)buf, IMD_FCOORDS, n);
  memcpy((void *)(buf+HEADERSIZE), (const void *)coords, 12*n);
  rc = (imd_writen(s, buf, size) != size);
  free(buf);
  return rc;
}

/* The IMD receive functions */

IMDType imd_recv_header_nolengthswap(void *s, int32_t *length) {
  IMDheader header;
  if (imd_readn(s, (char *)&header, HEADERSIZE) != HEADERSIZE)
    return IMD_IOERROR;
  *length = header.length;
  swap_header(&header);
  return (IMDType) header.type;
}

int imd_recv_handshake(void *s) {
  int32_t buf;
  IMDType type;

  /* Wait 5 seconds for the handshake to come */
  if (vmdsock_selread(s, 5) != 1) return -1;

  /* Check to see that a valid handshake was received */
  type = imd_recv_header_nolengthswap(s, &buf);
  if (type != IMD_HANDSHAKE) return -1;

  /* Check its endianness, as well as the IMD version. */
  if (buf == IMDVERSION) {
    if (!imd_go(s)) return 0;
    return -1;
  }
  swap4((char *)&buf, 4);
  if (buf == IMDVERSION) {
    if (!imd_go(s)) return 1;
  }
  
  /* We failed to determine endianness. */
  return -1; 
}

IMDType imd_recv_header(void *s, int32_t *length) {
  IMDheader header;
  if (imd_readn(s, (char *)&header, HEADERSIZE) != HEADERSIZE)
    return IMD_IOERROR;
  swap_header(&header);
  *length = header.length;
  return (IMDType)header.type; 
}

int imd_recv_mdcomm(void *s, int32_t n, int32_t *indices, float *forces) {
  if (imd_readn(s, (char *)indices, 4*n) != 4*n) return 1;
  if (imd_readn(s, (char *)forces, 12*n) != 12*n) return 1;
  return 0;
}

int imd_recv_energies(void *s, IMDEnergies *energies) {
  return (imd_readn(s, (char *)energies, sizeof(IMDEnergies))
          != sizeof(IMDEnergies));
}

int imd_recv_fcoords(void *s, int32_t n, float *coords) {
  return (imd_readn(s, (char *)coords, 12*n) != 12*n);
}

/***************************** Espresso stuff ************************/

int imd_drain_socket(Tcl_Interp *interp)
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

int imd_check_connect(Tcl_Interp *interp)
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

int imd_parse_pos(Tcl_Interp *interp, int argc, char **argv)
{
  enum flag {NONE, UNFOLDED, FOLD_CHAINS};
  double shift[3] = {0.0,0.0,0.0};

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
    if (imd_check_connect(interp) == TCL_ERROR)
      return (TCL_ERROR);
    
    /* no VMD is ok, but tell the user */
    if (!sock) {
      Tcl_AppendResult(interp, "no connection",
		       (char *) NULL);
      return (TCL_OK);
    }
  }

  if (imd_drain_socket(interp) == TCL_ERROR)
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
  coord = malloc(n_total_particles*3*sizeof(float));
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

int imd(ClientData data, Tcl_Interp *interp,
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
      if (imd_check_connect(interp) == TCL_ERROR)
	return (TCL_ERROR);
      sleep(1);
    }

    if (!sock)
      Tcl_AppendResult(interp, "no connection",
		       (char *) NULL);
    else {
      if (imd_drain_socket(interp) == TCL_ERROR)
	return (TCL_ERROR);
      Tcl_AppendResult(interp, "connected",
		       (char *) NULL);
    }
    return (TCL_OK);
  }

  if (ARG1_IS_S("positions")) 
    return imd_parse_pos(interp, argc, argv);
  
  if (ARG1_IS_S("energies")) {
    Tcl_AppendResult(interp, "Sorry. imd energies not yet implemented",
		     (char *) NULL);
    return (TCL_ERROR);      
  }

  Tcl_AppendResult(interp, "imd: unkown job.",
		   (char *) NULL);
  return (TCL_ERROR);      
}
