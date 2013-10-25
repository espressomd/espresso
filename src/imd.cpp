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
/** \file imd.cpp
    Implementation of \ref imd.hpp "imd.h".
 */
#include <climits>
#include <unistd.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cerrno>
#include <inttypes.h>
#include "utils.hpp"
#include "imd.hpp"
#include "communication.hpp"
#include "particle_data.hpp"
#include "statistics_chain.hpp"
#include "statistics_molecule.hpp"
#include "vmdsock.hpp"

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
int   imd_trate(void *, int32_t);

/* Send data */
int   imd_send_mdcomm(void *, int32_t, const int32_t *, const float *);
int   imd_send_energies(void *, const IMDEnergies *);

/* Receive header and data */
int imd_recv_handshake(void *);

int imd_recv_mdcomm(void *, int32_t, int32_t *, float *);
int imd_recv_energies(void *, IMDEnergies *);
int imd_recv_fcoords(void *, int32_t, float *);


typedef struct {
  int32_t type;
  int32_t length;
} IMDheader;

#define HEADERSIZE 8
#define IMDVERSION 2

int transfer_rate = 0;

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
  char *buf = (char*)malloc(sizeof(char)*size);
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
  char *buf = (char*)malloc(sizeof(char)*size);
  int rc;

  fill_header((IMDheader *)buf, IMD_ENERGIES, 1);
  memcpy((void *)(buf+HEADERSIZE), (const void *)energies, sizeof(IMDEnergies));
  rc = (imd_writen(s, buf, size) != size);
  free(buf);
  return rc;
}

int imd_send_fcoords(void *s, int32_t n, const float *coords) {
  int32_t size = HEADERSIZE+12*n;
  char *buf = (char*)malloc(sizeof(char)*size);
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

