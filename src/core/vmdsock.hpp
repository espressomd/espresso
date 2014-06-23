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
/** \file vmdsock.hpp

    DO NOT CHANGE !!!
***
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#if defined(VMDSOCKINTERNAL)

#ifndef WIN32
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/time.h>
#include <netinet/in.h>
#include <sys/file.h>
#endif

typedef struct {
  struct sockaddr_in addr; /* address of socket provided by bind() */
  int addrlen;             /* size of the addr struct */
  int sd;                  /* socket file descriptor */
} vmdsocket;

#endif /* VMDSOCKINTERNAL */

int   vmdsock_init(void);
void *vmdsock_create(void);
int   vmdsock_bind(void *, int);
int   vmdsock_listen(void *);
void *vmdsock_accept(void *);  /* return new socket */
int   vmdsock_connect(void *, const char *, int);
int   vmdsock_write(void *, const void *, int);
int   vmdsock_read(void *, void *, int);
int   vmdsock_selread(void *, int);
int   vmdsock_selwrite(void *, int);
void  vmdsock_destroy(void *);

