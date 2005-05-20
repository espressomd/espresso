// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2005; all rights reserved unless otherwise stated.
#ifndef UWERR_H
#define UWERR_H
/** \file uwerr.h
 *
 *  PLEASE INSERT DOCUMENTATION
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:holm@mpip-mainz.mpg.de">Christian</a>
*/


#include <tcl.h>

/** The C implementation of the tcl function uwerr \ref tcl_uwerr.
*/
int uwerr(ClientData data, Tcl_Interp *interp, int argc, char *argv[]);
#endif
