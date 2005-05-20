// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
#ifndef IMD_H__
#define IMD_H__
/** \file imd.h 
    The interface with VMD. This code just provides a wrapper for the IMD interface functions, which allow to send
    particle positions to VMD. Additionally, VMD can send back a single integer value, called transfer_rate, which
    is accessible both from c and from Tcl. The IMD force feedback is not implemented.

    <b>Responsible:</b>
    <a href="mailto:arnolda@mpip-mainz.mpg.de">Axel</a>
*/

#include <tcl.h>

int imd(ClientData data, Tcl_Interp *interp,
	int argc, char **argv);

extern int transfer_rate;

#endif
