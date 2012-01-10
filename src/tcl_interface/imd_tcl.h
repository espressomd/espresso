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
#ifndef IMD_TCL_H__
#define IMD_TCL_H__
/** \file imd.h 
    The interface with VMD. This code just provides a wrapper for the IMD interface functions, which allow to send
    particle positions to VMD. Additionally, VMD can send back a single integer value, called transfer_rate, which
    is accessible both from c and from Tcl. The IMD force feedback is not implemented.
*/

#include <tcl.h>

int tclcommand_imd(ClientData data, Tcl_Interp *interp,
	int argc, char **argv);

extern int transfer_rate;

#endif
