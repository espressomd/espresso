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
/** \file blockfile_tcl.h
    Contains only the tcl interface for block coded files.

    It is the header file for \ref blockfile_tcl.c "blockfile_tcl.c" and provides the
    tcl command blockfile.
*/
#ifndef BLOCKFILE_TCL_H
#define BLOCKFILE_TCL_H
#include "config.h"

/** Implementation of the Tcl command blockfile. Allows to read and write
    blockfile comfortably from Tcl. */
int tclcommand_blockfile(ClientData data, Tcl_Interp *interp,
	      int argc, char **argv);
#endif
