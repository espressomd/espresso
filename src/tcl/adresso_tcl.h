/*
  Copyright (C) 2010,2012,2013 The ESPResSo project
  Copyright (C) 2008,2009,2010 
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
#ifndef ADRESSO_TCL_H
#define ADRESSO_TCL_H
/** \file adresso.h
    This is the place for adaptive resolution scheme (adress)
    Implementation of adresso.h

    For more details about adress see:
    - M. Praprotnik, L. Delle Site and K. Kremer, JCP 123, 224106, 2005. 
    - M. Praprotnik, L. Delle Site and K. Kremer, Annu. Rev. Phys. Chem. 59, 545-571, 2008. 
    - S. Poblete, M. Praprotnik, K. Kremer and L. Delle Site, J. Chem. Phys. 132, 114101, 2010. 

    For more detail about the implementation here see:
    - C. Junghans and S. Poblete, Comp. Phys. Comm. 181, 1449, 2010.
*/
#include "parser.h"

/** \name Exported Functions */
/************************************************************/
/*@{*/
/** Implements the Tcl command "adress". This allows for seetings for adress
*/
int tclcommand_adress(ClientData data, Tcl_Interp *interp, int argc, char **argv);

int tclcommand_update_adress_weights(ClientData _data, Tcl_Interp * interp, int argc, char ** argv);

#ifdef ADRESS

/* #ifdef THERMODYNAMIC_FORCE */
int tclcommand_thermodynamic_force_parse_opt(Tcl_Interp * interp, int type, double prefactor, int argc, char ** argv);
int tclcommand_thermodynamic_force(ClientData _data, Tcl_Interp * interp, int argc, char ** argv);
/* #endif */

///
int adress_tab_parser(Tcl_Interp * interp,
		      int part_type_a, int part_type_b,
		      int argc, char ** argv);
#endif
/*@}*/
#endif
