/*
  Copyright (C) 2010,2011,2012 The ESPResSo project
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

/** \file reaction_tcl.c
 *
*/

#include "reaction.h"
#include "parser.h"
#include "utils.h"
#include "communication.h"
#include "integrate.h"

#ifdef REACTIONS
int tcl_command_reaction_print_usage(Tcl_Interp * interp){
  char buffer[256];
  sprintf(buffer, "Usage: reaction [off | reactant_type <rt> catalyzer_type <ct> product_type <pt> range <r> rate <k> [back_rate]]\n");
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return TCL_ERROR;
}

int tcl_command_reaction_print(Tcl_Interp * interp){
  char buffer[512];
  sprintf(buffer, "{reactant_type %d} {catalyzer_type %d} {product_type %d} {range %f} {rate %f} {back_rate %f}\n",
			reaction.reactant_type,
			reaction.catalyzer_type,
			reaction.product_type,
			reaction.range,
			reaction.rate,
			reaction.back_rate);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return TCL_OK;
}
#endif /* ifdef REACTIONS */

#define ARG_IS_S_EXACT(no, str) !strcmp(argv[(no)], (str))
int tclcommand_reaction(ClientData data, Tcl_Interp * interp, int argc, char ** argv){
#ifdef REACTIONS
  if (argc == 1  ) return tcl_command_reaction_print_usage(interp);
  if (argc == 2 ) { 
     if (ARG1_IS_S("off")) {
           reaction.rate=0.0;
           mpi_bcast_event(REACTION); 
           return TCL_OK;
     }
     if (ARG1_IS_S("print")) {
           return tcl_command_reaction_print(interp);
     }
  }
  if( argc!=11 && argc!=13) 
     return tcl_command_reaction_print_usage(interp);
     
  if(reaction.rate != 0.0) {
    Tcl_AppendResult(interp, "Currently a simulation can only contain a single reaction!", (char *) NULL);
    return (TCL_ERROR);
  }
     
  if(time_step < 0.0) {
    Tcl_AppendResult(interp, "Time step needs to be set before setting up a reaction!", (char *) NULL);
    return (TCL_ERROR);
  }

  argc--;
  argv++;
  while (argc>0){
      if (ARG_IS_S(0,"product_type")) {
          if (!ARG_IS_I(1,reaction.product_type)) 
            return tcl_command_reaction_print_usage(interp);
          argc-=2;
	  argv+=2;
      } else 
      if (ARG_IS_S(0,"reactant_type")) {
          if (!ARG_IS_I(1,reaction.reactant_type)) 
            return tcl_command_reaction_print_usage(interp);
          argc-=2;
	  argv+=2;
      } else 
      if (ARG_IS_S(0,"catalyzer_type")) {
          if (!ARG_IS_I(1,reaction.catalyzer_type)) 
            return tcl_command_reaction_print_usage(interp);
          argc-=2;
	  argv+=2;
      } else 
      if (ARG_IS_S_EXACT(0,"range")) {
          if (!ARG_IS_D(1,reaction.range)) 
            return tcl_command_reaction_print_usage(interp);
          argc-=2;
	  argv+=2;
      } else
      if (ARG_IS_S_EXACT(0,"rate")) {
          if (!ARG_IS_D(1,reaction.rate)) 
            return tcl_command_reaction_print_usage(interp);
          argc-=2;
	  argv+=2;
      } else
      if (ARG_IS_S_EXACT(0,"back_rate")) {
          if (!ARG_IS_D(1,reaction.back_rate)) 
            return tcl_command_reaction_print_usage(interp);
          argc-=2;
	  argv+=2;
      } else {
            return tcl_command_reaction_print_usage(interp);
      }
   }
   mpi_bcast_event(REACTION);
   return TCL_OK;
#else /* ifdef REACTIONS */
  Tcl_AppendResult(interp, "REACTIONS not compiled in!" ,(char *) NULL);
  return (TCL_ERROR);
#endif /* ifdef REACTIONS */
}
