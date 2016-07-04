/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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

/** \file reaction_tcl.cpp
 *
*/

#include "reaction.hpp"
#include "parser.hpp"
#include "utils.hpp"
#include "communication.hpp"
#include "integrate.hpp"

#ifdef CATALYTIC_REACTIONS
int tcl_command_reaction_print_usage(Tcl_Interp * interp){
  char buffer[256];
  sprintf(buffer, "Usage: reaction [off | print | reactant_type <rt> catalyzer_type <ct> product_type <pt> range <r> ct_rate <k_ct> [eq_rate <k_eq>] [react_once <on|off>] [swap <on|off>]\n");
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return TCL_ERROR;
}

int tcl_command_reaction_print(Tcl_Interp * interp){
  char buffer[512];
  sprintf(buffer, "{reactant_type %d} {catalyzer_type %d} {product_type %d} {range %f} {ct_rate %f} {eq_rate %f}",
    reaction.reactant_type,
    reaction.catalyzer_type,
    reaction.product_type,
    reaction.range,
    reaction.ct_rate,
    reaction.eq_rate);
  if ( reaction.sing_mult == 0 ) {
    sprintf(buffer + strlen(buffer)," {react_once off}");
  } else {
    sprintf(buffer + strlen(buffer)," {react_once on}");
  }
  if ( reaction.swap ) {
    sprintf(buffer + strlen(buffer)," {swap on}");
  } else {
    sprintf(buffer + strlen(buffer)," {swap off}");
  }
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return TCL_OK;
}
#endif

int tclcommand_reaction(ClientData data, Tcl_Interp * interp, int argc, char ** argv){
#ifdef CATALYTIC_REACTIONS

  /* Determine the currently set types, to block the user from 
     trying to set up multiple reactions */

  int react_type, prdct_type, catal_type;
  react_type = reaction.reactant_type;
  prdct_type = reaction.product_type;
  catal_type = reaction.catalyzer_type;

  if (argc == 1  ) return tcl_command_reaction_print_usage(interp);

  if (argc == 2 ) { 
     if (ARG1_IS_S("off")) {
       /* We only need to set ct_rate to zero and we
          do not enter the reaction integration loop */
       reaction.ct_rate=0.0;
       mpi_setup_reaction();
       return TCL_OK;
     }
     if (ARG1_IS_S("print")) {
           return tcl_command_reaction_print(interp);
     }
  }

  if( argc!=11 && argc!=13 && argc!=15 && argc!=17 && argc!=3) {
     return tcl_command_reaction_print_usage(interp);
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
    } 
    else if (ARG_IS_S(0,"reactant_type")) {
      if (!ARG_IS_I(1,reaction.reactant_type)) 
        return tcl_command_reaction_print_usage(interp);
      argc-=2;
	    argv+=2;
    } 
    else if (ARG_IS_S(0,"catalyzer_type")) {
      if (!ARG_IS_I(1,reaction.catalyzer_type)) 
        return tcl_command_reaction_print_usage(interp);
    argc-=2;
	  argv+=2;
    } 
    else if (ARG_IS_S_EXACT(0,"range")) {
      if (!ARG_IS_D(1,reaction.range)) 
        return tcl_command_reaction_print_usage(interp);
      argc-=2;
      argv+=2;
    }
    else if (ARG_IS_S_EXACT(0,"ct_rate")) {
      if (!ARG_IS_D(1,reaction.ct_rate)) 
        return tcl_command_reaction_print_usage(interp);
      argc-=2;
      argv+=2;
    } 
    else if (ARG_IS_S_EXACT(0,"eq_rate")) {
      if (!ARG_IS_D(1,reaction.eq_rate)) 
        return tcl_command_reaction_print_usage(interp);
      argc-=2;
	    argv+=2;
    } 
    else if (ARG_IS_S_EXACT(0,"react_once")) {
      if (!ARG_IS_S(1,"on")&&!ARG_IS_S(1,"off")) {
        return tcl_command_reaction_print_usage(interp);}
      if (ARG_IS_S(1,"on")) reaction.sing_mult = 1;
      if (ARG_IS_S(1,"off")) reaction.sing_mult = 0;
    argc-=2;
	  argv+=2;
    }
    else if (ARG_IS_S_EXACT(0,"swap")) {
#ifndef ROTATION
      char buffer[80];
      sprintf(buffer, "WARNING: Parameter \"swap\" has no effect when ROTATION is not compiled in.");
      Tcl_AppendResult(interp, buffer, (char *)NULL);
#endif
      if (ARG_IS_S(1,"on"))
        reaction.swap = 1;
      else if (ARG_IS_S(1,"off"))
        reaction.swap = 0;
      else
        return tcl_command_reaction_print_usage(interp);
      argc-=2;
      argv+=2;
    } else {
      return tcl_command_reaction_print_usage(interp);
    }
  }

  if( reaction.ct_rate < 0.0 ) {
    Tcl_AppendResult(interp, "Negative catalytic reaction rate constant is not allowed!", (char *) NULL);
    return (TCL_ERROR);
  }

  if( reaction.eq_rate < 0.0 && fabs(reaction.eq_rate + 1.0) > 0.001 ) {
    Tcl_AppendResult(interp, "Negative equilibrium reaction rate contstant is not allowed!", (char *) NULL);
    return (TCL_ERROR);
  }

  if( (reaction.product_type == reaction.reactant_type) || (reaction.product_type == reaction.catalyzer_type) || (reaction.catalyzer_type == reaction.reactant_type) ) {
    Tcl_AppendResult(interp, "One particle type cannot be a part more than one reaction species!", (char *) NULL);
    return (TCL_ERROR);
  }

  if ( ((react_type != reaction.reactant_type) || (prdct_type != reaction.product_type) || (catal_type != reaction.catalyzer_type)) && (react_type != prdct_type) ) {
    Tcl_AppendResult(interp, "A simulation can only contain a single reaction!", (char *) NULL);
    return (TCL_ERROR); 
  }

  mpi_setup_reaction();
  return TCL_OK;
#else
  Tcl_AppendResult(interp, "CATALYTIC_REACTIONS not compiled in!" ,(char *) NULL);
  return (TCL_ERROR);
#endif
}
