/*
  Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
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

#include "rotate_system.hpp"
#include "parser.hpp"

/* ############### */
int tclcommand_rotate_system(ClientData data, Tcl_Interp * interp, int argc, char ** argv){
  double alpha,theta,phi;
  if (argc != 4) { 
    fprintf(stderr,"needs 3 angles\n");
    return ES_ERROR;
  }
  if (! (ARG_IS_D(1,phi)))
  {
      fprintf(stderr,"Expects 3 floats\n");
      return ES_ERROR;
  }

  if (!(ARG_IS_D(2,theta)))
  {
      fprintf(stderr,"Expects 3 floats\n");
      return ES_ERROR;
  }
  if (! (ARG_IS_D(3,alpha)))
  {
      fprintf(stderr,"Expects 3 floats\n");
      return ES_ERROR;
  }
  

  rotate_system(phi,theta,alpha);

  return ES_OK;
}
