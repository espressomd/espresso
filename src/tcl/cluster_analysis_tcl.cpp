/*
  Copyright (C) 2011,2012,2013,2014,2015,2016 The ESPResSo project
  
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
#include "particle_data.hpp"
#include "interaction_data.hpp"
#include "cluster_analysis.hpp"
#include "parser.hpp"
#include <sstream>

int tclcommand_cluster_analysis(ClientData data, Tcl_Interp *interp, int argc, char **argv) 
{
  // If no argumens are given, print status
  if (argc==1) {
      if (cluster_analysis().get_criterion()) {
        Tcl_AppendResult(interp, cluster_analysis().get_criterion()->name().c_str(), (char*) NULL);
      } else {
        Tcl_AppendResult(interp, "off", (char*) NULL);
      }        
      return TCL_OK;
    }

  argc--; argv++;

  // Otherwise, we set parameters
  if (ARG0_IS_S("off")) {
    cluster_analysis().set_criterion(NULL);
    return TCL_OK;
  }
  if (ARG0_IS_S("distance")) {
      if (argc != 2) {
      	Tcl_AppendResult(interp, "The distnace criterion needs a distance as argument.", (char*) NULL);
      	return TCL_ERROR;
      }
      double d;
      if (!ARG_IS_D(1,d)) {
        	Tcl_AppendResult(interp, "Need a distance as 1st arg.", (char*) NULL);
        	return TCL_ERROR;
      }
      cluster_analysis().set_criterion(new DistanceCriterion(d));
      argc -= 2; argv += 2;
    }
    else if (ARG0_IS_S("analyze_pair")) {
      cluster_analysis().analyze_pair();
      argc -= 1; argv += 1;
      return TCL_OK;
    }
    else if (ARG0_IS_S("print")) {
      std::stringstream res;
      for (auto it : cluster_analysis().clusters) {
        res << "{ "<<it.first<<" {";
        Cluster cluster = it.second;
        for (int pid : cluster.particles) {
          res << " "<<pid<<" ";
        }
        res << "} } ";
      }
      argc -= 1; argv += 1;
    	Tcl_AppendResult(interp, res.str().c_str(), (char*) NULL);
      return TCL_OK;
    }
    else {
    	Tcl_AppendResult(interp, "Unknown argument.", (char*) NULL);
	    return TCL_ERROR;
      }
}
