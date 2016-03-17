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

#include <boost/mpi.hpp>

#include "ParallelScriptObject.hpp"
#include "MpiCallbacks.hpp"

namespace ScriptInterface {

void ParallelScriptObject::broadcast_parameters() {
  Parameters master_parameters = get_parameters();
  
  set_parameters_all_nodes(master_parameters);
}

void ParallelScriptObject::set_parameters_all_nodes(Parameters &parameters) {
  call_slaves(SET_PARAMETERS, 0);

  boost::mpi::communicator comm;
  boost::mpi::broadcast(comm, parameters, 0);
    
  set_parameters(parameters);
}

ParallelScriptObject::~ParallelScriptObject() {  
  if(mpiCallbacks().mpi_comm.rank() == 0) {
    call_slaves(DELETE, 0);
  }
}

void ParallelScriptObject::callback(int par, int) {
  switch(par) {
    case SET_PARAMETERS:
      {
        Parameters param;
        boost::mpi::communicator comm;
        boost::mpi::broadcast(comm, param, 0);
        set_parameters(param);          
      }
      break;
    case DELETE:
      delete this;
      return;
    default:
      break;
  }
}

}
