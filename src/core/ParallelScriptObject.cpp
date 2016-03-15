#include "ParallelScriptObject.hpp"

#include <iostream>

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
