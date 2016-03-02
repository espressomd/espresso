/*
  Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
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
/** \file errorhandling.cpp
    Implementation of \ref errorhandling.hpp.
*/
#include <cstring>
#include <cstdlib>
#include <csignal>
#include <iostream>
#include "utils.hpp"
#include "errorhandling.hpp"
#include "communication.hpp"
#include "RuntimeErrorCollector.hpp"

using namespace std;

static void sigint_handler(int sig) {
  /* without this exit handler the nodes might exit asynchronously
   * without calling MPI_Finalize, which may cause MPI to hang
   * (e.g. this handler makes CTRL-c work properly with poe)
   *
   * NOTE: mpirun installs its own handler for SIGINT/SIGTERM
   *       and takes care of proper cleanup and exit */


  static int numcalls = 0;
  if (numcalls++ > 0) exit(sig); // catch sig only once

  /* we use runtime_error to indicate that sig was called;
   * upon next call of mpi_gather_runtime_errors all nodes 
   * will clean up and exit. */
  ostringstream msg;
  msg <<"caught signal "<< sig ;
  runtimeError(msg);
}

void register_sigint_handler() {
  signal(SIGINT, sigint_handler);
}

/* NEW RUNTIME ERROR HANDLING. */

/** list that contains the runtime error messages */
RuntimeErrorCollector *runtimeErrorCollector = NULL;

void
initRuntimeErrorCollector() {
  runtimeErrorCollector = new RuntimeErrorCollector(MPI_COMM_WORLD);
}

void _runtimeWarning(const std::string &msg, 
                     const char* function, const char* file, const int line) {
  runtimeErrorCollector->warning(msg, function, file, line);
}

void _runtimeWarning(const char* msg, 
                     const char* function, const char* file, const int line) {
  runtimeErrorCollector->warning(msg, function, file, line);
}

void _runtimeWarning(const std::ostringstream &msg, 
                     const char* function, const char* file, const int line) {
  runtimeErrorCollector->warning(msg, function, file, line);
}

void _runtimeError(const std::string &msg, 
                     const char* function, const char* file, const int line) {
  runtimeErrorCollector->error(msg, function, file, line);
}

void _runtimeError(const char* msg, 
                     const char* function, const char* file, const int line) {
  runtimeErrorCollector->error(msg, function, file, line);
}

void _runtimeError(const std::ostringstream &msg, 
                     const char* function, const char* file, const int line) {
  runtimeErrorCollector->error(msg, function, file, line);
}

ErrorHandling::RuntimeErrorStream _runtimeErrorStream(const std::string &file, const int line, const std::string &function) {
  return ErrorHandling::RuntimeErrorStream(*runtimeErrorCollector, file, line, function);
}

int check_runtime_errors() {
  return runtimeErrorCollector->count();
}

list<string>
mpiRuntimeErrorCollectorGather() {
  // Tell other processors to send their erros
  mpi_call(mpiRuntimeErrorCollectorGatherSlave, -1, 0);
  return runtimeErrorCollector->gather();
}

void
mpiRuntimeErrorCollectorGatherSlave(int node, int parm) {
  runtimeErrorCollector->gatherSlave();
}

void errexit() {
#ifdef FORCE_CORE
  core();
#endif
  mpi_abort();
  exit(1);
}
