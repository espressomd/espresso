/*
  Copyright (C) 2010,2012,2013,2014 The ESPResSo project
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
    Implementation of \ref errorhandling.hpp "errorhandling.h".
*/
#include <cstring>
#include <cstdlib>
#include <csignal>
#include <iostream>
#include "utils.hpp"
#include "errorhandling.hpp"
#include "communication.hpp"

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

/* New runtime error handling. */
/** list that contains the runtime error messages */
ParallelRuntimeErrors *runtimeErrors = NULL;

enum RuntimeErrorType { ERROR, WARNING };

static const string
createRuntimeErrorMessage(const string &_message,
                          const char* function, const char* file, const int line, 
                          RuntimeErrorType type = ERROR) {
  ostringstream ostr;
  ostr << "{ ";
  switch (type) {
  case ERROR:
    ostr << "ERROR: "; break;
  case WARNING:
    ostr << "WARNING: "; break;
  }
  ostr << _message;
  ostr << " in function " << function << " (" << file << ":" << line 
       << ") on node " << this_node;
  ostr << " } ";
  return ostr.str();
}

ParallelRuntimeErrors::
ParallelRuntimeErrors(MPI_Comm comm) {
  this->comm = comm;
}

void ParallelRuntimeErrors::
warning(const string &msg,
        const char* function, const char* file, const int line) {
  errors.push_back
    (createRuntimeErrorMessage(msg, function, file, line, WARNING));
}

void ParallelRuntimeErrors::
warning(const char *msg,
        const char* function, const char* file, const int line) {
  this->warning(string(msg), function, file, line);
}

void ParallelRuntimeErrors::
warning(ostringstream &mstr,
        const char* function, const char* file, const int line) {
  this->warning(mstr.str(), function, file, line);
}

void ParallelRuntimeErrors::
error(const string &msg,
      const char* function, const char* file, const int line) {
  errors.push_back
    (createRuntimeErrorMessage(msg, function, file, line, ERROR));
}

void ParallelRuntimeErrors::
error(const char *msg,
      const char* function, const char* file, const int line) {
  this->error(string(msg), function, file, line);
}

void ParallelRuntimeErrors::
error(ostringstream &mstr,
      const char* function, const char* file, const int line) {
  this->error(mstr.str(), function, file, line);
}

int ParallelRuntimeErrors::
count() {
  int numMessages = errors.size();
  MPI_Allreduce(MPI_IN_PLACE, &numMessages, 1, MPI_INT, MPI_SUM, this->comm);
  return numMessages;
}

list<string> &ParallelRuntimeErrors::
gather() {
  // Tell other processors to send their erros
  mpi_call(mpiParallelRuntimeErrorsGatherSlave, -1, 0);

  int numMessages = this->count();
  
  // If no processor encountered an error, return
  if (!numMessages) return errors;

  // subtract the number of messages on the master
  numMessages -= errors.size();

  MPI_Status status;
  int count;
  for (int i = 0; i < numMessages; ++i) {
    // get the next message
    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, this->comm, &status);
    MPI_Get_count(&status, MPI_CHAR, &count);

    char buffer[count];
    MPI_Recv(buffer, count, MPI_CHAR, 
             MPI_ANY_SOURCE, MPI_ANY_TAG, this->comm, MPI_STATUS_IGNORE);
    
    errors.push_back(string());
    string &s = errors.back();
    s.assign(buffer, count);
  }

  return errors;
}

void ParallelRuntimeErrors::
gatherSlave() {
  // If no processor encountered an error, return
  if (!this->count()) return;
  
  // send all messages
  for (list<string>::iterator it = errors.begin();
       it != errors.end(); ++it)
    MPI_Send(const_cast<char*>(it->data()), it->length(), MPI_CHAR, 0, 42, comm_cart);

  // finally empty the list
  this->clear();
}

void ParallelRuntimeErrors::
clear() {
  errors.clear();
}

int check_runtime_errors() {
  return runtimeErrors->count();
}

void
initRuntimeErrors() {
  runtimeErrors = new ParallelRuntimeErrors(comm_cart);
}

void
mpiParallelRuntimeErrorsGatherSlave(int node, int parm) {
  runtimeErrors->gatherSlave();
}


