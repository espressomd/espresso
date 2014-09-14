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
/** \file errorhandling.hpp 
    This file contains the errorhandling code for severe errors, like
    a broken bond or illegal parameter combinations. See section
    "Errorhandling for developers" for details on the error format and
    how to use this.
*/
#ifndef _ERRORHANDLING_HPP
#define _ERRORHANDLING_HPP

#include "config.hpp"
#include <string>
#include <sstream>
#include <list>
#include <mpi.h>

/** check for runtime errors on all nodes. This has to be called on all nodes synchronously.
    @return the number of characters in the error messages of all nodes together. */
int check_runtime_errors();

/** Gather all error messages from all nodes and return them

    @param errors contains the errors from all nodes. This has to point to an array
    of character pointers, one for each node.
    @return \ref ES_OK if no error occured, otherwise \ref ES_ERROR
*/
int mpi_gather_runtime_errors(char **errors);
void mpi_gather_runtime_errors_slave(int node, int parm);

/** exit ungracefully, core dump if switched on. */
void errexit();

/** register a handler for sigint that translates it into an runtime error. */
void register_sigint_handler();

#define runtimeWarning(msg)                                     \
 runtimeErrors->warning(msg, __PRETTYFUNC__, __FILE__, __LINE__)
#define runtimeError(msg) \
 runtimeErrors->error(msg, __PRETTYFUNC__, __FILE__, __LINE__)


class ParallelRuntimeErrors {
  std::list<std::string> errors;
  MPI_Comm comm;

public:
  ParallelRuntimeErrors(MPI_Comm comm);

  void warning(const std::string &msg,
               const char* function, const char* file, const int line);
  void warning(const char *msg,
               const char* function, const char* file, const int line);
  void warning(std::ostringstream &mstr,
               const char* function, const char* file, const int line);

  void error(const std::string &msg,
             const char* function, const char* file, const int line);
  void error(const char *msg,
             const char* function, const char* file, const int line);
  void error(std::ostringstream &mstr,
             const char* function, const char* file, const int line);

  int count();
  std::list<std::string> &gather();
  void gatherSlave();

  void clear();
};

// Function to initialize the global runtimeErrors object
void initRuntimeErrors();

// callback function for communicator
void mpiParallelRuntimeErrorsGatherSlave(int node, int parm);
extern ParallelRuntimeErrors *runtimeErrors;

#endif
