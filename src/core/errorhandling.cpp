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
#include <csignal>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>

#include "errorhandling.hpp"
#include "utils.hpp"

#include "MpiCallbacks.hpp"
#include "RuntimeErrorCollector.hpp"

using namespace std;

namespace ErrorHandling {

/* Forward declarations */

void mpi_gather_runtime_errors_slave(int node, int parm);

namespace {
/** RuntimeErrorCollector instance.
 *  This is a weak pointer so we don't
 *  leak on repeated calls of @f init_error_handling.
 */
unique_ptr<RuntimeErrorCollector> runtimeErrorCollector;

/** The callback loop we are on. */
Communication::MpiCallbacks *m_callbacks = nullptr;
}

/** Initialize the error collection system. */
void init_error_handling(Communication::MpiCallbacks &cb) {
  m_callbacks = &cb;

  m_callbacks->add(mpi_gather_runtime_errors_slave);

  runtimeErrorCollector = unique_ptr<RuntimeErrorCollector>(
      new RuntimeErrorCollector(m_callbacks->comm()));
}

void _runtimeWarning(const std::string &msg, const char *function,
                     const char *file, const int line) {
  runtimeErrorCollector->warning(msg, function, file, line);
}

void _runtimeWarning(const char *msg, const char *function, const char *file,
                     const int line) {
  runtimeErrorCollector->warning(msg, function, file, line);
}

void _runtimeWarning(const std::ostringstream &msg, const char *function,
                     const char *file, const int line) {
  runtimeErrorCollector->warning(msg, function, file, line);
}

void _runtimeError(const std::string &msg, const char *function,
                   const char *file, const int line) {
  runtimeErrorCollector->error(msg, function, file, line);
}

void _runtimeError(const char *msg, const char *function, const char *file,
                   const int line) {
  runtimeErrorCollector->error(msg, function, file, line);
}

void _runtimeError(const std::ostringstream &msg, const char *function,
                   const char *file, const int line) {
  runtimeErrorCollector->error(msg, function, file, line);
}

RuntimeErrorStream _runtimeMessageStream(RuntimeError::ErrorLevel level,
                                         const std::string &file,
                                         const int line,
                                         const std::string &function) {
  return RuntimeErrorStream(*runtimeErrorCollector, level, file, line,
                            function);
}

vector<RuntimeError> mpi_gather_runtime_errors() {
  // Tell other processors to send their erros
  m_callbacks->call(mpi_gather_runtime_errors_slave, -1, 0);
  return runtimeErrorCollector->gather();
}

void mpi_gather_runtime_errors_slave(int node, int parm) {
  runtimeErrorCollector->gatherSlave();
}

} /* ErrorHandling */

void errexit() {
  ErrorHandling::m_callbacks->comm().abort(1);
  exit(1);
}

int check_runtime_errors() {
  using namespace ErrorHandling;
  return runtimeErrorCollector->count(RuntimeError::ErrorLevel::ERROR);
}
