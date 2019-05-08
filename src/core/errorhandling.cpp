/*
  Copyright (C) 2010-2018 The ESPResSo project
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
/** \file
 *  Implementation of \ref errorhandling.hpp.
 */

#include "errorhandling.hpp"

#include "MpiCallbacks.hpp"
#include "RuntimeErrorCollector.hpp"

#include <boost/mpi/collectives.hpp>
#include <boost/mpi/communicator.hpp>

#include <memory>

namespace ErrorHandling {
namespace {
/** RuntimeErrorCollector instance.
 *  This is a weak pointer so we don't
 *  leak on repeated calls of init_error_handling.
 */
std::unique_ptr<RuntimeErrorCollector> runtimeErrorCollector;

/** The callback loop we are on. */
Communication::MpiCallbacks *m_callbacks = nullptr;
} // namespace

void init_error_handling(Communication::MpiCallbacks &cb) {
  m_callbacks = &cb;

  runtimeErrorCollector =
      std::make_unique<RuntimeErrorCollector>(m_callbacks->comm());
}

RuntimeErrorStream _runtimeMessageStream(RuntimeError::ErrorLevel level,
                                         const std::string &file,
                                         const int line,
                                         const std::string &function) {
  return RuntimeErrorStream(*runtimeErrorCollector, level, file, line,
                            function);
}

void mpi_gather_runtime_errors_slave() { runtimeErrorCollector->gatherSlave(); }
REGISTER_CALLBACK(mpi_gather_runtime_errors_slave)

std::vector<RuntimeError> mpi_gather_runtime_errors() {
  m_callbacks->call(mpi_gather_runtime_errors_slave);
  return runtimeErrorCollector->gather();
}
} // namespace ErrorHandling

void errexit() {
  ErrorHandling::m_callbacks->comm().abort(1);
  std::abort();
}

int check_runtime_errors_local() {
  using namespace ErrorHandling;
  return runtimeErrorCollector->count(RuntimeError::ErrorLevel::ERROR);
}

int check_runtime_errors(boost::mpi::communicator const &comm) {
  return boost::mpi::all_reduce(comm, check_runtime_errors_local(),
                                std::plus<int>());
}
