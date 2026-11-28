/*
 * Copyright (C) 2010-2026 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
/** \file
 *  Implementation of \ref errorhandling.hpp.
 */

#include "errorhandling.hpp"
#include "communication.hpp"

#include "MpiCallbacks.hpp"
#include "error_handling/RuntimeErrorCollector.hpp"

#include <boost/mpi/collectives.hpp>
#include <boost/mpi/communicator.hpp>

#include <cstdlib>
#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace ErrorHandling {
/** RuntimeErrorCollector instance.
 *  This is a unique pointer so we don't
 *  leak on repeated calls of @ref init_error_handling.
 */
static std::unique_ptr<RuntimeErrorCollector> runtimeErrorCollector;

void init_error_handling(boost::mpi::communicator const &comm) {
  runtimeErrorCollector = std::make_unique<RuntimeErrorCollector>(comm);
}

void deinit_error_handling() { runtimeErrorCollector.reset(); }

RuntimeErrorStream _runtimeMessageStream(RuntimeError::ErrorLevel level,
                                         const std::string &file,
                                         const int line,
                                         const std::string &function) {
  return {*runtimeErrorCollector, level, file, line, function};
}

static void mpi_gather_runtime_errors_local() {
  runtimeErrorCollector->gather_local();
}

REGISTER_CALLBACK(mpi_gather_runtime_errors_local)

std::vector<RuntimeError> mpi_gather_runtime_errors() {
  ::Communication::mpiCallbacks().call(mpi_gather_runtime_errors_local);
  return runtimeErrorCollector->gather();
}

std::vector<RuntimeError> mpi_gather_runtime_errors_all(bool is_head_node) {
  if (is_head_node) {
    return runtimeErrorCollector->gather();
  }
  runtimeErrorCollector->gather_local();
  return {};
}

static void mpi_poll_runtime_messages_local() {
  bool has_any_message =
      runtimeErrorCollector->count_local(
          ErrorHandling::RuntimeError::ErrorLevel::WARNING) != 0;
  boost::mpi::reduce(runtimeErrorCollector->comm(), has_any_message,
                     std::logical_or<bool>(), 0);
}

REGISTER_CALLBACK(mpi_poll_runtime_messages_local)

bool mpi_poll_runtime_messages() {
  bool has_any_message =
      runtimeErrorCollector->count_local(
          ErrorHandling::RuntimeError::ErrorLevel::WARNING) != 0;
  bool has_any_message_global = has_any_message;
  if (runtimeErrorCollector->comm().size() > 1) {
    ::Communication::mpiCallbacks().call(mpi_poll_runtime_messages_local);
    boost::mpi::reduce(runtimeErrorCollector->comm(), has_any_message,
                       has_any_message_global, std::logical_or<bool>(), 0);
  }
  return has_any_message_global;
}
} // namespace ErrorHandling

void errexit() {
  ErrorHandling::runtimeErrorCollector->comm().abort(1);

  std::abort();
}

int check_runtime_errors_local() {
  using namespace ErrorHandling;
  return runtimeErrorCollector->count_local(RuntimeError::ErrorLevel::ERROR);
}

int check_runtime_errors(boost::mpi::communicator const &comm) {
  return boost::mpi::all_reduce(comm, check_runtime_errors_local(),
                                std::plus<int>());
}

void flush_runtime_errors_local() {
  ErrorHandling::runtimeErrorCollector->flush();
}
