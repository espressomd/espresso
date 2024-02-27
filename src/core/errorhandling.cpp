/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

/** The callback loop we are on. */
static std::weak_ptr<Communication::MpiCallbacks> m_callbacks;

void init_error_handling(std::weak_ptr<Communication::MpiCallbacks> callbacks) {
  m_callbacks = std::move(callbacks);

  runtimeErrorCollector =
      std::make_unique<RuntimeErrorCollector>(m_callbacks.lock()->comm());
}

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
  m_callbacks.lock()->call(mpi_gather_runtime_errors_local);
  return runtimeErrorCollector->gather();
}

std::vector<RuntimeError> mpi_gather_runtime_errors_all(bool is_head_node) {
  if (is_head_node) {
    return runtimeErrorCollector->gather();
  }
  runtimeErrorCollector->gather_local();
  return {};
}
} // namespace ErrorHandling

void errexit() {
  ErrorHandling::m_callbacks.lock()->comm().abort(1);

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

void flush_runtime_errors_local() {
  ErrorHandling::runtimeErrorCollector->flush();
}
