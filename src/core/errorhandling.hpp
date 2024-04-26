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
 *  This file contains the errorhandling code for severe errors, like
 *  a broken bond or illegal parameter combinations.
 */

#ifndef ESPRESSO_SRC_CORE_ERROR_HANDLING_HPP
#define ESPRESSO_SRC_CORE_ERROR_HANDLING_HPP

#include "config/config.hpp"

#include "error_handling/RuntimeError.hpp"
#include "error_handling/RuntimeErrorStream.hpp"

#include <memory>
#include <string>
#include <vector>

/* Forward declaration of MpiCallbacks,
 * so we don't have to include the header.
 * It depends on mpi and cannot be in cuda
 * code.
 */
namespace Communication {
class MpiCallbacks;
}

namespace boost {
namespace mpi {
class communicator;
}
} // namespace boost

/**
 * @brief exit ungracefully,
 * core dump if switched on.
 */
[[noreturn]] void errexit();

/**
 * @brief Count runtime errors on all nodes.
 * This has to be called on all nodes synchronously.
 *
 * @return the number of error messages of all nodes together.
 */
int check_runtime_errors(boost::mpi::communicator const &comm);

/**
 * @brief Count runtime errors on the local node.
 * This has to be called on all nodes synchronously.
 *
 * @return the number of error messages on this node.
 */
int check_runtime_errors_local();

/**
 * @brief Flush runtime errors to standard error on the local node.
 * This is used to clear pending runtime error messages when the
 * call site is handling an exception that needs to be re-thrown
 * instead of being queued as an additional runtime error message.
 */
void flush_runtime_errors_local();

namespace ErrorHandling {
/**
 * @brief Initialize the error collection system.
 *
 * @param callbacks Callbacks system the error handler should be on.
 */
void init_error_handling(std::weak_ptr<Communication::MpiCallbacks> callbacks);

RuntimeErrorStream _runtimeMessageStream(RuntimeError::ErrorLevel level,
                                         const std::string &file, int line,
                                         const std::string &function);

#define runtimeErrorMsg()                                                      \
  ErrorHandling::_runtimeMessageStream(                                        \
      ErrorHandling::RuntimeError::ErrorLevel::ERROR, __FILE__, __LINE__,      \
      ESPRESSO_PRETTY_FUNCTION_EXTENSION)

#define runtimeWarningMsg()                                                    \
  ErrorHandling::_runtimeMessageStream(                                        \
      ErrorHandling::RuntimeError::ErrorLevel::WARNING, __FILE__, __LINE__,    \
      ESPRESSO_PRETTY_FUNCTION_EXTENSION)

/** @brief Gather messages on main rank. Only call from main rank. */
std::vector<RuntimeError> mpi_gather_runtime_errors();
/** @brief Gather messages on main rank. Call on all ranks. */
std::vector<RuntimeError> mpi_gather_runtime_errors_all(bool is_head_node);

} // namespace ErrorHandling

#endif
