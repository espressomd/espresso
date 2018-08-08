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
/** \file errorhandling.hpp
    This file contains the errorhandling code for severe errors, like
    a broken bond or illegal parameter combinations. See section
    "Errorhandling for developers" for details on the error format and
    how to use this.
*/
#ifndef _ERRORHANDLING_HPP
#define _ERRORHANDLING_HPP

#include <iosfwd>
#include <string>
#include <vector>

#include "config.hpp"

#include "RuntimeError.hpp"
#include "RuntimeErrorStream.hpp"

/* Forward declaration of MpiCallbacks,
 * so we don't have to include the header.
 * It depends on mpi and can not be in cuda
 * code.
 */
namespace Communication {
class MpiCallbacks;
}

/**
 * @brief exit ungracefully,
 * core dump if switched on.
 */
void errexit();

/**
 * @brief check for runtime errors on all nodes.
 * This has to be called on all nodes synchronously.
 *
 * @return the number of characters in the error messages of all nodes together.
 */
int check_runtime_errors();

namespace ErrorHandling {

/** register a handler for sigint that translates it into an runtime error. */
void register_sigint_handler();

/**
 * @brief Init error handling system.
 *
 * @param cb Callbacks system the error handler should be on.
 */
void init_error_handling(Communication::MpiCallbacks &callbacks);

void _runtimeMessage(RuntimeError::ErrorLevel level, const std::string &msg,
                     const char *function, const char *file, const int line);

void _runtimeWarning(const char *msg, const char *function, const char *file,
                     const int line);
void _runtimeWarning(const std::string &msg, const char *function,
                     const char *file, const int line);
void _runtimeWarning(const std::ostringstream &msg, const char *function,
                     const char *file, const int line);

void _runtimeError(const char *msg, const char *function, const char *file,
                   const int line);
void _runtimeError(const std::string &msg, const char *function,
                   const char *file, const int line);
void _runtimeError(const std::ostringstream &msg, const char *function,
                   const char *file, const int line);

RuntimeErrorStream _runtimeMessageStream(RuntimeError::ErrorLevel level,
                                         const std::string &file,
                                         const int line,
                                         const std::string &function);

#define runtimeWarning(msg)                                                    \
  ErrorHandling::_runtimeWarning(msg, __PRETTYFUNC__, __FILE__, __LINE__)
#define runtimeError(msg)                                                      \
  ErrorHandling::_runtimeError(msg, __PRETTYFUNC__, __FILE__, __LINE__)

#define runtimeErrorMsg()                                                      \
  ErrorHandling::_runtimeMessageStream(                                        \
      ErrorHandling::RuntimeError::ErrorLevel::ERROR, __FILE__, __LINE__,      \
      __PRETTYFUNC__)

#define runtimeWarningMsg()                                                      \
  ErrorHandling::_runtimeMessageStream(                                        \
      ErrorHandling::RuntimeError::ErrorLevel::WARNING, __FILE__, __LINE__,      \
      __PRETTYFUNC__)

#define debugMsg()                                                             \
  ErrorHandling::_runtimeMessageStream(                                        \
      ErrorHandling::RuntimeError::ErrorLevel::DEBUG, __FILE__, __LINE__,      \
      __PRETTYFUNC__)

std::vector<RuntimeError> mpi_gather_runtime_errors();

} /* ErrorHandling */

#endif
