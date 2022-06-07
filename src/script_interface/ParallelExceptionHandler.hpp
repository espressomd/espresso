/*
 * Copyright (C) 2022 The ESPResSo project
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
#ifndef ESPRESSO_SCRIPT_INTERFACE_PARALLEL_EXCEPTION_HANDLER_HPP
#define ESPRESSO_SCRIPT_INTERFACE_PARALLEL_EXCEPTION_HANDLER_HPP

#include "core/errorhandling.hpp"

#include <boost/mpi/communicator.hpp>

#include <stdexcept>
#include <string>
#include <utility>

namespace ScriptInterface {
/**
 * Handle exceptions thrown in MPI parallel code.
 *
 * Instantiate this class inside the catch block and after the catch block,
 * like so:
 * @code{.cpp}
 * boost::mpi::communicator world;
 * auto handler = ScriptInterface::ParallelExceptionHandler{world};
 * std::shared_ptr<MyClass> obj;
 * context()->parallel_try_catch([&obj]() {
 *   obj = std::make_shared<MyClass>(2., true);
 * });
 * @endcode
 *
 * Exceptions are handled as follows:
 * * the main rank throws: re-throw on main rank and throw @ref Exception
 *   on all other ranks
 * * one or more of the worker nodes throw: collect error messages from
 *   worker nodes and throw them on the main rank as a @c std::runtime_error,
 *   throw @ref Exception on all other ranks
 *
 * Throwing a @ref Exception guarantees that the partially initialized script
 * interface object won't be registered in the @ref GlobalContext dictionary;
 * this is the only side-effect on worker nodes, since the exception itself
 * is otherwise silently ignored. On the main rank, the thrown exception is
 * converted to a Python exception.
 */
class ParallelExceptionHandler {
public:
  ParallelExceptionHandler(boost::mpi::communicator comm)
      : m_comm(std::move(comm)) {}

  /**
   * @brief Handle exceptions in synchronous code.
   * Error messages queued in the runtime error collector are flushed
   * to standard error if the code throws on any rank.
   * @pre Must be called on all ranks.
   * @pre The @p callback cannot invoke remote functions from the
   *      @ref Communication::MpiCallbacks framework due to blocking
   *      communication (risk of MPI deadlock on worker nodes).
   * @param[in] callback  Callback to execute synchronously on all ranks.
   */
  template <typename T>
  void parallel_try_catch(std::function<void()> const &callback) const {
    try {
      callback();
    } catch (T const &error) {
      handle_impl(&error);
    }
    handle_impl(nullptr);
  }

private:
  void handle_impl(std::exception const *error) const;
  boost::mpi::communicator m_comm;
};
} // namespace ScriptInterface

#endif
