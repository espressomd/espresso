/*
 * Copyright (C) 2010,2012-2016,2019 The ESPResSo project
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
#ifndef SIGNAL_HANDLER_HPP
#define SIGNAL_HANDLER_HPP

#include <csignal>

#include "errorhandling.hpp"

/** @brief RAII guard for signal handling
 *
 * This object saves the current signal handler for @p signal,
 * replaces it with a custom handler @p handler and restores it on
 * destruction.
 */
class SignalHandler {
  void (*old_action)(int);

public:
  // Delete all copy and move constructors
  SignalHandler(SignalHandler &&) = delete;
  SignalHandler &operator=(SignalHandler &&) = delete;
  SignalHandler(SignalHandler const &) = delete;
  SignalHandler &operator=(SignalHandler const &) = delete;

  /** @brief Constructor
   *
   * @param[in] signal  Number of signal to replace
   * @param[in] handler Function to handle the signal
   */
  SignalHandler(int signal, void (*handler)(int)) {
    auto new_action = handler;

    auto res = std::signal(SIGINT, new_action);
    if (res == SIG_ERR) {
      throw std::runtime_error("Failed to replace signal handler!");
    }
    old_action = res;
  }

  /** @brief Destructor
   *
   * Restores the handler which was active at the time of
   * construction.
   */
  ~SignalHandler() {
    if (signal(SIGINT, old_action) == SIG_ERR) {
      runtimeErrorMsg() << "Failed to restore signal handler!";
    }
  }
};

#endif // SIGNAL_HANDLER_HPP
