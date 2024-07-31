/*
 * Copyright (C) 2021-2022 The ESPResSo project
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

#pragma once

#include "config/config.hpp"

#ifdef CUDA

#include <exception>
#include <stdexcept>
#include <string>

/**
 * @brief Wrapper for CUDA runtime exceptions.
 * When the exception cannot be recovered from,
 * prefer using @ref cuda_fatal_error instead.
 */
class cuda_runtime_error : public std::runtime_error {
public:
  cuda_runtime_error(std::string const &msg) : std::runtime_error(msg) {}
};

/**
 * @brief Fatal CUDA exception.
 * Best course of action is to terminate the program immediately.
 */

class cuda_fatal_error {
  std::string m_msg;
  std::terminate_handler m_terminate_handler;

public:
  explicit cuda_fatal_error(std::string msg);

  ~cuda_fatal_error() { terminate(); }

  auto get_terminate() noexcept { return m_terminate_handler; }

  auto set_terminate(std::terminate_handler callback) noexcept {
    auto old_handler = m_terminate_handler;
    m_terminate_handler = callback;
    return old_handler;
  }

  void terminate() noexcept;

  char const *what() const noexcept { return m_msg.c_str(); }
};

/**
 * @brief Invoke a function and silently ignore any thrown
 * @ref cuda_runtime_error error.
 * This is useful when querying the properties of CUDA devices that may
 * not have a suitable CUDA version, or when there is no compatible CUDA
 * device available.
 */
template <class F, class... Args>
void invoke_skip_cuda_exceptions(F &&f, Args &&...args) {
  try {
    return f(args...);
  } catch (cuda_runtime_error const &) { // NOLINT(bugprone-empty-catch)
    // pass
  }
}

#endif // CUDA
