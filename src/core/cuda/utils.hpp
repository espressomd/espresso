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

#include <stdexcept>
#include <string>

class cuda_runtime_error : public std::runtime_error {
public:
  cuda_runtime_error(std::string const &msg) : std::runtime_error(msg) {}
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
