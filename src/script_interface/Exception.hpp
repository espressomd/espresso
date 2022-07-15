/*
 * Copyright (C) 2020-2022 The ESPResSo project
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
#ifndef ESPRESSO_EXCEPTIONS_HPP
#define ESPRESSO_EXCEPTIONS_HPP

#include <stdexcept>
#include <string>
#include <utility>

namespace ScriptInterface {
struct Exception : public std::exception {
  explicit Exception(const char *msg) : message(msg) {}
  explicit Exception(std::string msg) : message(std::move(msg)) {}

  const char *what() const noexcept override { return message.c_str(); }

private:
  std::string message;
};
} // namespace ScriptInterface

#endif // ESPRESSO_EXCEPTIONS_HPP
