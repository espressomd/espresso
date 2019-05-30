/*
Copyright (C) 2010-2018 The ESPResSo project

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
#ifndef ERROR_HANDLING_RUNTIME_ERROR_STREAM_HPP
#define ERROR_HANDLING_RUNTIME_ERROR_STREAM_HPP

#include <sstream>
#include <string>

#include "RuntimeError.hpp"

namespace ErrorHandling {

class RuntimeErrorCollector;

/** \brief Allows creating a runtime error messing by using the streaming
 * operator */

class RuntimeErrorStream {
public:
  RuntimeErrorStream(const RuntimeErrorStream &rhs);
  RuntimeErrorStream(RuntimeErrorCollector &ec, RuntimeError::ErrorLevel level,
                     std::string file, int line, std::string function);
  ~RuntimeErrorStream();
  template <typename T> RuntimeErrorStream &operator<<(T const &value) {
    m_buff << value;

    return *this;
  }

private:
  RuntimeErrorCollector &m_ec;
  RuntimeError::ErrorLevel m_level;
  const int m_line;
  const std::string m_file, m_function;
  std::ostringstream m_buff;
};

} // namespace ErrorHandling

#endif
