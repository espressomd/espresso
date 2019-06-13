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
#include "RuntimeErrorStream.hpp"

#include "RuntimeErrorCollector.hpp"
#include <utility>

namespace ErrorHandling {
/** ostringstream is not copyable, but it is fine here to copy just the content.
 */
RuntimeErrorStream::RuntimeErrorStream(const RuntimeErrorStream &rhs)
    : m_ec(rhs.m_ec), m_line(rhs.m_line), m_file(rhs.m_file),
      m_function(rhs.m_function) {
  m_buff << rhs.m_buff.rdbuf();
}

RuntimeErrorStream::RuntimeErrorStream(RuntimeErrorCollector &ec,
                                       RuntimeError::ErrorLevel level,
                                       std::string file, const int line,
                                       std::string function)
    : m_ec(ec), m_level(level), m_line(line), m_file(std::move(file)),
      m_function(std::move(function)) {}

RuntimeErrorStream::~RuntimeErrorStream() {
  m_ec.message(m_level, m_buff.str(), m_function.c_str(), m_file.c_str(),
               m_line);
}

} // namespace ErrorHandling
