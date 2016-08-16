/*
  Copyright (C) 2014,2015,2016 The ESPResSo project

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

#include "RuntimeError.hpp"

#include <sstream>

namespace ErrorHandling {

std::string RuntimeError::format() const {
  std::string label;

  switch (m_level) {
  case ErrorLevel::DEBUG:
    label = "DEBUG";
    break;
  case ErrorLevel::WARNING:
    label = "WARNING";
    break;
  case ErrorLevel::ERROR:
    label = "ERROR";
    break;
  case ErrorLevel::INFO:
    label = "INFO";
    break;
  }

  std::ostringstream ostr;

  ostr << label << ": " << what();

#ifndef NDEBUG
  ostr << " in function " << function() << " (" << file() << ":" << line()
       << ") on node " << who();
#endif

  return ostr.str();
}

} /* ErrorHandling */
