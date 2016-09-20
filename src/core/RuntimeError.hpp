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

#ifndef ERROR_HANDLING_RUNTIME_ERROR_HPP
#define ERROR_HANDLING_RUNTIME_ERROR_HPP

#include <boost/serialization/access.hpp>
#include <string>

namespace ErrorHandling {

/** \brief A runtime error.
 * This class describes an runtime error,
 * including where it occured and its
 * severity.
 */
struct RuntimeError {
  /** The error level, warnings are only displayed to the user,
   *  errors are fatal.
   */
  enum class ErrorLevel { DEBUG, INFO, WARNING, ERROR };
  RuntimeError() {}
  RuntimeError(ErrorLevel level, int who, const std::string &what,
               const std::string &function, const std::string &file, int line)
      : m_level(level), m_who(who), m_what(what), m_function(function),
        m_file(file), m_line(line) {}

  /** The error level */
  ErrorLevel level() const { return m_level; }
  /** Which MPI node raised the error. */
  int who() const { return m_who; }
  /** The Error Message */
  std::string what() const { return m_what; }
  /** The function where the error occured. */
  std::string function() const { return m_function; }
  /** The file where the error occured. */
  std::string file() const { return m_file; }
  /** The line where the error occured. */
  int line() const { return m_line; }
  /** Get a string representation */
  std::string format() const;

private:
  /** Boost serialization */
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int) {
    ar &m_level;
    ar &m_who;
    ar &m_what;
    ar &m_function;
    ar &m_file;
    ar &m_line;
  }

  ErrorLevel m_level;
  int m_who;
  std::string m_what;
  std::string m_function;
  std::string m_file;
  int m_line;
};

} /* ErrorHandling */

#endif
