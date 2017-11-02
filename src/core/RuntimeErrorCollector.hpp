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

#ifndef ERROR_HANDLING_RUNTIMEERRORCOLLECTOR_HPP
#define ERROR_HANDLING_RUNTIMEERRORCOLLECTOR_HPP

#include <string>
#include <vector>

#include <boost/mpi/communicator.hpp>

#include "RuntimeError.hpp"

namespace ErrorHandling {

class RuntimeErrorCollector {
public:
  RuntimeErrorCollector(const boost::mpi::communicator &comm);

  void message(RuntimeError && message);
  void message(const RuntimeError &message);
  void message(RuntimeError::ErrorLevel level, const std::string &msg,
               const char *function, const char *file, const int line);

  void warning(const std::string &msg, const char *function, const char *file,
               const int line);
  void warning(const char *msg, const char *function, const char *file,
               const int line);
  void warning(const std::ostringstream &mstr, const char *function,
               const char *file, const int line);

  void error(const std::string &msg, const char *function, const char *file,
             const int line);
  void error(const char *msg, const char *function, const char *file,
             const int line);
  void error(const std::ostringstream &mstr, const char *function,
             const char *file, const int line);

  /**
   * \brief Return the number of all flying messages.
   *
   * @return Total number of messages.
   */
  int count() const;

  /**
   * \brief Number of Messages that have at least level level.
   *
   * @param level Severity filter.
   * @return Number of Messages that have at least level.
   */
  int count(RuntimeError::ErrorLevel level);

  /**
   * @brief Reset error messages.
   */
  void clear();

  std::vector<RuntimeError> gather();
  void gatherSlave();

private:
  std::vector<RuntimeError> m_errors;
  boost::mpi::communicator m_comm;
};

} /* ErrorHandling */

#endif
