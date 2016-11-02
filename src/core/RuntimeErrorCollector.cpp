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
#include "RuntimeErrorCollector.hpp"
#include <utility>

#include <boost/mpi/collectives.hpp>

#include "communication.hpp"

using namespace std;
using boost::mpi::communicator;
using boost::mpi::all_reduce;

namespace ErrorHandling {

RuntimeErrorCollector::RuntimeErrorCollector(const communicator &comm)
    : m_comm(comm) {}

void RuntimeErrorCollector::message(const RuntimeError &message) {
  m_errors.push_back(message);
}

void RuntimeErrorCollector::message(RuntimeError::ErrorLevel level,
                                    const std::string &msg,
                                    const char *function, const char *file,
                                    const int line) {
  m_errors.emplace_back(level, m_comm.rank(), msg, string(function),
                        string(file), line);
}

void RuntimeErrorCollector::warning(const string &msg, const char *function,
                                    const char *file, const int line) {
  m_errors.emplace_back(RuntimeError::ErrorLevel::WARNING, m_comm.rank(), msg,
                        string(function), string(file), line);
}

void RuntimeErrorCollector::warning(const char *msg, const char *function,
                                    const char *file, const int line) {
  warning(string(msg), function, file, line);
}

void RuntimeErrorCollector::warning(const ostringstream &mstr,
                                    const char *function, const char *file,
                                    const int line) {
  warning(mstr.str(), function, file, line);
}

void RuntimeErrorCollector::error(const string &msg, const char *function,
                                  const char *file, const int line) {
  m_errors.emplace_back(RuntimeError::ErrorLevel::ERROR, m_comm.rank(), msg,
                        string(function), string(file), line);
}

void RuntimeErrorCollector::error(const char *msg, const char *function,
                                  const char *file, const int line) {
  error(string(msg), function, file, line);
}

void RuntimeErrorCollector::error(const ostringstream &mstr,
                                  const char *function, const char *file,
                                  const int line) {
  error(mstr.str(), function, file, line);
}

int RuntimeErrorCollector::count() const {
  int totalMessages;
  const int numMessages = m_errors.size();

  all_reduce(m_comm, numMessages, totalMessages, std::plus<int>());

  return totalMessages;
}

int RuntimeErrorCollector::count(RuntimeError::ErrorLevel level) {
  const int numMessages = std::count_if(
      m_errors.begin(), m_errors.end(),
      [level](const RuntimeError &e) { return e.level() >= level; });
  int totalMessages;

  all_reduce(m_comm, numMessages, totalMessages, std::plus<int>());

  return totalMessages;
}

void RuntimeErrorCollector::clear() { m_errors.clear(); }

vector<RuntimeError> RuntimeErrorCollector::gather() {
  typedef vector<RuntimeError> return_type;

  if (count() == 0) {
    return return_type();
  }

  vector<return_type> all_error_vectors;
  return_type all_errors;

  // Gather the errors on the master
  boost::mpi::gather(m_comm, m_errors, all_error_vectors, 0);

  /** Flaten the vector of vectors */
  for (auto const &v : all_error_vectors) {
    all_errors.insert(all_errors.end(), v.begin(), v.end());
  }

  this->clear();

  return all_errors;
}

void RuntimeErrorCollector::gatherSlave() {
  // If no processor encountered an error, return
  if (count() == 0) {
    return;
  }

  // Gather the errors on the master
  boost::mpi::gather(m_comm, m_errors, 0);

  // finally empty the list
  this->clear();
}

} /* ErrorHandling */
