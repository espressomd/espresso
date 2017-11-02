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

#include "communication.hpp"
#include "utils/mpi/gather_buffer.hpp"

#include <boost/mpi/collectives.hpp>

#include <utility>

using namespace std;
using boost::mpi::communicator;
using boost::mpi::all_reduce;

namespace ErrorHandling {

RuntimeErrorCollector::RuntimeErrorCollector(const communicator &comm)
    : m_comm(comm) {}

void RuntimeErrorCollector::message(const RuntimeError &message) {
  m_errors.emplace_back(message);
}

void RuntimeErrorCollector::message(RuntimeError &&message) {
  m_errors.emplace_back(std::move(message));
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
  vector<RuntimeError> all_errors{};
  std::swap(all_errors, m_errors);

  Utils::Mpi::gather_buffer(all_errors, m_comm);

  return all_errors;
}

void RuntimeErrorCollector::gatherSlave() {
  Utils::Mpi::gather_buffer(m_errors, m_comm);

  this->clear();
}

} /* ErrorHandling */
