/*
 * Copyright (C) 2014-2026 The ESPResSo project
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

#include "error_handling/RuntimeErrorCollector.hpp"

#include <utils/mpi/gather_buffer.hpp>

#include <boost/mpi/collectives.hpp>

#include <algorithm>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace ErrorHandling {

RuntimeErrorCollector::RuntimeErrorCollector(boost::mpi::communicator comm)
    : m_comm(std::move(comm)) {}

RuntimeErrorCollector::~RuntimeErrorCollector() {
  if (!m_errors.empty()) {
    /* Print remaining error messages on destruction */
    std::cerr << "There were unhandled errors.\n";
    flush();
  }
}

void RuntimeErrorCollector::message(RuntimeError::ErrorLevel level,
                                    const std::string &msg,
                                    const char *function, const char *file,
                                    const int line) {
  std::lock_guard<std::mutex> lock(mutex);
  m_errors.emplace_back(level, m_comm.rank(), msg, std::string(function),
                        std::string(file), line);
}

void RuntimeErrorCollector::warning(const std::string &msg,
                                    const char *function, const char *file,
                                    const int line) {
  std::lock_guard<std::mutex> lock(mutex);
  m_errors.emplace_back(RuntimeError::ErrorLevel::WARNING, m_comm.rank(), msg,
                        std::string(function), std::string(file), line);
}

void RuntimeErrorCollector::error(const std::string &msg, const char *function,
                                  const char *file, const int line) {
  std::lock_guard<std::mutex> lock(mutex);
  m_errors.emplace_back(RuntimeError::ErrorLevel::ERROR, m_comm.rank(), msg,
                        std::string(function), std::string(file), line);
}

int RuntimeErrorCollector::count() const {
  std::lock_guard<std::mutex> lock(mutex);
  return boost::mpi::all_reduce(m_comm, static_cast<int>(m_errors.size()),
                                std::plus<>());
}

int RuntimeErrorCollector::count_local(RuntimeError::ErrorLevel level) {
  return static_cast<int>(std::ranges::count_if(
      m_errors, [level](auto const &e) { return e.level() >= level; }));
}

void RuntimeErrorCollector::clear() {
  std::lock_guard<std::mutex> lock(mutex);
  m_errors.clear();
}

void RuntimeErrorCollector::flush() {
  {
    std::lock_guard<std::mutex> lock(mutex);
    for (auto const &e : m_errors) {
      std::cerr << e.format() << std::endl;
    }
  }
  this->clear();
}

std::vector<RuntimeError> RuntimeErrorCollector::gather() {
  std::lock_guard<std::mutex> lock(mutex);
  std::vector<RuntimeError> all_errors{};
  std::swap(all_errors, m_errors);

  Utils::Mpi::gather_buffer(all_errors, m_comm);

  return all_errors;
}

void RuntimeErrorCollector::gather_local() {
  {
    std::lock_guard<std::mutex> lock(mutex);
    Utils::Mpi::gather_buffer(m_errors, m_comm);
  }
  this->clear();
}

} // namespace ErrorHandling
