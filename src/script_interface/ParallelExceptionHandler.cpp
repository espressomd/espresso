/*
 * Copyright (C) 2022 The ESPResSo project
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

#include "ParallelExceptionHandler.hpp"

#include "Exception.hpp"

#include "core/error_handling/RuntimeError.hpp"
#include "core/errorhandling.hpp"

#include <boost/mpi/collectives.hpp>
#include <boost/serialization/string.hpp>

#include <cassert>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

namespace ScriptInterface {

void ParallelExceptionHandler::handle_impl(std::exception const *error) const {
  auto const head_node = 0;
  auto const this_node = m_comm.rank();

  enum : unsigned char {
    NO_RANK_FAILED = 0u,
    SOME_RANK_FAILED = 1u,
    THIS_RANK_SUCCESS = 0u,
    THIS_RANK_FAILED = 1u,
    MAIN_RANK_FAILED = 2u,
  };
  auto const this_fail_flag =
      ((error)
           ? ((this_node == head_node) ? MAIN_RANK_FAILED : THIS_RANK_FAILED)
           : THIS_RANK_SUCCESS);
  auto const fail_flag = boost::mpi::all_reduce(
      m_comm, static_cast<unsigned char>(this_fail_flag), std::bit_or<>());
  auto const main_rank_failed = fail_flag & MAIN_RANK_FAILED;
  auto const some_rank_failed = fail_flag & SOME_RANK_FAILED;

  if (main_rank_failed) {
    flush_runtime_errors_local();
    if (this_node == head_node) {
      throw;
    }
    throw Exception("");
  }

  if (some_rank_failed) {
    flush_runtime_errors_local();
    std::vector<std::string> messages;
    std::string this_message{(error) ? error->what() : ""};
    boost::mpi::gather(m_comm, this_message, messages, head_node);
    if (this_node == head_node) {
      std::string error_message{"an error occurred on one or more MPI ranks:"};
      for (std::size_t i = 0; i < messages.size(); ++i) {
        error_message += "\n  rank " + std::to_string(i) + ": " + messages[i];
      }
      throw std::runtime_error(error_message.c_str());
    }
    throw Exception("");
  }
}

} // namespace ScriptInterface
