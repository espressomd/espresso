/*
 * Copyright (C) 2021-2023 The ESPResSo project
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

#pragma once

#include "config/config.hpp"

#ifdef WALBERLA

#include "script_interface/Context.hpp"

#include <utils/Vector.hpp>

#include <boost/mpi/collectives/broadcast.hpp>

#include <fstream>
#include <ios>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace ScriptInterface::walberla {

enum class CptMode : int {
  ascii = 0,
  binary = 1,
  unit_test_runtime_error = -1,
  unit_test_ios_failure = -2
};

/** Inject code for unit tests. */
inline void unit_test_handle(int mode) {
  switch (mode) {
  case static_cast<int>(CptMode::ascii):
  case static_cast<int>(CptMode::binary):
    return;
  case static_cast<int>(CptMode::unit_test_runtime_error):
    throw std::runtime_error("unit test error");
  case static_cast<int>(CptMode::unit_test_ios_failure):
    throw std::ios_base::failure("unit test error");
  default:
    throw std::domain_error("Unknown mode " + std::to_string(mode));
  }
}

/** Handle for a checkpoint file. */
class CheckpointFile {
private:
  bool m_binary;

public:
  std::fstream stream;

  CheckpointFile(std::string const &filename, std::ios_base::openmode mode,
                 bool binary) {
    m_binary = binary;
    auto flags = mode;
    if (m_binary)
      flags |= std::ios_base::binary;
    stream.open(filename, flags);
  }

  ~CheckpointFile() = default;

  template <typename T> void write(T const &value) {
    if (m_binary) {
      stream.write(reinterpret_cast<const char *>(&value), sizeof(T));
    } else {
      stream << value << "\n";
    }
  }

  template <typename T> void write(std::vector<T> const &vector) {
    if (m_binary) {
      stream.write(reinterpret_cast<const char *>(vector.data()),
                   vector.size() * sizeof(T));
    } else {
      for (auto const &value : vector) {
        stream << value << "\n";
      }
    }
  }

  template <typename T, std::size_t N>
  void write(Utils::Vector<T, N> const &vector) {
    if (m_binary) {
      stream.write(reinterpret_cast<const char *>(vector.data()),
                   N * sizeof(T));
    } else {
      stream << Utils::Vector<T, N>::formatter(" ") << vector << "\n";
    }
  }

  template <typename T> void read(T &value) {
    if (m_binary) {
      stream.read(reinterpret_cast<char *>(&value), sizeof(T));
    } else {
      stream >> value;
    }
  }

  template <typename T, std::size_t N> void read(Utils::Vector<T, N> &vector) {
    if (m_binary) {
      stream.read(reinterpret_cast<char *>(vector.data()), N * sizeof(T));
    } else {
      for (auto &value : vector) {
        stream >> value;
      }
    }
  }

  template <typename T> void read(std::vector<T> &vector) {
    if (m_binary) {
      stream.read(reinterpret_cast<char *>(vector.data()),
                  vector.size() * sizeof(T));
    } else {
      for (auto &value : vector) {
        stream >> value;
      }
    }
  }
};

template <typename F1, typename F2, typename F3>
void load_checkpoint_common(Context const &context, std::string const classname,
                            std::string const &filename, int mode,
                            F1 const read_metadata, F2 const read_data,
                            F3 const on_success) {
  auto const err_msg =
      std::string("Error while reading " + classname + " checkpoint: ");
  auto const binary = mode == static_cast<int>(CptMode::binary);
  auto const &comm = context.get_comm();
  auto const is_head_node = context.is_head_node();

  // open file and set exceptions
  CheckpointFile cpfile(filename, std::ios_base::in, binary);
  if (!cpfile.stream) {
    if (is_head_node) {
      throw std::runtime_error(err_msg + "could not open file " + filename);
    }
    return;
  }
  cpfile.stream.exceptions(std::ios_base::failbit | std::ios_base::badbit);

  try {
    read_metadata(cpfile);
    read_data(cpfile);
    comm.barrier();
    on_success();
    // check EOF
    if (!binary) {
      if (cpfile.stream.peek() == '\n') {
        static_cast<void>(cpfile.stream.get());
      }
    }
    if (cpfile.stream.peek() != EOF) {
      throw std::runtime_error(err_msg + "extra data found, expected EOF.");
    }
  } catch (std::ios_base::failure const &) {
    auto const eof_error = cpfile.stream.eof();
    cpfile.stream.close();
    if (eof_error) {
      if (is_head_node) {
        throw std::runtime_error(err_msg + "EOF found.");
      }
      return;
    }
    if (is_head_node) {
      throw std::runtime_error(err_msg + "incorrectly formatted data.");
    }
    return;
  } catch (std::runtime_error const &err) {
    cpfile.stream.close();
    if (is_head_node) {
      throw std::runtime_error(err_msg + err.what());
    }
    return;
  }
}

template <typename F1, typename F2, typename F3>
void save_checkpoint_common(Context const &context, std::string const classname,
                            std::string const &filename, int mode,
                            F1 const write_metadata, F2 const write_data,
                            F3 const on_failure) {
  auto const err_msg =
      std::string("Error while writing " + classname + " checkpoint: ");
  auto const binary = mode == static_cast<int>(CptMode::binary);
  auto const &comm = context.get_comm();
  auto const is_head_node = context.is_head_node();

  // open file and set exceptions
  auto failure = false;
  std::shared_ptr<CheckpointFile> cpfile;
  if (is_head_node) {
    cpfile =
        std::make_shared<CheckpointFile>(filename, std::ios_base::out, binary);
    failure = !cpfile->stream;
    boost::mpi::broadcast(comm, failure, 0);
    if (failure) {
      throw std::runtime_error(err_msg + "could not open file " + filename);
    }
    cpfile->stream.exceptions(std::ios_base::failbit | std::ios_base::badbit);
    if (!binary) {
      cpfile->stream.precision(16);
      cpfile->stream << std::fixed;
    }
  } else {
    boost::mpi::broadcast(comm, failure, 0);
    if (failure) {
      return;
    }
  }

  try {
    write_metadata(cpfile, context);
    write_data(cpfile, context);
  } catch (std::exception const &error) {
    on_failure(cpfile, context);
    if (is_head_node) {
      cpfile->stream.close();
      if (dynamic_cast<std::ios_base::failure const *>(&error)) {
        throw std::runtime_error(err_msg + "could not write to " + filename);
      }
      throw;
    }
  }
}

} // namespace ScriptInterface::walberla

#endif // WALBERLA
