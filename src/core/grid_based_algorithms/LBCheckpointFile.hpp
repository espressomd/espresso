/*
 * Copyright (C) 2021 The ESPResSo project
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
#ifndef GRID_BASED_ALGORITHMS_LBCHECKPOINTFILE_HPP
#define GRID_BASED_ALGORITHMS_LBCHECKPOINTFILE_HPP

#include <utils/Vector.hpp>

#include <fstream>
#include <ios>
#include <string>
#include <vector>

/** Handle for a LB checkpoint file. */
class LBCheckpointFile {
private:
  bool m_binary;

public:
  std::fstream stream;

  LBCheckpointFile(std::string const &filename, std::ios_base::openmode mode,
                   bool binary) {
    m_binary = binary;
    auto flags = mode;
    if (m_binary)
      flags |= std::ios_base::binary;
    stream.open(filename, flags);
  }

  ~LBCheckpointFile() = default;

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

#endif
