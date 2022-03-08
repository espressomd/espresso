/*
 * Copyright (C) 2010-2019 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

#ifndef ESPRESSO_H5MD_CORE_HPP
#define ESPRESSO_H5MD_CORE_HPP

#include "BoxGeometry.hpp"
#include "ParticleRange.hpp"

#include <utils/Vector.hpp>

#include <boost/filesystem.hpp>
#include <boost/mpi/communicator.hpp>

#include <h5xx/h5xx.hpp>

#include <cstddef>
#include <exception>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>

namespace h5xx {
template <typename T, std::size_t size>
struct is_array<Utils::Vector<T, size>> : std::true_type {};
} // namespace h5xx

namespace Writer {
namespace H5md {

/**
 * @brief Class for writing H5MD files.
 */
class File {
public:
  /**
   * @brief Constructor of the File class.
   * @param file_path Name for the hdf5 file on disk.
   * @param script_path Path to the simulation script.
   * @param mass_unit The unit for mass.
   * @param length_unit The unit for length.
   * @param time_unit The unit for time.
   * @param force_unit The unit for force.
   * @param velocity_unit The unit for velocity.
   * @param charge_unit The unit for charge.
   * @param comm The MPI communicator.
   */
  File(std::string file_path, std::string script_path, std::string mass_unit,
       std::string length_unit, std::string time_unit, std::string force_unit,
       std::string velocity_unit, std::string charge_unit,
       boost::mpi::communicator comm = boost::mpi::communicator())
      : m_script_path(std::move(script_path)),
        m_mass_unit(std::move(mass_unit)),
        m_length_unit(std::move(length_unit)),
        m_time_unit(std::move(time_unit)), m_force_unit(std::move(force_unit)),
        m_velocity_unit(std::move(velocity_unit)),
        m_charge_unit(std::move(charge_unit)), m_comm(std::move(comm)) {
    init_file(file_path);
  }
  ~File() = default;

  /**
   * @brief Method to perform the renaming of the temporary file from
   * "filename" + ".bak" to "filename".
   */
  void close();

  /**
   * @brief Write data to the hdf5 file.
   * @param particles Particle range for which to write data.
   * @param time Simulation time.
   * @param step Simulation step (monotonically increasing).
   * @param geometry The box dimensions.
   */
  void write(const ParticleRange &particles, double time, int step,
             BoxGeometry const &geometry);

  /**
   * @brief Retrieve the path to the hdf5 file.
   * @return The path as a string.
   */
  std::string file_path() const { return m_h5md_file.name(); }

  /**
   * @brief Retrieve the path to the simulation script.
   * @return The path as a string.
   */
  std::string &script_path() { return m_script_path; }

  /**
   * @brief Retrieve the set mass unit.
   * @return The unit as a string.
   */
  std::string &mass_unit() { return m_mass_unit; }

  /**
   * @brief Retrieve the set length unit.
   * @return The unit as a string.
   */
  std::string &length_unit() { return m_length_unit; }

  /**
   * @brief Retrieve the set time unit.
   * @return The unit as a string.
   */
  std::string &time_unit() { return m_time_unit; }

  /**
   * @brief Retrieve the set force unit.
   * @return The unit as a string.
   */
  std::string &force_unit() { return m_force_unit; }

  /**
   * @brief Retrieve the set velocity unit.
   * @return The unit as a string.
   */
  std::string &velocity_unit() { return m_velocity_unit; }

  /**
   * @brief Retrieve the set charge unit.
   * @return The unit as a string.
   */
  std::string &charge_unit() { return m_charge_unit; }

  /**
   * @brief Method to enforce flushing the buffer to disk.
   */
  void flush();

private:
  /**
   * @brief Initialize the File object.
   */
  void init_file(std::string const &file_path);

  /**
   * @brief Creates a new H5MD file.
   * @param file_path The filename.
   */
  void create_file(const std::string &file_path);

  /**
   * @brief Loads an existing H5MD file.
   * @param file_path The filename.
   */
  void load_file(const std::string &file_path);

  /**
   * @brief Create the HDF5 groups according to the H5MD specification.
   */
  void create_groups();

  /**
   * @brief Creates the necessary HDF5 datasets according to the H5MD
   * specification.
   */
  void create_datasets();

  /**
   * @brief Load datasets of the file.
   */
  void load_datasets();

  /**
   * @brief Write the particle bonds (currently only pairs).
   * @param particles Particle range for which to write bonds.
   */
  void write_connectivity(const ParticleRange &particles);
  /**
   * @brief Write the unit attributes.
   */
  void write_units();
  /**
   * @brief Create hard links for the time and step entries of time-dependent
   * datasets.
   */
  void create_hard_links();
  std::string m_script_path;
  std::string m_mass_unit;
  std::string m_length_unit;
  std::string m_time_unit;
  std::string m_force_unit;
  std::string m_velocity_unit;
  std::string m_charge_unit;
  boost::mpi::communicator m_comm;
  std::string m_backup_filename;
  boost::filesystem::path m_absolute_script_path;
  h5xx::file m_h5md_file;
  std::unordered_map<std::string, h5xx::dataset> datasets;
};

struct incompatible_h5mdfile : public std::exception {
  const char *what() const noexcept override {
    return "The given hdf5 file does not have a valid h5md structure!";
  }
};

struct left_backupfile : public std::exception {
  const char *what() const noexcept override {
    return "A backup of the .h5 file exists. This usually means \
that either you forgot to call the 'close' method or your simulation \
crashed.";
  }
};

} /* namespace H5md */
} /* namespace Writer */
#endif /* ESPRESSO_H5MD_CORE_HPP */
