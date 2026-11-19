/*
 * Copyright (C) 2010-2026 The ESPResSo project
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

#pragma once

#include "BoxGeometry.hpp"
#include "ParticleRange.hpp"
#include "h5md_specification.hpp"

#include <utils/Vector.hpp>

#include <boost/mpi/communicator.hpp>

#include <cstddef>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace HighFive {
class File;
class DataSet;
} // namespace HighFive

namespace Writer {
namespace H5md {

/**
 * @brief Constants which indicate what to output.
 * To indicate the output of multiple fields, OR the
 * corresponding values.
 */
enum H5MDOutputFields : unsigned int {
  H5MD_OUT_NONE = 0u,
  H5MD_OUT_TYPE = 1u,
  H5MD_OUT_POS = 2u,
  H5MD_OUT_IMG = 4u,
  H5MD_OUT_VEL = 8u,
  H5MD_OUT_FORCE = 16u,
  H5MD_OUT_MASS = 32u,
  H5MD_OUT_CHARGE = 64u,
  H5MD_OUT_BONDS = 128u,
  H5MD_OUT_BOX_L = 256u,
  H5MD_OUT_LE_OFF = 512u,
  H5MD_OUT_LE_DIR = 1024u,
  H5MD_OUT_LE_NORMAL = 2048u,
  H5MD_OUT_ALL = 0b1111111111111111u,
};

/**
 * @brief Class for writing H5MD files.
 */
class File {
public:
  /**
   * @brief Constructor.
   * @param file_path Name for the hdf5 file on disk.
   * @param script_path Path to the simulation script.
   * @param output_fields Properties to write to disk.
   * @param mass_unit The unit for mass.
   * @param length_unit The unit for length.
   * @param time_unit The unit for time.
   * @param force_unit The unit for force.
   * @param velocity_unit The unit for velocity.
   * @param charge_unit The unit for charge.
   * @param chunk_size The chunk size for DataSet in hdf5 file
   */
  File(std::filesystem::path file_path, std::filesystem::path script_path,
       std::vector<std::string> const &output_fields, std::string mass_unit,
       std::string length_unit, std::string time_unit, std::string force_unit,
       std::string velocity_unit, std::string charge_unit, int chunk_size);
  ~File();

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
   * @return The path as a file system object.
   */
  auto const &file_path() const { return m_file_path; }

  /**
   * @brief Retrieve the path to the simulation script.
   * @return The path as a file system object.
   */
  auto const &script_path() const { return m_script_path; }

  /**
   * @brief Retrieve the set mass unit.
   * @return The unit as a string.
   */
  auto const &mass_unit() const { return m_mass_unit; }

  /**
   * @brief Retrieve the set length unit.
   * @return The unit as a string.
   */
  auto const &length_unit() const { return m_length_unit; }

  /**
   * @brief Retrieve the set time unit.
   * @return The unit as a string.
   */
  auto const &time_unit() const { return m_time_unit; }

  /**
   * @brief Retrieve the set force unit.
   * @return The unit as a string.
   */
  auto const &force_unit() const { return m_force_unit; }

  /**
   * @brief Retrieve the set velocity unit.
   * @return The unit as a string.
   */
  auto const &velocity_unit() const { return m_velocity_unit; }

  /**
   * @brief Retrieve the set charge unit.
   * @return The unit as a string.
   */
  auto const &charge_unit() const { return m_charge_unit; }

  /**
   * @brief Retrieve the set chunk size.
   */
  auto const &chunk_size() const { return m_chunk_size; }

  /**
   * @brief Build the list of valid output fields.
   * @return The list as a vector of strings.
   */
  std::vector<std::string> valid_fields() const;

  /**
   * @brief Method to enforce flushing the buffer to disk.
   * This pushes all buffers to the Virtual File Driver (VFD) layer.
   * However, the file superblock, which contains the address space information
   * (End of Allocation, EOA), won't be flushed. Therefore, the file won't
   * necessarily be in a guaranteed consistent state for a read operation
   * by another HDF5 parser. Only the @ref close() method can flush the EOA.
   */
  void flush();

private:
  /**
   * @brief Initialize the File object.
   */
  void init_file();

  /**
   * @brief Creates a new H5MD file.
   */
  void create_file();

  /**
   * @brief Loads an existing H5MD file.
   */
  void load_file();

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
   * @brief Create hard links for the simulation time and simulation step
   * entries of time-dependent datasets.
   */
  void create_hard_links();

  std::filesystem::path m_file_path;
  std::filesystem::path m_backup_path;
  std::filesystem::path m_script_path;
  std::filesystem::path m_absolute_script_path;
  std::string m_mass_unit;
  std::string m_length_unit;
  std::string m_time_unit;
  std::string m_force_unit;
  std::string m_velocity_unit;
  std::string m_charge_unit;
  int m_chunk_size;
  boost::mpi::communicator m_comm;
  unsigned int m_fields;
  std::unique_ptr<HighFive::File> m_h5md_file;
  std::unordered_map<std::string, HighFive::DataSet> m_datasets;
  Specification m_h5md_specification;
};

struct incompatible_h5mdfile : public std::exception {
  const char *what() const noexcept override {
    return "The given .h5 file does not match the specifications in 'fields'.";
  }
};

struct left_backupfile : public std::exception {
  const char *what() const noexcept override {
    return "A backup of the .h5 file exists. This usually means that either "
           "you forgot to call the 'close' method or your simulation crashed.";
  }
};

} /* namespace H5md */
} /* namespace Writer */
