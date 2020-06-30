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

#include <boost/filesystem.hpp>
#include <boost/mpi/communicator.hpp>
#include <h5xx/h5xx.hpp>

#include <BoxGeometry.hpp>
#include <algorithm>
#include <string>
#include <unordered_map>
#include <utility>

#include "PartCfg.hpp"
#include "ParticleRange.hpp"

namespace h5xx {
template <typename T, size_t size>
struct is_array<Utils::Vector<T, size>> : std::true_type {};
} // namespace h5xx

namespace Writer {
namespace H5md {

/**
 * @brief Class for writing H5MD files.
 **/
class File {
public:
  /**
   * Constructor/destructor without arguments (due to script_interface).
   * @brief Constructor of the File class.
   */
  File(std::string filename, std::string scriptname, std::string mass_unit,
       std::string length_unit, std::string time_unit, std::string force_unit,
       std::string velocity_unit, std::string charge_unit,
       boost::mpi::communicator comm = boost::mpi::communicator())
      : m_filename(std::move(filename)), m_scriptname(std::move(scriptname)),
        m_mass_unit(std::move(mass_unit)),
        m_length_unit(std::move(length_unit)),
        m_time_unit(std::move(time_unit)), m_force_unit(std::move(force_unit)),
        m_velocity_unit(std::move(velocity_unit)),
        m_charge_unit(std::move(charge_unit)), m_comm(std::move(comm)) {
    init_file();
  };
  ~File() = default;

  /**
   * @brief Method to perform the renaming of the temporary file from
   * "filename" + ".bak" to "filename".
   */
  void close();

  /**
   * @brief General method to write to the datasets which calls more specific
   * write methods.
   * Boolean values for position, velocity, force and mass.
   */
  void write(const ParticleRange &particles, double time, int step);

  std::string &filename() { return m_filename; };
  std::string &scriptname() { return m_scriptname; };
  std::string &mass_unit() { return m_mass_unit; };
  std::string &length_unit() { return m_length_unit; };
  std::string &time_unit() { return m_time_unit; };
  std::string &force_unit() { return m_force_unit; };
  std::string &velocity_unit() { return m_velocity_unit; };
  std::string &charge_unit() { return m_charge_unit; };

  /**
   * @brief Method to force flush to h5md file.
   */
  void flush();

private:
  /**
   * @brief Initialize the File object.
   */
  void init_file();

  /**
   * @brief Creates a new H5MD file.
   * @param filename The filename
   */
  void create_new_file(const std::string &filename, BoxGeometry &geometry);

  /**
   * @brief Loads an existing H5MD file.
   * @param filename The filename
   */
  void load_file(const std::string &filename);

  void create_groups();
  /**
   * @brief Creates the necessary HDF5 datasets.
   */
  void create_datasets();
  void load_datasets();

  std::string m_filename;
  std::string m_scriptname;
  std::string m_mass_unit;
  std::string m_length_unit;
  std::string m_time_unit;
  std::string m_force_unit;
  std::string m_velocity_unit;
  std::string m_charge_unit;
  boost::mpi::communicator m_comm;

  std::string m_backup_filename;
  boost::filesystem::path m_absolute_script_path = "nullptr";
  h5xx::file m_h5md_file;

  std::unordered_map<std::string, h5xx::dataset> datasets;
  void write_connectivity(const ParticleRange &particles);
  void write_units();
  void create_hard_links();
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
