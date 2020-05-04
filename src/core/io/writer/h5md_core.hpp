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

#include <algorithm>
#include <string>
#include <unordered_map>

#include "PartCfg.hpp"
#include "ParticleRange.hpp"

namespace Writer {
namespace H5md {

typedef boost::multi_array<double, 3> double_array_3d;
typedef boost::multi_array<int, 3> int_array_3d;

/**
 * @brief Class for writing H5MD files.
 **/
class File {
public:
  /**
   * Constructor/destructor without arguments (due to script_interface).
   * @brief Constructor of the File class.
   */
  File() = default;
  ~File() = default;
  /**
   * @brief Initialize the File object.
   */
  void InitFile();
  /**
   * @brief Method to perform the renaming of the temporary file from
   * "filename" + ".bak" to "filename".
   */
  void Close();

  enum WriteData {
    W_POS = 1 << 0,
    W_V = 1 << 1,
    W_F = 1 << 2,
    W_TYPE = 1 << 3,
    W_MASS = 1 << 4,
    W_CHARGE = 1 << 5
  };
  /**
   * @brief General method to write to the datasets which calls more specific
   * write methods.
   * Boolean values for position, velocity, force and mass.
   */
  void Write(int write_dat, PartCfg &partCfg, const ParticleRange &particles);

  std::string &filename() { return m_filename; };
  std::string &scriptname() { return m_scriptname; };
  // Returns the int that describes which data should be written to the dataset.
  int &what() { return m_what; };
  // Returns the boolean value that describes whether data should be written to
  // the dataset in the order of ids (possibly slower on output for many
  // particles).
  bool &write_ordered() { return m_write_ordered; };
  /**
   * @brief Method to force flush to h5md file.
   */
  void Flush();

private:
  boost::mpi::communicator m_hdf5_comm;
  bool m_already_wrote_bonds = false;

  /**
   * @brief Method to check if the H5MD structure is present in the file.
   * Only call this on valid HDF5 files.
   * @param filename The Name of the hdf5-file to check.
   * @return TRUE if H5MD structure is present, FALSE else.
   */
  bool check_for_H5MD_structure(std::string const &filename);
  /**
   * @brief Method that performs all the low-level stuff for writing the
   * particle
   * positions to the dataset.
   */
  template <typename T>
  void WriteDataset(T &data, const std::string &path,
                    const std::vector<int> &change_extent, hsize_t *offset,
                    hsize_t *count);

  /**
   * @brief Method that extends datasets by the given extent.
   */
  void ExtendDataset(const std::string &path,
                     const std::vector<int> &change_extent);

  /**
   * @brief Method that returns chunk dimensions.
   */
  std::vector<hsize_t> create_chunk_dims(hsize_t dim, hsize_t size,
                                         hsize_t chunk_size);

  /*
   * @brief Method to fill the arrays that are used by WriteDataset particle by
   * particle.
   */
  void fill_arrays_for_h5md_write_with_particle_property(
      int particle_index, int_array_3d &id, int_array_3d &typ,
      double_array_3d &mass, double_array_3d &pos, int_array_3d &image,
      double_array_3d &vel, double_array_3d &f, double_array_3d &charge,
      Particle const &current_particle, int write_dat, int_array_3d &bond);
  /*
   * @brief Method to write the simulation script to the dataset.
   */
  void WriteScript(std::string const &filename);

  /**
   * @brief Creates a new H5MD file.
   * @param filename The filename
   */
  void create_new_file(const std::string &filename);

  /**
   * @brief Loads an existing H5MD file.
   * @param filename The filename
   */
  void load_file(const std::string &filename);

  /**
   * @brief Initializes the necessary data to create the datasets and
   * attributes.
   */
  void init_filestructure();

  /**
   * @brief Creates the necessary HDF5 datasets.
   * @param only_load Set this to true if you want to append to an existing
   * file.
   */
  void create_datasets(bool only_load);

  /**
   * @brief Links the time and step datasets to of all properties to the time
   * and step dataset of the id property. All properties are written at the same
   * time.
   */
  void create_links_for_time_and_step_datasets();

  /**
   * Member variables.
   */
  int m_max_n_part = 0;
  std::string m_filename;
  std::string m_scriptname;
  int m_what;
  bool m_write_ordered;
  std::string m_backup_filename;
  boost::filesystem::path m_absolute_script_path = "nullptr";
  h5xx::file m_h5md_file;

  struct DatasetDescriptor {
    std::string path;
    hsize_t dim;
    h5xx::datatype type;
  };
  std::vector<std::string> group_names;
  std::vector<DatasetDescriptor> dataset_descriptors;
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
