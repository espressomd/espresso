/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
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


#ifndef ESPRESSO_H5MD_CORE_HPP
#define ESPRESSO_H5MD_CORE_HPP
#include <iostream>
#include <fstream>
#include <string>
#include <mpi.h>
#include "cells.hpp"
#include "global.hpp"
#include <h5xx/h5xx.hpp>
extern double sim_time;
extern double time_step;


namespace writer {
namespace h5md {


typedef boost::multi_array<double,3> double_array_3d;
typedef boost::multi_array<double,1> double_array_1d;
typedef boost::multi_array<int,1> int_array_1d;
typedef boost::multi_array<char, 1> char_array_1d;


/**
 * Class for writing H5MD files.
**/    
class File
{
    public:
    /**
     * Constructor of the "File" class.
     * In the constructor the H5MD structure is initialized.
     * @param filename Path to the .h5 file.
     * @param python_script_path Path to the python simulation script.
     */
    File(std::string const &filename, std::string const &script_path);

    /**
     * Method to write to the datasets.
     * @param An integer denoting what will be written to the dataset. 0:
     * positions, 1: velocities, 2: forces.
     */
    void Write(int what);

    /**
     * Method to write the energy contributions to the H5MD file.
     * @param names An optional vector of strings can be given to decide which
     * energies should be written to the file.
     */
    void WriteEnergy(std::vector<std::string> names);
    private:
    /**
     * Method which uses the low level hdf5 library functions to write to a
     * dataset.
     * @param An integer denoting what shall be written to the dataset. 0:
     * positions, 1: velocities, 2: forces.
     * @param A boost::multiarray containing the three-dimensional data.
     */
    void _write_time_dependent(int what);
    /**
     * Method to check if the H5MD structure is present in the file.
     * @param filename The Name of the hdf5-file to check.
     * @return TRUE if H5MD structure is present, FALSE else.
     */
    bool _check_for_H5MD_structure(std::string const &filename);

    bool _has_H5MD_structure = false;

    int _write_species();

    int _dump_script(std::string const &);
    hsize_t max_dims[3] = {H5S_UNLIMITED, H5S_UNLIMITED, H5S_UNLIMITED};
    hsize_t max_dims_single[1] = {H5S_UNLIMITED};
    int _max_n_part;
    h5xx::file _file;
    // datatypes
    h5xx::datatype type_double;
    h5xx::datatype type_int;
    // particles -- atoms -- box -- edges
    h5xx::group _group_particles;
    h5xx::group _group_particles_atoms;
    h5xx::group _group_particles_atoms_box;
    h5xx::attribute _attribute_particles_atoms_box_dimension;
    h5xx::attribute _attribute_particles_atoms_box_boundary;
    h5xx::dataset _dataset_particles_atoms_box_edges;
    // particles -- atoms -- mass
    h5xx::group _group_particles_atoms_mass;
    h5xx::dataspace _dataspace_particles_atoms_mass_value;
    h5xx::dataset _dataset_particles_atoms_mass_value;
    h5xx::dataspace _dataspace_particles_atoms_mass_time;
    h5xx::dataset _dataset_particles_atoms_mass_time;
    h5xx::dataspace _dataspace_particles_atoms_mass_step;
    h5xx::dataset _dataset_particles_atoms_mass_step;
    // particles -- atoms -- position
    h5xx::group _group_particles_atoms_position;
    h5xx::dataspace _dataspace_particles_atoms_position_value;
    h5xx::dataset _dataset_particles_atoms_position_value;
    h5xx::dataspace _dataspace_particles_atoms_position_time;
    h5xx::dataset _dataset_particles_atoms_position_time;
    h5xx::dataspace _dataspace_particles_atoms_position_step;
    h5xx::dataset _dataset_particles_atoms_position_step;
    // particles -- atoms -- velocity
    h5xx::group _group_particles_atoms_velocity;
    h5xx::dataspace _dataspace_particles_atoms_velocity_value;
    h5xx::dataset _dataset_particles_atoms_velocity_value;
    h5xx::dataspace _dataspace_particles_atoms_velocity_time;
    h5xx::dataset _dataset_particles_atoms_velocity_time;
    h5xx::dataspace _dataspace_particles_atoms_velocity_step;
    h5xx::dataset _dataset_particles_atoms_velocity_step;
    // particles -- atoms -- force
    h5xx::group _group_particles_atoms_force;
    h5xx::dataspace _dataspace_particles_atoms_force_value;
    h5xx::dataset _dataset_particles_atoms_force_value;
    h5xx::dataspace _dataspace_particles_atoms_force_time;
    h5xx::dataset _dataset_particles_atoms_force_time;
    h5xx::dataspace _dataspace_particles_atoms_force_step;
    h5xx::dataset _dataset_particles_atoms_force_step;
    // particles -- atoms -- image
    h5xx::group _group_particles_atoms_image;
    h5xx::dataspace _dataspace_particles_atoms_image_value;
    h5xx::dataset _dataset_particles_atoms_image_value;
    h5xx::dataspace _dataspace_particles_atoms_image_time;
    h5xx::dataset _dataset_particles_atoms_image_time;
    h5xx::dataspace _dataspace_particles_atoms_image_step;
    h5xx::dataset _dataset_particles_atoms_image_step;
    // particles -- atoms -- species
    h5xx::dataspace _dataspace_particles_atoms_species;
    h5xx::dataset _dataset_particles_atoms_species;
    // particles -- atoms -- parameters
    h5xx::group _group_parameters;
    h5xx::group _group_parameters_vmd_structure;
    h5xx::group _group_parameters_files;
    h5xx::dataspace _dataspace_parameters_files_script;
    h5xx::dataset _dataset_parameters_files_script;
};


inline bool file_exists(const std::string& name)
{
    std::ifstream f(name.c_str());
    return f.good();
};

struct incompatible_h5mdfile : public std::exception {
    const char* what () const throw ()
    {
        return "The given hdf5 file does not have a valid h5md structure!";
    }
};
} // namespace h5md
} // namespace writer
#endif //ESPRESSO_H5MD_CORE_HPP
