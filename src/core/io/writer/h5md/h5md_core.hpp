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
#include <h5xx/h5xx.hpp>

namespace writer {
namespace h5md {


typedef boost::multi_array<double,2> particle_data_3d;
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
        File(std::string const& filename, std::string const& python_script_path);
        /**
         * Destructor of the "File" class.
         * Deletes the member variables and closes the H5MD file.
         */
        ~File();
        /**
         * Method to write to the datasets
         */
        int Write();
        /**
         * Method to write the energy contributions to the H5MD file.
         * @param names An optional vector of strings can be given to decide which
         * energies should be written to the file.
         */
        int WriteEnergy(std::vector<std::string> names);
    private:
        int write_species();
        int dump_script(std::string const&);
        int max_n_part;
        h5xx::file file;
        h5xx::group* group_particles;
        h5xx::group* group_particles_atoms;
        h5xx::group* group_particles_atoms_box;
        h5xx::group* group_particles_atoms_mass;
        h5xx::group* group_particles_atoms_position;
        h5xx::group* group_particles_atoms_species;
        h5xx::group* group_particles_atoms_image;
        h5xx::group* group_particles_atoms_force;
        h5xx::group* group_parameters;
        h5xx::group* group_parameters_vmd_structure;
        h5xx::group* group_parameters_files;
        h5xx::dataset* dataset_parameters_files_script;
};
} // namespace h5md
} // namespace writer
#endif //ESPRESSO_H5MD_CORE_HPP
