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

#include "h5md_core.hpp"

// Constructor for the H5mdCore class
H5mdCore::H5mdCore(std::string const& filename, std::string const& python_script_path)
{
    // Create a new h5xx file object.
    file = new h5xx::file(filename, h5xx::file::out);
    // Ensure the H5MD structure is present in the file
    group_particles = new h5xx::group(*file, "particles");
    group_particles_atoms = new h5xx::group(*group_particles, "atoms");
    group_particles_atoms_box = new h5xx::group(*group_particles_atoms, "box");
    group_particles_atoms_mass = new h5xx::group(*group_particles_atoms, "mass");
    group_particles_atoms_position = new h5xx::group(*group_particles_atoms, "position");
    group_particles_atoms_force = new h5xx::group(*group_particles_atoms, "force");
    group_parameters = new h5xx::group(*file, "parameters");
    group_parameters_vmd_structure = new h5xx::group(*group_parameters, "vmd_structure");
    group_parameters_files = new h5xx::group(*group_parameters, "files");
    // Write the python script to the H5MD file
    H5mdCore::dump_script(python_script_path);
}


// Destructor of the H5mdCore class
H5mdCore::~H5mdCore()
{
    // Delete all member objects
    delete file;
    delete group_particles;
    delete group_particles_atoms;
    delete group_particles_atoms_box;
    delete group_particles_atoms_mass;
    delete group_particles_atoms_position;
    delete group_particles_atoms_force;
    delete group_parameters;
    delete group_parameters_vmd_structure;
    delete group_parameters_files;
}


int H5mdCore::dump_script(std::string const& python_script_path)
{
    // String vector for storing the python script
    std::vector<std::string> script_string_vector;
    // File stream for the python script
    std::fstream python_file;
    // String to temporarily store the current line of the python script
    std::string line;
    // Open the file
    python_file.open(python_script_path, std::fstream::in);
    if (python_file.is_open())
    {
        while (std::getline(python_file, line))
        {
            std::cout << line << std::endl;
            script_string_vector.push_back(line);
        }
        python_file.close();
    }
    // Create dataset for the dump
    // TODO: Make this running, currently the h5xx::create_dataset method wants fundamental types in the array
    dataset_parameters_files_script = h5xx::create_dataset(*file, "parameters/files/script", script_string_vector)
}
