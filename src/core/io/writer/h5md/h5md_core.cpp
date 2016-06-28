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


namespace writer {
namespace h5md {

// Constructor for the H5mdCore class
File::File(std::string const& filename, std::string const& python_script_path)
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
    // TODO: This method does not work yet.
    //File::dump_script(python_script_path);
}


// Destructor of the H5mdCore class
File::~File()
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
    this->file->close();
}


// Method to write particle positions
int File::WritePositions(std::vector<int> ids)
{
    // Get number of local particles
    int n_local_part = cells_get_n_particles();
    // Keep static buffers in order not having to allocate them on every
    // function call
    particle_data_3d::extent_gen extent_particle_data_3d;
    particle_data_3d pos(extent_particle_data_3d[1][1]);
    particle_data_3d vel(extent_particle_data_3d[1][1]);
    std::vector<int> id;
    std::vector<int> type;
    Cell *local_cell;

    // Realloc static buffers if necessary
    if (n_local_part > pos.size()) pos.resize(extent_particle_data_3d[n_local_part][3]);
    if (n_local_part > vel.size()) vel.resize(extent_particle_data_3d[n_local_part][3]);
    if (n_local_part > id.size()) id.resize(n_local_part);
    if (n_local_part > type.size()) type.resize(n_local_part);

    // Prepare data for writing
    // loop over all local cells 
    int particle_index = 0;
    for (int cell_id = 0; cell_id < local_cells.n; ++cell_id)
    {
        local_cell = local_cells.cell[cell_id];
        for (int local_part_id = 0; local_part_id < local_cell->n; ++local_part_id)
        {
            auto current_particle = local_cell->part[local_part_id];
            id[particle_index] = current_particle.p.identity;
            type[particle_index] = current_particle.p.type;
            pos[particle_index][0] = current_particle.r.p[0];
            pos[particle_index][1] = current_particle.r.p[1];
            pos[particle_index][2] = current_particle.r.p[2];
            vel[particle_index][0] = current_particle.m.v[0];
            vel[particle_index][1] = current_particle.m.v[1];
            vel[particle_index][2] = current_particle.m.v[2];
        }
        particle_index++;
    }
}


// Method to write particle velocities
int File::WriteVelocities(std::vector<int> ids)
{
}


// Method to write particle velocities
int File::WriteForces(std::vector<int> ids)
{
}


// Method to write energy of the system
//int File::WriteEnergy(std::vector<std::string> names = {"total", "nonbonded", "bonded"})
//{
//}

/**
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
*/
} // namespace h5md
} // namespace writer
