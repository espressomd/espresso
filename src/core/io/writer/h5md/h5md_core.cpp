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
File::File(std::string const& filename, std::string const& script_path)
{
    // Check if a file with given filename exists.
    bool _file_exists = file_exists(filename);
    // If it exists, check for the H5MD structure.
    if (_file_exists)
    {
        this->_has_H5MD_structure = _check_for_H5MD_structure(filename);
    }
    if(_file_exists && this->_has_H5MD_structure)
    {
        // If the H5MD structure is present in the file, just open it
        this->_file = h5xx::file(filename, MPI_COMM_WORLD, MPI_INFO_NULL, h5xx::file::out);
        this->_group_particles =  h5xx::group(this->_file, "particles");
        this->_group_particles_atoms = h5xx::group(this->_group_particles, "atoms");
        this->_group_particles_atoms_box = h5xx::group(this->_group_particles_atoms, "box");
        h5xx::write_attribute(_group_particles_atoms_box, "dimension", 3);
        h5xx::write_attribute(_group_particles_atoms_box, "boundary", "periodic");
        this->_dataset_particles_atoms_box_edges = h5xx::dataset(this->_group_particles_atoms_box, "edges");
        this->_group_particles_atoms_mass = h5xx::group(this->_group_particles_atoms, "mass");
        this->_dataset_particles_atoms_mass_value = h5xx::dataset(this->_group_particles_atoms_mass, "value");
        this->_dataset_particles_atoms_mass_time = h5xx::dataset(this->_group_particles_atoms_mass, "time");
        this->_dataset_particles_atoms_mass_step = h5xx::dataset(this->_group_particles_atoms_mass, "step");
        this->_group_particles_atoms_position = h5xx::group(this->_group_particles_atoms, "position");
        this->_dataset_particles_atoms_position_value = h5xx::dataset(this->_group_particles_atoms_position,
                                                                             "value");
        this->_dataset_particles_atoms_position_time = h5xx::dataset(this->_group_particles_atoms_position,
                                                                            "time");
        this->_dataset_particles_atoms_position_step = h5xx::dataset(this->_group_particles_atoms_position,
                                                                            "step");
        this->_group_particles_atoms_velocity = h5xx::group(this->_group_particles_atoms, "velocity");
        this->_dataset_particles_atoms_velocity_value = h5xx::dataset(_group_particles_atoms_velocity, "value");
        this->_dataset_particles_atoms_velocity_time = h5xx::dataset(_group_particles_atoms_velocity, "time");
        this->_dataset_particles_atoms_velocity_step = h5xx::dataset(_group_particles_atoms_velocity, "step");
        this->_group_particles_atoms_force = h5xx::group(_group_particles_atoms, "force");
        this->_dataset_particles_atoms_force_value = h5xx::dataset(_group_particles_atoms_force, "value");
        this->_dataset_particles_atoms_force_time = h5xx::dataset(_group_particles_atoms_force, "time");
        this->_dataset_particles_atoms_force_step = h5xx::dataset(_group_particles_atoms_force, "step");
        this->_group_particles_atoms_image = h5xx::group(_group_particles_atoms, "image");
        this->_dataset_particles_atoms_image_value = h5xx::dataset(_group_particles_atoms_image, "value");
        this->_dataset_particles_atoms_image_time = h5xx::dataset(_group_particles_atoms_image, "time");
        this->_dataset_particles_atoms_image_step = h5xx::dataset(_group_particles_atoms_image, "step");
        this->_dataset_particles_atoms_species = h5xx::dataset(_group_particles_atoms, "species");
        this->_group_parameters = h5xx::group(this->_file, "parameters");
        this->_group_parameters_vmd_structure = h5xx::group(_group_parameters, "vmd_structure");
        this->_group_parameters_files = h5xx::group(_group_parameters, "files");
        this->_dataset_parameters_files_script = h5xx::dataset(_group_parameters_files, "script");
        return;
    } else if (_file_exists && !this->_has_H5MD_structure)
    {
        throw incompatible_h5mdfile();
    } else
    {
        // Create a new h5xx file object.
        this->_file = h5xx::file(filename, MPI_COMM_WORLD, MPI_INFO_NULL, h5xx::file::out);
        // Sample multi_array for dataset creations
        double_array_3d tmp_3d(boost::extents[1][n_part][3]);
        double_array_1d tmp_1d(boost::extents[n_part]);
        int_array_1d int_tmp_1d(boost::extents[n_part]);
        char_array_1d char_tmp_1d(boost::extents[100]);
        hsize_t chunk_dims_3d[3] = {1, static_cast<hsize_t>(n_part), 3};
        hsize_t chunk_dims_1d[1] = {static_cast<hsize_t>(n_part)};

        // Ensure the H5MD structure is present in the file
        this->_group_particles = h5xx::group(this->_file, "particles");
        this->_group_particles_atoms = h5xx::group(this->_group_particles, "atoms");
        this->_group_particles_atoms_box = h5xx::group(this->_group_particles_atoms, "box");
        h5xx::write_attribute(_group_particles_atoms_box, "dimension", 3);
        h5xx::write_attribute(_group_particles_atoms_box, "boundary", "periodic");
        this->_dataset_particles_atoms_box_edges = h5xx::create_dataset(this->_group_particles_atoms_box, "edges", tmp_3d);
        this->_group_particles_atoms_mass = h5xx::group(this->_group_particles_atoms, "mass");
        this->_dataset_particles_atoms_mass_value = h5xx::create_dataset(this->_group_particles_atoms_mass, "value",
                                                                         tmp_1d, h5xx::policy::storage::chunked(1,
                                                                                                                chunk_dims_1d));
        this->_dataset_particles_atoms_mass_time = h5xx::create_dataset(this->_group_particles_atoms_mass, "time",
                                                                        tmp_1d, h5xx::policy::storage::chunked(1,
                                                                                                               chunk_dims_1d));
        this->_dataset_particles_atoms_mass_step = h5xx::create_dataset(this->_group_particles_atoms_mass, "step",
                                                                        int_tmp_1d,
                                                                        h5xx::policy::storage::chunked(1,
                                                                                                       chunk_dims_1d));
        this->_group_particles_atoms_position = h5xx::group(this->_group_particles_atoms, "position");
        this->_dataset_particles_atoms_position_value = h5xx::create_dataset(this->_group_particles_atoms_position,
                                                                             "value", tmp_3d,
                                                                             h5xx::policy::storage::chunked(3,
                                                                                                            chunk_dims_3d));
        this->_dataset_particles_atoms_position_time = h5xx::create_dataset(this->_group_particles_atoms_position,
                                                                            "time", tmp_1d,
                                                                            h5xx::policy::storage::chunked(1,
                                                                                                           chunk_dims_1d));
        this->_dataset_particles_atoms_position_step = h5xx::create_dataset(this->_group_particles_atoms_position,
                                                                            "step", int_tmp_1d,
                                                                            h5xx::policy::storage::chunked(1,
                                                                                                           chunk_dims_1d));
        this->_group_particles_atoms_velocity = h5xx::group(this->_group_particles_atoms, "velocity");
        this->_dataset_particles_atoms_velocity_value = h5xx::create_dataset(_group_particles_atoms_velocity, "value", tmp_3d,
                             h5xx::policy::storage::chunked(3, chunk_dims_3d));
        this->_dataset_particles_atoms_velocity_time = h5xx::create_dataset(_group_particles_atoms_velocity, "time", tmp_1d,
                             h5xx::policy::storage::chunked(1, chunk_dims_1d));
        this->_dataset_particles_atoms_velocity_step = h5xx::create_dataset(_group_particles_atoms_velocity, "step", int_tmp_1d,
                             h5xx::policy::storage::chunked(1, chunk_dims_1d));
        this->_group_particles_atoms_force = h5xx::group(_group_particles_atoms, "force");
        this->_dataset_particles_atoms_force_value = h5xx::create_dataset(_group_particles_atoms_force, "value", tmp_3d,
                             h5xx::policy::storage::chunked(3, chunk_dims_3d));
        this->_dataset_particles_atoms_force_time = h5xx::create_dataset(_group_particles_atoms_force, "time", tmp_1d,
                             h5xx::policy::storage::chunked(1, chunk_dims_1d));
        this->_dataset_particles_atoms_force_step = h5xx::create_dataset(_group_particles_atoms_force, "step", int_tmp_1d,
                             h5xx::policy::storage::chunked(1, chunk_dims_1d));
        this->_group_particles_atoms_image = h5xx::group(_group_particles_atoms, "image");
        this->_dataset_particles_atoms_image_value = h5xx::create_dataset(_group_particles_atoms_image, "value", int_tmp_1d,
                                                                          h5xx::policy::storage::chunked(1, chunk_dims_1d) );
        this->_dataset_particles_atoms_image_time = h5xx::create_dataset(_group_particles_atoms_image, "time", tmp_1d,
                                                                         h5xx::policy::storage::chunked(1, chunk_dims_1d));
        this->_dataset_particles_atoms_image_step = h5xx::create_dataset(_group_particles_atoms_image, "step", int_tmp_1d,
                                                                         h5xx::policy::storage::chunked(1, chunk_dims_1d));
        this->_dataset_particles_atoms_species = h5xx::create_dataset(_group_particles_atoms, "species", int_tmp_1d);
        this->_group_parameters = h5xx::group(this->_file, "parameters");
        this->_group_parameters_vmd_structure = h5xx::group(_group_parameters, "vmd_structure");
        this->_group_parameters_files = h5xx::group(_group_parameters, "files");
        this->_dataset_parameters_files_script = h5xx::create_dataset(_group_parameters_files, "script", char_tmp_1d);
    }
}


// Method to write particle positions
void File::Write()
{
    // Get number of local particles
    int n_local_part = cells_get_n_particles();
    printf("n_local_part: %d.\n", n_local_part);
    // Get the number of particles on all other nodes
    int pref = 1;
    MPI_Exscan(&n_local_part, &pref, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    // Keep static buffers in order not having to allocate them on every
    // function call
    double_array_3d pos(boost::extents[1][1][1]);
    double_array_3d vel(boost::extents[1][1][1]);
    std::vector<int> id;
    std::vector<int> type;
    Cell *local_cell;

    // Realloc static buffers if necessary
    pos.resize(boost::extents[1][n_local_part][3]);
    vel.resize(boost::extents[1][n_local_part][3]);
    id.reserve(n_local_part);
    type.reserve(n_local_part);

    // Prepare data for writing
    // loop over all local cells 
    int particle_index = 0;
    for (int cell_id = 0; cell_id < local_cells.n; ++cell_id)
    {
        local_cell = local_cells.cell[cell_id];
        for (int local_part_id = 0; local_part_id < local_cell->n; ++local_part_id)
        {
            auto current_particle = local_cell->part[local_part_id];
            // store the particle ids
            id.push_back(current_particle.p.identity);
            type.push_back(current_particle.p.type);
            // store folded particle positions
            pos[0][particle_index][0] = current_particle.r.p[0];
            pos[0][particle_index][1] = current_particle.r.p[1];
            pos[0][particle_index][2] = current_particle.r.p[2];
            // store the particles velocities
            vel[0][particle_index][0] = current_particle.m.v[0]/time_step;
            vel[0][particle_index][1] = current_particle.m.v[1]/time_step;
            vel[0][particle_index][2] = current_particle.m.v[2]/time_step;
            particle_index++;
        }
    }
    std::cout << "Generating dataspace." << std::endl;
    h5xx::dataspace dataspace_pos = h5xx::create_dataspace(pos);
    std::cout << "Generating filespace." << std::endl;
    h5xx::dataspace filespace_pos(this->_dataset_particles_atoms_position_value);
    std::cout << "Writing to the dataset." << std::endl;
    h5xx::write_dataset(this->_dataset_particles_atoms_position_value, pos, dataspace_pos, filespace_pos);
}




// Method to write energy of the system
//int File::WriteEnergy(std::vector<std::string> names = {"total", "nonbonded", "bonded"})
//{
//}

/**
int File::_dump_script(std::string const& script_path)
{
    // String vector for storing the python script
    std::vector<std::string> script_string_vector;
    // File stream for the python script
    std::fstream python_file;
    // String to temporarily store the current line of the python script
    std::string line;
    // Open the file
    python_file.open(script_path, std::fstream::in);
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


bool File::_check_for_H5MD_structure(std::string const &filename)
{
    if (h5xx::is_hdf5_file(filename))
    {
        h5xx::file h5mdfile(filename, h5xx::file::in);
        // Check if all groups are present in the file.
        bool groups_exist[11] = {h5xx::exists_group(h5mdfile, "particles"),
                                 h5xx::exists_group(h5mdfile, "particles/atoms"),
                                 h5xx::exists_group(h5mdfile, "particles/atoms/box"),
                                 h5xx::exists_group(h5mdfile, "particles/atoms/mass"),
                                 h5xx::exists_group(h5mdfile, "particles/atoms/position"),
                                 h5xx::exists_group(h5mdfile, "particles/atoms/velocity"),
                                 h5xx::exists_group(h5mdfile, "particles/atoms/force"),
                                 h5xx::exists_group(h5mdfile, "particles/atoms/image"),
                                 h5xx::exists_group(h5mdfile, "parameters"),
                                 h5xx::exists_group(h5mdfile, "parameters/vmd_structure"),
                                 h5xx::exists_group(h5mdfile, "parameters/files")
        };
        // Only if all boolean are TRUE, groups_all_exist will be TRUE.
        bool groups_all_exist = std::all_of(std::begin(groups_exist), std::end(groups_exist),
                           [](bool i)
                           {
                               return i;
                           });
        // Check if all datasets are present in the file.
        bool datasets_exist[18] = {
                h5xx::exists_dataset(h5mdfile, "particles/atoms/box/edges"),
                h5xx::exists_dataset(h5mdfile, "particles/atoms/mass/value"),
                h5xx::exists_dataset(h5mdfile, "particles/atoms/mass/time"),
                h5xx::exists_dataset(h5mdfile, "particles/atoms/mass/step"),
                h5xx::exists_dataset(h5mdfile, "particles/atoms/position/value"),
                h5xx::exists_dataset(h5mdfile, "particles/atoms/position/time"),
                h5xx::exists_dataset(h5mdfile, "particles/atoms/position/step"),
                h5xx::exists_dataset(h5mdfile, "particles/atoms/velocity/value"),
                h5xx::exists_dataset(h5mdfile, "particles/atoms/velocity/time"),
                h5xx::exists_dataset(h5mdfile, "particles/atoms/velocity/step"),
                h5xx::exists_dataset(h5mdfile, "particles/atoms/force/value"),
                h5xx::exists_dataset(h5mdfile, "particles/atoms/force/time"),
                h5xx::exists_dataset(h5mdfile, "particles/atoms/force/step"),
                h5xx::exists_dataset(h5mdfile, "particles/atoms/species"),
                h5xx::exists_dataset(h5mdfile, "particles/atoms/image/value"),
                h5xx::exists_dataset(h5mdfile, "particles/atoms/image/time"),
                h5xx::exists_dataset(h5mdfile, "particles/atoms/image/step"),
                h5xx::exists_dataset(h5mdfile, "parameters/files/script")
        };
        // Only if all boolean are TRUE, datasets_all_exist will be TRUE.
        bool datasets_all_exist = std::all_of(std::begin(datasets_exist), std::end(datasets_exist),
                                            [](bool i)
                                            {
                                                return i;
                                            });
        // Return the logical AND of the two boolean.
        return (groups_all_exist && datasets_all_exist);
    } else {
        return false;
    }
}
} // namespace h5md
} // namespace writer
