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


/* Constructor for the File class. */
File::File(std::string const &filename, std::string const &script_name)
{
    boost::filesystem::path script_path(script_name);
    boost::filesystem::path absolute_script_path = boost::filesystem::canonical(script_path);
    /* Store the filename in a member variable. */
    this->user_filename = filename;
    /* Get number of local particles. */
    this->n_local_part = cells_get_n_particles();
    /* Check if a file with given filename exists. */
    bool file_exists = this->check_file_exists(filename);
    /* If it exists, check for the H5MD structure. */
    if (file_exists)
    {
        this->has_H5MD_structure = check_for_H5MD_structure(filename);
    }
    if (file_exists && this->has_H5MD_structure)
    {
        /*
         * If the file exists and has a valid H5MD structure, lets create a new file with links to the old file.
         * This has the advantage, that the new file can just be deleted if the simulation crashes at some point and we
         * still have a valid trajectory.
        */
        /* If the H5MD structure is present in the file, just open it. */
        std::string temp_name = filename + "_tmp";
        this->lastfile = h5xx::file(filename, MPI_COMM_WORLD, MPI_INFO_NULL,
                                     h5xx::file::out);
        this->h5md_file = h5xx::file(temp_name, MPI_COMM_WORLD, MPI_INFO_NULL,
                                     h5xx::file::out);
        /* particles -- atoms -- box -- edges */
        /* Create a link to all the groups and datasets. */
        H5Ocopy(lastfile.hid(), "/particles", this->h5md_file.hid(), "/particles", H5P_DEFAULT,
                           H5P_DEFAULT);
        H5Ocopy(lastfile.hid(), "/particles/atoms", this->h5md_file.hid(), "/particles/atoms", H5P_DEFAULT,
                           H5P_DEFAULT);
        H5Ocopy(lastfile.hid(), "/particles/atoms/box", this->h5md_file.hid(), "/particles/atoms/box",
                           H5P_DEFAULT, H5P_DEFAULT);
        H5Ocopy(lastfile.hid(), "/particles/atoms/box/edges", this->h5md_file.hid(),
                           "/particles/atoms/box/edges", H5P_DEFAULT, H5P_DEFAULT);
        this->group_particles = h5xx::group(this->h5md_file, "particles");
        this->group_particles_atoms = h5xx::group(this->group_particles, "atoms");
        this->group_particles_atoms_box = h5xx::group(this->group_particles_atoms, "box");
        this->dataset_particles_atoms_box_edges = h5xx::dataset(this->group_particles_atoms_box, "edges");
        /* particles -- atoms -- mass */
        H5Ocopy(lastfile.hid(), "/particles/atoms/mass", this->h5md_file.hid(), "/particles/atoms/mass",
                           H5P_DEFAULT, H5P_DEFAULT);
        H5Ocopy(lastfile.hid(), "/particles/atoms/mass/value", this->h5md_file.hid(),
                           "/particles/atoms/mass/value",
                           H5P_DEFAULT, H5P_DEFAULT);
        H5Ocopy(lastfile.hid(), "/particles/atoms/mass/time", this->h5md_file.hid(),
                           "/particles/atoms/mass/time",
                           H5P_DEFAULT, H5P_DEFAULT);
        H5Ocopy(lastfile.hid(), "/particles/atoms/mass/step", this->h5md_file.hid(),
                           "/particles/atoms/mass/step",
                           H5P_DEFAULT, H5P_DEFAULT);
        this->group_particles_atoms_mass = h5xx::group(
                this->group_particles_atoms, "mass");
        this->dataset_particles_atoms_mass_value = h5xx::dataset(
                this->group_particles_atoms_mass, "value");
        this->dataset_particles_atoms_mass_time = h5xx::dataset(
                this->group_particles_atoms_mass, "time");
        this->dataset_particles_atoms_mass_step = h5xx::dataset(
                this->group_particles_atoms_mass, "step");
        /* particles -- atoms -- position */
        H5Ocopy(lastfile.hid(), "/particles/atoms/position", this->h5md_file.hid(),
                           "/particles/atoms/position",
                           H5P_DEFAULT, H5P_DEFAULT);
        H5Ocopy(lastfile.hid(), "/particles/atoms/position/value", this->h5md_file.hid(),
                           "/particles/atoms/position/value",
                           H5P_DEFAULT, H5P_DEFAULT);
        H5Ocopy(lastfile.hid(), "/particles/atoms/position/time", this->h5md_file.hid(),
                           "/particles/atoms/position/time",
                           H5P_DEFAULT, H5P_DEFAULT);
        H5Ocopy(lastfile.hid(), "/particles/atoms/position/step", this->h5md_file.hid(),
                           "/particles/atoms/position/step",
                           H5P_DEFAULT, H5P_DEFAULT);
        this->group_particles_atoms_position = h5xx::group(
                this->group_particles_atoms, "position");
        this->dataset_particles_atoms_position_value = h5xx::dataset(
                this->group_particles_atoms_position,
                "value");
        this->dataset_particles_atoms_position_time = h5xx::dataset(
                this->group_particles_atoms_position,
                "time");
        this->dataset_particles_atoms_position_step = h5xx::dataset(
                this->group_particles_atoms_position,
                "step");
        /* particles -- atoms -- velocity */
        H5Ocopy(lastfile.hid(), "/particles/atoms/velocity", this->h5md_file.hid(),
                           "/particles/atoms/velocity",
                           H5P_DEFAULT, H5P_DEFAULT);
        H5Ocopy(lastfile.hid(), "/particles/atoms/velocity/value", this->h5md_file.hid(),
                           "/particles/atoms/velocity/value",
                           H5P_DEFAULT, H5P_DEFAULT);
        H5Ocopy(lastfile.hid(), "/particles/atoms/velocity/time", this->h5md_file.hid(),
                           "/particles/atoms/velocity/time",
                           H5P_DEFAULT, H5P_DEFAULT);
        H5Ocopy(lastfile.hid(), "/particles/atoms/velocity/step", this->h5md_file.hid(),
                           "/particles/atoms/velocity/step",
                           H5P_DEFAULT, H5P_DEFAULT);
        this->group_particles_atoms_velocity = h5xx::group(
                this->group_particles_atoms, "velocity");
        this->dataset_particles_atoms_velocity_value = h5xx::dataset(
                this->group_particles_atoms_velocity, "value");
        this->dataset_particles_atoms_velocity_time = h5xx::dataset(
                this->group_particles_atoms_velocity, "time");
        this->dataset_particles_atoms_velocity_step = h5xx::dataset(
                this->group_particles_atoms_velocity, "step");
        /* particles -- atoms -- force */
        H5Ocopy(lastfile.hid(), "/particles/atoms/force", this->h5md_file.hid(), "/particles/atoms/force",
                           H5P_DEFAULT, H5P_DEFAULT);
        H5Ocopy(lastfile.hid(), "/particles/atoms/force/value", this->h5md_file.hid(),
                           "/particles/atoms/force/value",
                           H5P_DEFAULT, H5P_DEFAULT);
        H5Ocopy(lastfile.hid(), "/particles/atoms/force/time", this->h5md_file.hid(),
                           "/particles/atoms/force/time",
                           H5P_DEFAULT, H5P_DEFAULT);
        H5Ocopy(lastfile.hid(), "/particles/atoms/force/step", this->h5md_file.hid(),
                           "/particles/atoms/force/step",
                           H5P_DEFAULT, H5P_DEFAULT);
        this->group_particles_atoms_force = h5xx::group(
                this->group_particles_atoms, "force");
        this->dataset_particles_atoms_force_value = h5xx::dataset(
                this->group_particles_atoms_force, "value");
        this->dataset_particles_atoms_force_time = h5xx::dataset(
                this->group_particles_atoms_force, "time");
        this->dataset_particles_atoms_force_step = h5xx::dataset(
                this->group_particles_atoms_force, "step");
        /* particles -- atoms -- image */
        H5Ocopy(lastfile.hid(), "/particles/atoms/image", this->h5md_file.hid(), "/particles/atoms/image",
                           H5P_DEFAULT, H5P_DEFAULT);
        H5Ocopy(lastfile.hid(), "/particles/atoms/image/value", this->h5md_file.hid(),
                           "/particles/atoms/image/value",
                           H5P_DEFAULT, H5P_DEFAULT);
        H5Ocopy(lastfile.hid(), "/particles/atoms/image/time", this->h5md_file.hid(),
                           "/particles/atoms/image/time",
                           H5P_DEFAULT, H5P_DEFAULT);
        H5Ocopy(lastfile.hid(), "/particles/atoms/image/step", this->h5md_file.hid(),
                           "/particles/atoms/image/step",
                           H5P_DEFAULT, H5P_DEFAULT);
        this->group_particles_atoms_image = h5xx::group(
                this->group_particles_atoms, "image");
        this->dataset_particles_atoms_image_value = h5xx::dataset(
                this->group_particles_atoms_image, "value");
        this->dataset_particles_atoms_image_time = h5xx::dataset(
                this->group_particles_atoms_image, "time");
        this->dataset_particles_atoms_image_step = h5xx::dataset(
                this->group_particles_atoms_image, "step");
        /* particles -- atoms -- species */
        H5Ocopy(lastfile.hid(), "/particles/atoms/species", this->h5md_file.hid(),
                           "/particles/atoms/species",
                           H5P_DEFAULT, H5P_DEFAULT);
        this->dataset_particles_atoms_species = h5xx::dataset(
                this->group_particles_atoms, "species");
        /* particles -- atoms -- parameters */
        H5Ocopy(lastfile.hid(), "/parameters", this->h5md_file.hid(), "/parameters",
                           H5P_DEFAULT, H5P_DEFAULT);
        H5Ocopy(lastfile.hid(), "/parameters/vmd_structure", this->h5md_file.hid(),
                           "/parameters/vmd_structure",
                           H5P_DEFAULT, H5P_DEFAULT);
        H5Ocopy(lastfile.hid(), "/parameters/files", this->h5md_file.hid(), "/parameters/vmd_structure",
                           H5P_DEFAULT, H5P_DEFAULT);
        this->group_parameters = h5xx::group(this->h5md_file, "parameters");
        this->group_parameters_vmd_structure = h5xx::group(this->group_parameters, "vmd_structure");
        this->group_parameters_files = h5xx::group(this->group_parameters, "files");
        return;
    } else if (file_exists && !this->has_H5MD_structure) {
        throw incompatible_h5mdfile();
    } else
    {
        /* Create a new h5xx file object. */
        this->h5md_file = h5xx::file(filename, MPI_COMM_WORLD, MPI_INFO_NULL,
                                 h5xx::file::out);
        /* Array for max. dimensions of the dataspaces. */
        const std::vector<hsize_t> maxdims = {H5S_UNLIMITED, H5S_UNLIMITED,
                                            H5S_UNLIMITED};
        const std::vector<hsize_t> maxdims_single = {H5S_UNLIMITED};
        /* Sample multi_array for dataset creations. */
        const std::vector<hsize_t> chunk_dims_3d = {1, static_cast<hsize_t>
        (n_part), 3};
        const std::vector<hsize_t> chunk_dims_1d = {1, static_cast<hsize_t>
        (n_part), 1};
        hsize_t chunk_dims_1d_single = 1;
        const std::vector<hsize_t> dims_3d = {0, static_cast<hsize_t>
        (n_part), 3};
        const std::vector<hsize_t> dims_1d_single = {0};
        boost::multi_array<double, 2> box_vec(boost::extents[3][3]);

        box_vec[0][0]=box_l[0]; box_vec[0][1]=0.0; box_vec[0][2]=0.0;
        box_vec[1][0]=0.0; box_vec[1][1]=box_l[1]; box_vec[1][2]=0.0;
        box_vec[2][0]=0.0; box_vec[2][1]=0.0; box_vec[2][2]=box_l[2];

        /* Ensure the H5MD structure is present in the file. */
        /* particles -- atoms -- box -- edges */
        this->group_particles = h5xx::group(this->h5md_file, "particles");
        this->group_particles_atoms = h5xx::group(this->group_particles,
                                                   "atoms");
        this->group_particles_atoms_box = h5xx::group(
                this->group_particles_atoms, "box");
        h5xx::write_attribute(this->group_particles_atoms_box, "dimension", 3);
        h5xx::write_attribute(this->group_particles_atoms_box, "boundary",
                              "periodic");
        this->dataset_particles_atoms_box_edges =
                h5xx::create_dataset(this->group_particles_atoms_box,
                                     "edges", box_vec);
        h5xx::write_dataset(this->dataset_particles_atoms_box_edges,
                            box_vec);
        /* particles -- atoms -- mass */
        this->group_particles_atoms_mass =
                h5xx::group(this->group_particles_atoms, "mass");
        this->dataspace_particles_atoms_mass_value = h5xx::dataspace(
                dims_3d, maxdims);
        this->dataset_particles_atoms_mass_value = h5xx::dataset(
                this->group_particles_atoms_mass,
                "value", this->type_int,
                this->dataspace_particles_atoms_mass_value,
                h5xx::policy::storage::chunked(3, chunk_dims_3d.data()));
        this->dataspace_particles_atoms_mass_value = h5xx::dataspace(
                dims_3d, maxdims);
        this->dataspace_particles_atoms_mass_time = h5xx::dataspace(
                dims_1d_single, maxdims_single);
        this->dataset_particles_atoms_mass_time =
                h5xx::dataset(this->group_particles_atoms_mass,
                              "time", this->type_double,
                              this->dataspace_particles_atoms_mass_time,
                              h5xx::policy::storage::chunked(1,
                                                             &chunk_dims_1d_single));
        this->dataspace_particles_atoms_mass_step = h5xx::dataspace(
                dims_1d_single, maxdims_single);
        this->dataset_particles_atoms_mass_step =
                h5xx::dataset(this->group_particles_atoms_mass,
                              "step", this->type_int,
                              this->dataspace_particles_atoms_mass_step,
                              h5xx::policy::storage::chunked(1,
                                                             &chunk_dims_1d_single));
        /* particles -- atoms -- position */
        this->group_particles_atoms_position =
                h5xx::group(this->group_particles_atoms, "position");
        this->dataspace_particles_atoms_position_value = h5xx::dataspace(
                dims_3d, maxdims);
        this->dataset_particles_atoms_position_value = h5xx::dataset(
                this->group_particles_atoms_position,
                "value", this->type_double,
                this->dataspace_particles_atoms_position_value,
                h5xx::policy::storage::chunked(3, chunk_dims_3d.data()));
        this->dataspace_particles_atoms_position_value = h5xx::dataspace(
                dims_3d, maxdims);
        this->dataspace_particles_atoms_position_time = h5xx::dataspace(
                dims_1d_single, maxdims_single);
        this->dataset_particles_atoms_position_time =
                h5xx::dataset(this->group_particles_atoms_position,
                              "time", this->type_double,
                              this->dataspace_particles_atoms_position_time,
                              h5xx::policy::storage::chunked(1,
                                                             &chunk_dims_1d_single));
        this->dataspace_particles_atoms_position_step = h5xx::dataspace(
                dims_1d_single, maxdims_single);
        this->dataset_particles_atoms_position_step =
                h5xx::dataset(this->group_particles_atoms_position,
                                     "step", this->type_int,
                              this->dataspace_particles_atoms_position_step,
                                     h5xx::policy::storage::chunked(1,
                                                                    &chunk_dims_1d_single));
        /* particles -- atoms -- velocity */
        this->group_particles_atoms_velocity =
                h5xx::group(this->group_particles_atoms, "velocity");
        this->dataspace_particles_atoms_velocity_value = h5xx::dataspace(
                dims_3d, maxdims);
        this->dataset_particles_atoms_velocity_value = h5xx::dataset(
                this->group_particles_atoms_velocity,
                "value", this->type_double,
                this->dataspace_particles_atoms_velocity_value,
                h5xx::policy::storage::chunked(3, chunk_dims_3d.data()));
        this->dataspace_particles_atoms_velocity_value = h5xx::dataspace(
                dims_3d, maxdims);
        this->dataspace_particles_atoms_velocity_time = h5xx::dataspace(
                dims_1d_single, maxdims_single);
        this->dataset_particles_atoms_velocity_time =
                h5xx::dataset(this->group_particles_atoms_velocity,
                              "time", this->type_double,
                              this->dataspace_particles_atoms_velocity_time,
                              h5xx::policy::storage::chunked(1,
                                                             &chunk_dims_1d_single));
        this->dataspace_particles_atoms_velocity_step = h5xx::dataspace(
                dims_1d_single, maxdims_single);
        this->dataset_particles_atoms_velocity_step =
                h5xx::dataset(this->group_particles_atoms_velocity,
                              "step", this->type_int,
                              this->dataspace_particles_atoms_velocity_step,
                              h5xx::policy::storage::chunked(1,
                                                             &chunk_dims_1d_single));
        /* particles -- atoms -- force */
        this->group_particles_atoms_force =
                h5xx::group(this->group_particles_atoms, "force");
        this->dataspace_particles_atoms_force_value = h5xx::dataspace(
                dims_3d, maxdims);
        this->dataset_particles_atoms_force_value = h5xx::dataset(
                this->group_particles_atoms_force,
                "value", this->type_double,
                this->dataspace_particles_atoms_force_value,
                h5xx::policy::storage::chunked(3, chunk_dims_3d.data()));
        this->dataspace_particles_atoms_force_value = h5xx::dataspace(
                dims_3d, maxdims);
        this->dataspace_particles_atoms_force_time = h5xx::dataspace(
                dims_1d_single, maxdims_single);
        this->dataset_particles_atoms_force_time =
                h5xx::dataset(this->group_particles_atoms_force,
                              "time", this->type_double,
                              this->dataspace_particles_atoms_force_time,
                              h5xx::policy::storage::chunked(1,
                                                             &chunk_dims_1d_single));
        this->dataspace_particles_atoms_force_step = h5xx::dataspace(
                dims_1d_single, maxdims_single);
        this->dataset_particles_atoms_force_step =
                h5xx::dataset(this->group_particles_atoms_force,
                              "step", this->type_int,
                              this->dataspace_particles_atoms_force_step,
                              h5xx::policy::storage::chunked(1,
                                                             &chunk_dims_1d_single));
        /* particles -- atoms -- image */
        this->group_particles_atoms_image =
                h5xx::group(this->group_particles_atoms, "image");
        this->dataspace_particles_atoms_image_value = h5xx::dataspace(
                dims_3d, maxdims);
        this->dataset_particles_atoms_image_value = h5xx::dataset(
                this->group_particles_atoms_image,
                "value", this->type_int,
                this->dataspace_particles_atoms_image_value,
                h5xx::policy::storage::chunked(3, chunk_dims_3d.data()));
        this->dataspace_particles_atoms_image_value = h5xx::dataspace(
                dims_3d, maxdims);
        this->dataspace_particles_atoms_image_time = h5xx::dataspace(
                dims_1d_single, maxdims_single);
        this->dataset_particles_atoms_image_time =
                h5xx::dataset(this->group_particles_atoms_image,
                              "time", this->type_double,
                              this->dataspace_particles_atoms_image_time,
                              h5xx::policy::storage::chunked(1,
                                                             &chunk_dims_1d_single));
        this->dataspace_particles_atoms_image_step = h5xx::dataspace(
                dims_1d_single, maxdims_single);
        this->dataset_particles_atoms_image_step =
                h5xx::dataset(this->group_particles_atoms_image,
                              "step", this->type_int,
                              this->dataspace_particles_atoms_image_step,
                              h5xx::policy::storage::chunked(1,
                                                             &chunk_dims_1d_single));

        /* particles -- atoms -- species */
        this->dataspace_particles_atoms_species = h5xx::dataspace
                (chunk_dims_1d);
        this->dataset_particles_atoms_species =
                h5xx::create_dataset(this->group_particles_atoms, "species",
                                     this->type_int, this->dataspace_particles_atoms_species);
        WriteSpecies();
        /* particles -- atoms -- parameters */
        this->group_parameters = h5xx::group(this->h5md_file, "parameters");
        this->group_parameters_vmd_structure = h5xx::group(
                this->group_parameters, "vmd_structure");
        this->group_parameters_files = h5xx::group(this->group_parameters,
                                                    "files");

    }
}


void File::Close()
{
    if (this->lastfile.valid())
    {
        std::string lastfile_name = this->lastfile.name();
        /* Close the h5xx::file object */
        this->lastfile.close();
        boost::filesystem::path lastfile_path(lastfile_name);
        /* If we arrived at this point, removing the last file is save. */
        boost::filesystem::remove(lastfile_path);
        std::string h5mdfile_name = this->h5md_file.name();
        boost::filesystem::path h5mdfile_path(h5mdfile_name);
        boost::filesystem::path user_filename_path(this->user_filename);
        boost::filesystem::rename(h5mdfile_path, user_filename_path);
    }
}


/* Method to write particle positions. */
void File::Write(bool position, bool velocity, bool force)
{
    this->WriteTimedependent3D(position, velocity, force);
}


void File::WriteSpecies()
{
    /* Get the number of particles on all other nodes. */
    int pref = 1;
    MPI_Exscan(&n_local_part, &pref, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    int_array_3d type_data(boost::extents[1][this->n_local_part][1]);
    Cell *local_cell;

    /* Prepare data for writing, loop over all local cells. */
    int particle_index = 0;
    for (int cell_id = 0; cell_id < local_cells.n; ++cell_id)
    {
        local_cell = local_cells.cell[cell_id];
        for (int local_part_id = 0;
             local_part_id < local_cell->n; ++local_part_id)
        {
            auto current_particle = local_cell->part[local_part_id];
            type_data[0][particle_index][0] = current_particle.p.type;
            particle_index++;
        }
    }
    h5xx::write_dataset(this->dataset_particles_atoms_species, type_data);
}


void File::WriteTimedependent3D(bool position, bool velocity, bool force)
{
    /* Get the number of particles on all other nodes. */
    /* TODO: Writing is not yet working in parallel. */
    int pref = 1;
    MPI_Exscan(&n_local_part, &pref, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    double_array_3d pos(boost::extents[1][n_local_part][3]);
    double_array_3d vel(boost::extents[1][n_local_part][3]);
    double_array_3d f(boost::extents[1][n_local_part][3]);
    int_array_3d image(boost::extents[1][n_local_part][3]);
    /* Prepare data for writing, loop over all local cells. */
    Cell *local_cell;
    int particle_index = 0;
    for (int cell_id = 0; cell_id < local_cells.n; ++cell_id)
    {
        local_cell = local_cells.cell[cell_id];
        for (int local_part_id = 0;
             local_part_id < local_cell->n; ++local_part_id)
        {
            auto current_particle = local_cell->part[local_part_id];
            /* store folded particle positions. */
            if (position)
                pos[0][particle_index][0] = current_particle.r.p[0],
                pos[0][particle_index][1] = current_particle.r.p[1],
                pos[0][particle_index][2] = current_particle.r.p[2];
            if (velocity)
                vel[0][particle_index][0] = current_particle.m.v[0] / time_step,
                vel[0][particle_index][1] = current_particle.m.v[1] / time_step,
                vel[0][particle_index][2] = current_particle.m.v[2] / time_step;
            if (force)
                f[0][particle_index][0] = current_particle.f.f[0],
                f[0][particle_index][1] = current_particle.f.f[1],
                f[0][particle_index][2] = current_particle.f.f[2];
            image[0][particle_index][0] = current_particle.l.i[0];
            image[0][particle_index][1] = current_particle.l.i[1];
            image[0][particle_index][2] = current_particle.l.i[2];
            particle_index++;
        }
    }
    if (position)
        WriteDataset(pos,
                     this->dataset_particles_atoms_position_value,
                     this->dataset_particles_atoms_position_time,
                     this->dataset_particles_atoms_position_step);
    if (velocity)
        WriteDataset(vel,
                     this->dataset_particles_atoms_velocity_value,
                     this->dataset_particles_atoms_velocity_time,
                     this->dataset_particles_atoms_velocity_step);
    if (force)
        WriteDataset(f,
                     this->dataset_particles_atoms_force_value,
                     this->dataset_particles_atoms_force_time,
                     this->dataset_particles_atoms_force_step);
    WriteDataset(image,
    this->dataset_particles_atoms_image_value,
                 this->dataset_particles_atoms_image_time,
                 this->dataset_particles_atoms_image_step);
}


template <typename T>
void File::WriteDataset(T &data, h5xx::dataset& dataset,
                        h5xx::dataset& time, h5xx::dataset& step)
{
    /* Until now the h5xx does not support dataset extending, so we
       have to use the lower level hdf5 library functions. */
    /* Get the ID of the dataspace. */
    hid_t dataspace_local_id = H5Dget_space(dataset.hid());
    hsize_t dims[3];
    /* Get the current dimensions of the dataspace. */
    H5Sget_simple_extent_dims(dataspace_local_id, dims, NULL);
    /* Close the dataspace. */
    H5Sclose(dataspace_local_id);
    // The offset is just set in the time dimension, which is the first one. */
    hsize_t offset[3] = {dims[0], 0, 0};
    hsize_t count[3] = {1, static_cast<hsize_t>(this->n_local_part), 3};
    dims[0] += 1;
    /* Extend the dataset for another timestep. */
    H5Dset_extent(dataset.hid(), dims);
    /* Refresh the dataset after extension. */
    dataspace_local_id = H5Dget_space(dataset.hid());
    /* Select the region in the dataspace. */
    H5Sselect_hyperslab(dataspace_local_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    /* Create a temporary dataspace for the current positions. */
    hid_t dataspace_simple_local_id = H5Screate_simple(3, count, this->max_dims);
    /* Finally write the data to the dataset. */
    H5Dwrite(dataset.hid(),
             dataset.get_type(),
             dataspace_simple_local_id, dataspace_local_id, H5P_DEFAULT,
             data.origin());
    H5Sclose(dataspace_simple_local_id);
    H5Sclose(dataspace_local_id);

    /* Write the md time to the position -- time dataset. */
    hsize_t dims_single;
    hid_t dataspace_pos_time_id = H5Dget_space
            (time.hid());
    H5Sget_simple_extent_dims(dataspace_pos_time_id, &dims_single, NULL);
    H5Sclose(dataspace_pos_time_id);
    hsize_t offset_single = dims_single;
    hsize_t count_single = 1;
    dims_single += 1;
    H5Dset_extent(time.hid(), &dims_single);
    dataspace_pos_time_id = H5Dget_space(time.hid());
    H5Sselect_hyperslab(dataspace_pos_time_id, H5S_SELECT_SET, &offset_single,
                        NULL,
                        &count_single,
                        NULL);
    hid_t dataspace_simple_pos_time_id = H5Screate_simple(1, &count_single, this->max_dims_single);
    H5Dwrite(time.hid(),
             time.get_type(),
             dataspace_simple_pos_time_id, dataspace_pos_time_id, H5P_DEFAULT,
             &sim_time);
    H5Sclose(dataspace_simple_pos_time_id);
    H5Sclose(dataspace_pos_time_id);

    /* Write the md step to the position -- step dataset. */
    hid_t dataspace_pos_step_id = H5Dget_space(step.hid());
    H5Sget_simple_extent_dims(dataspace_pos_step_id, &dims_single, NULL);
    H5Sclose(dataspace_pos_step_id);
    offset_single = dims_single;
    count_single = 1;
    dims_single += 1;
    H5Dset_extent(step.hid(), &dims_single);
    dataspace_pos_step_id = H5Dget_space(step.hid());
    H5Sselect_hyperslab(dataspace_pos_step_id, H5S_SELECT_SET, &offset_single, NULL, &count_single, NULL);
    hid_t dataspace_simple_pos_step_id = H5Screate_simple(1, &count_single, this->max_dims_single);
    int sim_step_data= (int)std::round(sim_time/time_step);
    H5Dwrite(step.hid(),
             step.get_type(),
             dataspace_simple_pos_step_id, dataspace_pos_step_id, H5P_DEFAULT,
             &sim_step_data);
    H5Sclose(dataspace_simple_pos_step_id);
    H5Sclose(dataspace_pos_step_id);
}


bool File::check_for_H5MD_structure(std::string const &filename)
{
    if (h5xx::is_hdf5_file(filename))
    {
        h5xx::file h5mdfile(filename, h5xx::file::in);
        /* Check if all groups are present in the file. */
        bool groups_exist[11] = {h5xx::exists_group(h5mdfile, "particles"),
                                 h5xx::exists_group(h5mdfile,
                                                    "particles/atoms"),
                                 h5xx::exists_group(h5mdfile,
                                                    "particles/atoms/box"),
                                 h5xx::exists_group(h5mdfile,
                                                    "particles/atoms/mass"),
                                 h5xx::exists_group(h5mdfile,
                                                    "particles/atoms/position"),
                                 h5xx::exists_group(h5mdfile,
                                                    "particles/atoms/velocity"),
                                 h5xx::exists_group(h5mdfile,
                                                    "particles/atoms/force"),
                                 h5xx::exists_group(h5mdfile,
                                                    "particles/atoms/image"),
                                 h5xx::exists_group(h5mdfile, "parameters"),
                                 h5xx::exists_group(h5mdfile,
                                                    "parameters/vmd_structure"),
                                 h5xx::exists_group(h5mdfile,
                                                    "parameters/files")
        };
        /* Only if all boolean are TRUE, groups_all_exist will be TRUE. */
        bool groups_all_exist = std::all_of(std::begin(groups_exist),
                                            std::end(groups_exist),
                                            [](bool i)
                                            {
                                                return i;
                                            });
        /* Check if all datasets are present in the file. */
        bool datasets_exist[17] = {
                h5xx::exists_dataset(h5mdfile, "particles/atoms/box/edges"),
                h5xx::exists_dataset(h5mdfile, "particles/atoms/mass/value"),
                h5xx::exists_dataset(h5mdfile, "particles/atoms/mass/time"),
                h5xx::exists_dataset(h5mdfile, "particles/atoms/mass/step"),
                h5xx::exists_dataset(h5mdfile,
                                     "particles/atoms/position/value"),
                h5xx::exists_dataset(h5mdfile, "particles/atoms/position/time"),
                h5xx::exists_dataset(h5mdfile, "particles/atoms/position/step"),
                h5xx::exists_dataset(h5mdfile,
                                     "particles/atoms/velocity/value"),
                h5xx::exists_dataset(h5mdfile, "particles/atoms/velocity/time"),
                h5xx::exists_dataset(h5mdfile, "particles/atoms/velocity/step"),
                h5xx::exists_dataset(h5mdfile, "particles/atoms/force/value"),
                h5xx::exists_dataset(h5mdfile, "particles/atoms/force/time"),
                h5xx::exists_dataset(h5mdfile, "particles/atoms/force/step"),
                h5xx::exists_dataset(h5mdfile, "particles/atoms/species"),
                h5xx::exists_dataset(h5mdfile, "particles/atoms/image/value"),
                h5xx::exists_dataset(h5mdfile, "particles/atoms/image/time"),
                h5xx::exists_dataset(h5mdfile, "particles/atoms/image/step")
        };
        /* Only if all boolean are TRUE, datasets_all_exist will be TRUE. */
        bool datasets_all_exist = std::all_of(std::begin(datasets_exist),
                                              std::end(datasets_exist),
                                              [](bool i)
                                              {
                                                  return i;
                                              });
        /* Return the logical AND of the two boolean. */
        return (groups_all_exist && datasets_all_exist);
    } else
    {
        return false;
    }
}



} /* namespace h5md */
} /* namespace writer */
