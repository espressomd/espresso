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


namespace Writer {
namespace H5md {

static void backup_file(const std::string &from, const std::string& to)
{
#ifdef H5MD_DEBUG
    std::cout << "Called " << __func__ << " on node " << this_node << std::endl;
#endif
    if (this_node == 0)
    {
        /*
         * If the file itself *and* a backup file exists something must
         * went wrong before.
        */
        boost::filesystem::path pfrom(from), pto(to);
        try
        {
            boost::filesystem::copy_file(pfrom, pto,
                boost::filesystem::copy_option::fail_if_exists);
        } catch(const boost::filesystem::filesystem_error& e)
        {
            throw left_backupfile();
        }
    }
}

static std::vector<hsize_t>
create_dims(hsize_t dim, hsize_t size, hsize_t chunk_size=0)
{
#ifdef H5MD_DEBUG
    std::cout << "Called " << __func__ << " on node " << this_node << std::endl;
#endif
    if (dim == 3)
        return std::vector<hsize_t>{chunk_size, size, dim};
    if (dim == 2)
        return std::vector<hsize_t>{size, dim};
    if (dim == 1)
        return std::vector<hsize_t>{size};
    else
        throw std::runtime_error("H5MD Error: datastets with this dimension are not implemented\n");
}

static std::vector<hsize_t> create_maxdims(hsize_t dim, hsize_t size)
{
#ifdef H5MD_DEBUG
    std::cout << "Called " << __func__ << " on node " << this_node << std::endl;
#endif
    if (dim > 0)
        return std::vector<hsize_t>{H5S_UNLIMITED, H5S_UNLIMITED, H5S_UNLIMITED};
    else
        return std::vector<hsize_t>{size};
}




/* Initialize the file related variables after parameters have been set. */
void File::InitFile()
{
#ifdef H5MD_DEBUG
    H5Eset_auto(H5E_DEFAULT, (H5E_auto_t) H5Eprint, stderr);
    std::cout << "Called " << __func__ << " on node " << this_node << std::endl;
#endif
    boost::filesystem::path script_path(m_scriptname);
    m_absolute_script_path = boost::filesystem::canonical(script_path);
    if(n_part <= 0) {
        throw std::runtime_error("Please first set up particles before initializing the H5md object.");
    }

    init_filestructure();
    bool fileexists = check_file_exists(m_filename);
    /* Perform a barrier synchronization. Otherwise one process might already
     * create the file while another still checks for its existence. */
    MPI_Barrier(MPI_COMM_WORLD);
    if (fileexists) {
        if (check_for_H5MD_structure(m_filename)) {
            /*
             * If the file exists and has a valid H5MD structure, lets create a
             * backup of it.  This has the advantage, that the new file can
             * just be deleted if the simulation crashes at some point and we
             * still have a valid trajectory, we can start from.
            */
            m_backup_filename = m_filename + ".bak";
            if (this_node == 0) backup_file(m_filename, m_backup_filename);
            load_file(m_filename);
        } else {
            throw incompatible_h5mdfile();
        }
    } else {
        create_new_file(m_filename);
    }
}



void File::init_filestructure()
{
#ifdef H5MD_DEBUG
    std::cout << "Called " << __func__ << " on node " << this_node << std::endl;
#endif
    group_names = {
        "particles",
        "particles/atoms",
        "particles/atoms/box",
        "particles/atoms/mass",
        "particles/atoms/id",
        "particles/atoms/type",
        "particles/atoms/position",
        "particles/atoms/velocity",
        "particles/atoms/force",
        "particles/atoms/image",
        "parameters",
        "parameters/vmd_structure",
        "parameters/files"
    };
    
    h5xx::datatype type_double = h5xx::datatype(H5T_NATIVE_DOUBLE);
    h5xx::datatype type_int = h5xx::datatype(H5T_NATIVE_INT);
    hsize_t npart = static_cast<hsize_t>(n_part);
    dataset_descriptors = {
            //path, dim, size, type
        { "particles/atoms/box/edges"     , 1, type_double },
        { "particles/atoms/mass/value"    , 2, type_double },
        { "particles/atoms/mass/time"     , 1, type_double },
        { "particles/atoms/mass/step"     , 1, type_int },
        { "particles/atoms/id/value"      , 2, type_int },
        { "particles/atoms/id/time"       , 1, type_double },
        { "particles/atoms/id/step"       , 1, type_int },
        { "particles/atoms/type/value"    , 2, type_int },
        { "particles/atoms/type/time"     , 1, type_double },
        { "particles/atoms/type/step"     , 1, type_int },
        { "particles/atoms/position/value", 3, type_double },
        { "particles/atoms/position/time" , 1, type_double },
        { "particles/atoms/position/step" , 1, type_int },
        { "particles/atoms/velocity/value", 3, type_double },
        { "particles/atoms/velocity/time" , 1, type_double },
        { "particles/atoms/velocity/step" , 1, type_int },
        { "particles/atoms/force/value"   , 3, type_double },
        { "particles/atoms/force/time"    , 1, type_double },
        { "particles/atoms/force/step"    , 1, type_int },
        { "particles/atoms/image/value"   , 3, type_int },
        { "particles/atoms/image/time"    , 1, type_double },
        { "particles/atoms/image/step"    , 1, type_int },
    };
}

void File::create_groups()
{
#ifdef H5MD_DEBUG
    std::cout << "Called " << __func__ << " on node " << this_node << std::endl;
#endif
    for (const auto& path: group_names) {
        auto i = path.find_last_of('/');

        if (i == std::string::npos) {
            groups[path] = h5xx::group(m_h5md_file, path);
        } else {
            auto basename = path.substr(i + 1);
            auto father = path.substr(0, i);
            groups[path] = h5xx::group(groups[father], basename);
        }
    }
}

void File::create_datasets(bool only_load)
{
#ifdef H5MD_DEBUG
    std::cout << "Called " << __func__ << " on node " << this_node << std::endl;
#endif
    for (const auto& descr: dataset_descriptors) {
        const std::string& path = descr.path;
        auto i = path.find_last_of('/');

        if (i != std::string::npos) {
            auto basename = path.substr(i + 1);
            auto father = path.substr(0, i);

            if (only_load) {
                datasets[path] = h5xx::dataset(groups[father],
                                               basename);
            } else {
                int creation_size_dataset=0;
                int chunk_size=1;
                if(descr.dim>1){
                    //we deal now with a particle based property
                    chunk_size=n_part;
                }
                auto dims = create_dims(descr.dim, creation_size_dataset);
                auto chunk_dims = create_dims(descr.dim, chunk_size, 1);
                auto maxdims = create_maxdims(descr.dim, creation_size_dataset);
                auto storage = h5xx::policy::storage::chunked(chunk_dims).set(
                        h5xx::policy::storage::fill_value(-10));
                auto dataspace = h5xx::dataspace(dims, maxdims);
                datasets[path] = h5xx::dataset(groups[father],basename, descr.type,dataspace,storage);
            }
        }
    }
}

void File::load_file(const std::string& filename)
{
#ifdef H5MD_DEBUG
    std::cout << "Called " << __func__ << " on node " << this_node << std::endl;
#endif
    m_h5md_file = h5xx::file(filename, MPI_COMM_WORLD, MPI_INFO_NULL,
                             h5xx::file::out);
#ifdef H5MD_DEBUG
    std::cout << "Finished opening the h5 file on node " << this_node << std::endl;
#endif
    create_groups();
    bool only_load = true;
    create_datasets(only_load);
}

void File::create_new_file(const std::string &filename)
{
#ifdef H5MD_DEBUG
    std::cout << "Called " << __func__ << " on node " << this_node << std::endl;
#endif
    this->WriteScript(filename);
    /* Create a new h5xx file object. */
    m_h5md_file = h5xx::file(filename, MPI_COMM_WORLD, MPI_INFO_NULL,
                             h5xx::file::out);

    create_groups();
    bool only_load = false;
    create_datasets(only_load);
    std::vector<double> boxvec = {box_l[0], box_l[1], box_l[2]};
    h5xx::write_attribute(groups["particles/atoms/box"], "dimension", 3);
    h5xx::write_attribute(groups["particles/atoms/box"], "boundary", "periodic");
    std::string path_edges="particles/atoms/box/edges";
    int extent_edges = 3;
    ExtendDataset(path_edges, extent_edges);
    h5xx::write_dataset(datasets[path_edges], boxvec);
}


void File::Close()
{
#ifdef H5MD_DEBUG
    std::cout << "Called " << __func__ << " on node " << this_node << std::endl;
#endif
    if (this_node == 0) boost::filesystem::remove(m_backup_filename);
}


void File::fill_arrays_for_h5md_write_with_particle_property(int particle_index, int_array_3d& id, int_array_3d& typ, double_array_3d& mass, double_array_3d& pos, int_array_3d& image, double_array_3d& vel, double_array_3d& f, Particle* current_particle,bool write_typ,bool write_mass,bool write_pos, bool write_vel, bool write_force ){
#ifdef H5MD_DEBUG
    /* Turn on hdf5 error messages */
    std::cout << "Called " << __func__ << " on node " << this_node << std::endl;
#endif
	        id[0][particle_index][0] = current_particle->p.identity;
            if (write_typ)
                typ[0][particle_index][0] =current_particle->p.type;
            if (write_mass)
                mass[0][particle_index][0] = current_particle->p.mass;
            /* store folded particle positions. */
            if (write_pos)
            {
                pos[0][particle_index][0] = current_particle->r.p[0];
                pos[0][particle_index][1] = current_particle->r.p[1];
                pos[0][particle_index][2] = current_particle->r.p[2];
                image[0][particle_index][0] = current_particle->l.i[0];
                image[0][particle_index][1] = current_particle->l.i[1];
                image[0][particle_index][2] = current_particle->l.i[2];
            }
            if (write_vel)
            {
                vel[0][particle_index][0] = current_particle->m.v[0] / time_step;
                vel[0][particle_index][1] = current_particle->m.v[1] / time_step;
                vel[0][particle_index][2] = current_particle->m.v[2] / time_step;
            }
            if (write_force)
            {
                /* Scale the stored force with m/(0.5*dt**2.0) to get a real 
                 * world force. */
                double fac = current_particle->p.mass / (0.5 * time_step * time_step);
                f[0][particle_index][0] = current_particle->f.f[0] * fac;
                f[0][particle_index][1] = current_particle->f.f[1] * fac;
                f[0][particle_index][2] = current_particle->f.f[2] * fac;
            }

}

void File::Write(int write_dat)
{
#ifdef H5MD_DEBUG
    std::cout << "Called " << __func__ << " on node " << this_node << std::endl;
#endif
    bool write_typ = write_dat & W_TYPE;
    bool write_pos = write_dat & W_POS;
    bool write_vel = write_dat & W_V;
    bool write_force = write_dat & W_F;
    bool write_mass = write_dat & W_MASS;
    int num_particles_to_be_written;
    num_particles_to_be_written=cells_get_n_particles();

    double_array_3d pos(boost::extents[1][num_particles_to_be_written][3]);
    double_array_3d vel(boost::extents[1][num_particles_to_be_written][3]);
    double_array_3d f(boost::extents[1][num_particles_to_be_written][3]);
    int_array_3d image(boost::extents[1][num_particles_to_be_written][3]);
    int_array_3d id(boost::extents[1][num_particles_to_be_written][1]);
    int_array_3d typ(boost::extents[1][num_particles_to_be_written][1]);
    double_array_3d mass(boost::extents[1][num_particles_to_be_written][1]);

    /* loop over all local cells. */
    Cell *local_cell;
    int particle_index = 0;
    for (int cell_id = 0; cell_id < local_cells.n; ++cell_id)
    {
        local_cell = local_cells.cell[cell_id];
        for (int local_part_id = 0;local_part_id < local_cell->n; ++local_part_id)
        {
            auto & current_particle = local_cell->part[local_part_id];
            fill_arrays_for_h5md_write_with_particle_property(particle_index, id, typ, mass, pos, image, vel, f, &current_particle, write_typ, write_mass, write_pos, write_vel, write_force );
            particle_index++;
        }
    }

    WriteDataset(id, "particles/atoms/id", id);

    if (write_typ)
    {
        WriteDataset(typ, "particles/atoms/type", id);
    }
    if (write_mass)
    {
        WriteDataset(mass, "particles/atoms/mass", id);
    }
    if (write_pos)
    {
        WriteDataset(pos, "particles/atoms/position", id);
        WriteDataset(image, "particles/atoms/image", id);
    }
    if (write_vel)
    {
        WriteDataset(vel, "particles/atoms/velocity", id);
    }
    if (write_force)
    {
        WriteDataset(f, "particles/atoms/force", id);
    }

}


void File::ExtendDataset(std::string path, int extent=1){
    /* Until now the h5xx does not support dataset extending, so we
       have to use the lower level hdf5 library functions. */
    auto& dataset = datasets[path];
    /* Get the current dimensions of the dataspace. */
    hid_t ds = H5Dget_space(dataset.hid());
    hsize_t rank = H5Sget_simple_extent_ndims(ds);
    hsize_t dims[rank], maxdims[rank];
    H5Sget_simple_extent_dims(ds, dims, maxdims);
    H5Sclose(ds);
    /* Extend the dataset for another timestep. */
    dims[0] += extent;
    
    // Extend the dataset for more particles if the particle number increased
    if(rank>1){
        /* Get the number of particles on all other nodes. */
        int n_part = max_seen_particle+1;
        if (n_part > m_max_n_part) {
        	m_max_n_part = n_part;
        }
        dims[1] = m_max_n_part;
    }
    H5Dset_extent(dataset.hid(), dims); //extend all dims is collective
}

/* data is assumed to be three dimensional */
template <typename T>
void File::WriteDataset(T &data, const std::string& path, int_array_3d id )
{
#ifdef H5MD_DEBUG
    /* Turn on hdf5 error messages */
    H5Eset_auto(H5E_DEFAULT, (H5E_auto_t) H5Eprint, stderr);
    std::cout << "Called " << __func__ << " on node " << this_node << std::endl;
    std::cout << "Dataset: " << path << std::endl;
#endif
    H5Eset_auto(H5E_DEFAULT, (H5E_auto_t) H5Eprint, stderr);
    auto& dataset = datasets[path + "/value"];
    ExtendDataset(path+ "/value");
    auto& time = datasets[path + "/time"];
    ExtendDataset(path+ "/time");
    auto& step = datasets[path + "/step"];
    ExtendDataset(path+ "/step");

    int num_particles_to_be_written = data.shape()[1];
    int prefix = 0;
    hid_t ds = H5Dget_space(dataset.hid());
    hsize_t rank = H5Sget_simple_extent_ndims(ds);
    /* Get the current dimensions of the dataspace. */
    hsize_t dims[rank], maxdims[rank];
    H5Sget_simple_extent_dims(ds, dims, maxdims);
    /* Select the region in the dataspace. */
    hsize_t offset[rank];
    hsize_t count[rank];
    hid_t ds_new;
    if(m_write_ordered==false){
        ds = H5Dget_space(dataset.hid());
        //calculate offset (prefix) based on local particle number
        MPI_Exscan(&num_particles_to_be_written, &prefix, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        offset[0]= dims[0]-1;
        count [0] = 1;
        if(rank>1){
            offset[1]= static_cast<hsize_t>(prefix);
            count [1] = static_cast<hsize_t>(num_particles_to_be_written);
        }
        if(rank>2){
            offset[2]= 0;
            count [2] = data.shape()[2];
        }
            
        H5Sselect_hyperslab(ds, H5S_SELECT_SET, offset, NULL, count, NULL);
        /* Create a temporary dataspace. */
        ds_new = H5Screate_simple(rank, count, maxdims);
        /* Finally write the data to the dataset. */
        H5Dwrite(dataset.hid(),
                 dataset.get_type(),
                 ds_new,
                 ds, H5P_DEFAULT,
                 data.origin());
        H5Sclose(ds_new);
        H5Sclose(ds);
    }else {
        for(int i=0;i<num_particles_to_be_written;i++){
            ds = H5Dget_space(dataset.hid());
            prefix=id[0][i][0];
            offset[0]= dims[0]-1;
            count[0]= 1; 
            if(rank>1){
                offset[1]= static_cast<hsize_t>(prefix);
                count[1]= 1; 
            }
            if(rank>2){
                offset[2]= 0;
                count[2]= data.shape()[2];
            }
            H5Sselect_hyperslab(ds, H5S_SELECT_SET, offset, NULL, count, NULL);
            /* Create a temporary dataspace. */
            ds_new = H5Screate_simple(rank, count, maxdims);
            /* Finally write the data to the dataset. */
            T data_single_particle(boost::extents[1][1][3]);
            data_single_particle[0][0][0]=data[0][i][0];
            data_single_particle[0][0][1]=data[0][i][1];
            data_single_particle[0][0][2]=data[0][i][2];
            H5Dwrite(dataset.hid(), dataset.get_type(), ds_new, ds, H5P_DEFAULT, data_single_particle.origin());
            H5Sclose(ds_new);
            H5Sclose(ds);
        }
    }

    /* Write the md time to the position -- time dataset. */
    ds = H5Dget_space(time.hid());
    rank = H5Sget_simple_extent_ndims(ds);
    hsize_t dims_time[rank];
    hsize_t maxdims_time[rank];
    H5Sget_simple_extent_dims(ds, dims_time, maxdims_time);

    hsize_t time_offset[rank];
    time_offset[0] = dims_time[0]-1;
    hsize_t time_count[rank];
    time_count[0] = 0;
    /* Only master node writes time and timestep */
    if (this_node == 0)
        time_count[0] = 1; //write 1 timestep
    H5Sselect_hyperslab(ds, H5S_SELECT_SET, time_offset, NULL, time_count, NULL);
    ds_new = H5Screate_simple(rank, time_count, maxdims_time);
    H5Dwrite(time.hid(),
             time.get_type(),
             ds_new,
             ds, H5P_DEFAULT,
             &sim_time);
    H5Sclose(ds);
    /* Write the md step to the position -- step dataset. */
    /* Same offset, count and dims as the time dataset */
    ds = H5Dget_space(step.hid());
    H5Sget_simple_extent_dims(ds, dims_time, maxdims_time);
    H5Sselect_hyperslab(ds, H5S_SELECT_SET, time_offset, NULL, time_count, NULL);
    int sim_step_data = (int)std::round(sim_time/time_step);
    H5Dwrite(step.hid(),
             step.get_type(),
             ds_new,
             ds,
             H5P_DEFAULT,
             &sim_step_data);
    H5Sclose(ds_new);
    H5Sclose(ds);
}


void File::WriteScript(std::string const &filename)
{
#ifdef H5MD_DEBUG
    std::cout << "Called " << __func__ << " on node " << this_node << std::endl;
#endif
    /* First get the number of lines of the script. */
    hsize_t dims[1] = {1};
    std::string tmp;
    std::ifstream scriptfile(m_absolute_script_path.string());
    /* Read the whole script into a buffer. */
    scriptfile.seekg(0, std::ios::end);
    auto filelen = scriptfile.tellg();
    scriptfile.seekg(0);
    std::vector<char> buffer;
    buffer.reserve(filelen);
    buffer.assign(std::istreambuf_iterator<char>(scriptfile),
                  std::istreambuf_iterator<char>());

    hid_t filetype, dtype, space, dset, file_id, group1_id, group2_id;
    file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    group1_id = H5Gcreate2(file_id, "parameters", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    group2_id = H5Gcreate2(group1_id, "files", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    dtype = H5Tcopy(H5T_C_S1);
    H5Tset_size(dtype, filelen * sizeof(char));

    space = H5Screate_simple(1, dims, NULL);
    /* Create the dataset. */
    dset = H5Dcreate(group2_id, "script", dtype, space, H5P_DEFAULT, H5P_DEFAULT,
                      H5P_DEFAULT);
    /* Write data from buffer to dataset. */
    if (this_node == 0)
        H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data());
    /* Clean up. */
    H5Dclose(dset);
    H5Sclose(space);
    H5Tclose(dtype);
    H5Gclose(group1_id);
    H5Gclose(group2_id);
    H5Fclose(file_id);
}


bool File::check_for_H5MD_structure(std::string const &filename)
{
#ifdef H5MD_DEBUG
    std::cout << "Called " << __func__ << " on node " << this_node << std::endl;
#endif
    h5xx::file h5mdfile(filename, h5xx::file::in);

    for (const auto& gnam: group_names)
        if (!h5xx::exists_group(h5mdfile, gnam))
            return false;

    for (const auto& ddesc: dataset_descriptors)
        if (!h5xx::exists_dataset(h5mdfile, ddesc.path))
            return false;

    return true;
}

/* Constructor */
File::File ()
{
#ifdef H5MD_DEBUG
    std::cout << "Called " << __func__ << " on node " << this_node << std::endl;
#endif
}

/* Desctructor */
File::~File()
{
#ifdef H5MD_DEBUG
    std::cout << "Called " << __func__ << " on node " << this_node << std::endl;
#endif
}
} /* namespace H5md */
} /* namespace Writer */
