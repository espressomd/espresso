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
#include <h5xx/h5xx.hpp>


class H5mdCore
{
    public:
        H5mdCore(std::string const&, std::string const&);
        ~H5mdCore() {};
        int write_positions();
        int write_velocities();
        int write_forces();
    private:
        int write_species();
        int dump_script(std::string const&);
        h5xx::file* file;
        h5xx::group* group_particles;
        h5xx::group* group_particles_atoms;
        h5xx::group* group_particles_atoms_box;
        h5xx::group* group_particles_atoms_mass;
        h5xx::group* group_particles_atoms_position;
        h5xx::group* group_particles_atoms_force;
        h5xx::group* group_parameters;
        h5xx::group* group_parameters_vmd_structure;
        h5xx::group* group_parameters_files;
        h5xx::dataset* dataset_parameters_files_script;
};
#endif //ESPRESSO_H5MD_CORE_HPP
