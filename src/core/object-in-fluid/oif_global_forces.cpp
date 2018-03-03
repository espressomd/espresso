/*
  Copyright (C) 2012,2013,2014,2015,2016 The ESPResSo project

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

/** \file oif_global_forces.hpp
 *  Routines to calculate the OIF_GLOBAL_FORCES energy or/and and force
 *  for a particle triple (triangle from mesh). (Dupin2007)
 *  \ref forces.cpp
*/

#include "oif_global_forces.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"
#include "interaction_data.hpp"
#include "lb.hpp"
#include "particle_data.hpp"
#include "utils.hpp"
#include "utils/make_unique.hpp" //for creating a unique ptr to a bond class object

/** set parameters for the OIF_GLOBAL_FORCES potential.
*/
int oif_global_forces_set_params(int bond_type, double A0_g, double ka_g,
                                 double V0, double kv) {
  if (bond_type < 0)
    return ES_ERROR;

  //create bond class
  bond_container.set_bond_by_type(bond_type, Utils::make_unique<Bond::OifGlobalForces>
				  (A0_g, ka_g, V0, kv));

  return ES_OK;
}

/** called in force_calc() from within forces.cpp
 *  calculates the global area and global volume for a cell before the forces
 * are handled
 *  sums up parts for area with mpi_reduce from local triangles
 *  synchronization with allreduce
 *
 *  !!! loop over particles from domain_decomposition !!!
 */
void calc_oif_global(double *area_volume, int molType)
{
  double partArea, VOL_partVol;
  partArea = VOL_partVol = 0.0;
  double part_area_volume[2];

  Particle *p, *p1;
  
  for (auto &p : local_cells.particles()) {
    p1 = &p;
    if(p1->p.mol_id == molType){
      bond_container.oif_global_loop(p1, &partArea, &VOL_partVol);
    };
  }//for

  part_area_volume[0] = partArea;
  part_area_volume[1] = VOL_partVol;

  MPI_Allreduce(part_area_volume, area_volume, 2, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
}

void add_oif_global_forces(double *area_volume, int molType) 
{

  double area = area_volume[0];
  double VOL_volume = area_volume[1];

  Particle *p, *p1;

  /*set area and volume for oif_global_forces- Bonds
  which were calculated in calc_oif_global*/
  Bond::OifGlobalForces::set_area_VOL(area, VOL_volume);

  for (auto &p : local_cells.particles()) {
    p1 = &p;
    if(p1->p.mol_id == molType){
      bond_container.oif_global_force_loop(p1);
    };
  }//for
}
