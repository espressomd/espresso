/*
  Copyright (C) 2010,2012,2013 The ESPResSo project
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
/** \file energy.cpp
    Implementation of \ref energy.hpp "energy.hpp".
*/
#include "energy.hpp"
#include "cells.hpp"
#include "integrate.hpp"
#include "initialize.hpp"
#include "domain_decomposition.hpp"
#include "nsquare.hpp"
#include "layered.hpp"
#include "elc.hpp"
#include "magnetic_non_p3m_methods.hpp"
#include "mdlc_correction.hpp"
#include "external_potential.hpp"

Observable_stat energy = {0, {NULL,0,0}, 0,0,0};
Observable_stat total_energy = {0, {NULL,0,0}, 0,0,0};

/************************************************************/
/* local prototypes                                         */
/************************************************************/

/** Calculate long range energies (P3M, MMM2d...). */
void calc_long_range_energies();

/************************************************************/

void energy_calc(double *result)
{
  if (!interactions_sanity_checks())
    return;

  init_energies(&energy);

  on_observable_calc();
  
  switch (cell_structure.type) {
  case CELL_STRUCTURE_LAYERED:
    layered_calculate_energies();
    break;
  case CELL_STRUCTURE_DOMDEC: 
    if(dd.use_vList) {
      if (rebuild_verletlist)  
	build_verlet_lists();
      calculate_verlet_energies();
    }
    else
      calculate_link_cell_energies();
    break;
  case CELL_STRUCTURE_NSQUARE:
    nsq_calculate_energies();
  }
  /* rescale kinetic energy */
  energy.data.e[0] /= (2.0*time_step*time_step);

  calc_long_range_energies();
  
  /* gather data */
  MPI_Reduce(energy.data.e, result, energy.data.n, MPI_DOUBLE, MPI_SUM, 0, comm_cart);

  if (n_external_potentials > 0) {
    double* energies = (double*) malloc(n_external_potentials*sizeof(double));
    for (int i=0; i<n_external_potentials; i++) {
      energies[i]=external_potentials[i].energy;
    }
    double* energies_sum =  (double*) malloc(n_external_potentials*sizeof(double)); 
    MPI_Reduce(energies, energies_sum, n_external_potentials, MPI_DOUBLE, MPI_SUM, 0, comm_cart); 
    for (int i=0; i<n_external_potentials; i++) {
      external_potentials[i].energy=energies_sum[i];
    }
    free(energies);
    free(energies_sum);
  }

}

/************************************************************/

void calc_long_range_energies()
{
#ifdef ELECTROSTATICS  
  /* calculate k-space part of electrostatic interaction. */
  switch (coulomb.method) {
#ifdef P3M
  case COULOMB_P3M_GPU:
    printf("long range energy calculation not implemented for GPU P3M\n"); //TODO make right
    break;
  case COULOMB_P3M:
    p3m_charge_assign(); 
    energy.coulomb[1] = p3m_calc_kspace_forces(0,1);
    break;
  case COULOMB_ELC_P3M:
    // assign the original charges first
    // they may not have been assigned yet
    p3m_charge_assign(); 
    if(!elc_params.dielectric_contrast_on)
      energy.coulomb[1] = p3m_calc_kspace_forces(0,1);
    else {
      energy.coulomb[1] = 0.5*p3m_calc_kspace_forces(0,1); 
      energy.coulomb[1]+= 0.5*ELC_P3M_dielectric_layers_energy_self(); 

      //  assign both original and image charges now
      ELC_p3m_charge_assign_both();
      ELC_P3M_modify_p3m_sums_both();

      energy.coulomb[1] += 0.5*p3m_calc_kspace_forces(0,1); 

      //assign only the image charges now
      ELC_p3m_charge_assign_image();
      ELC_P3M_modify_p3m_sums_image();

      energy.coulomb[1]-= 0.5*p3m_calc_kspace_forces(0,1); 
    }
    energy.coulomb[2] = ELC_energy();
    break;
#endif
  case COULOMB_MMM2D:
    *energy.coulomb += MMM2D_far_energy();
    *energy.coulomb += MMM2D_dielectric_layers_energy_contribution();
    break;
    /* calculate electric part of energy (only for MAGGS) */
  case COULOMB_MAGGS:
    *energy.coulomb += maggs_electric_energy();
    break;
  }
#endif  /* ifdef ELECTROSTATICS */

#ifdef DIPOLES
  switch (coulomb.Dmethod) {
#ifdef DP3M
  case DIPOLAR_P3M:
    dp3m_dipole_assign(); 
    energy.dipolar[1] = dp3m_calc_kspace_forces(0,1);
    break;
  case DIPOLAR_MDLC_P3M:
    dp3m_dipole_assign(); 
    energy.dipolar[1] = dp3m_calc_kspace_forces(0,1);
    energy.dipolar[2] = add_mdlc_energy_corrections();
    break;
#endif
  case DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA:
    energy.dipolar[1] = dawaanr_calculations(0,1);
    break;
  case DIPOLAR_MDLC_DS:
    energy.dipolar[1] = magnetic_dipolar_direct_sum_calculations(0,1);
    energy.dipolar[2] = add_mdlc_energy_corrections();
    break;
  case DIPOLAR_DS:
    energy.dipolar[1] = magnetic_dipolar_direct_sum_calculations(0,1);
    break;
  
  } 
#endif /* ifdef DIPOLES */

}

/************************************************************/

void init_energies(Observable_stat *stat)
{
    int n_pre, n_non_bonded, n_coulomb, n_dipolar;

  n_pre        = 1;
  n_non_bonded = (n_particle_types*(n_particle_types+1))/2;

  n_coulomb    = 0;
#ifdef ELECTROSTATICS
  switch (coulomb.method) {
  case COULOMB_NONE:  n_coulomb = 0; break;
#ifdef P3M
  case COULOMB_ELC_P3M: n_coulomb = 3; break;
  case COULOMB_P3M_GPU:
  case COULOMB_P3M:   n_coulomb = 2; break;
#endif
  default: n_coulomb  = 1;
  }
#endif

  n_dipolar    = 0;
#ifdef DIPOLES

  switch (coulomb.Dmethod) {
  case DIPOLAR_NONE:  n_dipolar = 1; break;
#ifdef DP3M
  case DIPOLAR_MDLC_P3M: n_dipolar=3; break;
  case DIPOLAR_P3M:   n_dipolar = 2; break;
#endif
  case DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA:   n_dipolar = 2; break;
 case DIPOLAR_MDLC_DS: n_dipolar=3; break;
 case DIPOLAR_DS:   n_dipolar = 2; break;
  }

#endif
  
  obsstat_realloc_and_clear(stat, n_pre, n_bonded_ia, n_non_bonded, n_coulomb, n_dipolar, 0, 1);
  stat->init_status = 0;

  external_potential_init_energies();
}

/************************************************************/

void master_energy_calc() {
  mpi_gather_stats(1, total_energy.data.e, NULL, NULL, NULL);

  total_energy.init_status=1;
}
