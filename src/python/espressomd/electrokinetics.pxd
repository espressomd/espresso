# Copyright (C) 2010-2019 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
include "myconfig.pxi"
from libcpp cimport bool

IF ELECTROKINETICS and CUDA:
    cdef extern from "grid_based_algorithms/electrokinetics.hpp":

        DEF MAX_NUMBER_OF_SPECIES = 10

        # EK data struct
        IF EK_DEBUG:
            ctypedef struct EKParameters:
                float agrid
                float time_step
                float lb_density
                unsigned int dim_x
                unsigned int dim_y
                unsigned int dim_z
                unsigned int number_of_nodes
                float viscosity
                float bulk_viscosity
                float gamma_odd
                float gamma_even
                float friction
                float T
                float prefactor
                float lb_ext_force_density[3]
                unsigned int number_of_species
                int reaction_species[3]
                float rho_reactant_reservoir
                float rho_product0_reservoir
                float rho_product1_reservoir
                float reaction_ct_rate
                float reaction_fraction_0
                float reaction_fraction_1
                float mass_reactant
                float mass_product0
                float mass_product1
                int stencil
                int number_of_boundary_nodes
                float fluctuation_amplitude
                bool fluctuations
                bool advection
                bool fluidcoupling_ideal_contribution
                float * charge_potential
                float * j
                float * lb_force_density_previous
                float * j_fluc
                float * rho[MAX_NUMBER_OF_SPECIES]
                int species_index[MAX_NUMBER_OF_SPECIES]
                float density[MAX_NUMBER_OF_SPECIES]
                float D[MAX_NUMBER_OF_SPECIES]
                float d[MAX_NUMBER_OF_SPECIES]
                float valency[MAX_NUMBER_OF_SPECIES]
                float ext_force_density[3][MAX_NUMBER_OF_SPECIES]
                char * node_is_catalyst
                bool es_coupling
                float * charge_potential_buffer
                float * electric_field
        ELSE:
            ctypedef struct EKParameters:
                float agrid
                float time_step
                float lb_density
                unsigned int dim_x
                unsigned int dim_y
                unsigned int dim_z
                unsigned int number_of_nodes
                float viscosity
                float bulk_viscosity
                float gamma_odd
                float gamma_even
                float friction
                float T
                float prefactor
                float lb_ext_force_density[3]
                unsigned int number_of_species
                int reaction_species[3]
                float rho_reactant_reservoir
                float rho_product0_reservoir
                float rho_product1_reservoir
                float reaction_ct_rate
                float reaction_fraction_0
                float reaction_fraction_1
                float mass_reactant
                float mass_product0
                float mass_product1
                int stencil
                int number_of_boundary_nodes
                float fluctuation_amplitude
                bool fluctuations
                bool advection
                bool fluidcoupling_ideal_contribution
                float * charge_potential
                float * j
                float * lb_force_density_previous
                float * rho[MAX_NUMBER_OF_SPECIES]
                int species_index[MAX_NUMBER_OF_SPECIES]
                float density[MAX_NUMBER_OF_SPECIES]
                float D[MAX_NUMBER_OF_SPECIES]
                float d[MAX_NUMBER_OF_SPECIES]
                float valency[MAX_NUMBER_OF_SPECIES]
                float ext_force_density[3][MAX_NUMBER_OF_SPECIES]
                char * node_is_catalyst
                bool es_coupling
                float * charge_potential_buffer
                float * electric_field

        cdef extern EKParameters ek_parameters

        # EK functions
        void ek_print_parameters()
        void ek_print_lbpar()
        unsigned int ek_calculate_boundary_mass()
        int ek_print_vtk_density(int species, char * filename)
        int ek_print_vtk_flux(int species, char * filename)
        int ek_print_vtk_flux_fluc(int species, char * filename)
        int ek_print_vtk_flux_link(int species, char * filename)
        int ek_print_vtk_potential(char * filename)
        int ek_print_vtk_lbforce_density(char * filename)
        int ek_lb_print_vtk_density(char * filename)
        int ek_lb_print_vtk_velocity(char * filename)
        int ek_init()
        void ek_set_agrid(float agrid) except +
        void ek_set_lb_density(float lb_density) except +
        void ek_set_viscosity(float viscosity) except +
        void ek_set_friction(float friction) except +
        void ek_set_lb_ext_force_density(float lb_ext_force_dens_x, float lb_ext_force_dens_y, float lb_ext_force_dens_z) except +
        void ek_set_T(float T) except +
        void ek_set_prefactor(float prefactor) except +
        void ek_set_bulk_viscosity(float bulk_viscosity) except +
        void ek_set_gamma_odd(float gamma_odd) except +
        void ek_set_gamma_even(float gamma_even) except +
        void ek_set_density(int species, float density)
        void ek_set_D(int species, float D)
        void ek_set_valency(int species, float valency)
        void ek_set_ext_force_density(int species, float ext_force_density_x, float ext_force_density_y, float ext_force_density_z)
        void ek_set_stencil(int stencil) except +
        void ek_set_advection(bool advection) except +
        void ek_set_fluctuations(bool fluctuations) except +
        void ek_set_fluctuation_amplitude(float fluctuation_amplitude) except +
        void ek_set_fluidcoupling(bool ideal_contribution) except +
        void ek_set_electrostatics_coupling(bool electrostatics_coupling) except +
        int ek_node_get_density(int species, int x, int y, int z, double * density)
        int ek_node_get_flux(int species, int x, int y, int z, double * flux)
        int ek_node_get_potential(int x, int y, int z, double * potential)
        int ek_node_set_density(int species, int x, int y, int z, double density)
        float ek_calculate_net_charge()
        int ek_neutralize_system(int species)

        int ek_print_vtk_particle_potential(char * filename)
