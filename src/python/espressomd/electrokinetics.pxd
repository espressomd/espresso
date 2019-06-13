include "myconfig.pxi"
from libcpp cimport bool

from .utils cimport Vector3i

IF ELECTROKINETICS and CUDA:
    cdef extern from "grid_based_algorithms/electrokinetics.hpp":

        IF EK_DOUBLE_PREC:
            ctypedef double ekfloat
        ELSE:
            ctypedef float ekfloat

        DEF MAX_NUMBER_OF_SPECIES = 10

        #EK data struct
        ctypedef struct EK_parameters:
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
            float lb_force_density[3]
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
            ekfloat * j
            float * lb_force_density_previous
            ekfloat * rho[MAX_NUMBER_OF_SPECIES]
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

        cdef extern EK_parameters ek_parameters

        # EK functions
        void ek_integrate()
        void ek_print_parameters()
        void ek_print_lbpar()
        void lb_set_ek_pointer(EK_parameters * pointeradress)
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
        int ek_set_agrid(double agrid)
        int ek_set_lb_density(double lb_density)
        int ek_set_viscosity(double viscosity)
        int ek_set_friction(double friction)
        int ek_set_T(double T)
        int ek_set_prefactor(double prefactor)
        int ek_set_bulk_viscosity(double bulk_viscosity)
        int ek_set_gamma_odd(double gamma_odd)
        int ek_set_gamma_even(double gamma_even)
        int ek_set_lb_force_density(double * ext_force_density)
        int ek_set_density(int species, double density)
        int ek_set_D(int species, double D)
        int ek_set_valency(int species, double valency)
        int ek_set_ext_force_density(int species, double ext_force_density_x, double ext_force_density_y, double ext_force_density_z)
        int ek_set_stencil(int stencil)
        int ek_set_advection(bool advection)
        int ek_set_fluctuations(bool fluctuations)
        int ek_set_fluctuation_amplitude(float fluctuation_amplitude)
        int ek_set_fluidcoupling(bool ideal_contribution)
        int ek_node_print_velocity(int x, int y, int z, double * velocity)
        int ek_node_print_density(int species, int x, int y, int z, double * density)
        int ek_node_print_flux(int species, int x, int y, int z, double * flux)
        int ek_node_print_potential(int x, int y, int z, double * potential)
        int ek_node_set_density(int species, int x, int y, int z, double density)
        ekfloat ek_calculate_net_charge()
        int ek_neutralize_system(int species)
        int ek_save_checkpoint(char * filename, char * lb_filename)
        int ek_load_checkpoint(char * filename)

        int ek_set_electrostatics_coupling(bool electrostatics_coupling)
        void ek_calculate_electrostatic_coupling()
        int ek_print_vtk_particle_potential(char * filename)

        IF EK_BOUNDARIES:
            void ek_init_species_density_wallcharge(ekfloat * wallcharge_species_density, int wallcharge_species)
