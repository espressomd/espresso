include "myconfig.pxi"
from libcpp cimport bool

IF ELECTROKINETICS and CUDA:
    cdef extern from "electrokinetics.hpp":

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
            float lb_force[3]
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
            bool advection
            bool fluidcoupling_ideal_contribution
            float* charge_potential
            ekfloat* j
            float* lb_force_previous
            ekfloat* rho[MAX_NUMBER_OF_SPECIES]
            int species_index[MAX_NUMBER_OF_SPECIES]
            float density[MAX_NUMBER_OF_SPECIES]
            float D[MAX_NUMBER_OF_SPECIES]
            float d[MAX_NUMBER_OF_SPECIES]
            float valency[MAX_NUMBER_OF_SPECIES]
            float ext_force[3][MAX_NUMBER_OF_SPECIES]
            char* node_is_catalyst
            # IF EK_ELECTROSTATIC_COUPLING:
            #     bool es_coupling
            #     float *charge_potential_buffer
            #     float *electric_field

        cdef extern EK_parameters ek_parameters


        # EK functions
        void ek_integrate()
        void ek_print_parameters()
        void ek_print_lbpar()
        void lb_set_ek_pointer(EK_parameters* pointeradress)
        unsigned int ek_calculate_boundary_mass()
        int ek_print_vtk_density(int species, char* filename)
        int ek_print_vtk_flux(int species, char* filename)
        int ek_print_vtk_potential(char* filename)
        int ek_print_vtk_lbforce(char* filename)
        int ek_lb_print_vtk_density(char* filename)
        int ek_lb_print_vtk_velocity(char* filename)
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
        int ek_set_lb_force(double* ext_force)
        int ek_set_density(int species, double density)
        int ek_set_D(int species, double D)
        int ek_set_valency(int species, double valency)
        int ek_set_ext_force(int species, double ext_force_x, double ext_force_y, double ext_force_z)
        int ek_set_stencil(int stencil)
        int ek_set_advection(bool advection)
        int ek_set_fluidcoupling(bool ideal_contribution)
        int ek_node_print_velocity(int x, int y, int z, double* velocity)
        int ek_node_print_density(int species, int x, int y, int z, double* density)
        int ek_node_print_flux(int species, int x, int y, int z, double* flux)
        int ek_node_print_potential(int x, int y, int z, double* potential)
        int ek_node_set_density(int species, int x, int y, int z, double density)
        ekfloat ek_calculate_net_charge() 
        int ek_neutralize_system(int species) 
        int ek_save_checkpoint(char* filename)
        int ek_load_checkpoint(char* filename)

        IF EK_ELECTROSTATIC_COUPLING:
            int ek_set_electrostatics_coupling( bool electrostatics_coupling )
            void ek_calculate_electrostatic_coupling()
            int ek_print_vtk_particle_potential( char* filename )

        IF EK_BOUNDARIES:
            void ek_init_species_density_wallcharge(ekfloat* wallcharge_species_density, int wallcharge_species)


    cdef extern from "lb.hpp":
        int lb_lbfluid_print_vtk_boundary(char* filename)
        int lb_lbnode_get_pi(int* ind, double* p_pi)


    # def init_species(id, density, D, valency, ext_force):
    #     if ek_set_density(id, density):
    #         raise Exception('EK species init error', 'could not set density')

    #     if ek_set_D(id, D):
    #         raise Exception('EK species init error', 'could not set D')

    #     if ek_set_valency(id, valency):
    #         raise Exception('EK species init error', 'could not set valency')

    #     if ek_set_ext_force(id, ext_force[0], ext_force[1], ext_force[2]):
    #         raise Exception('EK species init error', 'could not set ext_force')

    #     ek_init_wrapper()

    # def neutralize_system(species_id):
    #     err = ek_neutralize_system(species_id)

    #     if err == 1:
    #         raise Exception(
    #             'EK neutralize_system error', 'Species used for neutralization must exist')
    #     elif err == 2:
    #         raise Exception(
    #             'EK neutralize_system error', 'Species used for neutralization must be charged')
    #     elif err == 3:
    #         raise Exception('EK neutralize_system error',
    #                         'Neutralization with specified species would result in negative density')
    #     elif err != 0:
    #         raise Exception(
    #             'EK neutralize_system error', 'Unknown error in EK neutralize_system')

    #     ek_init_wrapper()

    # def print_parameters():
    #     ek_print_parameters()

    # def print_lbpar():
    #     ek_print_lbpar()

    # def set_electrostatics_coupling(state):
    #     IF EK_ELECTROSTATIC_COUPLING:
    #         ek_set_electrostatics_coupling(state)
    #         ek_init_wrapper()
    #     ELSE:
    #         raise Exception(
    #             'missing feature', 'feature EK_ELECTROSTATICS_COUPLING needs to be enabled')

    # def set_lb_force(force):
    #     cdef double tmp[3]
    #     tmp[0] = force[0]
    #     tmp[1] = force[1]
    #     tmp[2] = force[2]
    #     ek_set_lb_force(tmp)

    #     ek_init_wrapper()

    # def net_charge():
    #     return ek_calculate_net_charge()

    # def save_checkpoint(path):
    #     if ek_save_checkpoint(path):
    #         raise Exception(
    #             'EK checkpointing error', 'could not save checkpoint')

    # def load_checkpoint(path):
    #     if ek_load_checkpoint(path):
    #         raise Exception(
    #             'EK checkpointing error', 'could not load checkpoint')

    # def print_density_vtk(species_id, path):
    #     if ek_print_vtk_density(species_id, path):
    #         raise Exception('EK output error', 'could not save density VTK')

    # def print_flux_vtk(species_id, path):
    #     if ek_print_vtk_flux(species_id, path):
    #         raise Exception('EK output error', 'could not save flux VTK')

    # def print_potential_vtk(path):
    #     if ek_print_vtk_potential(path):
    #         raise Exception('EK output error', 'could not save potential VTK')

    # def print_lb_force_vtk(path):
    #     if ek_print_vtk_lbforce(path):
    #         raise Exception('EK output error', 'could not save lbforce VTK')

    # def print_lb_density_vtk(path):
    #     if ek_lb_print_vtk_density(path):
    #         raise Exception('EK output error', 'could not save lbdensity VTK')

    # def print_velocity_vtk(path):
    #     if ek_lb_print_vtk_velocity(path):
    #         raise Exception('EK output error', 'could not save velocity VTK')
