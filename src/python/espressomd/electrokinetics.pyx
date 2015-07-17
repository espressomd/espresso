from libcpp cimport bool
import numpy as np

include "myconfig.pxi"

IF ELECTROKINETICS == 1:
    cdef extern from "electrokinetics.hpp":
        void ek_print_parameters()
        void ek_print_lbpar()
        unsigned int ek_calculate_boundary_mass()

        int ek_print_vtk_density(int species, char * filename)
        int ek_print_vtk_flux(int species, char * filename)
        int ek_print_vtk_potential(char * filename)

        int ek_print_vtk_lbforce(char * filename)
        int ek_print_vtk_reaction_tags(char * filename)
        int ek_lb_print_vtk_density(char * filename)
        int ek_lb_print_vtk_velocity(char * filename)

        int ek_init()

        int ek_set_agrid(double agrid)
        int ek_set_lb_density(double lb_density)
        int ek_set_viscosity(double viscosity)
        int ek_set_friction(double friction)
        int ek_set_T(double T)
        int ek_set_bjerrumlength(double bjerrumlength)

        IF EK_ELECTROSTATIC_COUPLING == 1:
            int ek_print_vtk_particle_potential(char * filename)
            int ek_set_electrostatics_coupling(bool electrostatics_coupling)

    #    int ek_set_bulk_viscosity(double bulk_viscosity)
    #    int ek_set_gamma_odd(double gamma_odd)
    #    int ek_set_gamma_even(double gamma_even)
        int ek_set_lb_force(double * ext_force)
        int ek_set_density(int species, double density)
        int ek_set_D(int species, double D)
        int ek_set_valency(int species, double valency)
        int ek_set_ext_force(int species, double ext_force_x, double ext_force_y, double ext_force_z)
        int ek_set_stencil(int stencil)

        int ek_node_print_velocity(int x, int y, int z, double * velocity)
        int ek_node_print_density(int species, int x, int y, int z, double * density)
        int ek_node_print_flux(int species, int x, int y, int z, double * flux)
        int ek_node_set_density(int species, int x, int y, int z, double density)

        float ek_calculate_net_charge()

        int ek_neutralize_system(int species)
        int ek_save_checkpoint(char * filename)
        int ek_load_checkpoint(char * filename)

    def ek_init_wrapper():
        err = ek_init()

        if err == 2:
            raise Exception(
                'EK init failed', 'agrid incompatible with box size')
        elif err != 0:
            raise Exception('EK init failed', 'unknown error')

    def init(agrid=1.0, lb_density=1.0, viscosity=1.0, friction=1.0, bjerrum_length=0.7095, T=1.0, stencil='linkcentered'):
        if ek_set_agrid(agrid):
            raise Exception('EK init error', 'could not set agrid')

        if ek_set_lb_density(lb_density):
            raise Exception('EK init error', 'could not set lb_density')

        if ek_set_viscosity(viscosity):
            raise Exception('EK init error', 'could not set viscosity')

        if ek_set_friction(friction):
            raise Exception('EK init error', 'could not set friction')

        if ek_set_bjerrumlength(bjerrum_length):
            raise Exception('EK init error', 'could not set bjerrum_length')

        if ek_set_T(T):
            raise Exception('EK init error', 'could not set T')

        if stencil == 'linkcentered':
            if ek_set_stencil(0):
                raise Exception(
                    'EK init error', 'could not set linkcentered stencil')
        elif stencil == 'nonlinear':
            if ek_set_stencil(1):
                raise Exception(
                    'EK init error', 'could not set nonlinear stencil')
        elif stencil == 'nodecentered':
            if ek_set_stencil(2):
                raise Exception(
                    'EK init error', 'could not set nodecentered stencil')
        else:
            raise Exception('EK init error', 'unknown stencil')

        ek_init_wrapper()

    def init_species(id, density, D, valency, ext_force):
        if ek_set_density(id, density):
            raise Exception('EK species init error', 'could not set density')

        if ek_set_D(id, D):
            raise Exception('EK species init error', 'could not set D')

        if ek_set_valency(id, valency):
            raise Exception('EK species init error', 'could not set valency')

        if ek_set_ext_force(id, ext_force[0], ext_force[1], ext_force[2]):
            raise Exception('EK species init error', 'could not set ext_force')

        ek_init_wrapper()

    def neutralize_system(species_id):
        err = ek_neutralize_system(species_id)

        if err == 1:
            raise Exception(
                'EK neutralize_system error', 'Species used for neutralization must exist')
        elif err == 2:
            raise Exception(
                'EK neutralize_system error', 'Species used for neutralization must be charged')
        elif err == 3:
            raise Exception('EK neutralize_system error',
                            'Neutralization with specified species would result in negative density')
        elif err != 0:
            raise Exception(
                'EK neutralize_system error', 'Unknown error in EK neutralize_system')

        ek_init_wrapper()

    def print_parameters():
        ek_print_parameters()

    def print_lbpar():
        ek_print_lbpar()

    def setElectrostaticsCoupling(state):
        IF EK_ELECTROSTATIC_COUPLING == 1:
            ek_set_electrostatics_coupling(state)
            ek_init_wrapper()
        ELSE:
            raise Exception(
                'missing feature', 'feature EK_ELECTROSTATICS_COUPLING needs to be enabled')

    def setLbForce(force):
        cdef double tmp[3]
        tmp[0] = force[0]
        tmp[1] = force[1]
        tmp[2] = force[2]
        ek_set_lb_force(tmp)

        ek_init_wrapper()

    def netCharge():
        return ek_calculate_net_charge()

    def saveCheckpoint(path):
        if ek_save_checkpoint(path):
            raise Exception(
                'EK checkpointing error', 'could not save checkpoint')

    def loadCheckpoint(path):
        if ek_load_checkpoint(path):
            raise Exception(
                'EK checkpointing error', 'could not load checkpoint')

    def printDensityVTK(species_id, path):
        if ek_print_vtk_density(species_id, path):
            raise Exception('EK output error', 'could not save density VTK')

    def printFluxVTK(species_id, path):
        if ek_print_vtk_flux(species_id, path):
            raise Exception('EK output error', 'could not save flux VTK')

    def printPotentialVTK(path):
        if ek_print_vtk_potential(path):
            raise Exception('EK output error', 'could not save potential VTK')

    def printLbForceVTK(path):
        if ek_print_vtk_lbforce(path):
            raise Exception('EK output error', 'could not save lbforce VTK')

    def printReactionTagsVTK(path):
        if ek_print_vtk_reaction_tags(path):
            raise Exception(
                'EK output error', 'could not save reaction tags VTK')

    def printLbDensityVTK(path):
        if ek_lb_print_vtk_density(path):
            raise Exception('EK output error', 'could not save lbdensity VTK')

    def printVelocityVTK(path):
        if ek_lb_print_vtk_velocity(path):
            raise Exception('EK output error', 'could not save velocity VTK')
