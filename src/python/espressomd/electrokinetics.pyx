from __future__ import print_function, absolute_import
include "myconfig.pxi"
from .lb cimport HydrodynamicInteraction
from .ekboundaries import EKBoundary
from . import utils
import numpy as np
from espressomd.utils import is_valid_type

IF ELECTROKINETICS:
    cdef class Electrokinetics(HydrodynamicInteraction):
        species_list = []

        def __getitem__(self, key):
            if isinstance(key, tuple) or isinstance(key, list) or isinstance(key, np.ndarray):
                if len(key) == 3:
                    return ElectrokineticsRoutines(np.array(key))
            else: 
                raise Exception("%s is not a valid key. Should be a point on the nodegrid e.g. ek[0,0,0]," %key)


        def validate_params(self):
            default_params = self.default_params()

            if not (self._params["stencil"] in ["linkcentered", "nonlinear", "nodecentered"]):
                raise ValueError("stencil has to be 'linkcentered', 'nonlinear' or 'nodecentered'.")

            if not (self._params["fluid_coupling"] in ["friction", "estatics"]):
                raise ValueError("fluid_coupling has to be 'friction' or 'estatics'.")


        def valid_keys(self):
            return "agrid", "lb_density", "viscosity", "friction", "bulk_viscosity", "gamma_even", "gamma_odd", "T", "prefactor", "stencil", "advection", "fluid_coupling"

        def required_keys(self):
            return ["agrid", "lb_density", "viscosity", "friction", "T", "prefactor"]

        def default_params(self):
            return {"agrid": -1,
                    "lb_density": -1,
                    "viscosity": -1,
                    "bulk_viscosity": -1,
                    "gamma_odd": 0.0,
                    "gamma_even": 0.0,
                    "friction": 0.0,
                    "T": -1,
                    "prefactor": -1,
                    "stencil": "linkcentered",
                    "advection": True,
                    "fluid_coupling": "friction"}

        def _get_params_from_es_core(self):
            if ek_parameters.stencil == 0:
                stencil = "linkcentered"
            elif ek_parameters.stencil == 1:
                stencil = "nonlinear"
            elif ek_parameters.stencil == 2:
                stencil = "nodecentered"
            else: 
                raise Exception("Value of stencil could not be identified.")

            if ek_parameters.fluidcoupling_ideal_contribution == True:
                fluid_coupling = "friction"
            else:
                fluid_coupling = "estatics"

            return {"agrid": ek_parameters.agrid,
                    "lb_density": ek_parameters.lb_density,
                    "viscosity": ek_parameters.viscosity,
                    "bulk_viscosity": ek_parameters.bulk_viscosity,
                    "gamma_odd": ek_parameters.gamma_odd,
                    "gamma_even": ek_parameters.gamma_even,
                    "friction": ek_parameters.friction,
                    "T": ek_parameters.T,
                    "prefactor":ek_parameters.prefactor,
                    "stencil": stencil,
                    "advection": ek_parameters.advection,
                    "fluid_coupling": fluid_coupling}


        def _set_params_in_es_core(self):
            if self._params["stencil"] == "linkcentered":
                ek_set_stencil(0)
            elif self._params["stencil"] == "nonlinear":
                ek_set_stencil(1)
            elif self._params["stencil"] == "nodecentered":
                ek_set_stencil(2)

            if self._params["fluid_coupling"] == "friction":
                ek_set_fluidcoupling(True)
            elif self._params["fluid_coupling"] == "estatics":
                ek_set_fluidcoupling(False)

            ek_set_agrid(self._params["agrid"])
            ek_set_lb_density(self._params["lb_density"])
            ek_set_viscosity(self._params["viscosity"])
            ek_set_friction(self._params["friction"])
            ek_set_T(self._params["T"])
            ek_set_prefactor(self._params["prefactor"])
            ek_set_bulk_viscosity(self._params["bulk_viscosity"])
            ek_set_gamma_odd(self._params["gamma_odd"])
            ek_set_gamma_even(self._params["gamma_even"])
            ek_set_advection(self._params["advection"])


        def set_density(self, species=None, density=None, node=None):
            if species == None or density == None:
                raise ValueError("species and density has to be set.")
            if not is_valid_type(species, int):
                raise ValueError("species needs to be an integer.")
            if node == None:
                ek_set_density(species, density)
            else:
                if not (isinstance(node, list) or isinstance(node, np.ndarray)):
                    if len(node) != 3:
                        raise ValueError("node has to be an array of length 3 of integers.")
                ek_node_set_density(species, node[0], node[1], node[2], density)


        def _activate_method(self):
            self._set_params_in_es_core()
            for species in self.species_list:
                species._activate_method()
            self.ek_init()


        def neutralize_system(self, species):
            err = ek_neutralize_system(species.id)

            if err == 1:
                raise Exception('Species used for neutralization must be added to electrokinetics')
            elif err == 2:
                raise Exception('Species used for neutralization must be charged')
            elif err == 3:
                raise Exception('Neutralization with specified species would result in negative density')
            elif err != 0: 
                raise Exception('Unknown error')

            self.ek_init()


        def ek_init(self):
            err = ek_init()
            if err == 2:
                raise Exception('EK init failed', 'agrid incompatible with box size')
            elif err != 0:
                raise Exception('EK init failed', 'unknown error')


        def add_species(self, species):
            self.species_list.append(species)

        def get_params(self):
            self._params.update(self._get_params_from_es_core())
            return self._params

        def print_vtk_boundary(self, path):
            lb_lbfluid_print_vtk_boundary(utils.to_char_pointer(path))

        def print_vtk_velocity(self, path):
            ek_lb_print_vtk_velocity(utils.to_char_pointer(path))

        def print_vtk_density(self, path):
            ek_lb_print_vtk_density(utils.to_char_pointer(path))

        def print_vtk_potential(self, path):
            ek_print_vtk_potential(utils.to_char_pointer(path))

        def print_vtk_lbforce(self, path):
            ek_print_vtk_lbforce(utils.to_char_pointer(path))

        def print_vtk_particle_potential(self, path):
            IF EK_ELECTROSTATIC_COUPLING:
                ek_print_vtk_particle_potential(utils.to_char_pointer(path))
            ELSE:
                raise Exception("'EK_ELECTROSTATIC_COUPLING' ist not active.")


        # TODO:
        def checkpoint(self):
            raise Exception("Please implement this method in the pickle routine.")

        def add_reaction(self, shape):
            raise Exception("This method is not implemented yet.")

        def add_boundary(self, shape):
            raise Exception("This method is not implemented yet.")


    cdef class ElectrokineticsRoutines(object):
        cdef int node[3]

        def __init__(self, key):
            self.node[0] = key[0]
            self.node[1] = key[1]
            self.node[2] = key[2]

        property potential:
            def __get__(self):
                cdef double potential
                ek_node_print_potential(self.node[0], self.node[1], self.node[2], &potential)
                return potential

            def __set__(self, value):
                raise Exception("Potential can not be set.")

        property velocity:
            def __get__(self):
                cdef double velocity[3]
                ek_node_print_velocity(self.node[0], self.node[1], self.node[2], velocity)
                return [velocity[0], velocity[1], velocity[2]]

            def __set__(self, value):
                raise Exception("Not implemented.")

        property pressure:
            def __get__(self):
                cdef double pi[6]
                lb_lbnode_get_pi(self.node, pi)
                return np.array([[pi[0],pi[1],pi[3]],
                                 [pi[1],pi[2],pi[4]],
                                 [pi[3],pi[4],pi[5]]])

            def __set__(self, value):
                raise Exception("Not implemented.")

    class Species(object):
        """Creates a species object that is passed to the ek instance"""
        py_number_of_species = 0
        id = -1
        _params = {}


        def __getitem__(self, key):
            if isinstance(key, tuple) or isinstance(key, list) or isinstance(key, np.ndarray):
                if len(key) == 3:
                    return SpecieRoutines(np.array(key), self.id)
            else: 
                raise Exception("%s is not a valid key. Should be a point on the nodegrid e.g. species[0,0,0]," %key)

        def __init__(self, **kwargs):
            Species.py_number_of_species += 1
            self.id = Species.py_number_of_species
            self._params = self.default_params()

            # Check if all required keys are given
            for k in self.required_keys():
                if k not in kwargs:
                    raise ValueError(
                        "At least the following keys have to be given as keyword arguments: " + self.required_keys().__str__() + " got " + kwargs.__str__())
                self._params[k] = kwargs[k]

            for k in kwargs:
                if k in self.valid_keys():
                    self._params[k] = kwargs[k]
                else:
                    raise KeyError("%s is not a vaild key" % k)

        def valid_keys(self):
            return "density", "D", "valency", "ext_force"

        def required_keys(self):
            return ["density", "D", "valency"]

        def default_params(self):
            return {"ext_force": [0, 0, 0]}

        def _get_params_from_es_core(self):
            return {"density": ek_parameters.density[ek_parameters.species_index[self.id]],
                    "D": ek_parameters.D[ek_parameters.species_index[self.id]],
                    "valency": ek_parameters.valency[ek_parameters.species_index[self.id]],
                    "ext_force": [ek_parameters.ext_force[0][ek_parameters.species_index[self.id]],
                                  ek_parameters.ext_force[1][ek_parameters.species_index[self.id]],
                                  ek_parameters.ext_force[2][ek_parameters.species_index[self.id]]]}

        def _set_params_in_es_core(self):
            ek_set_D(self.id, self._params["D"])
            ek_set_valency(self.id, self._params["valency"])
            ek_set_density(self.id, self._params["density"])
            ek_set_ext_force(self.id, self._params["ext_force"][0], self._params["ext_force"][1], self._params["ext_force"][2])

        def _activate_method(self):
            self._set_params_in_es_core()

        def get_params(self):
            self._params.update(self._get_params_from_es_core())
            return self._params

        def print_vtk_density(self, path):
            ek_print_vtk_density(self.id, utils.to_char_pointer(path))

        def print_vtk_flux(self, path):
            ek_print_vtk_flux(self.id, utils.to_char_pointer(path))



    cdef class SpecieRoutines(object):
        cdef int node[3]
        cdef int id

        def __init__(self, key, id):
            self.node = key
            self.id = id

        property density:
            def __set__(self, value):
                if is_valid_type(value, float) or is_valid_type(value, int):
                    if ek_node_set_density(self.id, self.node[0], self.node[1], self.node[2], value) != 0:
                        raise Exception("Species has not been added to EK.")

                else:
                    raise ValueError("Type of property is wrong. Expected: float.")

            def __get__(self):
                cdef double density
                if ek_node_print_density(self.id, self.node[0], self.node[1], self.node[2], &density) != 0:
                    raise Exception("Species has not been added to EK.")
                return density

        property flux:
            def __set__(self, value):
                raise ValueError("Node flux is not settable.")

            def __get__(self):
                cdef double flux[3]
                if ek_node_print_flux(self.id, self.node[0], self.node[1], self.node[2], flux) != 0:
                    raise Exception("Species has not been added to EK.")

                return np.array(flux[0], flux[1], flux[2])
