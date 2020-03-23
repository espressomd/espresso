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
IF CUDA:
    from .lb cimport HydrodynamicInteraction
    from .lb cimport lb_lbfluid_print_vtk_boundary
    from .lb cimport python_lbnode_get_stress
    from .lb cimport lb_lbnode_is_index_valid
    from .lb cimport lb_lbfluid_set_lattice_switch
    from .lb cimport GPU
from . import utils
import tempfile
import shutil
from .utils import is_valid_type
from .utils cimport Vector3i, Vector6d, handle_errors
import numpy as np

IF ELECTROKINETICS:
    cdef class Electrokinetics(HydrodynamicInteraction):
        """
        Creates the electrokinetic method using the GPU unit.

        """

        def __getitem__(self, key):
            if isinstance(key, tuple) or isinstance(
                    key, list) or isinstance(key, np.ndarray):
                if len(key) == 3:
                    return ElectrokineticsRoutines(np.array(key))
            else:
                raise Exception(
                    "%s is not a valid key. Should be a point on the nodegrid e.g. ek[0,0,0]," % key)

        def validate_params(self):
            """
            Checks if the parameters for "stencil" and "fluid_coupling" are valid.

            """
            default_params = self.default_params()

            if self._params["stencil"] not in ["linkcentered", "nodecentered"]:
                raise ValueError(
                    "stencil has to be 'linkcentered' or 'nodecentered'.")

            if self._params["fluid_coupling"] not in ["friction", "estatics"]:
                raise ValueError(
                    "fluid_coupling has to be 'friction' or 'estatics'.")

        def valid_keys(self):
            """
            Returns the valid options used for the electrokinetic method.
            """

            return ["agrid", "lb_density", "viscosity", "friction",
                    "bulk_viscosity", "gamma_even", "gamma_odd", "T",
                    "prefactor", "stencil", "advection", "fluid_coupling",
                    "fluctuations", "fluctuation_amplitude", "es_coupling",
                    "species"]

        def required_keys(self):
            """
            Returns the necessary options to initialize the electrokinetic method.

            """
            return ["agrid", "lb_density", "viscosity",
                    "friction", "T", "prefactor"]

        def default_params(self):
            """
            Returns the default parameters.

            """
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
                    "fluid_coupling": "friction",
                    "fluctuations": False,
                    "fluctuation_amplitude": 0.0,
                    "es_coupling": False,
                    "species": []}

        def _get_params_from_es_core(self):
            if ek_parameters.stencil == 0:
                stencil = "linkcentered"
            elif ek_parameters.stencil == 1:
                stencil = "nodecentered"
            else:
                raise Exception("Value of stencil could not be identified.")

            if ek_parameters.fluidcoupling_ideal_contribution:
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
                    "prefactor": ek_parameters.prefactor,
                    "stencil": stencil,
                    "advection": ek_parameters.advection,
                    "fluid_coupling": fluid_coupling,
                    "fluctuations": ek_parameters.fluctuations,
                    "fluctuation_amplitude":
                        ek_parameters.fluctuation_amplitude,
                    "es_coupling": ek_parameters.es_coupling}

        def _set_params_in_es_core(self):
            if self._params["stencil"] == "linkcentered":
                ek_set_stencil(0)
            elif self._params["stencil"] == "nodecentered":
                ek_set_stencil(1)

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
            ek_set_fluctuations(self._params["fluctuations"])
            ek_set_fluctuation_amplitude(self._params["fluctuation_amplitude"])
            ek_set_electrostatics_coupling(self._params["es_coupling"])

        def set_density(self, species=None, density=None, node=None):
            """
            Sets the density of a species at a specific node.
            If no node is given the density will be set global for the species.

            Parameters
            ----------
            species : :obj:`int`
                species for which the density will apply.
            density : :obj:`float`
                The value to which the density will be set to.
            node : numpy-array of type :obj:`int` of length (3)
                If set the density will be only applied on this specific node.

            """

            if species is None or density is None:
                raise ValueError("species and density has to be set.")
            if not is_valid_type(species, int):
                raise ValueError("species needs to be an integer.")
            if node is None:
                ek_set_density(species, density)
            else:
                if not isinstance(node, (list, np.ndarray)):
                    if len(node) != 3:
                        raise ValueError(
                            "node has to be an array of length 3 of integers.")
                ek_node_set_density(
                    species, node[0], node[1], node[2], density)

        def _activate_method(self):
            self._set_params_in_es_core()
            for species in self._params["species"]:
                species._activate_method()
            lb_lbfluid_set_lattice_switch(GPU)
            self.ek_init()

        def neutralize_system(self, species):
            """
            Sets the global density of a species to a specific value
            for which the whole system will have no net charge.

            .. note :: The previous density of the species will be ignored and
                       it will be homogeneous distributed over the whole system
                       The species must be charged to begin with. If the
                       neutralization would lead to a negative species density
                       an exception will be raised.

            Parameters
            ----------
            species : :obj:`int`
                The species which will be changed to neutralize the system.

            """
            err = ek_neutralize_system(species.id)

            if err == 1:
                raise Exception(
                    'Species used for neutralization must be added to electrokinetics')
            elif err == 2:
                raise Exception(
                    'Species used for neutralization must be charged')
            elif err == 3:
                raise Exception(
                    'Neutralization with specified species would result in negative density')
            elif err != 0:
                raise Exception('Unknown error')

            self.ek_init()

        def ek_init(self):
            """
            Initializes the electrokinetic system.
            This automatically initializes the lattice-Boltzmann method on the GPU.

            """
            err = ek_init()
            if err == 2:
                raise Exception(
                    'EK init failed', 'agrid incompatible with box size')
            elif err != 0:
                raise Exception('EK init failed, unknown error', err)

        def add_species(self, species):
            """
            Initializes a new species for the electrokinetic method.

            Parameters
            ----------
            species : :obj:`int`
                Species to be initialized.

            """
            self._params["species"].append(species)

        def get_params(self):
            """
            Prints out the parameters of the electrokinetic system.

            """
            self._params.update(self._get_params_from_es_core())
            return self._params

        def print_vtk_boundary(self, path):
            """
            Writes the boundary information into a vtk-file.

            Parameters
            ----------
            path : :obj:`str`
                Path of the .vtk file the boundary is written to.

            """
            lb_lbfluid_print_vtk_boundary(utils.to_char_pointer(path))

        def print_vtk_velocity(self, path):
            """
            Writes the lattice-Boltzmann velocity information into a vtk-file.

            Parameters
            ----------
            path : :obj:`str`
                Path of the .vtk file the velocity is written to.

            """
            ek_lb_print_vtk_velocity(utils.to_char_pointer(path))

        def print_vtk_density(self, path):
            """
            Writes the LB density information into a vtk-file.

            Parameters
            ----------
            path : :obj:`str`
                Path of the .vtk file the LB density is written to.

            """
            ek_lb_print_vtk_density(utils.to_char_pointer(path))

        def print_vtk_potential(self, path):
            """
            Writes the electrostatic potential into a vtk-file.

            Parameters
            ----------
            path : :obj:`str`
                Path of the .vtk file the electrostatic potential is written to.

            """
            ek_print_vtk_potential(utils.to_char_pointer(path))

        def print_vtk_lbforce(self, path):
            """
            Writes the LB force information into a vtk-file.

            Parameters
            ----------
            path : :obj:`str`
                Path of the .vtk file the LB force is written to.

            """
            ek_print_vtk_lbforce_density(utils.to_char_pointer(path))

        def print_vtk_particle_potential(self, path):
            """
            Writes the electrostatic particle potential into a vtk-file.

            .. note :: This only works if 'es_coupling' is active.

            Parameters
            ----------
            path : :obj:`str`
                Path of the .vtk file the electrostatic potential is written to.

            """

            if(self._params["es_coupling"]):
                ek_print_vtk_particle_potential(utils.to_char_pointer(path))
            else:
                raise Exception("'es_coupling' is not active.")

        def save_checkpoint(self, path):
            tmp_path = path + ".__tmp__"
            tmpfile_ek = tempfile.NamedTemporaryFile(mode='w+b')
            tmpfile_lb = tempfile.NamedTemporaryFile(mode='w+b')
            ek_save_checkpoint(utils.to_char_pointer(tmpfile_ek.name),
                               utils.to_char_pointer(tmpfile_lb.name))
            ek_path = tmp_path + ".ek"
            lb_path = tmp_path + ".lb"
            shutil.copy(tmpfile_ek.name, path + ".ek")
            shutil.copy(tmpfile_lb.name, path + ".lb")

        def load_checkpoint(self, path):
            self._activate_method()
            ek_load_checkpoint(utils.to_char_pointer(path))

        def add_reaction(self, shape):
            raise Exception("This method is not implemented yet.")

        def add_boundary(self, shape):
            raise Exception("This method is not implemented yet.")

    cdef class ElectrokineticsRoutines:
        cdef Vector3i node

        def __init__(self, key):
            self.node[0] = key[0]
            self.node[1] = key[1]
            self.node[2] = key[2]
            if not lb_lbnode_is_index_valid(self.node):
                raise ValueError("LB node index out of bounds")

        property potential:
            def __get__(self):
                cdef double potential
                ek_node_print_potential(self.node[0], self.node[1], self.node[2], & potential)
                return potential

            def __set__(self, value):
                raise Exception("Potential can not be set.")

        property velocity:
            def __get__(self):
                cdef double velocity[3]
                ek_node_print_velocity(
                    self.node[0], self.node[1], self.node[2], velocity)
                return [velocity[0], velocity[1], velocity[2]]

            def __set__(self, value):
                raise Exception("Not implemented.")

        property pressure:
            def __get__(self):
                cdef Vector6d stress
                stress = python_lbnode_get_stress(self.node)
                return np.array([[stress[0], stress[1], stress[3]],
                                 [stress[1], stress[2], stress[4]],
                                 [stress[3], stress[4], stress[5]]])

            def __set__(self, value):
                raise Exception("Not implemented.")

    class Species:

        """
        Creates a species object that is passed to the ek instance.

        """

        py_number_of_species = 0
        id = -1
        _params = {}

        # __getstate__ and __setstate__ define the pickle interaction
        def __getstate__(self):
            odict = {}
            odict["id"] = self.id
            odict.update(self._params.copy())
            return odict

        def __setstate__(self, params):
            self.id = params["id"]
            params.pop("id")
            self._params = params
            self._set_params_in_es_core()

        def __str__(self):
            return self.__class__.__name__ + "(" + str(self.get_params()) + ")"

        def __getitem__(self, key):
            if isinstance(key, (tuple, list, np.ndarray)):
                if len(key) == 3:
                    return SpecieRoutines(np.array(key), self.id)
            else:
                raise Exception(
                    "%s is not a valid key. Should be a point on the nodegrid e.g. species[0,0,0]," % key)

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
                    raise KeyError("%s is not a valid key" % k)

        def valid_keys(self):
            """
            Returns the valid keys for the species.

            """
            return "density", "D", "valency", "ext_force_density"

        def required_keys(self):
            """
            Returns the required keys for the species.

            """
            return ["density", "D", "valency"]

        def default_params(self):
            """
            Returns the default parameters for the species.

            """
            return {"ext_force_density": [0, 0, 0]}

        def _get_params_from_es_core(self):
            return {
                "density": ek_parameters.density[
                    ek_parameters.species_index[self.id]],
                "D": ek_parameters.D[ek_parameters.species_index[self.id]],
                "valency": ek_parameters.valency[
                    ek_parameters.species_index[self.id]],
                "ext_force_density":
                    [ek_parameters.ext_force_density[0][ek_parameters.species_index[self.id]],
                     ek_parameters.ext_force_density[1][ek_parameters.species_index[self.id]],
                     ek_parameters.ext_force_density[2][ek_parameters.species_index[self.id]]]}

        def _set_params_in_es_core(self):
            ek_set_D(self.id, self._params["D"])
            ek_set_valency(self.id, self._params["valency"])
            ek_set_density(self.id, self._params["density"])
            ek_set_ext_force_density(self.id,
                                     self._params["ext_force_density"][0],
                                     self._params["ext_force_density"][1],
                                     self._params["ext_force_density"][2])

        def _activate_method(self):
            self._set_params_in_es_core()

        def get_params(self):
            """
            Returns the parameters of the species.

            """
            self._params.update(self._get_params_from_es_core())
            return self._params

        def print_vtk_density(self, path):
            """
            Writes the species density into a vtk-file.

            Parameters
            ----------
            path : :obj:`str`
                Path of the .vtk file the species density is written to.

            """
            ek_print_vtk_density(self.id, utils.to_char_pointer(path))

        def print_vtk_flux(self, path):
            """
            Writes the species flux into a vtk-file.

            Parameters
            ----------
            path : :obj:`str`
                Path of the .vtk file the species flux is written to.

            """
            ek_print_vtk_flux(self.id, utils.to_char_pointer(path))

        def print_vtk_flux_fluc(self, path):
            ek_print_vtk_flux_fluc(self.id, utils.to_char_pointer(path))

        def print_vtk_flux_link(self, path):
            ek_print_vtk_flux_link(self.id, utils.to_char_pointer(path))

    cdef class SpecieRoutines:
        cdef Vector3i node
        cdef int id

        def __init__(self, key, id):
            self.node[0] = key[0]
            self.node[1] = key[1]
            self.node[2] = key[2]
            self.id = id
            if not lb_lbnode_is_index_valid(self.node):
                raise ValueError("LB node index out of bounds")

        property density:
            def __set__(self, value):
                if is_valid_type(value, float) or is_valid_type(value, int):
                    if ek_node_set_density(
                            self.id, self.node[0], self.node[1], self.node[2], value) != 0:
                        raise Exception("Species has not been added to EK.")
                else:
                    raise ValueError(
                        "Type of property is wrong. Expected: float.")

            def __get__(self):
                cdef double density
                if ek_node_print_density(self.id, self.node[0], self.node[1], self.node[2], & density) != 0:
                    raise Exception("Species has not been added to EK.")
                return density

        property flux:
            def __set__(self, value):
                raise ValueError("Node flux is not settable.")

            def __get__(self):
                cdef double flux[3]
                if ek_node_print_flux(
                        self.id, self.node[0], self.node[1], self.node[2], flux) != 0:
                    raise Exception("Species has not been added to EK.")

                return np.array(flux[0], flux[1], flux[2])
