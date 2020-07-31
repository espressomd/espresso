#
# Copyright (C) 2013-2019 The ESPResSo project
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
#
include "myconfig.pxi"
import os
import cython
import itertools
import numpy as np
cimport numpy as np
from libc cimport stdint
from .actors cimport Actor
from . cimport cuda_init
from . import cuda_init
from . import utils
from .utils import array_locked, is_valid_type, to_char_pointer
from .utils cimport Vector3i, Vector3d, Vector6d, Vector19d, make_array_locked
from .globals cimport time_step


IF LB_WALBERLA:

    vtk_registry = {}

    class VTKOutput:
        """
        VTK callback.

        Parameters
        ----------
        vtk_uid: :obj:`str`
            Unique identifier for the VTK callback.
        delta_N: :obj:`int`
            Frequency of writing in LB steps (0 if manual callback).
        """

        def __init__(self, vtk_uid, delta_N):
            self.vtk_uid = vtk_uid
            self.delta_N = delta_N

        def write(self):
            raise NotImplementedError()

        def disable(self):
            raise NotImplementedError()

        def enable(self):
            raise NotImplementedError()

    class VTKOutputAutomatic(VTKOutput):
        """
        Automatic VTK callback. Can be disabled at any time and re-enabled later.

        Please note that the internal VTK counter is no longer incremented
        when the callback is disabled, which means the number of LB steps
        between two frames will not always be an integer multiple of delta_N.
        """

        def __init__(self, *args):
            super().__init__(*args)
            self.enabled = True

        def __repr__(self):
            return "<{}.{}: writes to '{}' every {} LB steps{}>".format(
                self.__class__.__module__, self.__class__.__name__, self.vtk_uid, self.delta_N, '' if self.enabled else ' (disabled)')

        def write(self):
            raise RuntimeError(
                'This is an automatic VTK callback that writes every {} LB time steps, cannot be triggered manually'.format(
                    self.delta_N))

        def enable(self):
            lb_lbfluid_switch_vtk(to_char_pointer(self.vtk_uid), 1)
            self.enabled = True

        def disable(self):
            lb_lbfluid_switch_vtk(to_char_pointer(self.vtk_uid), 0)
            self.enabled = False

    class VTKOutputManual(VTKOutput):
        """
        Manual VTK callback. Can be called at any time to take a snapshot
        of the current state of the LB fluid.
        """

        def __repr__(self):
            return "<{}.{}: writes to '{}' on demand>".format(
                self.__class__.__module__, self.__class__.__name__, self.vtk_uid)

        def write(self):
            lb_lbfluid_write_vtk(to_char_pointer(self.vtk_uid))

        def disable(self):
            raise RuntimeError('Manual VTK callbacks cannot be disabled')

        def enable(self):
            raise RuntimeError('Manual VTK callbacks cannot be enabled')


def _construct(cls, params):
    obj = cls(**params)
    obj._params = params
    return obj


cdef class HydrodynamicInteraction(Actor):
    def _lb_init(self):
        raise Exception(
            "Subclasses of HydrodynamicInteraction must define the _lb_init() method.")

    def __reduce__(self):
        return _construct, (self.__class__, self._params), None

    def __getitem__(self, key):
        if isinstance(key, (tuple, list, np.ndarray)):
            if len(key) == 3:
                return LBFluidRoutines(np.array(key))
        else:
            raise Exception(
                "%s is not a valid key. Should be a point on the nodegrid e.g. lbf[0,0,0]," % key)

    # validate the given parameters on actor initialization
    #
    def validate_params(self):
        default_params = self.default_params()

        utils.check_type_or_throw_except(
            self._params["kT"], 1, float, "kT must be a number")
        if self._params["kT"] > 0. and not self._params["seed"]:
            raise ValueError(
                "seed has to be given if temperature is not 0.")

        if self._params["dens"] == default_params["dens"]:
            raise Exception("LB_FLUID density not set")
        elif not (self._params["dens"] > 0.0 and (is_valid_type(self._params["dens"], float) or is_valid_type(self._params["dens"], int))):
            raise ValueError("Density must be a positive double")

        if self._params["tau"] <= 0.:
            raise ValueError("tau has to be a positive double")

    def valid_keys(self):
        return "agrid", "dens", "ext_force_density", "visc", "tau", "bulk_visc", "gamma_odd", "gamma_even", "kT", "seed"

    def required_keys(self):
        return ["dens", "agrid", "visc", "tau"]

    def default_params(self):
        return {"agrid": -1.0,
                "dens": -1.0,
                "ext_force_density": [0.0, 0.0, 0.0],
                "visc": -1.0,
                "bulk_visc": -1.0,
                "tau": -1.0,
                "seed": 0,
                "kT": 0.}

    def _set_lattice_switch(self):
        raise Exception(
            "Subclasses of HydrodynamicInteraction must define the _set_lattice_switch() method.")

    def _set_params_in_es_core(self):
        pass

    def _get_params_from_es_core(self):
        default_params = self.default_params()
        self._params['agrid'] = self.agrid
        self._params["tau"] = self.tau
        #self._params['dens'] = self.density
        self._params["kT"] = self.kT
        if self._params['kT'] > 0.0:
            self._params['seed'] = self.seed
        self._params['visc'] = self.viscosity
        if not self._params["bulk_visc"] == default_params["bulk_visc"]:
            self._params['bulk_visc'] = self.bulk_viscosity
        self._params['ext_force_density'] = self.ext_force_density

        return self._params

    def get_interpolated_velocity(self, pos):
        """Get LB fluid velocity at specified position.

        Parameters
        ----------
        pos : (3,) array_like of :obj:`float`
              The position at which velocity is requested.

        Returns
        -------
        v : (3,) array_like :obj:`float`
            The LB fluid velocity at ``pos``.

        """
        cdef Vector3d p

        for i in range(3):
            p[i] = pos[i]
        cdef Vector3d v = lb_lbfluid_get_interpolated_velocity(p) * lb_lbfluid_get_lattice_speed()
        return make_array_locked(v)

    def add_force_at_pos(self, pos, force):
        """Adds a force to the fluid at given position

        Parameters
        ----------
        pos : (3,) array_like of :obj:`float`
              The position at which the force will be added.
        force : (3,) array_like of :obj:`float`
              The force vector which will be dirtibuted at the position.

        """
        cdef Vector3d p
        cdef Vector3d f

        for i in range(3):
            p[i] = pos[i]
            f[i] = force[i]
        lb_lbfluid_add_force_at_pos(p, f)

    def save_checkpoint(self, path, binary):
        tmp_path = path + ".__tmp__"
        lb_lbfluid_save_checkpoint(utils.to_char_pointer(tmp_path), binary)
        os.rename(tmp_path, path)

    def load_checkpoint(self, path, binary):
        lb_lbfluid_load_checkpoint(utils.to_char_pointer(path), binary)

    def _activate_method(self):
        raise Exception(
            "Subclasses of HydrodynamicInteraction have to implement _activate_method.")

    def _deactivate_method(self):
        lb_lbfluid_set_lattice_switch(NONE)

    property shape:
        def __get__(self):
            cdef Vector3i shape = lb_lbfluid_get_shape()
            return (shape[0], shape[1], shape[2])

    property kT:
        def __get__(self):
            return lb_lbfluid_get_kT()

    property seed:
        def __get__(self):
            return lb_lbfluid_get_rng_state()

        def __set__(self, seed):
            cdef stdint.uint64_t _seed = seed
            lb_lbfluid_set_rng_state(seed)

    property pressure_tensor:
        def __get__(self):
            cdef Vector6d tensor = python_lbfluid_get_pressure_tensor(self.agrid, self.tau)
            return array_locked(np.array([[tensor[0], tensor[1], tensor[3]],
                                          [tensor[1], tensor[2], tensor[4]],
                                          [tensor[3], tensor[4], tensor[5]]]))

        def __set__(self, value):
            raise NotImplementedError

    property ext_force_density:
        def __get__(self):
            cdef Vector3d res
            res = python_lbfluid_get_ext_force_density(
                self.agrid, self.tau)
            return make_array_locked(res)

        def __set__(self, ext_force_density):
            python_lbfluid_set_ext_force_density(
                ext_force_density, self.agrid, self.tau)

    property viscosity:
        def __get__(self):
            return python_lbfluid_get_viscosity(self.agrid, self.tau)

    property tau:
        def __get__(self):
            return lb_lbfluid_get_tau()

    property agrid:
        def __get__(self):
            return lb_lbfluid_get_agrid()

    def nodes(self):
        """Provides a generator for iterating over all lb nodes"""

        shape = self.shape
        for i, j, k in itertools.product(
                range(shape[0]), range(shape[1]), range(shape[2])):
            yield self[i, j, k]


IF LB_WALBERLA:
    cdef class LBFluidWalberla(HydrodynamicInteraction):
        """
        Initialize the lattice-Boltzmann method for hydrodynamic flow using Walberla

        """

        def _set_params_in_es_core(self):
            pass

        def default_params(self):
            return {"agrid": -1.0,
                    "dens": -1.0,
                    "ext_force_density": [0.0, 0.0, 0.0],
                    "visc": -1.0,
                    "bulk_visc": -1.0,
                    "tau": -1.0,
                    "kT": 0.0,
                    "seed": 0}

        def _set_lattice_switch(self):
            raise Exception("This may not be called")

        def _activate_method(self):
            self.validate_params()
            mpi_init_lb_walberla(
                self._params["visc"] * self._params['tau'] /
                self._params['agrid']**2, self._params["dens"], self._params["agrid"], self._params["tau"],
                self._params['kT'] / (self._params['agrid']
                                      ** 2 / self._params['tau']**2),
                self._params['seed'])
            utils.handle_errors("LB fluid activation")
            self.ext_force_density = self._params["ext_force_density"]

        def _deactivate_method(self):
            mpi_destruct_lb_walberla()
            super()._deactivate_method()

        def add_vtk_writer(self, identifier, observables, delta_N=0,
                           base_folder='vtk_out', prefix='simulation_step'):
            """
            Create a VTK observable.

            Files are written to ``<base_folder>/<identifier>/<prefix>_*.vtu``.
            Summary is written to ``<base_folder>/<identifier>.pvd``.

            Parameters
            ----------
            identifier : :obj:`str`
                Name of the VTK dataset.
            observables : :obj:`list`, \{'density', 'velocity_vector', 'pressure_tensor'\}
                List of observables to write to the VTK files.
            delta_N : :obj:`int`
                Write frequency, if 0 write a single frame (default),
                otherwise add a callback to write every ``delta_N`` LB steps
                to a new file.
            base_folder : :obj:`str`
                Path to the output VTK folder.
            prefix : :obj:`str`
                Prefix for VTK files.
            """
            # sanity checks
            utils.check_type_or_throw_except(
                delta_N, 1, int, "delta_N must be 1 integer")
            assert delta_N >= 0, 'delta_N must be a positive integer'
            utils.check_type_or_throw_except(
                identifier, 1, str, "identifier must be 1 string")
            assert os.path.sep not in identifier, 'identifier must be a string, not a filepath'
            vtk_uid = base_folder + '/' + identifier
            vtk_path = os.path.abspath(vtk_uid)
            if vtk_path in vtk_registry:
                raise RuntimeError(
                    'VTK identifier "{}" would overwrite files written by VTK identifier "{}"'.format(
                        vtk_path, vtk_registry[vtk_path]))
            vtk_registry[vtk_path] = vtk_uid
            if isinstance(observables, str):
                observables = [observables]
            # construct VTK callback
            observable2enum = {
                'density': < int > output_vtk_density,
                'velocity_vector': < int > output_vtk_velocity_vector,
                'pressure_tensor': < int > output_vtk_pressure_tensor}
            flag = 0
            for obs in set(observables):
                if obs not in observable2enum:
                    raise ValueError('Unknown VTK observable ' + obs)
                flag += observable2enum[obs]
            lb_lbfluid_create_vtk(delta_N, flag, to_char_pointer(identifier),
                                  to_char_pointer(base_folder),
                                  to_char_pointer(prefix))
            if delta_N:
                return VTKOutputAutomatic(vtk_uid, delta_N)
            else:
                return VTKOutputManual(vtk_uid, delta_N)


cdef class LBFluidRoutines:
    cdef Vector3i node

    def __init__(self, key):
        utils.check_type_or_throw_except(
            key, 3, int, "The index of an lb fluid node consists of three integers.")
        self.node[0] = key[0]
        self.node[1] = key[1]
        self.node[2] = key[2]
        if not lb_lbnode_is_index_valid(self.node):
            raise ValueError("LB node index out of bounds")

    property index:
        def __get__(self):
            return (self.node[0], self.node[1], self.node[2])

    property velocity:
        def __get__(self):
            return make_array_locked(python_lbnode_get_velocity(self.node))

        def __set__(self, value):
            cdef Vector3d c_velocity
            if all(is_valid_type(v, float) for v in value) and len(value) == 3:
                c_velocity[0] = value[0]
                c_velocity[1] = value[1]
                c_velocity[2] = value[2]
                python_lbnode_set_velocity(self.node, c_velocity)
            else:
                raise ValueError(
                    "Velocity has to be of shape 3 and type float.")
    property density:
        def __get__(self):
            return python_lbnode_get_density(self.node)

        def __set__(self, value):
            python_lbnode_set_density(self.node, value)

    property pressure_tensor:
        def __get__(self):
            cdef Vector6d tensor = python_lbnode_get_pressure_tensor(self.node)
            return array_locked(np.array([[tensor[0], tensor[1], tensor[3]],
                                          [tensor[1], tensor[2], tensor[4]],
                                          [tensor[3], tensor[4], tensor[5]]]))

        def __set__(self, value):
            raise NotImplementedError

    property population:
        def __get__(self):
            cdef Vector19d double_return
            double_return = lb_lbnode_get_pop(self.node)
            return array_locked(np.array([double_return[0],
                                          double_return[1],
                                          double_return[2],
                                          double_return[3],
                                          double_return[4],
                                          double_return[5],
                                          double_return[6],
                                          double_return[7],
                                          double_return[8],
                                          double_return[9],
                                          double_return[10],
                                          double_return[11],
                                          double_return[12],
                                          double_return[13],
                                          double_return[14],
                                          double_return[15],
                                          double_return[16],
                                          double_return[17],
                                          double_return[18]]
                                         ))

        def __set__(self, population):
            cdef Vector19d _population
            for i in range(19):
                _population[i] = population[i]
            lb_lbnode_set_pop(self.node, _population)

    property is_boundary:
        def __get__(self):
            return lb_lbnode_is_boundary(self.node)

        def __set__(self, value):
            raise NotImplementedError

    property last_applied_force:
        def __get__(self):
            return make_array_locked(
                python_lbnode_get_last_applied_force(self.node))
