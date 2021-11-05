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
import functools
import numpy as np
cimport numpy as np
from .script_interface import ScriptInterfaceHelper, script_interface_register
from .shapes import Shape
from . import lbboundaries
from . cimport cuda_init
from . import cuda_init
from . import utils
from .utils import is_valid_type, to_char_pointer
from .utils cimport Vector3i
from .utils cimport Vector3d
from .utils cimport Vector6d
from .utils cimport make_array_locked
from .utils cimport make_Vector3d
from .utils cimport create_nparray_from_double_array


IF LB_WALBERLA:

    import lxml.etree

    class VTKRegistry:

        def __init__(self):
            self.map = {}
            self.collisions = {}

        def _register_vtk_object(self, vtk_uid, vtk_obj):
            self.map[vtk_uid] = vtk_obj
            self.collisions[os.path.abspath(vtk_uid)] = vtk_uid

        def __getstate__(self):
            return self.map

        def __setstate__(self, active_vtk_objects):
            self.map = {}
            self.collisions = {}
            for vtk_uid, vtk_obj in active_vtk_objects.items():
                self.map[vtk_uid] = vtk_obj
                self.collisions[os.path.abspath(vtk_uid)] = vtk_uid

    _vtk_registry = VTKRegistry()

    class VTKOutput:
        """
        VTK callback.
        """
        observable2enum = {
            'density': < int > output_vtk_density,
            'velocity_vector': < int > output_vtk_velocity_vector,
            'pressure_tensor': < int > output_vtk_pressure_tensor}

        def __init__(self, vtk_uid, identifier, observables, delta_N,
                     base_folder, prefix):
            observables = set(observables)
            flag = sum(self.observable2enum[obs] for obs in observables)
            self._params = {
                'vtk_uid': vtk_uid, 'identifier': identifier, 'prefix': prefix,
                'delta_N': delta_N, 'base_folder': base_folder, 'flag': flag,
                'enabled': True, 'initial_count': 0, 'observables': observables}
            self._set_params_in_es_core()

        def __getstate__(self):
            odict = self._params.copy()
            # get initial execution counter
            pvd = os.path.abspath(self._params['vtk_uid']) + '.pvd'
            if os.path.isfile(pvd):
                tree = lxml.etree.parse(pvd)
                nodes = tree.xpath('/VTKFile/Collection/DataSet')
                if nodes:
                    odict['initial_count'] = int(
                        nodes[-1].attrib['timestep']) + 1
            return odict

        def __setstate__(self, params):
            self._params = params
            self._set_params_in_es_core()

        def _set_params_in_es_core(self):
            _vtk_registry._register_vtk_object(self._params['vtk_uid'], self)
            lb_lbfluid_create_vtk(self._params['delta_N'],
                                  self._params['initial_count'],
                                  self._params['flag'],
                                  to_char_pointer(self._params['identifier']),
                                  to_char_pointer(self._params['base_folder']),
                                  to_char_pointer(self._params['prefix']))
            if not self._params['enabled']:
                self.disable()

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

        def __repr__(self):
            return "<{}.{}: writes to '{}' every {} LB step{}{}>".format(
                self.__class__.__module__, self.__class__.__name__,
                self._params['vtk_uid'], self._params['delta_N'],
                '' if self._params['delta_N'] == 1 else 's',
                '' if self._params['enabled'] else ' (disabled)')

        def write(self):
            raise RuntimeError(
                'Automatic VTK callbacks writes automaticall, cannot be triggered manually')

        def enable(self):
            lb_lbfluid_switch_vtk(to_char_pointer(self._params['vtk_uid']), 1)
            self._params['enabled'] = True

        def disable(self):
            lb_lbfluid_switch_vtk(to_char_pointer(self._params['vtk_uid']), 0)
            self._params['enabled'] = False

    class VTKOutputManual(VTKOutput):
        """
        Manual VTK callback. Can be called at any time to take a snapshot
        of the current state of the LB fluid.
        """

        def __repr__(self):
            return "<{}.{}: writes to '{}' on demand>".format(
                self.__class__.__module__, self.__class__.__name__,
                self._params['vtk_uid'])

        def write(self):
            lb_lbfluid_write_vtk(to_char_pointer(self._params['vtk_uid']))

        def disable(self):
            raise RuntimeError('Manual VTK callbacks cannot be disabled')

        def enable(self):
            raise RuntimeError('Manual VTK callbacks cannot be enabled')


def _construct(cls, params):
    obj = cls(**params)
    obj._params = params
    return obj


class HydrodynamicInteraction(ScriptInterfaceHelper):
    """
    Base class for LB implementations.

    Parameters
    ----------
    agrid : :obj:`float`
        Lattice constant. The box size in every direction must be an integer
        multiple of ``agrid``.
    tau : :obj:`float`
        LB time step, must be an integer multiple of the MD time step.
    density : :obj:`float`
        Fluid density.
    viscosity : :obj:`float`
        Fluid kinematic viscosity.
    ext_force_density : (3,) array_like of :obj:`float`, optional
        Force density applied on the fluid.
    kT : :obj:`float`, optional
        Thermal energy of the simulated heat bath (for thermalized fluids).
        Set it to 0 for an unthermalized fluid.
    seed : :obj:`int`, optional
        Initial counter value (or seed) of the philox RNG.
        Required for a thermalized fluid. Must be positive.
    """

    def __init__(self, *args, **kwargs):
        if 'sip' not in kwargs:
            if not self.required_keys().issubset(kwargs):
                raise ValueError(
                    f"At least the following keys have to be given as keyword arguments: {self.required_keys()}")
            params = self.default_params()
            params.update(kwargs)
            self.validate_params(params)
            super().__init__(*args, **params)
            utils.handle_errors(
                "HydrodynamicInteraction initialization failed")
            self._params = {k: params.get(k) for k in self.valid_keys()}
        else:
            super().__init__(**kwargs)

    def _lb_init(self):
        raise Exception(
            "Subclasses of HydrodynamicInteraction must define the _lb_init() method.")

    def __reduce__(self):
        return _construct, (self.__class__, self._params), None

    def __getitem__(self, key):
        if isinstance(key, (tuple, list, np.ndarray)):
            if len(key) == 3:
                if any(isinstance(item, slice) for item in key):
                    return LBSlice(key, self.shape)
                else:
                    return LBFluidRoutines(np.array(key))
        else:
            raise Exception(
                "%s is not a valid key. Should be a point on the nodegrid e.g. lbf[0,0,0], or a slice" % key)

    def __str__(self):
        return f"{self.__class__.__name__}({self.get_params()})"

    def _activate(self):
        self._activate_method()

    def _deactivate(self):
        self._deactivate_method()

    def _activate_method(self):
        self._lb_init()
        self.call_method('activate')
        utils.handle_errors("HydrodynamicInteraction activation failed")

    def _deactivate_method(self):
        self.call_method('deactivate')
        utils.handle_errors("HydrodynamicInteraction deactivation failed")

    def get_params(self):
        """Get interaction parameters"""
        # If this instance refers to an actual interaction defined in the es
        # core, load current parameters from there
        if self.is_active:
            update = self._get_params_from_es_core()
            self._params.update(update)
        return self._params

    def validate_params(self, params):
        default_params = self.default_params()

        for key in ('agrid', 'tau', 'density', 'viscosity'):
            utils.check_type_or_throw_except(
                params[key], 1, float, f'{key} must be a number')
            if params[key] <= 0.:
                raise ValueError(f'{key} must be a strictly positive number')

        utils.check_type_or_throw_except(
            params['kT'], 1, float, 'kT must be a number')
        if params['kT'] > 0.:
            utils.check_type_or_throw_except(
                params['seed'], 1, int, 'seed must be a number')
            if params['seed'] < 0:
                raise ValueError(f'seed must be a positive number')
        utils.check_type_or_throw_except(
            params['ext_force_density'], 3, float, 'ext_force_density must be 3 floats')

    def valid_keys(self):
        return {"agrid", "density", "ext_force_density", "viscosity", "tau", "lattice",
                "kT", "seed"}

    def required_keys(self):
        return {"density", "agrid", "viscosity", "tau"}

    def default_params(self):
        return {"agrid": -1.0,
                "density": -1.0,
                "ext_force_density": [0.0, 0.0, 0.0],
                "viscosity": -1.0,
                "tau": -1.0,
                "seed": 0,
                "lattice": None,
                "kT": 0.}

    def mach_limit(self):
        """
        The fluid velocity is limited to :math:`v_{\\mathrm{max}} = 0.20`
        (see *quasi-incompressible limit* in :cite:`kruger17a`,
        chapter 7, page 272), which corresponds to Mach 0.35.

        The relative error in the fluid density between a compressible fluid
        and an incompressible fluid at Mach 0.30 is less than 5% (see
        *constant density assumption* in :cite:`kundu01a` chapter 16, page
        663). Since the speed of sound is :math:`c_s = 1 / \\sqrt{3}` in LB
        velocity units in a D3Q19 lattice, the velocity limit at Mach 0.30
        is :math:`v_{\\mathrm{max}} = 0.30 / \\sqrt{3} \\approx 0.17`.
        At Mach 0.35 the relative error is around 6% and
        :math:`v_{\\mathrm{max}} = 0.35 / \\sqrt{3} \\approx 0.20`.

        Returns
        -------
        v_max : :obj:`float`
            The Mach limit expressed in LB velocity units.

        """
        return 0.20

    @classmethod
    def _check_mach_limit(cls, velocities):
        lattice_vel = lb_lbfluid_get_lattice_speed()
        vel_max = cls.mach_limit(cls)
        velocities = np.reshape(velocities, (-1, 3))
        if np.any(np.linalg.norm(velocities, axis=1) > vel_max * lattice_vel):
            speed_of_sound = 1. / np.sqrt(3.)
            mach_number = vel_max / speed_of_sound
            raise ValueError(f'Slip velocity exceeds Mach {mach_number:.2f}')

    def _set_params_in_es_core(self):
        pass

    def _get_params_from_es_core(self):
        default_params = self.default_params()
        self._params['agrid'] = self.agrid
        self._params["tau"] = self.tau
        self._params['dens'] = self.density
        self._params["kT"] = self.kT
        if self._params['kT'] > 0.0:
            self._params['seed'] = self.seed
        self._params['viscosity'] = self.viscosity
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
        return python_lbnode_get_interpolated_velocity(make_Vector3d(pos))

    def add_force_at_pos(self, pos, force):
        """Adds a force to the fluid at given position

        Parameters
        ----------
        pos : (3,) array_like of :obj:`float`
              The position at which the force will be added.
        force : (3,) array_like of :obj:`float`
              The force vector which will be distributed at the position.

        """
        lb_lbfluid_add_force_at_pos(make_Vector3d(pos), make_Vector3d(force))

    def save_checkpoint(self, path, binary):
        '''
        Write LB node populations to a file.
        :class:`~espressomd.lbboundaries.LBBoundaries`
        information is not written to the file.
        '''
        tmp_path = path + ".__tmp__"
        lb_lbfluid_save_checkpoint(utils.to_char_pointer(tmp_path), binary)
        os.rename(tmp_path, path)

    def load_checkpoint(self, path, binary):
        '''
        Load LB node populations from a file.
        :class:`~espressomd.lbboundaries.LBBoundaries`
        information is not available in the file. The boundary
        information of the grid will be set to zero,
        even if :class:`~espressomd.lbboundaries.LBBoundaries`
        contains :class:`~espressomd.lbboundaries.LBBoundary`
        objects (they are ignored).
        '''
        lb_lbfluid_load_checkpoint(utils.to_char_pointer(path), binary)

    @property
    def pressure_tensor(self):
        tensor = python_lbfluid_get_pressure_tensor(self.agrid, self.tau)
        return utils.array_locked(tensor)

    def nodes(self):
        """Provides a generator for iterating over all lb nodes"""

        shape = self.shape
        for i, j, k in itertools.product(
                range(shape[0]), range(shape[1]), range(shape[2])):
            yield self[i, j, k]


IF LB_WALBERLA:
    @script_interface_register
    class LatticeWalberla(ScriptInterfaceHelper):
        """Interface to the Walberla lattice

        """
        _so_name = "walberla::LatticeWalberla"
        _so_creation_policy = "GLOBAL"

        def __init__(self, *args, **kwargs):
            if 'sip' not in kwargs:
                if not self.required_keys().issubset(kwargs):
                    raise ValueError(
                        f"At least the following keys have to be given as keyword arguments: {self.required_keys()}")
                params = self.default_params()
                params.update(kwargs)
                self.validate_params(params)
                super().__init__(*args, **params)
                utils.handle_errors("LatticeWalberla initialization failed")
                self._params = {k: getattr(self, k) for k in self.valid_keys()}
            else:
                super().__init__(**kwargs)

        def validate_params(self, params):
            utils.check_type_or_throw_except(
                params['agrid'], 1, float, 'agrid must be a real number')
            utils.check_type_or_throw_except(
                params['n_ghost_layers'], 1, int, 'n_ghost_layers must be an integer')
            if params['agrid'] <= 0.:
                raise ValueError('agrid has to be a positive double')
            if params['n_ghost_layers'] < 0:
                raise ValueError('n_ghost_layers has to be a positive integer')

        def valid_keys(self):
            return {"agrid", "n_ghost_layers"}

        def required_keys(self):
            return self.valid_keys()

        def default_params(self):
            return {}

    @script_interface_register
    class LBFluidWalberla(HydrodynamicInteraction):
        """
        Initialize the lattice-Boltzmann method for hydrodynamic flow using waLBerla.
        See :class:`HydrodynamicInteraction` for the list of parameters.
        """
        _so_name = "walberla::FluidWalberla"
        _so_creation_policy = "GLOBAL"

        def _set_params_in_es_core(self):
            pass

        def _lb_init(self):
            if not self.is_initialized:
                if self._params.get('lattice') is None:
                    self._params['lattice'] = LatticeWalberla(
                        agrid=self._params['agrid'], n_ghost_layers=1)
                self.call_method('instantiate', **self._params)
                utils.handle_errors(
                    "HydrodynamicInteraction instantiation failed")

        def clear_boundaries(self):
            """
            Remove boundary conditions.
            """
            lb_lbfluid_clear_boundaries()

        def add_boundary_from_shape(
                self, shape, velocity=np.zeros(3, dtype=float)):
            """
            Set boundary conditions from a shape.

            Parameters
            ----------
            shape : :obj:`espressomd.shapes.Shape`
                Shape to rasterize.
            velocity : (3,) or (L, M, N, 3) array_like of :obj:`float`, optional
                Slip velocity. By default no-slip boundary conditions are used.
                If a vector of 3 values, a uniform slip velocity is used,
                otherwise ``L, M, N`` must be equal to the LB grid dimensions.
            """
            utils.check_type_or_throw_except(
                shape, 1, Shape, "expected an espressomd.shapes.Shape")
            if np.shape(velocity) not in [(3,), tuple(self.shape) + (3,)]:
                raise ValueError(
                    f'Cannot process velocity grid of shape {np.shape(velocity)}')

            # range checks
            velocity = np.array(velocity, dtype=float).reshape((-1, 3))
            self._check_mach_limit(velocity)

            velocity *= 1. / lb_lbfluid_get_lattice_speed()
            velocity_flat = velocity.reshape((-1,))
            mask_flat = shape.call_method('rasterize', grid_size=self.shape,
                                          grid_spacing=self._params['agrid'],
                                          grid_offset=0.5)

            cdef double[::1] velocity_view = np.ascontiguousarray(
                velocity_flat, dtype=np.double)
            cdef int[::1] raster_view = np.ascontiguousarray(
                mask_flat, dtype=np.int32)
            python_lb_lbfluid_update_boundary_from_shape(
                raster_view, velocity_view)

        def add_boundary_from_list(
                self, nodes, velocity=np.zeros(3, dtype=float)):
            """
            Set boundary conditions from a list of node indices.

            Parameters
            ----------
            nodes : (N, 3) array_like of :obj:`int`
                List of node indices to update. If they were originally not
                boundary nodes, they will become boundary nodes.
            velocity : (3,) or (N, 3) array_like of :obj:`float`, optional
                Slip velocity. By default no-slip boundary conditions are used.
                If a vector of 3 values, a uniform slip velocity is used,
                otherwise ``N`` must be identical to the ``N`` of ``nodes``.
            """

            nodes = np.array(nodes, dtype=int)
            velocity = np.array(velocity, dtype=float)
            if len(nodes.shape) != 2 or nodes.shape[1] != 3:
                raise ValueError(
                    f'Cannot process node list of shape {nodes.shape}')
            if velocity.shape == (3,):
                velocity = np.tile(velocity, (nodes.shape[0], 1))
            elif nodes.shape != velocity.shape:
                raise ValueError(
                    f'Node indices and velocities must have the same shape, got {nodes.shape} and {velocity.shape}')

            # range checks
            self._check_mach_limit(velocity)
            if np.any(np.logical_or(np.max(nodes, axis=0) >= self.shape,
                                    np.min(nodes, axis=0) < 0)):
                raise ValueError(
                    f'Node indices must be in the range (0, 0, 0) to {self.shape}')

            velocity *= 1. / lb_lbfluid_get_lattice_speed()
            cdef double[::1] velocity_view = np.ascontiguousarray(
                velocity.reshape((-1,)), dtype=np.double)
            cdef int[::1] nodes_view = np.ascontiguousarray(
                nodes.reshape((-1,)), dtype=np.int32)
            python_lb_lbfluid_update_boundary_from_list(
                nodes_view, velocity_view)

        def get_nodes_in_shape(self, shape):
            """Provides a generator for iterating over all lb nodes inside the given shape"""
            utils.check_type_or_throw_except(
                shape, 1, Shape, "expected a espressomd.shapes.Shape")
            lb_shape = self.shape
            idxs = itertools.product(
                range(lb_shape[0]), range(lb_shape[1]), range(lb_shape[2]))
            for idx in idxs:
                pos = (np.asarray(idx) + 0.5) * self._params['agrid']
                if shape.is_inside(position=pos):
                    yield self[idx]

        def get_shape_bitmask(self, shape):
            """Create a bitmask for the given shape."""
            utils.check_type_or_throw_except(
                shape, 1, Shape, "expected a espressomd.shapes.Shape")
            mask_flat = shape.call_method('rasterize', grid_size=self.shape,
                                          grid_spacing=self._params['agrid'],
                                          grid_offset=0.5)
            return np.reshape(mask_flat, self.shape).astype(type(True))

        # TODO WALBERLA: maybe split this method in 2 methods with clear names
        # like add_vtk_writer_auto_update() and add_vtk_writer_manual()
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
            if isinstance(observables, str):
                observables = [observables]
            for obs in observables:
                if obs not in VTKOutput.observable2enum:
                    raise ValueError('Unknown VTK observable ' + obs)
            vtk_uid = base_folder + '/' + identifier
            vtk_path = os.path.abspath(vtk_uid)
            if vtk_path in _vtk_registry.collisions:
                raise RuntimeError(
                    'VTK identifier "{}" would overwrite files written by VTK identifier "{}"'.format(
                        vtk_uid, _vtk_registry.collisions[vtk_path]))
            args = (
                vtk_uid,
                identifier,
                observables,
                delta_N,
                base_folder,
                prefix)
            if delta_N:
                obj = VTKOutputAutomatic(*args)
            else:
                obj = VTKOutputManual(*args)
            return obj

    IF CUDA:
        class LBFluidWalberlaGPU(HydrodynamicInteraction):
            """
            Initialize the lattice-Boltzmann method for hydrodynamic flow using
            waLBerla for the GPU. See :class:`HydrodynamicInteraction` for the
            list of parameters.

            """

            def __init__(self, *args, **kwargs):
                raise NotImplementedError("Not implemented yet")


cdef class LBFluidRoutines:

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
            return python_lbnode_get_velocity(self.node)

        def __set__(self, value):
            utils.check_type_or_throw_except(
                value, 3, float, "velocity has to be 3 floats")
            python_lbnode_set_velocity(self.node, make_Vector3d(value))

    property boundary:
        def __get__(self):
            """
            Returns
            -------
            :ref:`espressomd.lbboundaries.VelocityBounceBack`
                If the node is a boundary node
            None
                If the node is not a boundary node
            """

            is_boundary = lb_lbnode_is_boundary(self.node)
            if is_boundary:
                vel = python_lbnode_get_velocity_at_boundary(self.node)
                return lbboundaries.VelocityBounceBack(vel)
            return None

        def __set__(self, value):
            """
            Parameters
            ----------
            value : :ref:`espressomd.lbboundaries.VelocityBounceBack` or None
                If value is :ref:`espressomd.lbboundaries.VelocityBounceBack`,
                set the node to be a boundary node with the specified velocity.
                If value is ``None``, the node will become a fluid node.
            """

            if isinstance(value, lbboundaries.VelocityBounceBack):
                HydrodynamicInteraction._check_mach_limit(value.velocity)
                python_lbnode_set_velocity_at_boundary(
                    self.node, make_Vector3d(value.velocity))
            elif value is None:
                lb_lbnode_remove_from_boundary(self.node)
            else:
                raise ValueError(
                    "value must be an instance of lbboundaries.VelocityBounceBack or None")

    property boundary_force:
        def __get__(self):
            return python_lbnode_get_boundary_force(self.node)

        def __set__(self, val):
            raise NotImplementedError(
                "The boundary force can only be read, never set.")

    property density:
        def __get__(self):
            return python_lbnode_get_density(self.node)

        def __set__(self, value):
            python_lbnode_set_density(self.node, value)

    property pressure_tensor:
        def __get__(self):
            tensor = python_lbnode_get_pressure_tensor(self.node)
            return utils.array_locked(tensor)

        def __set__(self, value):
            raise NotImplementedError

    property population:
        def __get__(self):
            cdef vector[double] pop
            pop = lb_lbnode_get_pop(self.node)
            return utils.array_locked(
                create_nparray_from_double_array(pop.data(), pop.size()))

        def __set__(self, population):
            cdef vector[double] _population
            _population.resize(len(population))
            for i in range(len(population)):
                _population[i] = population[i]
            lb_lbnode_set_pop(self.node, _population)

    property is_boundary:
        def __get__(self):
            return lb_lbnode_is_boundary(self.node)

        def __set__(self, value):
            raise NotImplementedError

    property last_applied_force:
        def __get__(self):
            return python_lbnode_get_last_applied_force(self.node)

        def __set__(self, force):
            python_lbnode_set_last_applied_force(
                self.node, make_Vector3d(force))

    def __eq__(self, obj1):
        index_1 = np.array(self.index)
        index_2 = np.array(obj1.index)
        return all(index_1 == index_2)

    def __hash__(self):
        return hash(self.index)


class LBSlice:

    def __init__(self, key, shape):
        self.x_indices, self.y_indices, self.z_indices = self.get_indices(
            key, shape[0], shape[1], shape[2])

    def get_indices(self, key, shape_x, shape_y, shape_z):
        x_indices = np.atleast_1d(np.arange(shape_x)[key[0]])
        y_indices = np.atleast_1d(np.arange(shape_y)[key[1]])
        z_indices = np.atleast_1d(np.arange(shape_z)[key[2]])
        return x_indices, y_indices, z_indices

    def get_values(self, x_indices, y_indices, z_indices, prop_name):
        shape_res = np.shape(
            getattr(LBFluidRoutines(np.array([0, 0, 0])), prop_name))
        res = np.zeros(
            (x_indices.size,
             y_indices.size,
             z_indices.size,
             *shape_res))
        for i, x in enumerate(x_indices):
            for j, y in enumerate(y_indices):
                for k, z in enumerate(z_indices):
                    res[i, j, k] = getattr(LBFluidRoutines(
                        np.array([x, y, z])), prop_name)
        if shape_res == (1,):
            res = np.squeeze(res, axis=-1)
        return utils.array_locked(res)

    def set_values(self, x_indices, y_indices, z_indices, prop_name, value):
        for i, x in enumerate(x_indices):
            for j, y in enumerate(y_indices):
                for k, z in enumerate(z_indices):
                    setattr(LBFluidRoutines(
                        np.array([x, y, z])), prop_name, value[i, j, k])

    def __iter__(self):
        indices = [(x, y, z) for (x, y, z) in itertools.product(
            self.x_indices, self.y_indices, self.z_indices)]
        return (LBFluidRoutines(np.array(index)) for index in indices)


def _add_lb_slice_properties():
    """
    Automatically add all of LBFluidRoutines's properties to LBSlice.

    """

    def set_attribute(lb_slice, value, attribute):
        """
        Setter function that sets attribute on every member of lb_slice.
        If values contains only one element, all members are set to it.

        """

        indices = [lb_slice.x_indices, lb_slice.y_indices, lb_slice.z_indices]
        N = [len(x) for x in indices]

        if N[0] * N[1] * N[2] == 0:
            raise AttributeError("Cannot set properties of an empty LBSlice")

        value = np.copy(value)
        attribute_shape = lb_slice.get_values(
            *np.zeros((3, 1), dtype=int), attribute).shape[3:]
        target_shape = (*N, *attribute_shape)

        # broadcast if only one element was provided
        if value.shape == attribute_shape:
            value = np.ones(target_shape) * value

        if value.shape != target_shape:
            raise ValueError(
                f"Input-dimensions of {attribute} array {value.shape} does not match slice dimensions {target_shape}.")

        lb_slice.set_values(*indices, attribute, value)

    def get_attribute(lb_slice, attribute):
        """
        Getter function that copies attribute from every member of
        lb_slice into an array (if possible).

        """

        indices = [lb_slice.x_indices, lb_slice.y_indices, lb_slice.z_indices]
        N = [len(x) for x in indices]

        if N[0] * N[1] * N[2] == 0:
            return np.empty(0, dtype=type(None))

        return lb_slice.get_values(*indices, attribute)

    for attribute_name in dir(LBFluidRoutines):
        if attribute_name in dir(LBSlice) or not isinstance(
                getattr(LBFluidRoutines, attribute_name), type(LBFluidRoutines.density)):
            continue

        # synthesize a new property
        new_property = property(
            functools.partial(get_attribute, attribute=attribute_name),
            functools.partial(set_attribute, attribute=attribute_name),
            doc=getattr(LBFluidRoutines, attribute_name).__doc__ or f'{attribute_name} for a slice')
        # attach the property to LBSlice
        setattr(LBSlice, attribute_name, new_property)


_add_lb_slice_properties()
