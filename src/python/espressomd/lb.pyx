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
from . cimport cuda_init
from . import cuda_init
from . import utils
from .utils import is_valid_type, to_char_pointer, check_type_or_throw_except
from .utils cimport Vector3i
from .utils cimport Vector3d
from .utils cimport Vector6d
from .utils cimport make_array_locked
from .utils cimport make_Vector3d
from .utils cimport create_nparray_from_double_array


class VelocityBounceBack:
    """
    Holds velocity information for the velocity bounce back boundary
    condition at a single node.
    """

    def __init__(self, velocity):
        check_type_or_throw_except(
            velocity, 3, float, "VelocityBounceBack velocity must be three floats")
        self.velocity = velocity


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
                    return LBFluidNode(index=np.array(key))
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
        return {"agrid", "tau", "density", "ext_force_density", "viscosity",
                "lattice", "kT", "seed"}

    def required_keys(self):
        return {"density", "agrid", "viscosity", "tau"}

    def default_params(self):
        return {"agrid": -1.0,
                "density": -1.0,
                "ext_force_density": [0.0, 0.0, 0.0],
                "viscosity": -1.0,
                "tau": -1.0,
                "lattice": None,
                "seed": 0,
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
        return self.call_method('get_interpolated_velocity', pos)

    def add_force_at_pos(self, pos, force):
        """Adds a force to the fluid at given position

        Parameters
        ----------
        pos : (3,) array_like of :obj:`float`
              The position at which the force will be added.
        force : (3,) array_like of :obj:`float`
              The force vector which will be distributed at the position.

        """
        return self.call_method('add_force_at_pos', pos=pos, force=force)

    def save_checkpoint(self, path, binary):
        '''
        Write LB node populations to a file.
        Boundary information is not written to the file.
        '''
        tmp_path = path + ".__tmp__"
        lb_lbfluid_save_checkpoint(utils.to_char_pointer(tmp_path), binary)
        os.rename(tmp_path, path)

    def load_checkpoint(self, path, binary):
        '''
        Load LB node populations from a file.
        Boundary information is not available in the file. The boundary
        information of the grid will be set to zero.
        '''
        lb_lbfluid_load_checkpoint(utils.to_char_pointer(path), binary)

    @property
    def pressure_tensor(self):
        tensor = self.call_method('get_pressure_tensor')
        return utils.array_locked(tensor)

    @pressure_tensor.setter
    def pressure_tensor(self, value):
        raise RuntimeError(f"Property 'pressure_tensor' is read-only")

    def nodes(self):
        """Provides a generator for iterating over all lb nodes"""

        shape = self.shape
        for i, j, k in itertools.product(
                range(shape[0]), range(shape[1]), range(shape[2])):
            yield self[i, j, k]


IF LB_WALBERLA:

    class VTKRegistry:

        def __init__(self):
            self.map = {}
            self.collisions = {}

        def _register_vtk_object(self, vtk_obj):
            vtk_uid = vtk_obj.vtk_uid
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

    @script_interface_register
    class VTKOutput(ScriptInterfaceHelper):
        """
        Create a VTK writer.

        Files are written to ``<base_folder>/<identifier>/<prefix>_*.vtu``.
        Summary is written to ``<base_folder>/<identifier>.pvd``.

        Manual VTK callbacks can be called at any time to take a snapshot
        of the current state of the LB fluid.

        Automatic VTK callbacks can be disabled at any time and re-enabled later.
        Please note that the internal VTK counter is no longer incremented when
        an automatic callback is disabled, which means the number of LB steps
        between two frames will not always be an integer multiple of ``delta_N``.

        Parameters
        ----------
        identifier : :obj:`str`
            Name of the VTK writer.
        observables : :obj:`list`, \{'density', 'velocity_vector', 'pressure_tensor'\}
            List of observables to write to the VTK files.
        delta_N : :obj:`int`
            Write frequency. If this value is 0 (default), the object is a
            manual VTK callback that must be triggered manually. Otherwise,
            it is an automatic callback that is added to the time loop and
            writes every ``delta_N`` LB steps.
        base_folder : :obj:`str` (optional), default is 'vtk_out'
            Path to the output VTK folder.
        prefix : :obj:`str` (optional), default is 'simulation_step'
            Prefix for VTK files.

        """
        _so_name = "walberla::VTKHandle"
        _so_creation_policy = "GLOBAL"
        _so_bind_methods = ("enable", "disable", "write")

        def __init__(self, *args, **kwargs):
            if 'sip' not in kwargs:
                params = self.default_params()
                params.update(kwargs)
                if isinstance(params['observables'], str):
                    params['observables'] = [params['observables']]
                self.validate_params(params)
                super().__init__(*args, **params)
                utils.handle_errors(
                    f"{self.__class__.__name__} initialization failed")
            else:
                super().__init__(**kwargs)
            _vtk_registry._register_vtk_object(self)

        def validate_params(self, params):
            if not self.required_keys().issubset(params):
                raise ValueError(
                    f"At least the following keys have to be given as keyword arguments: {sorted(self.required_keys())}")
            if not self.valid_observables().issuperset(params['observables']):
                raise ValueError(
                    f"Only the following VTK observables are supported: {sorted(self.valid_observables())}, got {params['observables']}")
            if not isinstance(params['lb_fluid'], HydrodynamicInteraction):
                raise ValueError("'lb_fluid' must be an LB actor")
            utils.check_type_or_throw_except(
                params['delta_N'], 1, int, "'delta_N' must be 1 integer")
            if params['delta_N'] < 0:
                raise ValueError("'delta_N' must be a positive integer")
            utils.check_type_or_throw_except(
                params['base_folder'], 1, str, "'base_folder' must be a string")
            utils.check_type_or_throw_except(
                params['identifier'], 1, str, "'identifier' must be a string")
            if os.path.sep in params['identifier']:
                raise ValueError(
                    "'identifier' must be a string, not a filepath")
            vtk_uid = os.path.join(params['base_folder'],
                                   params['identifier'])
            vtk_path = os.path.abspath(vtk_uid)
            if vtk_path in _vtk_registry.collisions:
                raise RuntimeError(
                    f"VTK object '{vtk_uid}' would overwrite files written "
                    f"by VTK object '{_vtk_registry.collisions[vtk_path]}'")
            params['vtk_uid'] = vtk_uid

        def valid_observables(self):
            return {'density', 'velocity_vector', 'pressure_tensor'}

        def valid_keys(self):
            return {'lb_fluid', 'delta_N', 'execution_count', 'observables',
                    'identifier', 'base_folder', 'prefix', 'enabled'}

        def required_keys(self):
            return self.valid_keys() - self.default_params().keys()

        def default_params(self):
            return {'delta_N': 0, 'enabled': True, 'execution_count': 0,
                    'base_folder': 'vtk_out', 'prefix': 'simulation_step'}

        def __repr__(self):
            class_id = f"{self.__class__.__module__}.{self.__class__.__name__}"
            if self.delta_N:
                write_when = f"every {self.delta_N} LB steps"
                if not self.enabled:
                    write_when += " (disabled)"
            else:
                write_when = "on demand"
            return f"<{class_id}: write to '{self.vtk_uid}' {write_when}>"

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
        # TODO WALBERLA: here we cannot use _so_bind_methods without lb_vtk.py
        # failing: the walberla::FluidWalberla script interface object doesn't
        # expire even when the LBFluidWalberla is removed from the actors list
        # and the last python variable holding a reference to it is deleted.
        #_so_bind_methods = (
        #    "add_force_at_pos",
        #    "get_interpolated_velocity",
        #    "get_pressure_tensor")

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

        def default_params(self):
            return {"single_precision": False, **super().default_params()}

        def valid_keys(self):
            return {"single_precision", *super().valid_keys()}

        def clear_boundaries(self):
            """
            Remove boundary conditions.
            """
            lb_lbfluid_clear_boundaries()

        def add_boundary_from_shape(self, shape,
                                    velocity=np.zeros(3, dtype=float),
                                    boundary_type=VelocityBounceBack):
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
            boundary_type : Union[:class:`~espressomd.lb.VelocityBounceBack`] (optional)
                Type of the boundary condition.
            """
            if not issubclass(boundary_type, VelocityBounceBack):
                raise ValueError(
                    "boundary_type must be a subclass of VelocityBounceBack")

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

        def add_boundary_from_list(self, nodes,
                                   velocity=np.zeros(3, dtype=float),
                                   boundary_type=VelocityBounceBack):
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
            boundary_type : Union[:class:`~espressomd.lb.VelocityBounceBack`] (optional)
                Type of the boundary condition.
            """
            if not issubclass(boundary_type, VelocityBounceBack):
                raise ValueError(
                    "boundary_type must be a subclass of VelocityBounceBack")

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

        # TODO WALBERLA: deprecated function, consider removing it
        def add_vtk_writer(self, identifier, observables, **kwargs):
            """
            Create a VTK observable.

            """
            return VTKOutput(lb_fluid=self, identifier=identifier,
                             observables=observables, **kwargs)

    IF CUDA:
        class LBFluidWalberlaGPU(HydrodynamicInteraction):
            """
            Initialize the lattice-Boltzmann method for hydrodynamic flow using
            waLBerla for the GPU. See :class:`HydrodynamicInteraction` for the
            list of parameters.

            """

            def __init__(self, *args, **kwargs):
                raise NotImplementedError("Not implemented yet")
    ELSE:
        class LBFluidWalberlaGPU(HydrodynamicInteraction):
            def __init__(self, *args, **kwargs):
                raise NotImplementedError("Feature not compiled in")

ELSE:
    class LatticeWalberla(ScriptInterfaceHelper):
        def __init__(self, *args, **kwargs):
            raise NotImplementedError("Feature not compiled in")

    class LBFluidWalberla(HydrodynamicInteraction):
        def __init__(self, *args, **kwargs):
            raise NotImplementedError("Feature not compiled in")

    class LBFluidWalberlaGPU(HydrodynamicInteraction):
        def __init__(self, *args, **kwargs):
            raise NotImplementedError("Feature not compiled in")


@script_interface_register
class LBFluidNode(ScriptInterfaceHelper):
    _so_name = "walberla::FluidNodeWalberla"
    _so_creation_policy = "GLOBAL"

    def required_keys(self):
        return {'index'}

    def validate_params(self, params):
        utils.check_type_or_throw_except(
            params['index'], 3, int, "The index of an lb fluid node consists of three integers.")

    def __init__(self, *args, **kwargs):
        if 'sip' not in kwargs:
            if not self.required_keys().issubset(kwargs):
                raise ValueError(
                    f"At least the following keys have to be given as keyword arguments: {self.required_keys()}")
            self.validate_params(kwargs)
            super().__init__(*args, **kwargs)
            utils.handle_errors("LBFluidNode instantiation failed")
        else:
            super().__init__(**kwargs)

    @property
    def index(self):
        return tuple(self._index)

    @index.setter
    def index(self, value):
        raise RuntimeError("Property 'index' is read-only.")

    @property
    def population(self):
        return utils.array_locked(self._population)

    @population.setter
    def population(self, value):
        self._population = value

    @property
    def pressure_tensor(self):
        p_tensor = self._pressure_tensor
        p_tensor = [[p_tensor[0], p_tensor[1], p_tensor[3]],
                    [p_tensor[1], p_tensor[2], p_tensor[4]],
                    [p_tensor[3], p_tensor[4], p_tensor[5]]]
        return utils.array_locked(p_tensor)

    @pressure_tensor.setter
    def pressure_tensor(self, value):
        raise RuntimeError("Property 'pressure_tensor' is read-only.")

    @property
    def boundary(self):
        """
        Returns
        -------
        :class:`~espressomd.lb.VelocityBounceBack`
            If the node is a boundary node
        None
            If the node is not a boundary node
        """

        velocity = self.velocity_at_boundary
        if velocity is not None:
            return VelocityBounceBack(velocity)
        return None

    @boundary.setter
    def boundary(self, value):
        """
        Parameters
        ----------
        value : :class:`~espressomd.lb.VelocityBounceBack` or None
            If value is :class:`~espressomd.lb.VelocityBounceBack`,
            set the node to be a boundary node with the specified velocity.
            If value is ``None``, the node will become a fluid node.
        """

        if isinstance(value, VelocityBounceBack):
            HydrodynamicInteraction._check_mach_limit(value.velocity)
            self.velocity_at_boundary = value.velocity
        elif value is None:
            self.velocity_at_boundary = None
        else:
            raise ValueError(
                "value must be an instance of VelocityBounceBack or None")

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

    def __getattr__(self, attr):
        node = LBFluidNode(index=np.array([0, 0, 0]))
        if not hasattr(node, attr):
            if attr in self.__dict__:
                return self.__dict__[attr]
            raise AttributeError(
                f"Object '{self.__class__.__name__}' has no attribute '{attr}'")

        dimensions = [self.x_indices.size,
                      self.y_indices.size,
                      self.z_indices.size]
        if 0 in dimensions:
            return np.empty(0, dtype=type(None))

        res = getattr(node, attr)
        shape_res = np.shape(res)
        dtype = res.dtype if isinstance(res, np.ndarray) else type(res)
        res = np.zeros((*dimensions, *shape_res), dtype=dtype)
        for i, x in enumerate(self.x_indices):
            for j, y in enumerate(self.y_indices):
                for k, z in enumerate(self.z_indices):
                    err = node.call_method("override_index", index=[x, y, z])
                    assert err == 0
                    res[i, j, k] = getattr(node, attr)
        if shape_res == (1,):
            res = np.squeeze(res, axis=-1)
        return utils.array_locked(res)

    def __setattr__(self, attr, value):
        node = LBFluidNode(index=np.array([0, 0, 0]))
        if not hasattr(node, attr):
            self.__dict__[attr] = value
            return

        dimensions = [self.x_indices.size,
                      self.y_indices.size,
                      self.z_indices.size]
        if 0 in dimensions:
            raise AttributeError("Cannot set properties of an empty LBSlice")

        value = np.copy(value)
        shape_res = np.shape(getattr(node, attr))
        target_shape = (*dimensions, *shape_res)

        # broadcast if only one element was provided
        if value.shape == shape_res:
            value = np.ones(target_shape) * value

        if value.shape != target_shape:
            raise ValueError(
                f"Input-dimensions of '{attr}' array {value.shape} does not match slice dimensions {target_shape}.")

        for i, x in enumerate(self.x_indices):
            for j, y in enumerate(self.y_indices):
                for k, z in enumerate(self.z_indices):
                    err = node.call_method("override_index", index=[x, y, z])
                    assert err == 0
                    setattr(node, attr, value[i, j, k])

    def __iter__(self):
        indices = [(x, y, z) for (x, y, z) in itertools.product(
            self.x_indices, self.y_indices, self.z_indices)]
        return (LBFluidNode(index=np.array(index)) for index in indices)


def edge_detection(boundary_mask, periodicity):
    """
    Find boundary nodes in contact with the fluid. Relies on a convolution
    kernel constructed from the D3Q19 stencil.

    Parameters
    ----------
    boundary_mask : (N, M, L) array_like of :obj:`bool`
        Bitmask for the rasterized boundary geometry.
    periodicity : (3,) array_like of :obj:`bool`
        Bitmask for the box periodicity.

    Returns
    -------
    (N, 3) array_like of :obj:`int`
        The indices of the boundary nodes at the interface with the fluid.
    """
    import scipy.signal
    import itertools

    fluid_mask = np.logical_not(boundary_mask)

    # edge kernel
    edge = -np.ones((3, 3, 3))
    for i, j, k in itertools.product((0, 2), (0, 2), (0, 2)):
        edge[i, j, k] = 0
    edge[1, 1, 1] = -np.sum(edge)

    # periodic convolution
    wrapped_mask = np.pad(fluid_mask.astype(int), 3 * [(2, 2)], mode='wrap')
    if not periodicity[0]:
        wrapped_mask[:2, :, :] = 0
        wrapped_mask[-2:, :, :] = 0
    if not periodicity[1]:
        wrapped_mask[:, :2, :] = 0
        wrapped_mask[:, -2:, :] = 0
    if not periodicity[2]:
        wrapped_mask[:, :, :2] = 0
        wrapped_mask[:, :, -2:] = 0
    convolution = scipy.signal.convolve(
        wrapped_mask, edge, mode='same', method='direct')[2:-2, 2:-2, 2:-2]
    convolution = np.multiply(convolution, boundary_mask)

    return np.array(np.nonzero(convolution < 0)).T


def calc_cylinder_tangential_vectors(center, agrid, offset, node_indices):
    """
    Utility function to calculate a constant slip velocity tangential to the
    surface of a cylinder.

    Parameters
    ----------
    center : (3,) array_like of :obj:`float`
        Center of the cylinder.
    agrid : :obj:`float`
        LB agrid.
    offset : :obj:`float`
        LB offset.
    node_indices : (N, 3) array_like of :obj:`int`
        Indices of the boundary surface nodes.

    Returns
    -------
    (N, 3) array_like of :obj:`float`
        The unit vectors tangential to the surface of a cylinder.
    """
    velocities = []
    for ijk in node_indices:
        p = (ijk + offset) * agrid
        r = center - p
        norm = np.linalg.norm(r[:2])
        if norm < 1e-10:
            velocities.append(np.zeros(3))
            continue
        angle_r = np.arccos(np.dot(r[:2] / norm, [1, 0]))
        angle_v = angle_r - np.pi / 2
        flip = np.sign(r[1])
        slip_velocity = np.array([flip * np.cos(angle_v), np.sin(angle_v), 0.])
        velocities.append(slip_velocity)
    return np.array(velocities)
