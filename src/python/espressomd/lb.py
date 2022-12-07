#
# Copyright (C) 2013-2022 The ESPResSo project
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

import os
import itertools
import numpy as np
from .script_interface import ScriptInterfaceHelper, script_interface_register, array_variant
from .shapes import Shape
from . import utils
from .code_features import has_features


class VelocityBounceBack:
    """
    Holds velocity information for the velocity bounce back boundary
    condition at a single node.
    """

    def __init__(self, velocity):
        utils.check_type_or_throw_except(
            velocity, 3, float, "VelocityBounceBack velocity must be three floats")
        self.velocity = velocity


class HydrodynamicInteraction(ScriptInterfaceHelper):
    """
    Base class for LB implementations.

    Parameters
    ----------
    lattice :
        Lattice object. If not provided, a default one can be constructed
        using the ``agrid`` parameter.
    agrid : :obj:`float`
        Lattice constant. The box size in every direction must be an integer
        multiple of ``agrid``. Cannot be provided together with ``lattice``
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

    def __getitem__(self, key):
        raise NotImplementedError("Derived classes must implement this method")

    def __str__(self):
        return f"{self.__class__.__name__}({self.get_params()})"

    def _activate(self):
        self._activate_method()

    def _deactivate(self):
        self._deactivate_method()

    def _activate_method(self):
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
        for key in ('agrid', 'tau', 'density', 'viscosity'):
            if key not in params or key not in self.valid_keys():
                continue
            utils.check_type_or_throw_except(
                params[key], 1, float, f'{key} must be a number')
            if params[key] <= 0.:
                raise ValueError(f'{key} must be a strictly positive number')
        utils.check_type_or_throw_except(
            params['kT'], 1, float, 'kT must be a number')
        if params['kT'] < 0.:
            raise ValueError('kT must be a positive number')
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
        return {"lattice", "density", "viscosity", "tau"}

    def default_params(self):
        return {"lattice": None, "seed": 0, "kT": 0.,
                "ext_force_density": [0.0, 0.0, 0.0]}

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
        vel_max = cls.mach_limit(cls)
        velocities = np.reshape(velocities, (-1, 3))
        if np.any(np.linalg.norm(velocities, axis=1) > vel_max):
            speed_of_sound = 1. / np.sqrt(3.)
            mach_number = vel_max / speed_of_sound
            raise ValueError(f'Slip velocity exceeds Mach {mach_number:.2f}')

    def _set_params_in_es_core(self):
        pass

    def _get_params_from_es_core(self):
        self._params['agrid'] = self.agrid
        self._params['tau'] = self.tau
        self._params['dens'] = self.density
        self._params['kT'] = self.kT
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
        return self.call_method('get_interpolated_velocity', pos=pos)

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
        """
        Write LB node populations and boundary conditions to a file.
        """
        tmp_path = path + ".__tmp__"
        self.call_method(
            'save_checkpoint', path=tmp_path, mode=int(binary))
        os.rename(tmp_path, path)

    def load_checkpoint(self, path, binary):
        """
        Load LB node populations and boundary conditions from a file.
        """
        return self.call_method(
            'load_checkpoint', path=path, mode=int(binary))

    @property
    def pressure_tensor(self):
        tensor = self.call_method('get_pressure_tensor')
        return utils.array_locked(tensor)

    @pressure_tensor.setter
    def pressure_tensor(self, value):
        raise RuntimeError(f"Property 'pressure_tensor' is read-only")


class VTKRegistry:

    def __init__(self):
        if not has_features("WALBERLA"):
            raise NotImplementedError("Feature WALBERLA not compiled in")
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


if has_features("WALBERLA"):
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
    observables : :obj:`list`, {'density', 'velocity_vector', 'pressure_tensor'}
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
        if not has_features("WALBERLA"):
            raise NotImplementedError("Feature WALBERLA not compiled in")

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
        utils.check_required_keys(self.required_keys(), params.keys())
        utils.check_valid_keys(self.valid_keys(), params.keys())
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
        return set(self.call_method("get_valid_observable_names"))

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
        if not has_features("WALBERLA"):
            raise NotImplementedError("Feature WALBERLA not compiled in")

        if 'sip' not in kwargs:
            params = self.default_params()
            params.update(kwargs)
            self.validate_params(params)
            super().__init__(*args, **params)
            utils.handle_errors("LatticeWalberla initialization failed")
            self._params = {k: getattr(self, k) for k in self.valid_keys()}
        else:
            super().__init__(**kwargs)

    def validate_params(self, params):
        utils.check_required_keys(self.required_keys(), params.keys())
        utils.check_valid_keys(self.valid_keys(), params.keys())
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
    See :class:`HydrodynamicInteraction` for the list of parameters. If argument
    ``lattice`` is not provided, one will be default constructed if an argument
    ``agrid`` is provided.

    """
    _so_name = "walberla::FluidWalberla"
    _so_creation_policy = "GLOBAL"
    # TODO WALBERLA: here we cannot use _so_bind_methods without lb_vtk.py
    # failing: the walberla::FluidWalberla script interface object doesn't
    # expire even when the LBFluidWalberla is removed from the actors list
    # and the last python variable holding a reference to it is deleted.
    # _so_bind_methods = (
    #    "add_force_at_pos",
    #    "clear_boundaries",
    #    "get_interpolated_velocity",
    #    "get_pressure_tensor")

    def __init__(self, *args, **kwargs):
        if not has_features("WALBERLA"):
            raise NotImplementedError("Feature WALBERLA not compiled in")

        if 'sip' not in kwargs:
            params = self.default_params()
            params.update(kwargs)
            self.validate_params(params)
            super().__init__(*args, **params)
            utils.handle_errors(
                "HydrodynamicInteraction initialization failed")
            self._params = {k: params.get(k) for k in self.valid_keys()}
            self._params['agrid'] = self.agrid
            self._lattice_default_constructed = kwargs.get('lattice') is None
        else:
            super().__init__(**kwargs)

    @classmethod
    def _construct_from_cpt(cls, params):
        del params['agrid']
        lattice = params.pop('lattice')
        params['lattice'] = LatticeWalberla(
            agrid=lattice.agrid,
            n_ghost_layers=lattice.n_ghost_layers)
        obj = cls(**params)
        obj._params = params
        return obj

    def __reduce__(self):
        return self._construct_from_cpt, (self._params,), None

    def validate_params(self, params):
        super().validate_params(params)

        # construct default lattice if necessary
        if params.get('lattice') is None:
            if 'agrid' not in params:
                raise ValueError('missing argument "lattice" or "agrid"')
            params['lattice'] = LatticeWalberla(
                agrid=params.pop('agrid'), n_ghost_layers=1)
        elif 'agrid' in params:
            raise ValueError("cannot provide both 'lattice' and 'agrid'")

        utils.check_required_keys(self.required_keys(), params.keys())
        utils.check_valid_keys(self.valid_keys(), params.keys())

    def _set_params_in_es_core(self):
        pass

    def default_params(self):
        return {"single_precision": False, **super().default_params()}

    def valid_keys(self):
        return {"single_precision", *super().valid_keys()}

    def __getitem__(self, key):
        if isinstance(key, (tuple, list, np.ndarray)) and len(key) == 3:
            if any(isinstance(item, slice) for item in key):
                return LBFluidSliceWalberla(
                    parent_sip=self, slice_range=key, node_grid=self.shape)
            else:
                return LBFluidNodeWalberla(
                    parent_sip=self, index=np.array(key))

        raise TypeError(
            f"{key} is not a valid index. Should be a point on the "
            "nodegrid e.g. lbf[0,0,0], or a slice e.g. lbf[:,0,0]")

    def clear_boundaries(self):
        """
        Remove velocity bounce-back boundary conditions.
        """
        self.call_method("clear_boundaries")

    def add_boundary_from_shape(self, shape,
                                velocity=np.zeros(3, dtype=float),
                                boundary_type=VelocityBounceBack):
        """
        Set velocity bounce-back boundary conditions from a shape.

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
            raise TypeError(
                "Parameter 'boundary_type' must be a subclass of VelocityBounceBack")

        utils.check_type_or_throw_except(
            shape, 1, Shape, "expected an espressomd.shapes.Shape")
        if np.shape(velocity) not in [(3,), tuple(self.shape) + (3,)]:
            raise ValueError(
                f'Cannot process velocity value grid of shape {np.shape(velocity)}')

        # range checks
        lattice_speed = self.call_method('get_lattice_speed')
        velocity = np.array(velocity, dtype=float).reshape((-1, 3))
        velocity *= 1. / lattice_speed
        self._check_mach_limit(velocity)

        velocity_flat = velocity.reshape((-1,))
        mask_flat = shape.call_method('rasterize', grid_size=self.shape,
                                      grid_spacing=self._params['agrid'],
                                      grid_offset=0.5)

        self.call_method(
            'add_boundary_from_shape',
            raster=array_variant(mask_flat),
            velocity=array_variant(velocity_flat))

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
        lattice_speed = self.call_method('get_lattice_speed')
        velocity *= 1. / lattice_speed
        self._check_mach_limit(velocity)
        if np.any(np.logical_or(np.max(nodes, axis=0) >= self.shape,
                                np.min(nodes, axis=0) < 0)):
            raise ValueError(
                f'Node indices must be in the range (0, 0, 0) to {self.shape}')

        self.call_method('add_boundary_from_list', nodes=array_variant(
            nodes.reshape((-1,))), velocity=array_variant(velocity.reshape((-1,))))

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
        print("add_vtk_writer() is a deprecated function")
        return VTKOutput(lb_fluid=self, identifier=identifier,
                         observables=observables, **kwargs)


class LBFluidWalberlaGPU(HydrodynamicInteraction):
    """
    Initialize the lattice-Boltzmann method for hydrodynamic flow using
    waLBerla for the GPU. See :class:`HydrodynamicInteraction` for the
    list of parameters.

    """

    # pylint: disable=unused-argument
    def __init__(self, *args, **kwargs):
        if not has_features("CUDA"):
            raise NotImplementedError("Feature CUDA not compiled in")
        if not has_features("WALBERLA"):
            raise NotImplementedError("Feature WALBERLA not compiled in")
        raise NotImplementedError("Not implemented yet")


@script_interface_register
class LBFluidNodeWalberla(ScriptInterfaceHelper):
    _so_name = "walberla::FluidNodeWalberla"
    _so_creation_policy = "GLOBAL"

    def required_keys(self):
        return {"parent_sip", "index"}

    def validate_params(self, params):
        utils.check_required_keys(self.required_keys(), params.keys())
        utils.check_type_or_throw_except(
            params["index"], 3, int, "The index of a LB fluid node consists of three integers.")

    def __init__(self, *args, **kwargs):
        if "sip" not in kwargs:
            self.validate_params(kwargs)
            super().__init__(*args, **kwargs)
            utils.handle_errors("LBFluidNode instantiation failed")
        else:
            super().__init__(**kwargs)

    def __reduce__(self):
        raise NotImplementedError("Cannot serialize LB fluid node objects")

    def __eq__(self, obj):
        return isinstance(obj, LBFluidNodeWalberla) and self.index == obj.index

    def __hash__(self):
        return hash(self.index)

    @property
    def index(self):
        return tuple(self._index)

    @index.setter
    def index(self, value):
        raise RuntimeError("Parameter 'index' is read-only.")

    @property
    def density(self):
        return self.call_method("get_density")

    @density.setter
    def density(self, value):
        self.call_method("set_density", value=value)

    @property
    def population(self):
        return utils.array_locked(self.call_method("get_population"))

    @population.setter
    def population(self, value):
        self.call_method("set_population", value=value)

    @property
    def pressure_tensor(self):
        tensor = self.call_method("get_pressure_tensor")
        return utils.array_locked(tensor)

    @pressure_tensor.setter
    def pressure_tensor(self, value):
        raise RuntimeError("Property 'pressure_tensor' is read-only.")

    @property
    def is_boundary(self):
        return self.call_method("get_is_boundary")

    @is_boundary.setter
    def is_boundary(self, value):
        raise RuntimeError("Property 'is_boundary' is read-only.")

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

        velocity = self.call_method("get_velocity_at_boundary")
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
            lattice_speed = self.call_method("get_lattice_speed")
            HydrodynamicInteraction._check_mach_limit(
                np.array(value.velocity) / lattice_speed)
            self.call_method("set_velocity_at_boundary", value=value.velocity)
        elif value is None:
            self.call_method("set_velocity_at_boundary", value=None)
        else:
            raise TypeError(
                "Parameter 'value' must be an instance of VelocityBounceBack or None")

    @property
    def boundary_force(self):
        return self.call_method("get_boundary_force")

    @boundary_force.setter
    def boundary_force(self, value):
        raise RuntimeError("Property 'boundary_force' is read-only.")

    @property
    def velocity(self):
        return self.call_method("get_velocity")

    @velocity.setter
    def velocity(self, value):
        self.call_method("set_velocity", value=value)

    @property
    def last_applied_force(self):
        return self.call_method("get_last_applied_force")

    @last_applied_force.setter
    def last_applied_force(self, value):
        self.call_method("set_last_applied_force", value=value)


class LatticeSliceWalberla:

    def __init__(self, parent_sip, node_grid, slice_range):
        self.indices = [np.atleast_1d(np.arange(node_grid[i])[slice_range[i]])
                        for i in range(3)]
        self.dimensions = [ind.size for ind in self.indices]
        self._sip = parent_sip
        self._node = self._node_class(parent_sip=self._sip, index=[0, 0, 0])

    def __iter__(self):
        for index in itertools.product(*self.indices):
            yield self._node_class(parent_sip=self._sip, index=np.array(index))

    def setter_checks(self, attr, values):
        raise NotImplementedError("Derived classes must implement this method")

    def __getattr__(self, attr):
        node = self.__dict__.get("_node")
        if node is None or not hasattr(node, attr):
            if attr in self.__dict__:
                return self.__dict__[attr]
            raise AttributeError(
                f"Object '{self.__class__.__name__}' has no attribute '{attr}'")

        if 0 in self.dimensions:
            return np.empty(0, dtype=type(None))

        value = getattr(node, attr)
        shape_val = np.shape(value)
        dtype = value.dtype if isinstance(value, np.ndarray) else type(value)
        np_gen = np.zeros if dtype in (int, float, bool) else np.empty
        value_grid = np_gen((*self.dimensions, *shape_val), dtype=dtype)

        indices = itertools.product(*map(enumerate, self.indices))
        for (i, x), (j, y), (k, z) in indices:
            err = node.call_method("override_index", index=[x, y, z])
            assert err == 0
            value_grid[i, j, k] = getattr(node, attr)

        if shape_val == (1,):
            value_grid = np.squeeze(value_grid, axis=-1)
        return utils.array_locked(value_grid)

    def __setattr__(self, attr, values):
        node = self.__dict__.get("_node")
        if node is None or not hasattr(node, attr):
            self.__dict__[attr] = values
            return
        elif 0 in self.dimensions:
            raise AttributeError(
                f"Cannot set properties of an empty '{self.__class__.__name__}' object")

        values = np.copy(values)
        shape_val = np.shape(getattr(node, attr))
        target_shape = (*self.dimensions, *shape_val)

        # broadcast if only one element was provided
        if values.shape == shape_val:
            values = np.full(target_shape, values)

        if values.shape != target_shape:
            raise ValueError(
                f"Input-dimensions of '{attr}' array {values.shape} does not match slice dimensions {target_shape}")

        self.setter_checks(attr, values)

        indices = itertools.product(*map(enumerate, self.indices))
        for (i, x), (j, y), (k, z) in indices:
            err = node.call_method("override_index", index=[x, y, z])
            assert err == 0
            setattr(node, attr, values[i, j, k])


class LBFluidSliceWalberla(LatticeSliceWalberla):
    _node_class = LBFluidNodeWalberla

    def __reduce__(self):
        raise NotImplementedError("Cannot serialize LB fluid slice objects")

    def setter_checks(self, attr, values):
        node = self.__dict__.get("_node")
        if attr == "boundary":
            # Mach checks
            lattice_speed = node.call_method("get_lattice_speed")
            for value in values.flatten():
                if isinstance(value, VelocityBounceBack):
                    HydrodynamicInteraction._check_mach_limit(
                        np.array(value.velocity) / lattice_speed)
                elif value is not None:
                    raise TypeError(
                        "values must be instances of VelocityBounceBack or None")


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
