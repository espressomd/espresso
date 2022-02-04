from __future__ import print_function, absolute_import

import os

from .script_interface import ScriptObjectRegistry, ScriptInterfaceHelper, script_interface_register
import numpy as np
import itertools

from . import utils
from .__init__ import has_features

from .shapes import Shape
from .lb import VTKRegistry


@script_interface_register
class EKFFT(ScriptObjectRegistry):
    _so_name = "walberla::EKFFT"
    _so_creation_policy = "GLOBAL"


@script_interface_register
class EKNone(ScriptObjectRegistry):
    _so_name = "walberla::None"
    _so_creation_policy = "GLOBAL"


@script_interface_register
class EKContainer(ScriptObjectRegistry):
    _so_name = "walberla::EKContainer"

    @property
    def tau(self):
        return self.call_method("get_tau")

    @tau.setter
    def tau(self, tau):
        self.call_method("set_tau", tau=tau)

    @property
    def solver(self):
        raise NotImplementedError(
            "PoissonSolver object property is not implemented")

    @solver.setter
    def solver(self, solver):
        self.call_method("set_poissonsolver", object=solver)

    def add(self, ekspecies):
        if not isinstance(ekspecies, EKSpecies):
            raise TypeError("EKSpecies object is not of correct type.")

        self.call_method("add", object=ekspecies)

        return ekspecies

    def remove(self, ekspecies):
        self.call_method("remove", object=ekspecies)

    def clear(self):
        self.call_method("clear")


@script_interface_register
class EKSpecies(ScriptInterfaceHelper):
    """Interface to the Walberla EKSpecies
    """
    _so_name = "walberla::EKSpecies"
    _so_creation_policy = "GLOBAL"

    def __getitem__(self, key):
        if isinstance(key, (tuple, list, np.ndarray)):
            if len(key) == 3:
                if any(isinstance(typ, slice) for typ in key):
                    shape = self.shape
                    return EKSlice(self, key, (shape[0], shape[1], shape[2]))
                else:
                    return EKRoutines(self, np.array(key))
        else:
            raise Exception(
                f"{key} is not a valid key. Should be a point on the nodegrid e.g. ek[0,0,0], or a slice")

    def add_boundary_from_shape(self, shape,
                                value, boundary_type):
        """
        Set boundary conditions from a shape.

        Parameters
        ----------
        shape : :obj:`espressomd.shapes.Shape`
            Shape to rasterize.
        value :
        boundary_type :
            Type of the boundary condition.
        """
        if not (issubclass(boundary_type, FluxBoundary)
                or issubclass(boundary_type, DensityBoundary)):
            raise ValueError(
                "boundary_type must be a subclass of FluxBoundary or DensityBoundary")

        if not hasattr(value, "__iter__"):
            value = (value, )

        value = np.array(value, dtype=float)
        utils.check_type_or_throw_except(
            shape, 1, Shape, "expected an espressomd.shapes.Shape")
        if issubclass(boundary_type, FluxBoundary):
            if np.shape(value) not in [(3,), tuple(self.shape) + (3,)]:
                raise ValueError(
                    f'Cannot process flux value grid of shape {np.shape(value)}')
        if issubclass(boundary_type, DensityBoundary):
            if np.shape(value) not in [(1,), tuple(self.shape) + (1,)]:
                raise ValueError(
                    f'Cannot process density grid of shape {np.shape(value)}')

        # TODO unit conversion
        #        value *=
        value_flat = value.reshape((-1,))
        mask_flat = shape.call_method('rasterize', grid_size=self.shape,
                                      grid_spacing=self.lattice.agrid,
                                      grid_offset=0.5)

        value_view = np.ascontiguousarray(value_flat, dtype=np.double)
        raster_view = np.ascontiguousarray(mask_flat, dtype=np.int32)
        if issubclass(boundary_type, FluxBoundary):
            self.call_method(
                "update_flux_boundary_from_shape",
                raster_view=raster_view,
                value_view=value_view)
        if issubclass(boundary_type, DensityBoundary):
            self.call_method(
                "update_density_boundary_from_shape",
                raster_view=raster_view,
                value_view=value_view)

    def get_nodes_in_shape(self, shape):
        """Provides a generator for iterating over all ek nodes inside the given shape"""
        utils.check_type_or_throw_except(
            shape, 1, Shape, "expected a espressomd.shapes.Shape")
        ek_shape = self.shape
        idxs = itertools.product(
            range(ek_shape[0]), range(ek_shape[1]), range(ek_shape[2]))
        for idx in idxs:
            pos = (np.asarray(idx) + 0.5) * self._params['agrid']
            if shape.is_inside(position=pos):
                yield self[idx]

    def clear_boundaries(self):
        """
        Remove boundary conditions.
        """
        self.call_method("clear_boundary")

    def get_shape_bitmask(self, shape):
        """Create a bitmask for the given shape."""
        utils.check_type_or_throw_except(
            shape, 1, Shape, "expected a espressomd.shapes.Shape")
        mask_flat = shape.call_method('rasterize', grid_size=self.shape,
                                      grid_spacing=self.lattice.agrid,
                                      grid_offset=0.5)
        return np.reshape(mask_flat, self.shape).astype(type(True))


class FluxBoundary:
    """
    Holds flux information for the flux boundary
    condition at a single node.
    """

    def __init__(self, flux):
        utils.check_type_or_throw_except(
            flux, 3, float, "FluxBoundary flux must be three floats")
        self.flux = flux


class DensityBoundary:
    """
    Holds density information for the density boundary
    condition at a single node.
    """

    def __init__(self, density):
        utils.check_type_or_throw_except(
            density, 1, float, "DensityBoundary flux must be one float")
        self.density = density


if has_features("LB_WALBERLA"):
    _ek_vtk_registry = VTKRegistry()


@script_interface_register
class EKVTKOutput(ScriptInterfaceHelper):
    """
    Create a VTK writer.

    Files are written to ``<base_folder>/<identifier>/<prefix>_*.vtu``.
    Summary is written to ``<base_folder>/<identifier>.pvd``.

    Manual VTK callbacks can be called at any time to take a snapshot
    of the current state of the EK species.

    Automatic VTK callbacks can be disabled at any time and re-enabled later.
    Please note that the internal VTK counter is no longer incremented when
    an automatic callback is disabled, which means the number of EK steps
    between two frames will not always be an integer multiple of ``delta_N``.

    Parameters
    ----------
    identifier : :obj:`str`
        Name of the VTK writer.
    observables : :obj:`list`, \{'density',\}
        List of observables to write to the VTK files.
    delta_N : :obj:`int`
        Write frequency. If this value is 0 (default), the object is a
        manual VTK callback that must be triggered manually. Otherwise,
        it is an automatic callback that is added to the time loop and
        writes every ``delta_N`` EK steps.
    base_folder : :obj:`str` (optional), default is 'vtk_out'
        Path to the output VTK folder.
    prefix : :obj:`str` (optional), default is 'simulation_step'
        Prefix for VTK files.

    """
    _so_name = "walberla::EKVTKHandle"
    _so_creation_policy = "GLOBAL"
    _so_bind_methods = ("enable", "disable", "write")

    def __init__(self, *args, **kwargs):
        if not has_features("LB_WALBERLA"):
            raise NotImplementedError("Feature LB_WALBERLA not compiled in")
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
        _ek_vtk_registry._register_vtk_object(self)

    def validate_params(self, params):
        if not self.required_keys().issubset(params):
            raise ValueError(
                f"At least the following keys have to be given as keyword arguments: {sorted(self.required_keys())}")
        if not self.valid_observables().issuperset(params['observables']):
            raise ValueError(
                f"Only the following VTK observables are supported: {sorted(self.valid_observables())}, got {params['observables']}")
        if not isinstance(params['species'], EKSpecies):
            raise ValueError("'species' must be an EKSpecies")
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
        if vtk_path in _ek_vtk_registry.collisions:
            raise RuntimeError(
                f"VTK object '{vtk_uid}' would overwrite files written "
                f"by VTK object '{_ek_vtk_registry.collisions[vtk_path]}'")
        params['vtk_uid'] = vtk_uid

    def valid_observables(self):
        return {'density', }

    def valid_keys(self):
        return {'species', 'delta_N', 'execution_count', 'observables',
                'identifier', 'base_folder', 'prefix', 'enabled'}

    def required_keys(self):
        return self.valid_keys() - self.default_params().keys()

    def default_params(self):
        return {'delta_N': 0, 'enabled': True, 'execution_count': 0,
                'base_folder': 'vtk_out', 'prefix': 'simulation_step'}

    def __repr__(self):
        class_id = f"{self.__class__.__module__}.{self.__class__.__name__}"
        if self.delta_N:
            write_when = f"every {self.delta_N} EK steps"
            if not self.enabled:
                write_when += " (disabled)"
        else:
            write_when = "on demand"
        return f"<{class_id}: write to '{self.vtk_uid}' {write_when}>"


class EKRoutines:
    def __init__(self, species, node):
        self.node = node
        self.species = species

    @property
    def index(self):
        return self.node

    @property
    def density(self):
        return self.species.call_method("get_density", position=self.node)

    @density.setter
    def density(self, value):
        self.species.call_method(
            "set_density",
            position=self.node,
            value=value)

    @property
    def is_boundary(self):
        return self.species.call_method("is_boundary", position=self.node)

    @property
    def flux_boundary(self):
        flux = self.species.call_method(
            "get_node_flux_at_boundary", position=self.node)

        if flux is not None:
            return FluxBoundary(flux)
        return None

    @flux_boundary.setter
    def flux_boundary(self, flux):
        """
        Parameters
        ----------
        flux : :class:`~espressomd.EKSpecies.FluxBoundary` or None
            If flux is :class:`~espressomd.EkSpecies.FluxBoundary`,
            set the node to be a boundary node with the specified flux.
            If flux is ``None``, the node will become a domain node.
        """

        if isinstance(flux, FluxBoundary):
            self.species.call_method(
                'set_node_flux_boundary',
                position=self.node,
                flux=flux.flux)
        elif flux is None:
            self.species.call_method(
                'remove_node_flux_boundary',
                position=self.node)
        else:
            raise ValueError(
                "value must be an instance of FluxBoundary or None")


class EKSlice:
    def __init__(self, species, key, shape):
        self._species = species
        self.indices = [
            np.atleast_1d(
                np.arange(
                    shape[i])[
                    key[i]]) for i in range(3)]
        self.dimensions = [ind.size for ind in self.indices]

        self._node = EKRoutines(species=species, node=[0, 0, 0])

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
            node.node = [x, y, z]
            value_grid[i, j, k] = getattr(node, attr)

        if shape_val == (1,):
            value_grid = np.squeeze(value_grid, axis=-1)
        return utils.array_locked(value_grid)

    def __setattr__(self, attr, values):
        node = self.__dict__.get('_node')
        if node is None or not hasattr(node, attr):
            self.__dict__[attr] = values
            return
        elif 0 in self.dimensions:
            raise AttributeError("Cannot set properties of an empty LBSlice")

        values = np.copy(values)
        shape_val = np.shape(getattr(node, attr))
        target_shape = (*self.dimensions, *shape_val)

        # broadcast if only one element was provided
        if values.shape == shape_val:
            values = np.full(target_shape, values)

        if values.shape != target_shape:
            raise ValueError(
                f"Input-dimensions of '{attr}' array {values.shape} does not match slice dimensions {target_shape}.")

        indices = itertools.product(*map(enumerate, self.indices))
        for (i, x), (j, y), (k, z) in indices:
            node.node = [x, y, z]
            setattr(node, attr, values[i, j, k])

    def __iter__(self):
        return (EKRoutines(species=self._species, node=np.array(index))
                for index in itertools.product(*self.indices))


@script_interface_register
class EKReactant(ScriptObjectRegistry):
    _so_name = "walberla::EKReactant"
    _so_creation_policy = "GLOBAL"


@script_interface_register
class EKReaction(ScriptObjectRegistry):
    _so_name = "walberla::EKReaction"
    _so_creation_policy = "GLOBAL"
