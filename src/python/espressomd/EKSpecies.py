from __future__ import print_function, absolute_import

import os

from .script_interface import ScriptObjectRegistry, ScriptInterfaceHelper, script_interface_register
import numpy as np
import functools
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

    def add_boundary_from_list(self, nodes,
                               value, boundary_type):
        """
        Set boundary conditions from a list of node indices.

        Parameters
        ----------
        nodes : (N, 3) array_like of :obj:`int`
            List of node indices to update. If they were originally not
            boundary nodes, they will become boundary nodes.
        velocity : (
        boundary_type :
            Type of the boundary condition.
        """
        if not (issubclass(boundary_type, FluxBoundary)
                or issubclass(boundary_type, DensityBoundary)):
            raise ValueError(
                "boundary_type must be a subclass of FluxBoundary or DensityBoundary")

        nodes = np.array(nodes, dtype=int)
        value = np.array(value, dtype=float)
        if issubclass(boundary_type, FluxBoundary):
            array_length = 3
            if np.shape(value) not in [(3,), tuple(self.shape) + (3,)]:
                raise ValueError(
                    f'Cannot process flux value grid of shape {np.shape(value)}')
        if issubclass(boundary_type, DensityBoundary):
            array_length = 1
            if np.shape(value) not in [(1,), tuple(self.shape) + (1,)]:
                raise ValueError(
                    f'Cannot process density grid of shape {np.shape(value)}')

        if value.shape == (array_length,):
            value = np.tile(value, (nodes.shape[0], 1))
        elif nodes.shape != value.shape:
            raise ValueError(
                f'Node indices and velocities must have the same shape, got {nodes.shape} and {value.shape}')

        # range checks
        if np.any(np.logical_or(np.max(nodes, axis=0) >= self.shape,
                                np.min(nodes, axis=0) < 0)):
            raise ValueError(
                f'Node indices must be in the range (0, 0, 0) to {self.shape}')

        # TODO unit conversion
        #        value *=
        value_view = np.ascontiguousarray(
            value.reshape((-1,)), dtype=np.double)
        nodes_view = np.ascontiguousarray(
            nodes.reshape((-1,)), dtype=np.int32)
        if issubclass(boundary_type, FluxBoundary):
            self.call_method(
                "update_flux_boundary_from_list",
                nodes_view=nodes_view,
                value_view=value_view)
        if issubclass(boundary_type, DensityBoundary):
            self.call_method(
                "update_density_boundary_from_list",
                nodes_view=nodes_view,
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


class EKSlice:
    def __init__(self, species, key, shape):
        self.species = species
        self.x_indices, self.y_indices, self.z_indices = self.get_indices(
            key, shape[0], shape[1], shape[2])

    def get_indices(self, key, shape_x, shape_y, shape_z):
        x_indices = np.atleast_1d(np.arange(shape_x)[key[0]])
        y_indices = np.atleast_1d(np.arange(shape_y)[key[1]])
        z_indices = np.atleast_1d(np.arange(shape_z)[key[2]])
        return x_indices, y_indices, z_indices

    def get_values(self, x_indices, y_indices, z_indices, prop_name):
        shape_res = np.shape(
            getattr(EKRoutines(self.species, np.array([0, 0, 0])), prop_name))
        res = np.zeros(
            (x_indices.size,
             y_indices.size,
             z_indices.size,
             *shape_res))
        for i, x in enumerate(x_indices):
            for j, y in enumerate(y_indices):
                for k, z in enumerate(z_indices):
                    res[i, j, k] = getattr(EKRoutines(self.species,
                                                      np.array([x, y, z])), prop_name)
        if shape_res == (1,):
            return np.squeeze(res, axis=-1)
        return res

    def set_values(self, x_indices, y_indices, z_indices, prop_name, value):
        for i, x in enumerate(x_indices):
            for j, y in enumerate(y_indices):
                for k, z in enumerate(z_indices):
                    setattr(EKRoutines(self.species,
                                       np.array([x, y, z])), prop_name, value[i, j, k])


def _add_ek_slice_properties():
    """
    Automatically add all of EKRoutines's properties to EKSlice.

    """

    def set_attribute(ek_slice, value, attribute):
        """
        Setter function that sets attribute on every member of lb_slice.
        If values contains only one element, all members are set to it.

        """

        indices = [ek_slice.x_indices, ek_slice.y_indices, ek_slice.z_indices]
        N = [len(x) for x in indices]

        if N[0] * N[1] * N[2] == 0:
            raise AttributeError("Cannot set properties of an empty LBSlice")

        value = np.copy(value)
        attribute_shape = ek_slice.get_values(
            *np.zeros((3, 1), dtype=int), attribute).shape[3:]
        target_shape = (*N, *attribute_shape)

        # broadcast if only one element was provided
        if value.shape == attribute_shape:
            value = np.ones(target_shape) * value

        if value.shape != target_shape:
            raise ValueError(
                f"Input-dimensions of {attribute} array {value.shape} does not match slice dimensions {target_shape}.")

        ek_slice.set_values(*indices, attribute, value)

    def get_attribute(ek_slice, attribute):
        """
        Getter function that copies attribute from every member of
        ek_slice into an array (if possible).

        """

        indices = [ek_slice.x_indices, ek_slice.y_indices, ek_slice.z_indices]
        N = [len(x) for x in indices]

        if N[0] * N[1] * N[2] == 0:
            return np.empty(0, dtype=type(None))

        return ek_slice.get_values(*indices, attribute)

    for attribute_name in dir(EKRoutines):
        if attribute_name in dir(EKSlice) or not isinstance(
                getattr(EKRoutines, attribute_name), type(EKRoutines.density)):
            continue

        # synthesize a new property
        new_property = property(
            functools.partial(get_attribute, attribute=attribute_name),
            functools.partial(set_attribute, attribute=attribute_name),
            doc=getattr(EKRoutines, attribute_name).__doc__ or f'{attribute_name} for a slice')
        # attach the property to LBSlice
        setattr(EKSlice, attribute_name, new_property)


_add_ek_slice_properties()
