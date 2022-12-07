#
# Copyright (C) 2021-2022 The ESPResSo project
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

from .script_interface import ScriptObjectList, ScriptInterfaceHelper, script_interface_register
import numpy as np
import itertools

from . import utils
from .code_features import has_features

from .shapes import Shape
from .lb import VTKRegistry, LatticeSliceWalberla


class EKFFT(ScriptInterfaceHelper):
    _so_name = "walberla::EKFFT"
    _so_creation_policy = "GLOBAL"

    def __init__(self, *args, **kwargs):
        if not has_features("WALBERLA"):
            raise NotImplementedError("Feature WALBERLA not compiled in")

        super().__init__(*args, **kwargs)


class EKNone(ScriptInterfaceHelper):
    _so_name = "walberla::None"
    _so_creation_policy = "GLOBAL"

    def __init__(self, *args, **kwargs):
        if not has_features("WALBERLA"):
            raise NotImplementedError("Feature WALBERLA not compiled in")

        super().__init__(*args, **kwargs)


@script_interface_register
class EKContainer(ScriptObjectList):
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

    def __init__(self, *args, **kwargs):
        if not has_features("WALBERLA"):
            raise NotImplementedError("Feature WALBERLA not compiled in")

        super().__init__(*args, **kwargs)

    def __getitem__(self, key):
        if isinstance(key, (tuple, list, np.ndarray)) and len(key) == 3:
            if any(isinstance(item, slice) for item in key):
                return EKSpeciesSlice(
                    parent_sip=self, slice_range=key, node_grid=self.shape)
            else:
                return EKSpeciesNode(parent_sip=self, index=np.array(key))

        raise TypeError(
            f"{key} is not a valid index. Should be a point on the "
            "nodegrid e.g. ek[0,0,0], or a slice, e.g. ek[:,0,0]")

    def clear_density_boundaries(self):
        """
        Remove density boundary conditions.
        """
        self.call_method("clear_density_boundaries")

    def clear_flux_boundaries(self):
        """
        Remove flux boundary conditions.
        """
        self.call_method("clear_flux_boundaries")

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
        if not issubclass(boundary_type, (FluxBoundary, DensityBoundary)):
            raise TypeError(
                "Parameter 'boundary_type' must be a subclass of FluxBoundary or DensityBoundary")

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
                    f'Cannot process density value grid of shape {np.shape(value)}')

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


if has_features("WALBERLA"):
    _ek_vtk_registry = VTKRegistry()


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
    observables : :obj:`list`, {'density',}
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
        _ek_vtk_registry._register_vtk_object(self)

    def validate_params(self, params):
        if not self.required_keys().issubset(params):
            raise ValueError(
                f"At least the following keys have to be given as keyword arguments: {sorted(self.required_keys())}")
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
        return set(self.call_method("get_valid_observable_names"))

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


@script_interface_register
class EKSpeciesNode(ScriptInterfaceHelper):
    _so_name = "walberla::EKSpeciesNode"
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
            utils.handle_errors("EKSpeciesNode instantiation failed")
        else:
            super().__init__(**kwargs)

    def __reduce__(self):
        raise NotImplementedError("Cannot serialize EK species node objects")

    def __eq__(self, obj):
        return isinstance(obj, EKSpeciesNode) and self.index == obj.index

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
    def is_boundary(self):
        return self.call_method("get_is_boundary")

    @is_boundary.setter
    def is_boundary(self, value):
        raise RuntimeError("Property 'is_boundary' is read-only.")

    @property
    def density_boundary(self):
        """
        Returns
        -------
        :class:`~espressomd.EKSpecies.DensityBoundary`
            If the node is a boundary node
        None
            If the node is not a boundary node
        """
        density = self.call_method("get_node_density_at_boundary")
        if density is not None:
            return DensityBoundary(density)
        return None

    @density_boundary.setter
    def density_boundary(self, value):
        """
        Parameters
        ----------
        value : :class:`~espressomd.EKSpecies.DensityBoundary` or None
            If value is :class:`~espressomd.EkSpecies.DensityBoundary`,
            set the node to be a boundary node with the specified density.
            If value is ``None``, the node will become a domain node.
        """

        if isinstance(value, DensityBoundary):
            self.call_method(
                "set_node_density_at_boundary",
                value=value.density)
        elif value is None:
            self.call_method("set_node_density_at_boundary", value=None)
        else:
            raise TypeError(
                "Parameter 'value' must be an instance of DensityBoundary or None")

    @property
    def flux_boundary(self):
        """
        Returns
        -------
        :class:`~espressomd.EKSpecies.FluxBoundary`
            If the node is a boundary node
        None
            If the node is not a boundary node
        """
        flux = self.call_method("get_node_flux_at_boundary")
        if flux is not None:
            return FluxBoundary(flux)
        return None

    @flux_boundary.setter
    def flux_boundary(self, value):
        """
        Parameters
        ----------
        value : :class:`~espressomd.EKSpecies.FluxBoundary` or None
            If value is :class:`~espressomd.EkSpecies.FluxBoundary`,
            set the node to be a boundary node with the specified flux.
            If value is ``None``, the node will become a domain node.
        """

        if isinstance(value, FluxBoundary):
            self.call_method("set_node_flux_at_boundary", value=value.flux)
        elif value is None:
            self.call_method("set_node_flux_at_boundary", value=None)
        else:
            raise TypeError(
                "Parameter 'value' must be an instance of FluxBoundary or None")


class EKSpeciesSlice(LatticeSliceWalberla):
    _node_class = EKSpeciesNode

    def __reduce__(self):
        raise NotImplementedError("Cannot serialize EK species slice objects")

    def setter_checks(self, attr, values):
        pass


@script_interface_register
class EKReactant(ScriptInterfaceHelper):
    _so_name = "walberla::EKReactant"
    _so_creation_policy = "GLOBAL"


class EKBulkReaction(ScriptInterfaceHelper):
    _so_name = "walberla::EKBulkReaction"
    _so_creation_policy = "GLOBAL"


class EKIndexedReaction(ScriptInterfaceHelper):
    _so_name = "walberla::EKIndexedReaction"
    _so_creation_policy = "GLOBAL"

    def add_node_to_index(self, node):
        self.call_method("set_node_is_boundary", node=node, is_boundary=True)

    def remove_node_from_index(self, node):
        self.call_method("set_node_is_boundary", node=node, is_boundary=False)

    def __getitem__(self, key):
        if isinstance(key, (tuple, list, np.ndarray)) and len(key) == 3:
            if any(isinstance(typ, slice) for typ in key):
                shape = self.shape

                indices = [np.atleast_1d(np.arange(shape[i])[key[i]])
                           for i in range(3)]
                dimensions = [ind.size for ind in indices]

                value_grid = np.zeros((*dimensions,), dtype=bool)
                indices = itertools.product(*map(enumerate, indices))
                for (i, x), (j, y), (k, z) in indices:
                    value_grid[i, j, k] = self.call_method(
                        "get_node_is_boundary", node=(x, y, z))

                return utils.array_locked(value_grid)
            else:
                return self.call_method("get_node_is_boundary", node=key)
        raise TypeError(
            f"{key} is not a valid index. Should be a point on the nodegrid or a slice")

    def __setitem__(self, key, values):
        if isinstance(key, (tuple, list, np.ndarray)) and len(key) == 3:
            if any(isinstance(typ, slice) for typ in key):
                shape = self.shape

                indices = [np.atleast_1d(np.arange(shape[i])[key[i]])
                           for i in range(3)]
                dimensions = tuple(ind.size for ind in indices)

                values = np.copy(values)

                # broadcast if only one element was provided
                if values.shape == ():
                    values = np.full(dimensions, values)
                if values.shape != dimensions:
                    raise ValueError(
                        f"Input-dimensions of array {values.shape} does not match slice dimensions {dimensions}.")

                indices = itertools.product(*map(enumerate, indices))
                for (i, x), (j, y), (k, z) in indices:
                    self.call_method("set_node_is_boundary", node=(
                        x, y, z), is_boundary=bool(values[i, j, k]))
            else:
                return self.call_method(
                    "set_node_is_boundary", node=key, is_boundary=values)
        else:
            raise TypeError(
                f"{key} is not a valid index. Should be a point on the nodegrid or a slice")


@script_interface_register
class EKReactions(ScriptObjectList):
    _so_name = "walberla::EKReactions"
    _so_creation_policy = "GLOBAL"

    def add(self, reaction):
        if not isinstance(reaction, (EKBulkReaction, EKIndexedReaction)):
            raise TypeError("reaction object is not of correct type.")

        self.call_method("add", object=reaction)

        return reaction

    def remove(self, reaction):
        self.call_method("remove", object=reaction)

    def clear(self):
        self.call_method("clear")
