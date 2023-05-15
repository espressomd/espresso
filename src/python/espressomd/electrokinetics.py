#
# Copyright (C) 2021-2023 The ESPResSo project
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

import itertools
import numpy as np

from . import utils
from .detail.walberla import VTKOutputBase, LatticeWalberla  # pylint: disable=unused-import
from .script_interface import ScriptInterfaceHelper, script_interface_register, ScriptObjectList, array_variant
import espressomd.detail.walberla
import espressomd.shapes
import espressomd.code_features


@script_interface_register
class EKFFT(ScriptInterfaceHelper):
    """
    A FFT-based Poisson solver.

    """

    _so_name = "walberla::EKFFT"
    _so_creation_policy = "GLOBAL"

    def __init__(self, *args, **kwargs):
        if not espressomd.code_features.has_features("WALBERLA_FFT"):
            raise NotImplementedError("Feature WALBERLA not compiled in")

        super().__init__(*args, **kwargs)


@script_interface_register
class EKNone(ScriptInterfaceHelper):
    """
    The default Poisson solver.
    Imposes a null electrostatic potential everywhere.

    """
    _so_name = "walberla::EKNone"
    _so_creation_policy = "GLOBAL"

    def __init__(self, *args, **kwargs):
        if not espressomd.code_features.has_features("WALBERLA"):
            raise NotImplementedError("Feature WALBERLA not compiled in")

        super().__init__(*args, **kwargs)


@script_interface_register
class EKContainer(ScriptObjectList):
    _so_name = "walberla::EKContainer"

    def __init__(self, *args, **kwargs):
        if not espressomd.code_features.has_features("WALBERLA"):
            raise NotImplementedError("Feature WALBERLA not compiled in")

        super().__init__(*args, **kwargs)

    def add(self, ekspecies):
        self.call_method("add", object=ekspecies)

    def remove(self, ekspecies):
        self.call_method("remove", object=ekspecies)

    def clear(self):
        self.call_method("clear")


@script_interface_register
class EKSpecies(ScriptInterfaceHelper,
                espressomd.detail.walberla.LatticeModel):
    """
    The advection-diffusion-reaction method for chemical species using waLBerla.

    Parameters
    ----------
    lattice : :obj:`espressomd.electrokinetics.LatticeWalberla <espressomd.detail.walberla.LatticeWalberla>`
        Lattice object.
    tau : :obj:`float`
        EK time step, must be an integer multiple of the MD time step.
    density : :obj:`float`
        Species density.
    diffusion : :obj:`float`
        Species diffusion coefficient.
    valency : :obj:`float`
        Species valency.
    advection : :obj:`bool`
        Whether to enable advection.
    friction_coupling : :obj:`bool`
        Whether to enable friction coupling.
    ext_efield : (3,) array_like of :obj:`float`, optional
        External electrical field.
    kT : :obj:`float`, optional
        Thermal energy of the simulated heat bath (for thermalized species).
        Set it to 0 for an unthermalized species.
    single_precision : :obj:`bool`, optional
        Use single-precision floating-point arithmetic.

    Methods
    -------
    clear_density_boundaries()
        Remove density boundary conditions.

    clear_flux_boundaries()
        Remove flux boundary conditions.

    clear_boundaries()
        Remove all boundary conditions.

    save_checkpoint()
        Write EK densities and boundary conditions to a file.

        Parameters
        ----------
        path : :obj:`str`
            Destination file path.
        binary : :obj:`bool`
            Whether to write in binary or ASCII mode.

    load_checkpoint()
        Load EK densities and boundary conditions from a file.

        Parameters
        ----------
        path : :obj:`str`
            File path to read from.
        binary : :obj:`bool`
            Whether to read in binary or ASCII mode.

    add_vtk_writer()
        Attach a VTK writer.

        Parameters
        ----------
        vtk : :class:`espressomd.electrokinetics.VTKOutput`
            VTK writer.

    remove_vtk_writer()
        Detach a VTK writer.

        Parameters
        ----------
        vtk : :class:`espressomd.electrokinetics.VTKOutput`
            VTK writer.

    clear_vtk_writers()
        Detach all VTK writers.

    """

    _so_name = "walberla::EKSpecies"
    _so_creation_policy = "GLOBAL"
    _so_bind_methods = (
        "clear_density_boundaries",
        "clear_flux_boundaries",
        "clear_boundaries",
        "add_vtk_writer",
        "remove_vtk_writer",
        "clear_vtk_writers",
    )

    def __init__(self, *args, **kwargs):
        if not espressomd.code_features.has_features("WALBERLA"):
            raise NotImplementedError("Feature WALBERLA not compiled in")

        if "sip" not in kwargs:
            params = self.default_params()
            params.update(kwargs)
            super().__init__(*args, **params)
        else:
            super().__init__(**kwargs)

    def default_params(self):
        return {"single_precision": False,
                "kT": 0., "ext_efield": [0., 0., 0.]}

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

    def add_boundary_from_shape(self, shape, value, boundary_type):
        """
        Set boundary conditions from a shape.

        Parameters
        ----------
        shape : :obj:`espressomd.shapes.Shape`
            Shape to rasterize.
        value : (O,) or (L, M, N, O) array_like of :obj:`float`, optional
            Boundary numerical value. If a single value of shape ``(O,)``
            is given, it will be broadcast to all nodes inside the shape,
            otherwise ``L, M, N`` must be equal to the EK grid dimensions.
        boundary_type : Union[:class:`~espressomd.electrokinetics.DensityBoundary`,
                              :class:`~espressomd.electrokinetics.FluxBoundary`] (optional)
            Type of the boundary condition.

        """
        if not issubclass(boundary_type, (FluxBoundary, DensityBoundary)):
            raise TypeError(
                "Parameter 'boundary_type' must be a subclass of FluxBoundary or DensityBoundary")

        if not hasattr(value, "__iter__"):
            value = (value, )

        value = np.array(value, dtype=float)
        utils.check_type_or_throw_except(
            shape, 1, espressomd.shapes.Shape, "expected an espressomd.shapes.Shape")
        if issubclass(boundary_type, FluxBoundary):
            if np.shape(value) not in [(3,), tuple(self.shape) + (3,)]:
                raise ValueError(
                    f"Cannot process flux value grid of shape {np.shape(value)}")
        if issubclass(boundary_type, DensityBoundary):
            if np.shape(value) not in [(1,), tuple(self.shape) + (1,)]:
                raise ValueError(
                    f"Cannot process density value grid of shape {np.shape(value)}")

        mask = self.get_shape_bitmask(shape=shape).astype(int)
        if issubclass(boundary_type, FluxBoundary):
            boundaries_update_method = "update_flux_boundary_from_shape"
        else:
            boundaries_update_method = "update_density_boundary_from_shape"
        self.call_method(
            boundaries_update_method,
            raster=array_variant(mask.flatten()),
            values=array_variant(value.flatten()))


class FluxBoundary:
    """
    Hold flux information for the flux boundary
    condition at a single node.

    """

    def __init__(self, flux):
        utils.check_type_or_throw_except(
            flux, 3, float, "FluxBoundary flux must be three floats")
        self.flux = flux


class DensityBoundary:
    """
    Hold density information for the density boundary
    condition at a single node.

    """

    def __init__(self, density):
        utils.check_type_or_throw_except(
            density, 1, float, "DensityBoundary flux must be one float")
        self.density = density


@script_interface_register
class EKSpeciesNode(ScriptInterfaceHelper):
    _so_name = "walberla::EKSpeciesNode"
    _so_creation_policy = "GLOBAL"

    def required_keys(self):
        return {"parent_sip", "index"}

    def validate_params(self, params):
        utils.check_required_keys(self.required_keys(), params.keys())
        utils.check_type_or_throw_except(
            params["index"], 3, int, "The index of an EK species node consists of three integers.")

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
        :class:`~espressomd.electrokinetics.DensityBoundary`
            If the node is a boundary node
        ``None``
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
        value : :class:`~espressomd.electrokinetics.DensityBoundary` or ``None``
            If value is :class:`~espressomd.EkSpecies.DensityBoundary`,
            set the node to be a boundary node with the specified density.
            If value is ``None``, the node will become a domain node.

        """

        if isinstance(value, DensityBoundary):
            value = value.density
        elif value is not None:
            raise TypeError(
                "Parameter 'value' must be an instance of DensityBoundary or None")
        self.call_method("set_node_density_at_boundary", value=value)

    @property
    def flux_boundary(self):
        """
        Returns
        -------
        :class:`~espressomd.electrokinetics.FluxBoundary`
            If the node is a boundary node
        ``None``
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
        value : :class:`~espressomd.electrokinetics.FluxBoundary` or ``None``
            If value is :class:`~espressomd.EkSpecies.FluxBoundary`,
            set the node to be a boundary node with the specified flux.
            If value is ``None``, the node will become a domain node.

        """

        if isinstance(value, FluxBoundary):
            value = value.flux
        elif value is not None:
            raise TypeError(
                "Parameter 'value' must be an instance of FluxBoundary or None")
        self.call_method("set_node_flux_at_boundary", value=value)


@script_interface_register
class EKSpeciesSlice(ScriptInterfaceHelper):
    _so_name = "walberla::EKSpeciesSlice"
    _so_creation_policy = "GLOBAL"

    def required_keys(self):
        return {"parent_sip", "slice_range"}

    def validate_params(self, params):
        utils.check_required_keys(self.required_keys(), params.keys())

    def __init__(self, *args, **kwargs):
        if "sip" in kwargs:
            super().__init__(**kwargs)
        else:
            self.validate_params(kwargs)
            slice_range = kwargs.pop("slice_range")
            grid_size = kwargs["parent_sip"].shape
            extra_kwargs = espressomd.detail.walberla.get_slice_bounding_box(
                slice_range, grid_size)
            node = EKSpeciesNode(index=np.array([0, 0, 0]), **kwargs)
            super().__init__(*args, node_sip=node, **kwargs, **extra_kwargs)
            utils.handle_errors("EKSpeciesSlice instantiation failed")

    def __iter__(self):
        lower, upper = self.call_method("get_slice_ranges")
        indices = [list(range(lower[i], upper[i])) for i in range(3)]
        lb_sip = self.call_method("get_ek_sip")
        for index in itertools.product(*indices):
            yield EKSpeciesNode(parent_sip=lb_sip, index=np.array(index))

    def __reduce__(self):
        raise NotImplementedError("Cannot serialize EK species slice objects")

    def _getter(self, attr):
        value_grid, shape = self.call_method(f"get_{attr}")
        if attr == "flux_at_boundary":
            value_grid = [
                None if x is None else FluxBoundary(x) for x in value_grid]
        elif attr == "density_at_boundary":
            value_grid = [
                None if x is None else DensityBoundary(x) for x in value_grid]
        return utils.array_locked(np.reshape(value_grid, shape))

    def _setter(self, attr, values):
        dimensions = self.call_method("get_slice_size")
        if 0 in dimensions:
            raise AttributeError(
                f"Cannot set properties of an empty '{self.__class__.__name__}' object")

        values = np.copy(values)
        value_shape = tuple(self.call_method("get_value_shape", name=attr))
        target_shape = (*dimensions, *value_shape)

        # broadcast if only one element was provided
        if values.shape == value_shape or values.shape == () and value_shape == (1,):
            values = np.full(target_shape, values)

        def shape_squeeze(shape):
            return tuple(x for x in shape if x != 1)

        if shape_squeeze(values.shape) != shape_squeeze(target_shape):
            raise ValueError(
                f"Input-dimensions of '{attr}' array {values.shape} does not match slice dimensions {target_shape}")

        self.call_method(f"set_{attr}", values=values.flatten())

    @property
    def density(self):
        return self._getter("density",)

    @density.setter
    def density(self, value):
        self._setter("density", value)

    @property
    def is_boundary(self):
        return self._getter("is_boundary")

    @is_boundary.setter
    def is_boundary(self, value):
        raise RuntimeError("Property 'is_boundary' is read-only.")

    @property
    def density_boundary(self):
        """
        Returns
        -------
        (N, M, L) array_like of :class:`~espressomd.electrokinetics.DensityBoundary`
            If the nodes are boundary nodes
        (N, M, L) array_like of ``None``
            If the nodes are not boundary nodes

        """

        return self._getter("density_at_boundary")

    @density_boundary.setter
    def density_boundary(self, values):
        """
        Parameters
        ----------
        values : (N, M, L) array_like of :class:`~espressomd.electrokinetics.DensityBoundary` or obj:`None`
            If values are :class:`~espressomd.electrokinetics.DensityBoundary`,
            set the nodes to be boundary nodes with the specified density.
            If values are obj:`None`, the nodes will become domain nodes.

        """

        type_error_msg = "Parameter 'values' must be an array_like of DensityBoundary or None"
        values = np.copy(values)
        if values.dtype != np.dtype("O"):
            raise TypeError(type_error_msg)
        for index in np.ndindex(*values.shape):
            if values[index] is not None:
                if not isinstance(values[index], DensityBoundary):
                    raise TypeError(type_error_msg)
                values[index] = np.array(values[index].density)
        self._setter("density_at_boundary", values=values)

    @property
    def flux_boundary(self):
        """
        Returns
        -------
        (N, M, L) array_like of :class:`~espressomd.electrokinetics.FluxBoundary`
            If the nodes are boundary nodes
        (N, M, L) array_like of `None``
            If the nodes are not boundary nodes

        """

        return self._getter("flux_at_boundary")

    @flux_boundary.setter
    def flux_boundary(self, values):
        """
        Parameters
        ----------
        values : (N, M, L) array_like of :class:`~espressomd.electrokinetics.FluxBoundary` or obj:`None`
            If values are :class:`~espressomd.lb.FluxBoundary`,
            set the nodes to be boundary nodes with the specified flux.
            If values are obj:`None`, the nodes will become domain nodes.

        """

        type_error_msg = "Parameter 'values' must be an array_like of FluxBoundary or None"
        values = np.copy(values)
        if values.dtype != np.dtype("O"):
            raise TypeError(type_error_msg)
        for index in np.ndindex(*values.shape):
            if values[index] is not None:
                if not isinstance(values[index], FluxBoundary):
                    raise TypeError(type_error_msg)
                values[index] = np.array(values[index].flux)
        self._setter("flux_at_boundary", values=values)


@script_interface_register
class VTKOutput(VTKOutputBase):
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

    def required_keys(self):
        return self.valid_keys() - self.default_params().keys()

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
