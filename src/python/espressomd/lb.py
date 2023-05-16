#
# Copyright (C) 2013-2023 The ESPResSo project
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
from .detail.walberla import VTKOutputBase, LatticeWalberla
from .script_interface import ScriptInterfaceHelper, script_interface_register, array_variant
import espressomd.detail.walberla
import espressomd.shapes
import espressomd.code_features


class VelocityBounceBack:
    """
    Hold velocity information for the velocity bounce back boundary
    condition at a single node.

    """

    def __init__(self, velocity):
        utils.check_type_or_throw_except(
            velocity, 3, float, "VelocityBounceBack velocity must be three floats")
        self.velocity = velocity


class HydrodynamicInteraction(ScriptInterfaceHelper):
    """
    Base class for LB implementations.

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
        self.call_method("activate")
        utils.handle_errors("HydrodynamicInteraction activation failed")

    def _deactivate_method(self):
        self.call_method("deactivate")
        utils.handle_errors("HydrodynamicInteraction deactivation failed")

    def validate_params(self, params):
        pass

    def valid_keys(self):
        return {"agrid", "tau", "density", "ext_force_density",
                "kinematic_viscosity", "lattice", "kT", "seed"}

    def required_keys(self):
        return {"lattice", "density", "kinematic_viscosity", "tau"}

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
            raise ValueError(f"Slip velocity exceeds Mach {mach_number:.2f}")

    @property
    def pressure_tensor(self):
        tensor = self.call_method("get_pressure_tensor")
        return utils.array_locked(tensor)

    @pressure_tensor.setter
    def pressure_tensor(self, value):
        raise RuntimeError(f"Property 'pressure_tensor' is read-only")


@script_interface_register
class LBFluidWalberla(HydrodynamicInteraction,
                      espressomd.detail.walberla.LatticeModel):
    """
    The lattice-Boltzmann method for hydrodynamics using waLBerla.
    If argument ``lattice`` is not provided, one will be default
    constructed if an argument ``agrid`` is provided.

    Parameters
    ----------
    lattice : :obj:`espressomd.lb.LatticeWalberla <espressomd.detail.walberla.LatticeWalberla>`
        Lattice object. If not provided, a default one will be constructed
        using the ``agrid`` parameter.
    agrid : :obj:`float`
        Lattice constant. The box size in every direction must be an integer
        multiple of ``agrid``. Cannot be provided together with ``lattice``.
    tau : :obj:`float`
        LB time step, must be an integer multiple of the MD time step.
    density : :obj:`float`
        Fluid density.
    kinematic_viscosity : :obj:`float`
        Fluid kinematic viscosity.
    ext_force_density : (3,) array_like of :obj:`float`, optional
        Force density applied on the fluid.
    kT : :obj:`float`, optional
        Thermal energy of the simulated heat bath (for thermalized fluids).
        Set it to 0 for an unthermalized fluid.
    seed : :obj:`int`, optional
        Initial counter value (or seed) of the philox RNG.
        Required for a thermalized fluid. Must be positive.
    single_precision : :obj:`bool`, optional
        Use single-precision floating-point arithmetic.

    Methods
    -------
    get_interpolated_velocity()
        Get LB fluid velocity at specified position.

        Parameters
        ----------
        pos : (3,) array_like of :obj:`float`
            The position at which velocity is requested.

        Returns
        -------
        v : (3,) array_like :obj:`float`
            The LB fluid velocity at ``pos``.

    add_force_at_pos():
        Adds a force to the fluid at given position.

        Parameters
        ----------
        pos : (3,) array_like of :obj:`float`
            The position at which the force will be added.
        force : (3,) array_like of :obj:`float`
            The force vector which will be distributed at the position.

    clear_boundaries()
        Remove velocity bounce-back boundary conditions.

    save_checkpoint()
        Write LB node populations and boundary conditions to a file.

        Parameters
        ----------
        path : :obj:`str`
            Destination file path.
        binary : :obj:`bool`
            Whether to write in binary or ASCII mode.

    load_checkpoint()
        Load LB node populations and boundary conditions from a file.

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
        vtk : :class:`espressomd.lb.VTKOutput`
            VTK writer.

    remove_vtk_writer()
        Detach a VTK writer.

        Parameters
        ----------
        vtk : :class:`espressomd.lb.VTKOutput`
            VTK writer.

    clear_vtk_writers()
        Detach all VTK writers.

    """

    _so_name = "walberla::LBFluid"
    _so_creation_policy = "GLOBAL"
    _so_bind_methods = (
        "add_force_at_pos",
        "clear_boundaries",
        "get_interpolated_velocity",
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
            self.validate_params(params)
            super().__init__(*args, **params)
        else:
            super().__init__(**kwargs)

    def validate_params(self, params):
        super().validate_params(params)

        # construct default lattice if necessary
        if params.get("lattice") is None:
            if "agrid" not in params:
                raise ValueError("missing argument 'lattice' or 'agrid'")
            params["lattice"] = LatticeWalberla(
                agrid=params.pop("agrid"), n_ghost_layers=1)
        elif "agrid" in params:
            raise ValueError("cannot provide both 'lattice' and 'agrid'")

        utils.check_required_keys(self.required_keys(), params.keys())
        utils.check_valid_keys(self.valid_keys(), params.keys())

    def default_params(self):
        return {"single_precision": False, **super().default_params()}

    def valid_keys(self):
        return {"single_precision", *super().valid_keys()}

    def __getitem__(self, key):
        if isinstance(key, (tuple, list, np.ndarray)) and len(key) == 3:
            if any(isinstance(item, slice) for item in key):
                return LBFluidSliceWalberla(parent_sip=self, slice_range=key)
            else:
                return LBFluidNodeWalberla(
                    parent_sip=self, index=np.array(key))

        raise TypeError(
            f"{key} is not a valid index. Should be a point on the "
            "nodegrid e.g. lbf[0,0,0], or a slice e.g. lbf[:,0,0]")

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
            shape, 1, espressomd.shapes.Shape, "expected an espressomd.shapes.Shape")
        if np.shape(velocity) not in [(3,), tuple(self.shape) + (3,)]:
            raise ValueError(
                f'Cannot process velocity value grid of shape {np.shape(velocity)}')

        # range checks
        lattice_speed = self.call_method("get_lattice_speed")
        velocity = np.array(velocity, dtype=float).reshape((-1, 3))
        velocity *= 1. / lattice_speed
        self._check_mach_limit(velocity)

        mask = self.get_shape_bitmask(shape=shape).astype(int)
        self.call_method(
            "add_boundary_from_shape",
            raster=array_variant(mask.flatten()),
            values=array_variant(velocity.flatten()))


class LBFluidWalberlaGPU(HydrodynamicInteraction):
    """
    Initialize the lattice-Boltzmann method for hydrodynamic flow using
    waLBerla for the GPU. See :class:`HydrodynamicInteraction` for the
    list of parameters.

    """

    # pylint: disable=unused-argument
    def __init__(self, *args, **kwargs):
        if not espressomd.code_features.has_features("CUDA"):
            raise NotImplementedError("Feature CUDA not compiled in")
        if not espressomd.code_features.has_features("WALBERLA"):
            raise NotImplementedError("Feature WALBERLA not compiled in")
        raise NotImplementedError("Not implemented yet")


@script_interface_register
class LBFluidNodeWalberla(ScriptInterfaceHelper):
    _so_name = "walberla::LBFluidNode"
    _so_creation_policy = "GLOBAL"

    def required_keys(self):
        return {"parent_sip", "index"}

    def __init__(self, *args, **kwargs):
        if "sip" not in kwargs:
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
    def pressure_tensor_neq(self):
        tensor = self.call_method("get_pressure_tensor_neq")
        return utils.array_locked(tensor)

    @pressure_tensor_neq.setter
    def pressure_tensor_neq(self, value):
        raise RuntimeError("Property 'pressure_tensor_neq' is read-only.")

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
        value : :class:`~espressomd.lb.VelocityBounceBack` or ``None``
            If value is :class:`~espressomd.lb.VelocityBounceBack`,
            set the node to be a boundary node with the specified velocity.
            If value is ``None``, the node will become a fluid node.

        """

        if isinstance(value, VelocityBounceBack):
            value = value.velocity
            lattice_speed = self.call_method("get_lattice_speed")
            HydrodynamicInteraction._check_mach_limit(
                np.array(value) / lattice_speed)
        elif value is not None:
            raise TypeError(
                "Parameter 'value' must be an instance of VelocityBounceBack or None")
        self.call_method("set_velocity_at_boundary", value=value)

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


@script_interface_register
class LBFluidSliceWalberla(ScriptInterfaceHelper):
    _so_name = "walberla::LBFluidSlice"
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
            node = LBFluidNodeWalberla(index=np.array([0, 0, 0]), **kwargs)
            super().__init__(*args, node_sip=node, **kwargs, **extra_kwargs)
            utils.handle_errors("LBFluidSliceWalberla instantiation failed")

    def __iter__(self):
        lower, upper = self.call_method("get_slice_ranges")
        indices = [list(range(lower[i], upper[i])) for i in range(3)]
        lb_sip = self.call_method("get_lb_sip")
        for index in itertools.product(*indices):
            yield LBFluidNodeWalberla(parent_sip=lb_sip, index=np.array(index))

    def __reduce__(self):
        raise NotImplementedError("Cannot serialize LB fluid slice objects")

    def _getter(self, attr):
        value_grid, shape = self.call_method(f"get_{attr}")
        if attr == "velocity_at_boundary":
            value_grid = [
                None if x is None else VelocityBounceBack(x) for x in value_grid]
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
    def population(self):
        return self._getter("population")

    @population.setter
    def population(self, value):
        self._setter("population", value)

    @property
    def pressure_tensor(self):
        return self._getter("pressure_tensor")

    @pressure_tensor.setter
    def pressure_tensor(self, value):
        raise RuntimeError("Property 'pressure_tensor' is read-only.")

    @property
    def pressure_tensor_neq(self):
        return self._getter("pressure_tensor_neq")

    @pressure_tensor_neq.setter
    def pressure_tensor_neq(self, value):
        raise RuntimeError("Property 'pressure_tensor_neq' is read-only.")

    @property
    def is_boundary(self):
        return self._getter("is_boundary")

    @is_boundary.setter
    def is_boundary(self, value):
        raise RuntimeError("Property 'is_boundary' is read-only.")

    @property
    def boundary(self):
        """
        Returns
        -------
        (N, M, L) array_like of :class:`~espressomd.lb.VelocityBounceBack`
            If the nodes are boundary nodes
        (N, M, L) array_like of ``None``
            If the nodes are not boundary nodes

        """

        return self._getter("velocity_at_boundary")

    @boundary.setter
    def boundary(self, values):
        """
        Parameters
        ----------
        values : (N, M, L) array_like of :class:`~espressomd.lb.VelocityBounceBack` or ``None``
            If values are :class:`~espressomd.lb.VelocityBounceBack`,
            set the nodes to be boundary nodes with the specified velocity.
            If values are ``None``, the nodes will become fluid nodes.

        """

        type_error_msg = "Parameter 'values' must be an array_like of VelocityBounceBack or None"
        values = np.copy(values)
        lattice_speed = self.call_method("get_lattice_speed")
        if values.dtype != np.dtype("O"):
            raise TypeError(type_error_msg)
        for index in np.ndindex(*values.shape):
            if values[index] is not None:
                if not isinstance(values[index], VelocityBounceBack):
                    raise TypeError(type_error_msg)
                HydrodynamicInteraction._check_mach_limit(
                    np.array(values[index].velocity) / lattice_speed)
                values[index] = np.array(values[index].velocity)
        self._setter("velocity_at_boundary", values=values)

    @property
    def boundary_force(self):
        return self._getter("boundary_force")

    @boundary_force.setter
    def boundary_force(self, value):
        raise RuntimeError("Property 'boundary_force' is read-only.")

    @property
    def velocity(self):
        return self._getter("velocity")

    @velocity.setter
    def velocity(self, value):
        self._setter("velocity", value)

    @property
    def last_applied_force(self):
        return self._getter("last_applied_force")

    @last_applied_force.setter
    def last_applied_force(self, value):
        self._setter("last_applied_force", value)


@script_interface_register
class VTKOutput(VTKOutputBase):
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
    _so_name = "walberla::LBVTKHandle"
    _so_creation_policy = "GLOBAL"
    _so_bind_methods = ("enable", "disable", "write")

    def required_keys(self):
        return self.valid_keys() - self.default_params().keys()

    def __repr__(self):
        class_id = f"{self.__class__.__module__}.{self.__class__.__name__}"
        if self.delta_N:
            write_when = f"every {self.delta_N} LB steps"
            if not self.enabled:
                write_when += " (disabled)"
        else:
            write_when = "on demand"
        return f"<{class_id}: write to '{self.vtk_uid}' {write_when}>"


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
    wrapped_mask = np.pad(fluid_mask.astype(int), 3 * [(2, 2)], mode="wrap")
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
        wrapped_mask, edge, mode="same", method="direct")[2:-2, 2:-2, 2:-2]
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
