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
from .script_interface import ScriptObjectRegistry, ScriptInterfaceHelper, script_interface_register
from .utils import check_type_or_throw_except
from .__init__ import has_features
import numpy as np


if has_features(["LB_BOUNDARIES"]):
    @script_interface_register
    class LBBoundaries(ScriptObjectRegistry):

        """
        Creates a set of lattice-Boltzmann boundaries.

        """

        _so_name = "LBBoundaries::LBBoundaries"

        def add(self, *args, **kwargs):
            """
            Adds a boundary to the set of boundaries.
            Either pass a valid boundary as argument,
            or a valid set of parameters to create a boundary.

            """

            if len(args) == 1:
                if isinstance(args[0], LBBoundary):
                    lbboundary = args[0]
                else:
                    raise TypeError(
                        "Either a LBBoundary object or key-value pairs for the parameters of a LBBoundary object need to be passed.")
            else:
                lbboundary = LBBoundary(**kwargs)
            self.call_method("add", object=lbboundary)
            return lbboundary

        def remove(self, lbboundary):
            """
            Removes a boundary from the set.

            Parameters
            ----------
            lbboundary : :obj:`LBBoundary`
                The boundary to be removed from the set.

            """

            self.call_method("remove", object=lbboundary)

        def clear(self):
            """
            Removes all boundaries.

            """

            self.call_method("clear")

        def size(self):
            return self.call_method("size")

        def empty(self):

            return self.call_method("empty")

    @script_interface_register
    class LBBoundary(ScriptInterfaceHelper):

        """
        Creates a LB boundary.

        """

        _so_name = "LBBoundaries::LBBoundary"
        _so_bind_methods = ("get_force",)

    class VelocityBounceBack:
        """
        Holds velocity information for the velocity bounce back boundary condition at a single node.
        """

        def __init__(self, velocity):
            check_type_or_throw_except(
                velocity, 3, float, "VelocityBounceBack velocity must be three floats")
            self.velocity = velocity


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
