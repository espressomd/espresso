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
import numpy as np


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
