#
# Copyright (C) 2010-2022 The ESPResSo project
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

import pathlib
import itertools
import numpy as np


def is_float(value):
    float_types = (float, np.float16, np.float32, np.float64, np.longdouble)
    if hasattr(np, 'float128'):
        float_types = float_types + (np.float128,)
    return isinstance(value, float_types)


def is_array(value):
    return isinstance(value, (list, tuple, np.ndarray))


def assert_params_match(ut_obj, inParams, outParams, msg_long=""):
    """Check if the parameters set and gotten back match.
    Only check keys present in ``inParams``.
    """
    if msg_long:
        msg_long = "\n" + msg_long

    for k, value_in in inParams.items():
        ut_obj.assertIn(k, outParams)
        value_out = outParams[k]
        msg = f"Mismatching parameter {k!r}{msg_long}"
        if is_array(value_out) or is_array(value_in):
            value_out = np.copy(value_out)
            value_in = np.copy(value_in)
            if is_float(value_in[0]) or is_float(value_out[0]):
                np.testing.assert_allclose(
                    value_out, value_in, atol=1e-14, err_msg=msg)
            else:
                np.testing.assert_array_equal(value_out, value_in, err_msg=msg)
        elif is_float(value_in) or is_float(value_out):
            ut_obj.assertAlmostEqual(value_out, value_in, delta=1e-14, msg=msg)
        else:
            ut_obj.assertEqual(value_out, value_in, msg=msg)


def generate_test_for_actor_class(_system, _class_actor, _params):
    """
    Generate a test case for an actor to verify parameters in the interface
    and the core match.
    Only keys which are present in ``_params`` are checked.
    """
    params_in = _params
    class_actor = _class_actor
    system = _system

    def func(self):
        # This code is run at the execution of the generated function.
        # It will use the state of the variables in the outer function,
        # which was there, when the outer function was called

        # set parameters
        actor = class_actor(**params_in)
        system.actors.add(actor)
        # read them out again
        params_out = {key: getattr(actor, key) for key in params_in}
        system.actors.remove(actor)

        assert_params_match(self, params_in, params_out,
                            msg_long=f"Parameters set {params_in} vs. {params_out}")

    return func


def lj_force_vector(v_d, d, lj_params):
    """Returns LJ force for distance ``d`` and distance vector ``v_d``
    based on the given ``lj_params``.
    Supports epsilon and cutoff.
    """

    if d >= lj_params["cutoff"]:
        return np.zeros(3)

    return 4. * lj_params["epsilon"] * v_d * (-12.0 * d**-14 + 6.0 * d**-8)


def verify_lj_forces(system, tolerance, ids_to_skip=()):
    """Goes over all pairs of particles in system and compares the forces
    on them to what would be expected based on the systems LJ parameters.
    Particle ids listed in ``ids_to_skip`` are not checked.
    Do not run this with a thermostat enabled.
    """

    # Initialize dict with expected forces
    f_expected = {}
    for pid in system.part.all().id:
        f_expected[pid] = np.zeros(3)

    # Cache some stuff to speed up pair loop
    dist_vec = system.distance_vec
    norm = np.linalg.norm
    non_bonded_inter = system.non_bonded_inter
    # LJ parameters
    lj_params = {}
    all_types = np.unique(system.part.all().type)
    for i in all_types:
        for j in all_types:
            lj_params[i, j] = non_bonded_inter[i, j].lennard_jones.get_params()

    # Go over all pairs of particles
    for pair in system.part.pairs():
        p0 = pair[0]
        p1 = pair[1]
        if p0.id in ids_to_skip or p1.id in ids_to_skip:
            continue

        # Distance and distance vec
        v_d = dist_vec(p0, p1)
        d = norm(v_d)

        # Calculate and add expected LJ force
        f = lj_force_vector(v_d, d, lj_params[p0.type, p1.type])
        f_expected[p0.id] += f
        f_expected[p1.id] -= f

    # Check actual forces against expected forces
    for p in system.part:
        if p.id in ids_to_skip:
            continue
        np.testing.assert_allclose(
            np.copy(p.f), f_expected[p.id], rtol=0, atol=tolerance,
            err_msg=f"LJ force verification failed on particle {p.id}.")


def data_path(filename):
    """
    Resolve absolute path to a resource.
    """
    return pathlib.Path(__file__).resolve().parent / "data" / filename


def transform_pos_from_cartesian_to_polar_coordinates(pos):
    """Transform the given cartesian coordinates to cylindrical coordinates.

    Parameters
    ----------
    pos : array_like :obj:`float`
        ``x``, ``y``, and ``z``-component of the cartesian position.

    Returns
    -------
    array_like
        The given position in polar coordinates.

    """
    return np.array([np.sqrt(pos[0]**2.0 + pos[1]**2.0),
                     np.arctan2(pos[1], pos[0]), pos[2]])


def get_cylindrical_basis_vectors(pos):
    phi = transform_pos_from_cartesian_to_polar_coordinates(pos)[1]
    e_r = np.array([np.cos(phi), np.sin(phi), 0.])
    e_phi = np.array([-np.sin(phi), np.cos(phi), 0.])
    e_z = np.array([0., 0., 1.])
    return e_r, e_phi, e_z


def convert_vec_body_to_space(p, vec):
    A = rotation_matrix_quat(p)
    return np.dot(A.transpose(), vec)


def rodrigues_rot(vec, axis, angle):
    """
    https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula#Statement
    """
    axis /= np.linalg.norm(axis)
    return np.cos(angle) * vec + np.sin(angle) * np.cross(axis, vec) + \
        (1 - np.cos(angle)) * np.dot(axis, vec) * axis


def rotation_matrix_quat(p):
    """
    Return the rotation matrix associated with quaternion.

    Parameters
    ----------
    p : :obj:`ParticleHandle`
        Particle.

    """
    A = np.zeros((3, 3))
    quat = p.quat
    qq = np.power(quat, 2)

    A[0, 0] = qq[0] + qq[1] - qq[2] - qq[3]
    A[1, 1] = qq[0] - qq[1] + qq[2] - qq[3]
    A[2, 2] = qq[0] - qq[1] - qq[2] + qq[3]

    A[0, 1] = 2 * (quat[1] * quat[2] + quat[0] * quat[3])
    A[0, 2] = 2 * (quat[1] * quat[3] - quat[0] * quat[2])
    A[1, 0] = 2 * (quat[1] * quat[2] - quat[0] * quat[3])

    A[1, 2] = 2 * (quat[2] * quat[3] + quat[0] * quat[1])
    A[2, 0] = 2 * (quat[1] * quat[3] + quat[0] * quat[2])
    A[2, 1] = 2 * (quat[2] * quat[3] - quat[0] * quat[1])

    return A


def normalize_cylindrical_hist(histogram, cyl_obs_params):
    """
    normalize a histogram in cylindrical coordinates. Helper to test the output
    of cylindrical histogram observables

    Parameters
    ----------
    histogram : (N,3) array_like of :obj:`float`
        The histogram that needs to be normalized
    cyl_obs_params : :obj:`dict`
        A dictionary containing the common parameters of the cylindrical histogram observables.
        Needs to contain the information about number and range of bins.
    """

    n_r_bins = cyl_obs_params['n_r_bins']
    n_phi_bins = cyl_obs_params['n_phi_bins']
    n_z_bins = cyl_obs_params['n_z_bins']
    min_r = cyl_obs_params['min_r']
    max_r = cyl_obs_params['max_r']
    min_phi = cyl_obs_params['min_phi']
    max_phi = cyl_obs_params['max_phi']
    min_z = cyl_obs_params['min_z']
    max_z = cyl_obs_params['max_z']

    bin_volume = np.zeros(n_r_bins)
    r_bin_size = (max_r - min_r) / n_r_bins
    phi_bin_size = (max_phi - min_phi) / n_phi_bins
    z_bin_size = (max_z - min_z) / n_z_bins
    for i in range(n_r_bins):
        bin_volume = np.pi * ((min_r + r_bin_size * (i + 1))**2.0 -
                              (min_r + r_bin_size * i)**2.0) * \
            phi_bin_size / (2.0 * np.pi) * z_bin_size
        histogram[i, :, :] /= bin_volume

    return histogram


def get_histogram(pos, obs_params, coord_system, **kwargs):
    """
    Helper function for ``np.histogramdd()`` and observables.

    Parameters
    ----------
    pos : (N, 3) array_like of :obj:`float`
        Particle positions.
    obs_params : :obj:`dict`
        Parameters of the observable.
    coord_system : :obj:`str`, {'cartesian', 'cylindrical'}
        Coordinate system.
    **kwargs
        Optional parameters to ``np.histogramdd()``.

    Returns
    -------
    array_like
        Bins and bin edges.

    """
    assert coord_system in ('cartesian', 'cylindrical'), \
        f"Unknown coord system '{coord_system}'"
    if coord_system == 'cartesian':
        bins = (obs_params['n_x_bins'],
                obs_params['n_y_bins'],
                obs_params['n_z_bins'])
        extent = [(obs_params['min_x'], obs_params['max_x']),
                  (obs_params['min_y'], obs_params['max_y']),
                  (obs_params['min_z'], obs_params['max_z'])]
    elif coord_system == 'cylindrical':
        bins = (obs_params['n_r_bins'],
                obs_params['n_phi_bins'],
                obs_params['n_z_bins'])
        extent = [(obs_params['min_r'], obs_params['max_r']),
                  (obs_params['min_phi'], obs_params['max_phi']),
                  (obs_params['min_z'], obs_params['max_z'])]
    return np.histogramdd(pos, bins=bins, range=extent, **kwargs)


# Generic Lennard-Jones


def lj_generic_potential(r, epsilon, sigma, cutoff, offset=0., shift=0.,
                         e1=12., e2=6., b1=4., b2=4., delta=0., lam=1.):
    r = np.array(r)
    V = np.zeros_like(r)
    cutoffMask = (r <= cutoff + offset)
    # LJGEN_SOFTCORE transformations
    rroff = np.sqrt(
        np.power(r[cutoffMask] - offset, 2) + (1 - lam) * delta * sigma**2)
    V[cutoffMask] = epsilon * lam * \
        (b1 * np.power(sigma / rroff, e1) -
         b2 * np.power(sigma / rroff, e2) + shift)
    return V


def lj_generic_force(espressomd, r, epsilon, sigma, cutoff, offset=0., e1=12,
                     e2=6, b1=4., b2=4., delta=0., lam=1., generic=True):
    f = 1.
    if r >= offset + cutoff:
        f = 0.
    else:
        h = (r - offset)**2 + delta * (1. - lam) * sigma**2
        f = (r - offset) * epsilon * lam * (
            b1 * e1 * np.power(sigma / np.sqrt(h), e1) - b2 * e2 * np.power(sigma / np.sqrt(h), e2)) / h
        if (not espressomd.has_features("LJGEN_SOFTCORE")) and generic:
            f *= np.sign(r - offset)
    return f

# Lennard-Jones


def lj_potential(r, epsilon, sigma, cutoff, shift, offset=0.):
    V = lj_generic_potential(
        r, epsilon, sigma, cutoff, offset=offset, shift=shift * 4.)
    return V


def lj_force(espressomd, r, epsilon, sigma, cutoff, offset=0.):
    f = lj_generic_force(
        espressomd, r, epsilon, sigma, cutoff, offset=offset, generic=False)
    return f


def fold_index(idx, shape):
    """Fold index into the range 0<x<shape"""

    res = np.copy(idx)
    for i in range(len(idx)):
        while res[i] >= shape[i]:
            res[i] -= shape[i]
        while res[i] < 0:
            res[i] += shape[i]
    return res


def get_lb_nodes_around_pos(pos, lbf):
    """Returns LB node(s) relevant for interpolation around the given position"""

    pos_lb_units = pos / lbf.agrid - 0.5  # relative to node centers
    lower_left_index = np.array(np.floor(pos_lb_units), dtype=int)

    nodes = []
    for i, j, k in itertools.product([0, 1], repeat=3):
        index = lower_left_index + np.array((i, j, k), dtype=int)
        nodes.append(lbf[fold_index(index, lbf.shape)])
    return nodes


def random_dipoles(n_particles):
    """Generate random dipoles by sampling Euler angles uniformly at random."""
    cos_theta = 2 * np.random.random(n_particles) - 1
    sin_theta = np.sin(np.arcsin(cos_theta))
    phi = 2 * np.pi * np.random.random(n_particles)
    dip = np.array([sin_theta * np.cos(phi),
                    sin_theta * np.sin(phi),
                    cos_theta]).T
    return dip


def check_non_bonded_loop_trace(ut_obj, system, cutoff=None):
    """Validates that the distances used by the non-bonded loop
    match with the minimum image distance accessible by Python,
    checks that no pairs are lost or double-counted.
    """

    cs_pairs = system.cell_system.non_bonded_loop_trace()
    # format [id1, id2, pos1, pos2, vec2, mpi_node]

    distance_vec = system.distance_vec
    if cutoff is None:
        cutoff = system.cell_system.max_cut_nonbonded

    # Distance for all pairs of particles obtained by Python
    py_distances = {}
    for p1, p2 in system.part.pairs():
        py_distances[p1.id, p2.id] = np.copy(distance_vec(p1, p2))

    # Go through pairs found by the non-bonded loop and check distance
    for p in cs_pairs:
        # p is a tuple with (id1, id2, pos1, pos2, vec21, mpi_node)
        # Note that system.distance_vec uses the opposite sign convention
        # as the minimum image distance in the core
        ut_obj.assertTrue(
            (p[0], p[1]) in py_distances or
            (p[1], p[0]) in py_distances,
            msg=f"Extra pair from core {p}")
        if (p[0], p[1]) in py_distances:
            np.testing.assert_allclose(
                np.copy(p[4]), -py_distances[p[0], p[1]])
            del py_distances[p[0], p[1]]
        elif (p[1], p[0]) in py_distances:
            np.testing.assert_allclose(
                np.copy(p[4]), py_distances[p[1], p[0]])
            del py_distances[p[1], p[0]]

    for ids, dist in py_distances.items():
        ut_obj.assertGreaterEqual(np.linalg.norm(dist), cutoff,
                                  msg=f"Pair not found by the core {ids}")
