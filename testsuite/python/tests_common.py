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
import os
import numpy as np


def assert_params_match(ut_obj, inParams, outParams, msg_long=None):
    """Check if the parameters set and gotten back match.
    Only check keys present in ``inParams``.
    """
    if msg_long:
        msg_long = "\n" + msg_long
    else:
        msg_long = ""

    for k in inParams.keys():
        ut_obj.assertIn(k, outParams)
        if isinstance(inParams[k], float):
            ut_obj.assertAlmostEqual(
                outParams[k], inParams[k], delta=1E-14,
                msg=f"Mismatching parameter {k!r}{msg_long}")
        else:
            ut_obj.assertEqual(
                outParams[k], inParams[k],
                msg=f"Mismatching parameter {k!r}{msg_long}")


def generate_test_for_class(_system, _interClass, _params):
    """Generates test cases for checking interaction parameters set and gotten back
    from Es actually match. Only keys which are present in _params are checked
    1st: Interaction parameters as dictionary, i.e., ``{"k": 1., "r_0": 0}``
    2nd: Name of the interaction property to set (i.e. ``"P3M"``)
    """
    params = _params
    interClass = _interClass
    system = _system

    def func(self):
        # This code is run at the execution of the generated function.
        # It will use the state of the variables in the outer function,
        # which was there, when the outer function was called

        # set Parameter
        Inter = interClass(**params)
        Inter.validate_params()
        system.actors.add(Inter)
        # Read them out again
        outParams = Inter.get_params()
        del system.actors[0]

        assert_params_match(self, params, outParams,
                            f"Parameters set {params} vs. {outParams}")

    return func


def lj_force_vector(v_d, d, lj_params):
    """Returns lj force for distance d and distance vector v_d based on the given lj_params.
    Supports epsilon and cutoff."""

    if d >= lj_params["cutoff"]:
        return np.zeros(3)

    return 4. * lj_params["epsilon"] * v_d * (-12.0 * d**-14 + 6.0 * d**-8)


def verify_lj_forces(system, tolerance, ids_to_skip=()):
    """Goes over all pairs of particles in system and compares the forces on them
       to what would be expected based on the systems LJ parameters.
       Particle ids listed in ids_to_skip are not checked
       Do not run this with a thermostat enabled."""

    # Initialize dict with expected forces
    f_expected = {}
    for pid in system.part[:].id:
        f_expected[pid] = np.zeros(3)

    # Cache some stuff to speed up pair loop
    dist_vec = system.distance_vec
    norm = np.linalg.norm
    non_bonded_inter = system.non_bonded_inter
    # lj parameters
    lj_params = {}
    all_types = np.unique(system.part[:].type)
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

        # calc and add expected lj force
        f = lj_force_vector(v_d, d, lj_params[p0.type, p1.type])
        f_expected[p0.id] += f
        f_expected[p1.id] -= f

    # Check actual forces against expected
    for p in system.part:
        if p.id in ids_to_skip:
            continue
        if np.linalg.norm(p.f - f_expected[p.id]) >= tolerance:
            raise Exception(f"LJ force verification failed on particle "
                            f"{p.id}. Got {p.f}, expected {f_expected[p.id]}")


def abspath(path):
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), path)


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


def transform_vel_from_cartesian_to_polar_coordinates(pos, vel):
    """Transform the given cartesian velocities to polar velocities.

    Parameters
    ----------
    pos : array_like :obj:`float`
        ``x``, ``y``, and ``z``-component of the cartesian position.
    vel : array_like :obj:`float`
        ``x``, ``y``, and ``z``-component of the cartesian velocity.

    """
    return np.array([
        (pos[0] * vel[0] + pos[1] * vel[1]) / np.sqrt(pos[0]**2 + pos[1]**2),
        (pos[0] * vel[1] - pos[1] * vel[0]) / np.sqrt(pos[0]**2 + pos[1]**2), vel[2]])


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
    coord_system : :obj:`str`, \{'cartesian', 'cylindrical'\}
        Coordinate system.
    \*\*kwargs :
        Optional parameters to ``np.histogramdd()``.

    Returns
    -------
    array_like
        Bins and bin edges.

    """
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
    else:
        raise ValueError(f"Unknown coord system '{coord_system}'")
    return np.histogramdd(pos, bins=bins, range=extent, **kwargs)

#
# Analytical Expressions for interactions
#

# Harmonic bond


def harmonic_potential(scalar_r, k, r_0):
    return 0.5 * k * (scalar_r - r_0)**2


def harmonic_force(scalar_r, k, r_0):
    return -k * (scalar_r - r_0)

# FENE bond


def fene_potential(scalar_r, k, d_r_max, r_0):
    return -0.5 * k * d_r_max**2 * np.log(1 - ((scalar_r - r_0) / d_r_max)**2)


def fene_force(scalar_r, k, d_r_max, r_0):
    return k * (scalar_r - r_0) * d_r_max**2 / \
        ((scalar_r - r_0)**2 - d_r_max**2)


def fene_force2(bond_vector, k, d_r_max, r_0):
    r = np.linalg.norm(bond_vector)
    return k * (r - r_0) / (r * (1 - ((r - r_0) / d_r_max)**2)) * \
        np.array(bond_vector)

# Coulomb bond


def coulomb_potential(scalar_r, k, q1, q2):
    return k * q1 * q2 / scalar_r


def coulomb_force(scalar_r, k, q1, q2):
    return k * q1 * q2 / scalar_r**2

# QUARTIC bond


def quartic_force(k0, k1, r, r_cut, scalar_r):
    if scalar_r > r_cut:
        return 0.0
    return - k0 * (scalar_r - r) - k1 * (scalar_r - r)**3


def quartic_potential(k0, k1, r, r_cut, scalar_r):
    if scalar_r > r_cut:
        return 0.0
    return 0.5 * k0 * (scalar_r - r)**2 + 0.25 * k1 * (scalar_r - r)**4

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

# Lennard-Jones Cosine


def lj_cos_potential(r, epsilon, sigma, cutoff, offset):
    V = 0.
    r_min = offset + np.power(2., 1. / 6.) * sigma
    r_cut = cutoff + offset
    if r < r_min:
        V = lj_potential(r, epsilon=epsilon, sigma=sigma,
                         cutoff=cutoff, offset=offset, shift=0.)
    elif r < r_cut:
        alpha = np.pi / \
            (np.power(r_cut - offset, 2) - np.power(r_min - offset, 2))
        beta = np.pi - np.power(r_min - offset, 2) * alpha
        V = 0.5 * epsilon * \
            (np.cos(alpha * np.power(r - offset, 2) + beta) - 1.)
    return V


def lj_cos_force(espressomd, r, epsilon, sigma, cutoff, offset):
    f = 0.
    r_min = offset + np.power(2., 1. / 6.) * sigma
    r_cut = cutoff + offset
    if r < r_min:
        f = lj_force(espressomd, r, epsilon=epsilon, sigma=sigma,
                     cutoff=cutoff, offset=offset)
    elif r < r_cut:
        alpha = np.pi / \
            (np.power(r_cut - offset, 2) - np.power(r_min - offset, 2))
        beta = np.pi - np.power(r_min - offset, 2) * alpha
        f = (r - offset) * alpha * epsilon * \
            np.sin(alpha * np.power(r - offset, 2) + beta)
    return f

# Lennard-Jones Cosine^2


def lj_cos2_potential(r, epsilon, sigma, offset, width):
    V = 0.
    r_min = offset + np.power(2., 1. / 6.) * sigma
    r_cut = r_min + width
    if r < r_min:
        V = lj_potential(r, epsilon=epsilon, sigma=sigma,
                         offset=offset, cutoff=r_cut, shift=0.)
    elif r < r_cut:
        V = -epsilon * np.power(np.cos(np.pi /
                                       (2. * width) * (r - r_min)), 2)
    return V


def lj_cos2_force(espressomd, r, epsilon, sigma, offset, width):
    f = 0.
    r_min = offset + np.power(2., 1. / 6.) * sigma
    r_cut = r_min + width
    if r < r_min:
        f = lj_force(espressomd, r, epsilon=epsilon,
                     sigma=sigma, cutoff=r_cut, offset=offset)
    elif r < r_cut:
        f = - np.pi * epsilon * \
            np.sin(np.pi * (r - r_min) / width) / (2. * width)
    return f

# Smooth-Step


def smooth_step_potential(r, eps, sig, cutoff, d, n, k0):
    V = 0.
    if r < cutoff:
        V = np.power(d / r, n) + eps / \
            (1 + np.exp(2 * k0 * (r - sig)))
    return V


def smooth_step_force(r, eps, sig, cutoff, d, n, k0):
    f = 0.
    if r < cutoff:
        f = n * d / r**2 * np.power(d / r, n - 1) + 2 * k0 * eps * np.exp(
            2 * k0 * (r - sig)) / (1 + np.exp(2 * k0 * (r - sig))**2)
    return f

# BMHTF


def bmhtf_potential(r, a, b, c, d, sig, cutoff):
    V = 0.
    if r == cutoff:
        V = a * np.exp(b * (sig - r)) - c * np.power(
            r, -6) - d * np.power(r, -8)
    if r < cutoff:
        V = a * np.exp(b * (sig - r)) - c * np.power(
            r, -6) - d * np.power(r, -8)
        V -= bmhtf_potential(cutoff, a, b, c, d, sig, cutoff)
    return V


def bmhtf_force(r, a, b, c, d, sig, cutoff):
    f = 0.
    if r < cutoff:
        f = a * b * np.exp(b * (sig - r)) - 6 * c * np.power(
            r, -7) - 8 * d * np.power(r, -9)
    return f

# Morse


def morse_potential(r, eps, alpha, cutoff, rmin=0):
    V = 0.
    if r < cutoff:
        V = eps * (np.exp(-2. * alpha * (r - rmin)) -
                   2 * np.exp(-alpha * (r - rmin)))
        V -= eps * (np.exp(-2. * alpha * (cutoff - rmin)) -
                    2 * np.exp(-alpha * (cutoff - rmin)))
    return V


def morse_force(r, eps, alpha, cutoff, rmin=0):
    f = 0.
    if r < cutoff:
        f = 2. * np.exp((rmin - r) * alpha) * \
            (np.exp((rmin - r) * alpha) - 1) * alpha * eps
    return f

# Buckingham


def buckingham_potential(r, a, b, c, d, cutoff, discont, shift):
    V = 0.
    if r < discont:
        m = - buckingham_force(
            discont, a, b, c, d, cutoff, discont, shift)
        c = buckingham_potential(
            discont, a, b, c, d, cutoff, discont, shift) - m * discont
        V = m * r + c
    if (r >= discont) and (r < cutoff):
        V = a * np.exp(- b * r) - c * np.power(
            r, -6) - d * np.power(r, -4) + shift
    return V


def buckingham_force(r, a, b, c, d, cutoff, discont, shift):
    f = 0.
    if r < discont:
        f = buckingham_force(
            discont, a, b, c, d, cutoff, discont, shift)
    if (r >= discont) and (r < cutoff):
        f = a * b * np.exp(- b * r) - 6 * c * np.power(
            r, -7) - 4 * d * np.power(r, -5)
    return f

# Soft-sphere


def soft_sphere_potential(r, a, n, cutoff, offset=0):
    V = 0.
    if r < offset + cutoff:
        V = a * np.power(r - offset, -n)
    return V


def soft_sphere_force(r, a, n, cutoff, offset=0):
    f = 0.
    if (r > offset) and (r < offset + cutoff):
        f = n * a * np.power(r - offset, -(n + 1))
    return f

# Hertzian


def hertzian_potential(r, eps, sig):
    V = 0.
    if r < sig:
        V = eps * np.power(1 - r / sig, 5. / 2.)
    return V


def hertzian_force(r, eps, sig):
    f = 0.
    if r < sig:
        f = 5. / 2. * eps / sig * np.power(1 - r / sig, 3. / 2.)
    return f

# Gaussian


def gaussian_potential(r, eps, sig, cutoff):
    V = 0.
    if r < cutoff:
        V = eps * np.exp(-np.power(r / sig, 2) / 2)
    return V


def gaussian_force(r, eps, sig, cutoff):
    f = 0.
    if r < cutoff:
        f = eps * r / sig**2 * np.exp(-np.power(r / sig, 2) / 2)
    return f


def gay_berne_potential(r_ij, u_i, u_j, epsilon_0, sigma_0, mu, nu, k_1, k_2):
    r_normed = r_ij / np.linalg.norm(r_ij)
    r_u_i = np.dot(r_normed, u_i)
    r_u_j = np.dot(r_normed, u_j)
    u_i_u_j = np.dot(u_i, u_j)

    chi = (k_1**2 - 1.) / (k_1**2 + 1.)
    chi_d = (k_2**(1. / mu) - 1) / (k_2**(1. / mu) + 1)

    sigma = sigma_0 \
        / np.sqrt(
            (1 - 0.5 * chi * (
                (r_u_i + r_u_j)**2 / (1 + chi * u_i_u_j) +
                (r_u_i - r_u_j)**2 / (1 - chi * u_i_u_j))))

    epsilon = epsilon_0 *\
        (1 - chi**2 * u_i_u_j**2)**(-nu / 2.) *\
        (1 - chi_d / 2. * (
            (r_u_i + r_u_j)**2 / (1 + chi_d * u_i_u_j) +
            (r_u_i - r_u_j)**2 / (1 - chi_d * u_i_u_j)))**mu

    rr = np.linalg.norm((np.linalg.norm(r_ij) - sigma + sigma_0) / sigma_0)

    return 4. * epsilon * (rr**-12 - rr**-6)


def count_fluid_nodes(lbf):
    """Counts the non-boundary nodes in the passed lb fluid instance."""

    fluid_nodes = 0
    for n in lbf.nodes():
        if not n.boundary:
            fluid_nodes += 1

    return fluid_nodes


def random_dipoles(n_particles):
    """Generate random dipoles by sampling Euler angles uniformly at random."""
    cos_theta = 2 * np.random.random(n_particles) - 1
    sin_theta = np.sin(np.arcsin(cos_theta))
    phi = 2 * np.pi * np.random.random(n_particles)
    dip = np.array([sin_theta * np.cos(phi),
                    sin_theta * np.sin(phi),
                    cos_theta]).T
    return dip


def check_non_bonded_loop_trace(system):
    """Validates that the distances used by the non-bonded loop 
    match with the minimum image distance accessible by Python,
    checks that no pairs are lost or double-counted.
    """

    cs_pairs = system.cell_system.non_bonded_loop_trace()
    # format [id1, id2, pos1, pos2, vec2, mpi_node]

    distance_vec = system.distance_vec
    cutoff = system.cell_system.max_cut_nonbonded

    # Distance for all pairs of particles obtained by Python
    py_distances = {}
    for p1, p2 in system.part.pairs():
        py_distances[p1.id, p2.id] = np.copy(distance_vec(p1, p2))

    # Go through pairs found by the non-bonded loop and check distance
    for p in cs_pairs:
        # p is a tuple with (id1,id2,pos1,pos2,vec21)
        # Note that system.distance_vec uses the opposite sign convention
        # as the minimum image distance in the core

        if (p[0], p[1]) in py_distances:
            np.testing.assert_allclose(
                np.copy(p[4]), -py_distances[p[0], p[1]])
            del py_distances[p[0], p[1]]
        elif (p[1], p[0]) in py_distances:
            np.testing.assert_allclose(
                np.copy(p[4]), py_distances[p[1], p[0]])
            del py_distances[p[1], p[0]]
        else:
            raise Exception("Extra pair from core", p)

    for ids, dist in py_distances.items():
        if np.linalg.norm(dist) < cutoff:
            raise Exception("Pair not found by the core", ids)
