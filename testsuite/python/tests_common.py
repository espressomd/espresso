# Copyright (C) 2010-2018 The ESPResSo project
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
from __future__ import print_function
import os
import numpy as np
try:
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy
except BaseException:
    pass


def calculate_vtk_max_pointwise_difference(file1, file2, tol=1e-6):
    arrays = [0] * 2

    reader = vtk.vtkStructuredPointsReader()
    for i, fname in enumerate([file1, file2]):
        reader.SetFileName(fname)
        reader.Update()
        data = reader.GetOutput().GetPointData()
        arrays[i] = np.array([vtk_to_numpy(data.GetArray(n))
                              for n in range(data.GetNumberOfArrays())])

    try:
        return np.allclose(
            arrays[0], arrays[1], rtol=0, atol=tol), np.max(
            np.abs(
                arrays[0] - arrays[1]))
    except BaseException:
        return False, np.inf


try:
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy

    def calculate_vtk_max_pointwise_difference(file1, file2, tol=1e-6):
        arrays = [0] * 2

        reader = vtk.vtkStructuredPointsReader()
        for i, fname in enumerate([file1, file2]):
            reader.SetFileName(fname)
            reader.Update()
            data = reader.GetOutput().GetPointData()
            arrays[i] = np.array([vtk_to_numpy(data.GetArray(n))
                                  for n in range(data.GetNumberOfArrays())])

        try:
            return np.allclose(
                arrays[0], arrays[1], rtol=0, atol=tol), np.max(
                np.abs(
                    arrays[0] - arrays[1]))
        except BaseException:
            return False, np.inf
except BaseException:
    pass


def params_match(inParams, outParams):
    """Check, if the parameters set and gotten back match.
    Only check keys present in inParams.
    """

    for k in list(inParams.keys()):
        if k not in outParams:
            print(k, "missing from returned parameters")
            return False
        if isinstance(inParams[k], float):
            if abs(outParams[k] - inParams[k]) >= 1E-14:
                print("Mismatch in parameter ", k, inParams[k], outParams[k], type(
                    inParams[k]), type(outParams[k]), abs(inParams[k] - outParams[k]))
                return False
        else:
            if outParams[k] != inParams[k]:
                print("Mismatch in parameter ", k, inParams[
                      k], outParams[k], type(inParams[k]), type(outParams[k]))
                return False

    return True


def generate_test_for_class(_system, _interClass, _params):
    """Generates test cases for checking interaction parameters set and gotten back
    from Es actually match. Only keys which are present  in _params are checked
    1st: Interaction parameters as dictionary, i.e., {"k"=1.,"r_0"=0.
    2nd: Name of the interaction property to set (i.e. "P3M")
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

        self.assertTrue(
            params_match(
                params,
                outParams),
            "Missmatch of parameters.\nParameters set " +
            params.__str__() +
            " vs. output parameters " +
            outParams.__str__())

    return func


def lj_force_vector(v_d, d, lj_params):
    """Returns lj force for distance d and distance vecotr v_d based on the given lj_params.
    Supports epsilon and cutoff."""

    if d >= lj_params["cutoff"]:
        return np.zeros(3)

    return 4. * lj_params["epsilon"] * v_d * (-12.0 * d**-14 + 6.0 * d**-8)


def verify_lj_forces(system, tolerance, ids_to_skip=[]):
    """Goes over all pairs of particles in system and compares the forces on them
       to what would be expected based on the systems LJ parametes.
       Particle ids listed in ids_to_skip are not checked
       Do not run this with a thermostat enabled."""

    # Initialize dict with expected forces
    f_expected = {}
    for id in system.part[:].id:
        f_expected[id] = np.zeros(3)

    # Cache some stuff to speed up pair loop
    dist_vec = system.distance_vec
    norm = np.linalg.norm
    non_bonded_inter = system.non_bonded_inter
    # lj parameters
    lj_params = {}
    all_types = np.unique(system.part[:].type)
    for i in all_types:
        for j in all_types:
            lj_params[i, j] = non_bonded_inter[
                int(i), int(j)].lennard_jones.get_params()

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
    # Check actual forces agaisnt expected
    for id in system.part[:].id:
        if id in ids_to_skip:
            continue
        if np.linalg.norm(system.part[id].f - f_expected[id]) >= tolerance:
            raise Exception("LJ force verification failed on particle " +
                            str(id) +
                            ". Got " +
                            str(system.part[id].f) +
                            ", expected " +
                            str(f_expected[id]))


def abspath(path):
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), path)


def transform_pos_from_cartesian_to_polar_coordinates(pos):
    """Transform the given cartesian coordinates to polar coordinates.

    Parameters
    ----------
    pos : array_like :obj:`float`
          ``x``, ``y``, and ``z``-component of the cartesian position.

    Returns
    -------
    array_like
        The given position in polar coordinates.

    """
    return np.array([np.sqrt(pos[0]**2.0 + pos[1]**2.0), np.arctan2(pos[1], pos[0]), pos[2]])


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
        (pos[0] * vel[0] + pos[1] *
         vel[1]) / np.sqrt(pos[0]**2.0 + pos[1]**2.0),
        (pos[0] * vel[1] - pos[1] *
         vel[0]) / (pos[0]**2.0 + pos[1]**2.0),
        vel[2]])


def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.

    Parameters
    ----------
    axis : array_like :obj:`float`
           Axis to rotate around.
    theta : :obj:`float`
            Rotation angle.

    """
    axis = np.asarray(axis)
    axis = axis / np.sqrt(np.dot(axis, axis))
    a = np.cos(theta / 2.0)
    b, c, d = -axis * np.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])


def rotation_matrix_quat(system, part):
    """
    Return the rotation matrix associated with quaternion.

    Parameters
    ----------
    part : :obj:`int`
            Particle index.

    """
    A = np.zeros((3, 3))
    quat = system.part[part].quat
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


def get_cylindrical_bin_volume(
        n_r_bins,
        n_phi_bins,
        n_z_bins,
        min_r,
        max_r,
        min_phi,
        max_phi,
        min_z,
        max_z):
    """
    Return the bin volumes for a cylindrical histogram.

    Parameters
    ----------
    n_r_bins : :obj:`float`
               Number of bins in ``r`` direction.
    n_phi_bins : :obj:`float`
               Number of bins in ``phi`` direction.
    n_z_bins : :obj:`float`
               Number of bins in ``z`` direction.
    min_r : :obj:`float`
            Minimum considered value in ``r`` direction.
    max_r : :obj:`float`
            Maximum considered value in ``r`` direction.
    min_phi : :obj:`float`
              Minimum considered value in ``phi`` direction.
    max_phi : :obj:`float`
              Maximum considered value in ``phi`` direction.
    min_z : :obj:`float`
            Minimum considered value in ``z`` direction.
    max_z : :obj:`float`
            Maximum considered value in ``z`` direction.

    Returns
    -------
    array_like : Bin volumes.

    """
    bin_volume = np.zeros(n_r_bins)
    r_bin_size = (max_r - min_r) / n_r_bins
    phi_bin_size = (max_phi - min_phi) / n_phi_bins
    z_bin_size = (max_z - min_z) / n_z_bins
    for i in range(n_r_bins):
        bin_volume[i] = np.pi * ((min_r + r_bin_size * (i + 1))**2.0 -
                                 (min_r + r_bin_size * i)**2.0) * \
            phi_bin_size / (2.0 * np.pi) * z_bin_size
    return bin_volume

#
# Analytical Expressions for interactions
#

# Harmonic bond


def harmonic_potential(scalar_r, k, r_0, r_cut):
    return 0.5 * k * (scalar_r - r_0)**2


def harmonic_force(scalar_r, k, r_0, r_cut):
    return -k * (scalar_r - r_0)

# FENE bond


def fene_potential(scalar_r, k, d_r_max, r_0):
    return -0.5 * k * d_r_max**2 * np.log(1 - ((scalar_r - r_0) / d_r_max)**2)


def fene_force(scalar_r, k, d_r_max, r_0):
    return k * (scalar_r - r_0) * d_r_max**2 / ((scalar_r - r_0)**2 - d_r_max**2)


def fene_force2(bond_vector, k, d_r_max, r_0):
    r = np.linalg.norm(bond_vector)
    return k * (r - r_0) / (r * (1 - ((r - r_0) / d_r_max)**2)) * np.array(bond_vector)

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


def lj_generic_potential(r, eps, sig, cutoff, offset=0., shift=0., e1=12.,
                         e2=6., b1=4., b2=4., delta=0., lam=1.):
    r = np.array(r)
    V = np.zeros_like(r)
    cutoffMask = (r <= cutoff + offset)
    # LJGEN_SOFTCORE transformations
    rroff = np.sqrt(
        np.power(r[cutoffMask] - offset, 2) + (1 - lam) * delta * sig**2)
    V[cutoffMask] = eps * lam * \
        (b1 * np.power(sig / rroff, e1) -
         b2 * np.power(sig / rroff, e2) + shift)
    return V


def lj_generic_force(espressomd, r, eps, sig, cutoff, offset=0., e1=12, e2=6,
                     b1=4., b2=4., delta=0., lam=1., generic=True):
    f = 1.
    if (r >= offset + cutoff):
        f = 0.
    else:
        h = (r - offset)**2 + delta * (1. - lam) * sig**2
        f = (r - offset) * eps * lam * (
            b1 * e1 * np.power(sig / np.sqrt(h), e1) - b2 * e2 * np.power(sig / np.sqrt(h), e2)) / h
        if (not espressomd.has_features("LJGEN_SOFTCORE")) and generic:
            f *= np.sign(r - offset)
    return f

# Lennard-Jones


def lj_potential(r, eps, sig, cutoff, shift, offset=0.):
    V = lj_generic_potential(
        r, eps, sig, cutoff, offset=offset, shift=shift * 4.)
    return V


def lj_force(espressomd, r, eps, sig, cutoff, offset=0.):
    f = lj_generic_force(
        espressomd, r, eps, sig, cutoff, offset=offset, generic=False)
    return f

# Lennard-Jones Cosine


def lj_cos_potential(r, eps, sig, cutoff, offset):
    V = 0.
    r_min = offset + np.power(2., 1. / 6.) * sig
    r_cut = cutoff + offset
    if (r < r_min):
        V = lj_potential(r, eps=eps, sig=sig,
                         cutoff=cutoff, offset=offset, shift=0.)
    elif (r < r_cut):
        alpha = np.pi / \
            (np.power(r_cut - offset, 2) - np.power(r_min - offset, 2))
        beta = np.pi - np.power(r_min - offset, 2) * alpha
        V = 0.5 * eps * \
            (np.cos(alpha * np.power(r - offset, 2) + beta) - 1.)
    return V


def lj_cos_force(espressomd, r, eps, sig, cutoff, offset):
    f = 0.
    r_min = offset + np.power(2., 1. / 6.) * sig
    r_cut = cutoff + offset
    if (r < r_min):
        f = lj_force(espressomd, r, eps=eps, sig=sig,
                     cutoff=cutoff, offset=offset)
    elif (r < r_cut):
        alpha = np.pi / \
            (np.power(r_cut - offset, 2) - np.power(r_min - offset, 2))
        beta = np.pi - np.power(r_min - offset, 2) * alpha
        f = (r - offset) * alpha * eps * \
            np.sin(alpha * np.power(r - offset, 2) + beta)
    return f

# Lennard-Jones Cosine^2


def lj_cos2_potential(r, eps, sig, offset, width):
    V = 0.
    r_min = offset + np.power(2., 1. / 6.) * sig
    r_cut = r_min + width
    if (r < r_min):
        V = lj_potential(r, eps=eps, sig=sig,
                         offset=offset, cutoff=r_cut, shift=0.)
    elif (r < r_cut):
        V = -eps * np.power(np.cos(np.pi /
                                   (2. * width) * (r - r_min)), 2)
    return V


def lj_cos2_force(espressomd, r, eps, sig, offset, width):
    f = 0.
    r_min = offset + np.power(2., 1. / 6.) * sig
    r_cut = r_min + width
    if (r < r_min):
        f = lj_force(espressomd, r, eps=eps,
                     sig=sig, cutoff=r_cut, offset=offset)
    elif (r < r_cut):
        f = - np.pi * eps * \
            np.sin(np.pi * (r - r_min) / width) / (2. * width)
    return f

# Smooth-Step


def smooth_step_potential(r, eps, sig, cutoff, d, n, k0):
    V = 0.
    if (r < cutoff):
        V = np.power(d / r, n) + eps / \
            (1 + np.exp(2 * k0 * (r - sig)))
    return V


def smooth_step_force(r, eps, sig, cutoff, d, n, k0):
    f = 0.
    if (r < cutoff):
        f = n * d / r**2 * np.power(d / r, n - 1) + 2 * k0 * eps * np.exp(
            2 * k0 * (r - sig)) / (1 + np.exp(2 * k0 * (r - sig))**2)
    return f

# BMHTF


def bmhtf_potential(r, a, b, c, d, sig, cutoff):
    V = 0.
    if (r == cutoff):
        V = a * np.exp(b * (sig - r)) - c * np.power(
            r, -6) - d * np.power(r, -8)
    if (r < cutoff):
        V = a * np.exp(b * (sig - r)) - c * np.power(
            r, -6) - d * np.power(r, -8)
        V -= bmhtf_potential(cutoff, a, b, c, d, sig, cutoff)
    return V


def bmhtf_force(r, a, b, c, d, sig, cutoff):
    f = 0.
    if (r < cutoff):
        f = a * b * np.exp(b * (sig - r)) - 6 * c * np.power(
            r, -7) - 8 * d * np.power(r, -9)
    return f

# Morse


def morse_potential(r, eps, alpha, cutoff, rmin=0):
    V = 0.
    if (r < cutoff):
        V = eps * (np.exp(-2. * alpha * (r - rmin)) -
                   2 * np.exp(-alpha * (r - rmin)))
        V -= eps * (np.exp(-2. * alpha * (cutoff - rmin)
                           ) - 2 * np.exp(-alpha * (cutoff - rmin)))
    return V


def morse_force(r, eps, alpha, cutoff, rmin=0):
    f = 0.
    if (r < cutoff):
        f = 2. * np.exp((rmin - r) * alpha) * \
            (np.exp((rmin - r) * alpha) - 1) * alpha * eps
    return f

#  Buckingham


def buckingham_potential(r, a, b, c, d, cutoff, discont, shift):
    V = 0.
    if (r < discont):
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
    if (r < discont):
        f = buckingham_force(
            discont, a, b, c, d, cutoff, discont, shift)
    if (r >= discont) and (r < cutoff):
        f = a * b * np.exp(- b * r) - 6 * c * np.power(
            r, -7) - 4 * d * np.power(r, -5)
    return f

# Soft-sphere


def soft_sphere_potential(r, a, n, cutoff, offset=0):
    V = 0.
    if (r < offset + cutoff):
        V = a * np.power(r - offset, -n)
    return V


def soft_sphere_force(r, a, n, cutoff, offset=0):
    f = 0.
    if ((r > offset) and (r < offset + cutoff)):
        f = n * a * np.power(r - offset, -(n + 1))
    return f

# Hertzian


def hertzian_potential(r, eps, sig):
    V = 0.
    if (r < sig):
        V = eps * np.power(1 - r / sig, 5. / 2.)
    return V


def hertzian_force(r, eps, sig):
    f = 0.
    if (r < sig):
        f = 5. / 2. * eps / sig * np.power(1 - r / sig, 3. / 2.)
    return f

# Gaussian


def gaussian_potential(r, eps, sig, cutoff):
    V = 0.
    if (r < cutoff):
        V = eps * np.exp(-np.power(r / sig, 2) / 2)
    return V


def gaussian_force(r, eps, sig, cutoff):
    f = 0.
    if (r < cutoff):
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


class DynamicDict(dict):

    def __getitem__(self, key):
        value = super(DynamicDict, self).__getitem__(key)
        return eval(value, self) if isinstance(value, str) else value


def single_component_maxwell(x1, x2, kT):
    """Integrate the probability density from x1 to x2 using the trapezoidal rule"""
    x = np.linspace(x1, x2, 1000)
    return np.trapz(np.exp(-x**2 / (2. * kT)), x) / np.sqrt(2. * np.pi * kT)


def lists_contain_same_elements(list1, list2):
    return len(list1) == len(list2) and sorted(list1) == sorted(list2)
