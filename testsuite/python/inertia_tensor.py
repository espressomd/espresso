import numpy as np
import unittest as ut
try:
    import scipy.spatial.transform as sst
except BaseException:
    pass

import espressomd.inertia_tensor
import unittest_decorators as utx


def generate_cuboid_positions(rho, dx, dy, dz):
    """
    Generate a list of three dimensional mesh positions.

    Parameters
    ----------
    rho : :obj:`float`
        Samples per unit length.
    dx : :obj:`float`
        Range in dimension 0.
    dy : :obj:`float`
        Range in dimension 1.
    dz : :obj:`float`
        Range in dimension 2.

    Returns
    -------
    array_like of :obj:`float`
        Three dimensional mesh positions;

    """
    xs = np.linspace(-0.5 * dx, 0.5 * dx, int(rho * dx))
    ys = np.linspace(-0.5 * dy, 0.5 * dy, int(rho * dy))
    zs = np.linspace(-0.5 * dz, 0.5 * dz, int(rho * dz))
    return np.vstack(np.meshgrid(xs, ys, zs)).reshape(3, -1).T


def inertia_tensor_cuboid(mass, dx, dy, dz):
    """
    Reference values for the inertia tensor of a cuboid.

    Parameters
    ----------
    mass : :obj:`float`
        Mass of the cuboid.
    dx : :obj:`float`
        Extension in dimension 0.
    dy : :obj:`float`
        Extension in dimension 1.
    dz : :obj:`float`
        Extension in dimension 2.

    Returns
    -------
    array_like of :obj:`float`
        Inertia tensor of the cuboid.

    Notes
    -----
    See wikipedia_.

    .. _wikipedia: https://en.wikipedia.org/wiki/List_of_moments_of_inertia#List_of_3D_inertia_tensors

    """
    return 1. / 12. * mass * \
        np.diag([dy**2.0 + dz**2.0, dx**2.0 + dz**2.0, dx**2.0 + dy**2.0])


class TestInertiaTensor(ut.TestCase):
    """
    Tests for the inertia tensor utility functions.

    """
    @classmethod
    def setUpClass(cls):
        cls.dx = 1.32
        cls.dy = 2.12
        cls.dz = 3.23
        rho = 5
        cls.samples = generate_cuboid_positions(rho, cls.dx, cls.dy, cls.dz)
        cls.N_samples = cls.samples.shape[0]
        cls.m = 0.5
        cls.masses = np.ones(cls.N_samples) * cls.m / cls.N_samples

    def test_inertia_tensor(self):
        """
        Compare the calculated inertia tensor of a sampled cuboid with the
        respective literature values.

        """
        np.testing.assert_almost_equal(espressomd.inertia_tensor.inertia_tensor(
            self.samples, self.masses), inertia_tensor_cuboid(self.m, self.dx, self.dy, self.dz), decimal=1)

    def test_right_handedness_eigenvectormatrix(self):
        """
        Check that the eigenvectors form a right-handed basis.

        """
        _, eigenvectors = espressomd.inertia_tensor.diagonalized_inertia_tensor(
            self.samples, self.masses)
        for i in range(3):
            ev = np.roll(eigenvectors, axis=0, shift=i)
            np.testing.assert_allclose(
                np.cross(ev[0], ev[1]), ev[2], atol=1e-7)

    @utx.skipIfUnmetModuleVersionRequirement('scipy', '>1.2.0')
    def test_inertia_tensor_rotated_cuboid(self):
        """
        Rotate the samples and check that the principal axes return by the
        utility function corresponds to the rotation matrix.

        """
        angle = 2.0 * np.pi * np.random.random()
        quat = [np.sin(angle / 2.0), np.sin(angle / 2.0),
                np.sin(angle / 2.0), np.cos(angle / 2.0)]
        rotation = sst.Rotation.from_quat(quat)
        rotated_samples = rotation.apply(self.samples)
        _, eigenvectors = espressomd.inertia_tensor.diagonalized_inertia_tensor(
            rotated_samples, self.masses)
        rotated_basis = rotation.apply(np.identity(3))
        for i in range(3):
            # because there is no particular order in the eigenvalues
            # the corresponding eigenvectors are either (anti-) parallel or
            # perpendicular to the rotated basis
            self.assertAlmostEqual(
                abs(abs(np.dot(rotated_basis[i], eigenvectors[i])) - 0.5) - 0.5, 0.0)


if __name__ == "__main__":
    ut.main()
