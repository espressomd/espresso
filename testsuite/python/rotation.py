import numpy as np
import unittest as ut

try:
    import scipy.spatial.transform as sst
except BaseException:
    pass

import espressomd.rotation
import unittest_decorators as utx


@utx.skipIfUnmetModuleVersionRequirement('scipy', '>=1.4.0')
class TestRotation(ut.TestCase):
    """
    Tests for the rotation utility functions.

    """

    def setUp(self):
        angle = 2.0 * np.pi * np.random.random()
        quat = [np.sin(angle / 2.0), np.sin(angle / 2.0),
                np.sin(angle / 2.0), np.cos(angle / 2.0)]
        self.rotation = sst.Rotation.from_quat(quat)

    def test_quat_from_matrix(self):
        """
        Compare the calculated quaternion representation with scipy.

        """
        v_x = np.array([1.0, 0.0, 0.0])
        rotated_vector_ref = self.rotation.apply(v_x)
        quat_from_matrix = espressomd.rotation.matrix_to_quat(
            self.rotation.as_matrix())
        rotated_vector_matrix = sst.Rotation.from_quat(
            np.roll(quat_from_matrix, shift=-1)).apply(v_x)
        self.assertAlmostEqual(
            np.dot(rotated_vector_ref, rotated_vector_matrix), 1.0)

    def test_raise_if_improper(self):
        """
        Check that an improper rotation matrix as an argument to
        :meth:`espressomd.rotation.matrix_to_quat` raises an exception.

        """
        matrix = self.rotation.as_matrix()
        matrix[[0, 1], :] = matrix[[1, 0], :]
        with self.assertRaises(ValueError):
            espressomd.rotation.matrix_to_quat(matrix)


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
        np.testing.assert_almost_equal(espressomd.rotation.inertia_tensor(
            self.samples, self.masses), inertia_tensor_cuboid(self.m, self.dx, self.dy, self.dz), decimal=1)

    def test_right_handedness_eigenvectormatrix(self):
        """
        Check that the eigenvectors form a right-handed basis.

        """
        _, eigenvectors = espressomd.rotation.diagonalized_inertia_tensor(
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
        _, eigenvectors = espressomd.rotation.diagonalized_inertia_tensor(
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
