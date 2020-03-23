import math

import numpy as np


def matrix_to_quat(m):
    """
    Convert the (proper) rotation matrix to the corresponding quaternion representation based
    on pyquaternion_.

    .. _pyquaternion: https://github.com/KieranWynn/pyquaternion

    Parameters
    ----------
    m : (3,3) array_like of :obj:`float`
        Rotation matrix.

    Returns
    -------
    array_like of :obj:`float`
        Quaternion representation of the rotation.

    Raises
    ------
    ValueError
        If the matrix does not a correspond to a proper rotation (det(m) != 1).

    """
    if not math.isclose(np.linalg.det(m), 1.0, abs_tol=1e-7):
        raise ValueError("Only proper rotations are supported")
    m = m.copy().conj().transpose()
    if m[2, 2] < 0:
        if m[0, 0] > m[1, 1]:
            t = 1 + m[0, 0] - m[1, 1] - m[2, 2]
            return 0.5 / \
                math.sqrt(t) * np.array([m[1, 2] - m[2, 1],
                                         t, m[0, 1] + m[1, 0], m[2, 0] + m[0, 2]])
        else:
            t = 1 - m[0, 0] + m[1, 1] - m[2, 2]
            return 0.5 / \
                math.sqrt(t) * np.array([m[2, 0] - m[0, 2],
                                         m[0, 1] + m[1, 0], t, m[1, 2] + m[2, 1]])
    else:
        if m[0, 0] < -m[1, 1]:
            t = 1 - m[0, 0] - m[1, 1] + m[2, 2]
            return 0.5 / \
                math.sqrt(t) * np.array([m[0, 1] - m[1, 0],
                                         m[2, 0] + m[0, 2], m[1, 2] + m[2, 1], t])
        else:
            t = 1 + m[0, 0] + m[1, 1] + m[2, 2]
            return 0.5 / \
                math.sqrt(
                    t) * np.array([t, m[1, 2] - m[2, 1], m[2, 0] - m[0, 2], m[0, 1] - m[1, 0]])


def inertia_tensor(positions, masses):
    """
    Calculate the inertia tensor for given point masses
    at given positions.

    Parameters
    ----------
    positions : (N,3) array_like of :obj:`float`
        Point masses' positions.
    masses : (N,) array_like of :obj:`float`

    Returns
    -------
    (3,3) array_like of :obj:`float`
        Moment of inertia tensor.

    Notes
    -----
    See wikipedia_.

    .. _wikipedia: https://en.wikipedia.org/wiki/Moment_of_inertia#Inertia_tensor

    """
    inertia_tensor = np.zeros((3, 3))
    for ind, p in enumerate(positions):
        inertia_tensor += masses[ind] * \
            (np.dot(p, p) * np.identity(3)
             - np.outer(p, p))
    return inertia_tensor


def diagonalized_inertia_tensor(positions, masses):
    """
    Calculate the diagonalized inertia tensor
    with respect to the center of mass
    for given point masses at given positions.

    Parameters
    ----------
    positions : (N,3) array_like of :obj:`float`
        Positions of the masses.
    masses : (N,) array_like of :obj:`float`

    Returns
    -------
    (3,) array_like of :obj:`float`
        Principal moments of inertia.
    (3,3) array_like of :obj:`float`
        Principal axes of inertia.
        Note that the second axis is the coordinate axis (same as input).

    """
    assert np.array(masses).shape[0] == np.array(positions).shape[0]

    def center_of_mass(positions, masses): return np.average(
        positions, axis=0, weights=masses)
    positions = np.array(np.copy(positions)) - \
        center_of_mass(positions, masses)
    inertia = inertia_tensor(positions, masses)
    eig, eigv = np.linalg.eig(inertia)
    eigv = np.transpose(eigv)
    # check for right-handedness
    if not np.allclose(np.cross(eigv[0], eigv[1]), eigv[2]):
        eigv[[0, 1], :] = eigv[[1, 0], :]
    return eig, eigv
