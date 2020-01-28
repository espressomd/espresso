import numpy as np


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
        mass = masses[ind]
        inertia_tensor[0, 0] += mass * (p[1]**2.0 + p[2]**2.0)
        inertia_tensor[1, 1] += mass * (p[0]**2.0 + p[2]**2.0)
        inertia_tensor[2, 2] += mass * (p[0]**2.0 + p[1]**2.0)
        xy = -mass * p[0] * p[1]
        inertia_tensor[0, 1] += xy
        inertia_tensor[1, 0] += xy
        xz = -mass * p[0] * p[2]
        inertia_tensor[0, 2] += xz
        inertia_tensor[2, 0] += xz
        yz = -mass * p[1] * p[2]
        inertia_tensor[1, 2] += yz
        inertia_tensor[2, 1] += yz
    assert inertia_tensor[0, 1] == inertia_tensor[1, 0]
    assert inertia_tensor[0, 2] == inertia_tensor[2, 0]
    assert inertia_tensor[1, 2] == inertia_tensor[2, 1]
    return inertia_tensor


def center_of_mass(positions, masses):
    """
    Calculate the center of mass for given point masses at given positions.

    Parameters
    ----------
    positions : (N,3) array_like of :obj:`float`
        Point masses' positions.
    masses : (N,) array_like of :obj:`float`

    Returns
    -------
    :obj:`float`
        Center of mass.

    """
    assert np.array(masses).shape[0] == np.array(positions).shape[0]
    return np.average(positions, axis=0, weights=masses)


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
    positions = np.array(np.copy(positions)) - \
        center_of_mass(positions, masses)
    inertia = inertia_tensor(positions, masses)
    eig, eigv = np.linalg.eig(inertia)
    eigv = np.transpose(eigv)
    # check for right-handedness
    if not np.allclose(np.cross(eigv[0], eigv[1]), eigv[2]):
        eigv[[0, 1], :] = eigv[[1, 0], :]
    return eig, eigv
