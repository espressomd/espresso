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
        if type(inParams[k]) == float:
          if abs(outParams[k] -inParams[k])>=1E-14:
              print("Mismatch in parameter ", k, inParams[k], outParams[k],type(inParams[k]),type(outParams[k]),abs(inParams[k]-outParams[k]))
              return False
        else:
          if outParams[k] !=inParams[k]:
              print("Mismatch in parameter ", k, inParams[k], outParams[k],type(inParams[k]),type(outParams[k]))
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


def lj_force(v_d, d, lj_params):
    """Returns lj force for distnace d and distance vecotr v_d based on the given lj_params.
    Supports epsilon and cutoff."""

    if d >= lj_params["cutoff"]:
        return np.zeros(3)

    return 4. * lj_params["epsilon"] * v_d * (-12.0 * d**-14 + 6.0 * d**-8)


def verify_lj_forces(system, tolerance, ids_to_skip=[]):
    """Goes over all pairs of paritcles in system and compares the forces on them
       to what would be expected based on the systems lj parametes.
       Particle ids listed in ids_to_skip are not checked
       Do not run this with a thermostat enabled."""

    # Initialize dict with expected forces
    f_expected = {}
    for id in system.part[:].id:
        f_expected[id] = np.zeros(3)

    # Cache some stuff to speed up pair loop
    dist_vec=system.distance_vec
    norm=np.linalg.norm
    non_bonded_inter=system.non_bonded_inter
    # lj parameters
    lj_params={}
    all_types=np.unique(system.part[:].type)
    for i in all_types:
        for j in all_types:
            lj_params[i,j]=non_bonded_inter[int(i),int(j)].lennard_jones.get_params()
            


      

    # Go over all pairs of particles
    for pair in system.part.pairs():
        p0=pair[0]
        p1=pair[1]
        if p0.id in ids_to_skip or p1.id in ids_to_skip:
            continue

        # Distance and distance vec
        v_d = dist_vec(p0,p1)
        d = norm(v_d)

        # calc and add expected lj force
        f = lj_force(v_d, d, lj_params[p0.type,p1.type])
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
         vel[1]) / np.sqrt(pos[0]**2.0 + pos[1]**2.0),\
        (pos[0] * vel[1] - pos[1] *
         vel[0]) / (pos[0]**2.0 + pos[1]**2.0),\
        vel[2]])

