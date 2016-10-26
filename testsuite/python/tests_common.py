from __future__ import print_function
import numpy as np
try:
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy
except:
    pass

def calculate_vtk_max_pointwise_difference(file1,file2,tol=1e-6):
    arrays = [0]*2

    reader = vtk.vtkStructuredPointsReader()
    for i, fname in enumerate([file1, file2]):
        reader.SetFileName(fname)
        reader.Update()
        data = reader.GetOutput().GetPointData()
        arrays[i] = np.array([vtk_to_numpy(data.GetArray(n)) for n in range(data.GetNumberOfArrays())])

    try:
        return np.allclose(arrays[0],arrays[1],rtol=0,atol=tol), np.max(np.abs(arrays[0]-arrays[1]))
    except:
        return False, np.inf


try:
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy

    def calculate_vtk_max_pointwise_difference(file1,file2,tol=1e-6):
        arrays = [0]*2

        reader = vtk.vtkStructuredPointsReader()
        for i, fname in enumerate([file1, file2]):
            reader.SetFileName(fname)
            reader.Update()
            data = reader.GetOutput().GetPointData()
            arrays[i] = np.array([vtk_to_numpy(data.GetArray(n)) for n in range(data.GetNumberOfArrays())])

        try:
            return np.allclose(arrays[0],arrays[1],rtol=0,atol=tol), np.max(np.abs(arrays[0]-arrays[1]))
        except:
            return False, np.inf
except:
    pass

# Tests particle property setters/getters
import unittest as ut
import espressomd
import numpy as np


def params_match(inParams, outParams):
    """Check, if the parameters set and gotten back match.
    Only check keys present in inParams.
    """

    for k in list(inParams.keys()):
        if k not in outParams:
            print(k, "missing from returned parameters")
            return False
        if outParams[k] != inParams[k]:
            print("Mismatch in parameter ", k, inParams[k], outParams[k])
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

        self.assertTrue(params_match(params, outParams), "Missmatch of parameters.\nParameters set " +
                        params.__str__() + " vs. output parameters " + outParams.__str__())

    return func

