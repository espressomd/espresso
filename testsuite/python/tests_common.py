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

def lj_force(v_d,d,lj_params):
    """Returns lj force for distnace d and distance vecotr v_d based on the given lj_params.
    Supports epsilon and cutoff."""

    if d>=lj_params["cutoff"]:
        return np.array((0.,0.,0.))
    
    return 4.*lj_params["epsilon"] * v_d/d *(-12*d**-13+6*d**-7)

     
def verify_lj_forces(system,tolerance,ids_to_skip=[]):
    """Goes over all pairs of paritcles in system and compares the forces on them
       to what would be expected based on the systems lj parametes.
       Particle ids listed in ids_to_skip are not checked
       Do not run this with a thermostat enabled."""
    
    # Initialize dict with expected forces
    f_expected={}
    for id in system.part[:].id:
        f_expected[id]=np.zeros(3)
    
    # Go over all pairs of particles
    for pair in system.part.pairs():
        if pair[0].id in ids_to_skip or pair[1].id in ids_to_skip: continue

        # Distance and distance vec
        v_d=system.distance_vec(pair[0],pair[1])
        d=system.distance(pair[0],pair[1])
        
        # get lj params for the type combination from system
        lj_params=system.non_bonded_inter[pair[0].type,pair[1].type].lennard_jones.get_params()

        # calc and add expected lj force
        f=lj_force(v_d,d,lj_params)
        f_expected[pair[0].id]+=f
        f_expected[pair[1].id]-=f
    # Check actual forces agaisnt expected
    for id in system.part[:].id:
       if id in ids_to_skip: continue
       if np.linalg.norm(system.part[id].f-f_expected[id]) >=tolerance:
           raise Exception("LJ force verification failed on particle "+str(id)+". Got "+str(system.part[id].f)+", expected "+str(f_expected[id]))


   
