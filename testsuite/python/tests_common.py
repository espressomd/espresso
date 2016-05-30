import numpy as np
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
