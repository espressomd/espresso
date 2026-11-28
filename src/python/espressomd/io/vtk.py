#
# Copyright (C) 2023-2026 The ESPResSo project
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
#

import numpy as np
import vtk
import vtk.util.numpy_support

IMAGE_DATA_TYPE = "ImageData"
UNSTRUCTURED_GRID_TYPE = "UnstructuredGrid"


class VTKReader:
    """
    Reader for VTK multi-piece uniform grids written in XML format.
    """
    error_tolerance = 1e-5  # VTK data is written with 1e-7 precision
    type_map = {
        IMAGE_DATA_TYPE: vtk.vtkXMLImageDataReader(),
        UNSTRUCTURED_GRID_TYPE: vtk.vtkXMLUnstructuredGridReader(),
    }

    @classmethod
    def get_file_format(cls, filepath):
        with open(str(filepath), "rb") as f:
            line = f.readline().lstrip()
            while not line.startswith(b"<VTKFile"):
                line = f.readline()
                if not line:
                    break
                line = line.lstrip()
            if line.startswith(b"<VTKFile"):
                for token in line.split():
                    if token.startswith(b"type="):
                        return token.split(b'"')[1].decode()
            raise RuntimeError(f"File {filepath} is not a compliant XML file")

    def parse(self, filepath):
        type_name = self.get_file_format(filepath)
        if type_name not in self.type_map.keys():
            raise NotImplementedError(
                f"Unknown VTK file format '{type_name}' (supported "
                f"file formats: {', '.join(self.type_map.keys())})")
        reader = self.type_map[type_name]
        reader.SetFileName(str(filepath))
        reader.Update()
        grid = reader.GetOutput()
        result = {}

        for i in range(grid.GetCellData().GetNumberOfArrays()):
            arr = grid.GetCellData().GetArray(i)
            data = vtk.util.numpy_support.vtk_to_numpy(arr)
            n_comp = arr.GetNumberOfComponents()

            if type_name == UNSTRUCTURED_GRID_TYPE:
                positions = {}
                for j in range(grid.GetNumberOfCells()):
                    cell = grid.GetCell(j)
                    bounds = cell.GetBounds()
                    key = (round(bounds[0]), round(
                        bounds[2]), round(bounds[4]))
                    positions[key] = j

                if not positions:
                    continue

                xs = [p[0] for p in positions]
                ys = [p[1] for p in positions]
                zs = [p[2] for p in positions]
                x0, x1 = min(xs), max(xs)
                y0, y1 = min(ys), max(ys)
                z0, z1 = min(zs), max(zs)
                nx, ny, nz = x1 - x0 + 1, y1 - y0 + 1, z1 - z0 + 1

                grid_shape = (nx, ny, nz)
                if n_comp > 1:
                    grid_shape = grid_shape + (n_comp,)
                    data = data.reshape(data.shape[0], n_comp)

                arr_out = np.full(grid_shape, np.nan, dtype=data.dtype)
                for (x, y, z), cid in positions.items():
                    arr_out[x - x0, y - y0, z - z0] = data[cid]
            elif type_name == IMAGE_DATA_TYPE:
                dims = grid.GetDimensions()
                shape = (dims[0] - 1, dims[1] - 1, dims[2] - 1)
                if n_comp > 1:
                    shape = shape + (n_comp,)
                arr_out = data.reshape(shape, order="F")

            result[arr.GetName()] = arr_out
        return result
