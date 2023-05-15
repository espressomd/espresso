#
# Copyright (C) 2023 The ESPResSo project
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


class VTKReader:
    """
    Reader for VTK multi-piece uniform grids written in XML format.
    """
    error_tolerance = 1e-5  # VTK data is written with 1e-7 precision

    @classmethod
    def get_array_names(cls, reader):
        array_names = set()
        n_ghost_layers = reader.GetUpdateGhostLevel()
        n_pieces = reader.GetNumberOfPieces()
        for piece_index in range(n_pieces):
            reader.UpdatePiece(piece_index, n_pieces, n_ghost_layers)
            piece = reader.GetOutput()
            cell = piece.GetCellData()
            for i in range(cell.GetNumberOfArrays()):
                array_names.add(cell.GetArrayName(i))
        return array_names

    @classmethod
    def get_piece_topology(
            cls, piece, array, bounding_box_lower, bounding_box_upper):
        bounds = np.array(piece.GetBounds())
        box_l = bounds[1::2] - bounds[0:-1:2]
        n_grid_points = array.GetNumberOfTuples()
        shape_float = box_l / np.min(box_l)
        shape_float *= np.cbrt(n_grid_points / np.prod(shape_float))
        shape_int = np.around(shape_float).astype(int)
        assert np.linalg.norm(shape_int - shape_float) < cls.error_tolerance and np.prod(
            shape_int) == n_grid_points, "only cubic grids are supported"
        agrid = np.mean(box_l / shape_float)
        shape = tuple(shape_int.tolist())
        lower_corner = []
        for i in range(3):
            start = int(np.around(bounds[i * 2]))
            stop = start + shape[i]
            bounding_box_lower[i] = min(bounding_box_lower[i], start)
            bounding_box_upper[i] = max(bounding_box_upper[i], stop)
            lower_corner.append(start)
        return agrid, shape, lower_corner

    @classmethod
    def reconstruct_array(cls, reader, array_name):
        n_pieces = reader.GetNumberOfPieces()
        n_ghost_layers = reader.GetUpdateGhostLevel()
        # get bounding box
        info = []
        agrids = []
        bounding_box_lower = 3 * [float("inf")]
        bounding_box_upper = 3 * [-float("inf")]
        for piece_index in range(n_pieces):
            reader.UpdatePiece(piece_index, n_pieces, n_ghost_layers)
            piece = reader.GetOutput()
            cell = piece.GetCellData()
            array = cell.GetArray(array_name)
            if array is not None:
                agrid, shape, lower_corner = cls.get_piece_topology(
                    piece, array, bounding_box_lower, bounding_box_upper)
                agrids.append(agrid)
                info.append([piece_index, shape, lower_corner])

        if not info:
            return None

        # get array type and size
        assert float("inf") not in bounding_box_lower
        assert -float("inf") not in bounding_box_upper
        if np.std(agrids) / np.mean(agrids) > cls.error_tolerance:
            raise NotImplementedError(
                f"VTK non-uniform grids are not supported (got agrid = {agrids} when parsing array '{array_name}')")
        data_dims = np.array(bounding_box_upper) - np.array(bounding_box_lower)
        piece_index = info[0][0]
        reader.UpdatePiece(piece_index, n_pieces, n_ghost_layers)
        array = reader.GetOutput().GetCellData().GetArray(array_name)
        vector_length = array.GetNumberOfComponents()
        val_dims = [] if vector_length == 1 else [vector_length]
        data_type = array.GetDataTypeAsString()
        if data_type == "float":
            dtype = float
        elif data_type == "int":
            dtype = int
        else:
            raise NotImplementedError(
                f"Unknown VTK data type '{data_type}' (when parsing array '{array_name}')")

        # get data
        data = np.empty(data_dims.tolist() + val_dims, dtype=dtype)
        for piece_index, shape, lower_corner in info:
            reader.UpdatePiece(piece_index, n_pieces, n_ghost_layers)
            array = reader.GetOutput().GetCellData().GetArray(array_name)
            subset = []
            for i in range(3):
                start = lower_corner[i] - bounding_box_lower[i]
                stop = start + shape[i]
                subset.append(slice(start, stop))
            data[tuple(subset)] = vtk.util.numpy_support.vtk_to_numpy(
                array).reshape(list(shape) + val_dims, order='F')

        return data

    def parse(self, filepath):
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(str(filepath))
        reader.Update()

        arrays = {}
        array_names = self.get_array_names(reader)
        for array_name in sorted(array_names):
            arrays[array_name] = self.reconstruct_array(reader, array_name)

        return arrays
