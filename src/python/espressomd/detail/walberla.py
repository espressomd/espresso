#
# Copyright (C) 2020-2023 The ESPResSo project
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

import os
import numpy as np
import itertools

from espressomd import utils
import espressomd.code_features


def get_slice_bounding_box(slices, grid_size):
    shape = []
    slice_lower_corner = []
    slice_upper_corner = []
    for i in range(3):
        indices = np.arange(grid_size[i])
        if isinstance(slices[i], slice):
            if slices[i].step not in [None, 1]:
                raise NotImplementedError(
                    "Slices with step != 1 are not supported")
            indices = indices[slices[i]]
        else:
            if isinstance(slices[i], (int, np.integer)):
                indices = [indices[slices[i]]]
            else:
                raise NotImplementedError(
                    "Tuple-based indexing is not supported")
        if len(indices) == 0:
            slice_lower_corner.append(0)
            slice_upper_corner.append(0)
            shape.append(0)
        elif isinstance(slices[i], (int, np.integer)):
            slice_lower_corner.append(indices[0])
            slice_upper_corner.append(indices[0] + 1)
        else:
            slice_lower_corner.append(indices[0])
            slice_upper_corner.append(indices[-1] + 1)
            shape.append(len(indices))
    return {"slice_lower_corner": slice_lower_corner,
            "slice_upper_corner": slice_upper_corner,
            "shape": shape}


class VTKRegistry:

    def __init__(self):
        if not espressomd.code_features.has_features("WALBERLA"):
            raise NotImplementedError("Feature WALBERLA not compiled in")
        self.map = {}
        self.collisions = {}

    def _register_vtk_object(self, vtk_obj):
        vtk_uid = vtk_obj.vtk_uid
        self.map[vtk_uid] = vtk_obj
        self.collisions[os.path.abspath(vtk_uid)] = vtk_uid

    def __getstate__(self):
        return self.map

    def __setstate__(self, active_vtk_objects):
        self.map = {}
        self.collisions = {}
        for vtk_uid, vtk_obj in active_vtk_objects.items():
            self.map[vtk_uid] = vtk_obj
            self.collisions[os.path.abspath(vtk_uid)] = vtk_uid


class LatticeSliceWalberla:

    def __init__(self, parent_sip, node_grid, slice_range):
        self.indices = [np.atleast_1d(np.arange(node_grid[i])[slice_range[i]])
                        for i in range(3)]
        self.dimensions = [ind.size for ind in self.indices]
        self._sip = parent_sip
        self._node = self._node_class(parent_sip=self._sip, index=[0, 0, 0])

    def __iter__(self):
        for index in itertools.product(*self.indices):
            yield self._node_class(parent_sip=self._sip, index=np.array(index))

    def setter_checks(self, attr, values):
        raise NotImplementedError("Derived classes must implement this method")

    def __getattr__(self, attr):
        node = self.__dict__.get("_node")
        if node is None or not hasattr(node, attr):
            if attr in self.__dict__:
                return self.__dict__[attr]
            raise AttributeError(
                f"Object '{self.__class__.__name__}' has no attribute '{attr}'")

        if 0 in self.dimensions:
            return np.empty(0, dtype=type(None))

        value = getattr(node, attr)
        value_shape = np.shape(value)
        dtype = value.dtype if isinstance(value, np.ndarray) else type(value)
        np_gen = np.zeros if dtype in (int, float, bool) else np.empty
        value_grid = np_gen((*self.dimensions, *value_shape), dtype=dtype)

        indices = itertools.product(*map(enumerate, self.indices))
        for (i, x), (j, y), (k, z) in indices:
            err = node.call_method("override_index", index=[x, y, z])
            assert err == 0
            value_grid[i, j, k] = getattr(node, attr)

        if value_shape == (1,):
            value_grid = np.squeeze(value_grid, axis=-1)
        return utils.array_locked(value_grid)

    def __setattr__(self, attr, values):
        node = self.__dict__.get("_node")
        if node is None or not hasattr(node, attr):
            self.__dict__[attr] = values
            return
        elif 0 in self.dimensions:
            raise AttributeError(
                f"Cannot set properties of an empty '{self.__class__.__name__}' object")

        values = np.copy(values)
        value_shape = np.shape(getattr(node, attr))
        target_shape = (*self.dimensions, *value_shape)

        # broadcast if only one element was provided
        if values.shape == value_shape:
            values = np.full(target_shape, values)

        if values.shape != target_shape:
            raise ValueError(
                f"Input-dimensions of '{attr}' array {values.shape} does not match slice dimensions {target_shape}")

        self.setter_checks(attr, values)

        indices = itertools.product(*map(enumerate, self.indices))
        for (i, x), (j, y), (k, z) in indices:
            err = node.call_method("override_index", index=[x, y, z])
            assert err == 0
            setattr(node, attr, values[i, j, k])
