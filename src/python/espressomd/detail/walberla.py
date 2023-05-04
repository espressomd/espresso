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

import numpy as np

import espressomd.code_features
from espressomd.script_interface import ScriptInterfaceHelper


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

    def __getitem__(self, vtk_uid):
        return self.map[vtk_uid]

    def __contains__(self, vtk_uid):
        return vtk_uid in self.map

    def _register_vtk_object(self, vtk_obj):
        vtk_uid = vtk_obj.vtk_uid
        self.map[vtk_uid] = vtk_obj

    def __getstate__(self):
        return self.map

    def __setstate__(self, active_vtk_objects):
        self.map = {}
        for vtk_uid, vtk_obj in active_vtk_objects.items():
            self.map[vtk_uid] = vtk_obj


class VTKOutputBase(ScriptInterfaceHelper):

    def __init__(self, *args, **kwargs):
        if not espressomd.code_features.has_features("WALBERLA"):
            raise NotImplementedError("Feature WALBERLA not compiled in")
        if "sip" not in kwargs:
            params = self.default_params()
            params.update(kwargs)
            if isinstance(params["observables"], str):
                params["observables"] = [params["observables"]]
            super().__init__(*args, **params)
        else:
            super().__init__(**kwargs)
        self._add_to_registry()

    def valid_observables(self):
        return set(self.call_method("get_valid_observable_names"))

    def default_params(self):
        return {"delta_N": 0, "enabled": True, "execution_count": 0,
                "base_folder": "vtk_out", "prefix": "simulation_step"}
