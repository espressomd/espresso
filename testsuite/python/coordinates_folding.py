#
# Copyright (C) 2024 The ESPResSo project
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

import sys
import pathlib
import tempfile
import contextlib
import unittest as ut
import numpy as np
import espressomd
import espressomd.io.writer
import espressomd.accumulators
import espressomd.observables
import espressomd.io.writer
with contextlib.suppress(ImportError):
    import h5py  # h5py has to be imported *after* espressomd (MPI)


class Test(ut.TestCase):

    system = espressomd.System(box_l=[10., 10., 10.])
    system.cell_system.skin = 0.4
    system.time_step = 0.01

    def tearDown(self):
        self.system.part.clear()
        self.system.auto_update_accumulators.clear()
        self.system.cell_system.skin = 0.4
        self.system.time_step = 0.01
        self.system.time = 0.

    def check_trajectory(self, pos_unfolded_ref, **kwargs):
        tol = {"atol": 1e-9, "rtol": 0.}
        box_l = np.copy(self.system.box_l)
        pos_unfolded_ref = np.copy(pos_unfolded_ref)
        if "image_box" in kwargs:
            image_box = np.copy(kwargs["image_box"])
            np.testing.assert_equal(image_box[:, 0], 5 * [+0] + 5 * [+1])
            np.testing.assert_equal(image_box[:, 1], 5 * [+3] + 5 * [+2])
            np.testing.assert_equal(image_box[:, 2], 5 * [-2] + 5 * [-1])
        if "pos_unfolded" in kwargs:
            pos_unfolded = np.copy(kwargs["pos_unfolded"])
            np.testing.assert_allclose(pos_unfolded, pos_unfolded_ref, **tol)
        if "pos_folded" in kwargs:
            pos_folded = np.copy(kwargs["pos_folded"])
            image_box = np.copy(kwargs["image_box"])
            pos_folded_ref = pos_unfolded_ref - image_box * box_l
            np.testing.assert_allclose(pos_folded, pos_folded_ref, **tol)

    def check_folding(self, description, skin, le_flag, temp_path):
        has_hdf5 = "h5py" in sys.modules and espressomd.has_features(["H5MD"])
        system = self.system
        temp_file = temp_path / f"skin_{skin:.5f}_le_{le_flag}.h5"
        system.part.clear()
        system.auto_update_accumulators.clear()
        system.lees_edwards.protocol = None
        system.cell_system.skin = skin
        system.time_step = 0.01
        system.time = 0.
        v0 = np.array([0.1, -0.1, 0.1])
        x0 = np.array([9.95 - 1e-9, 30.05 + 1e-9, -10.05 - 1e-9])
        p = system.part.add(pos=x0, v=v0)
        if le_flag:
            protocol = espressomd.lees_edwards.LinearShear(
                shear_velocity=0., initial_pos_offset=0.5, time_0=0.)
            system.lees_edwards.set_boundary_conditions(
                shear_direction="y", shear_plane_normal="x", protocol=protocol)
        obs = espressomd.observables.ParticlePositions(ids=[0])
        acc = espressomd.accumulators.TimeSeries(obs=obs, delta_N=10)
        system.auto_update_accumulators.add(acc)
        if has_hdf5:
            h5 = espressomd.io.writer.h5md.H5md(file_path=str(temp_file))
        pos_unfolded_ref = []
        pos_unfolded = []
        pos_folded = []
        image_box = []
        for _ in range(10):
            system.integrator.run(10)
            pos = x0 + v0 * system.time
            if le_flag and pos[0] > system.box_l[1]:
                pos[1] -= protocol.initial_pos_offset
            pos_unfolded_ref.append(pos)
            pos_unfolded.append(p.pos)
            pos_folded.append(p.pos_folded)
            image_box.append(p.image_box)
            if has_hdf5:
                h5.write()
        if has_hdf5:
            h5.flush()
            h5.close()
        with self.subTest(msg=f"{description}; trajectory from particle handle"):
            self.check_trajectory(pos_unfolded_ref, image_box=image_box,
                                  pos_folded=pos_folded, pos_unfolded=pos_unfolded)
        with self.subTest(msg=f"{description}; trajectory from observable"):
            obs_lag = (acc.delta_N - 1) * system.time_step * v0
            obs_unfolded = acc.time_series().reshape((-1, 3)) + obs_lag
            self.check_trajectory(pos_unfolded_ref, pos_unfolded=obs_unfolded)
        if has_hdf5:
            with self.subTest(msg=f"{description}; trajectory from hdf5 file"):
                with h5py.File(temp_file, "r") as h5:
                    prop = "particles/atoms/{}/value"
                    h5py_pos = h5[prop.format("position")][:].reshape((-1, 3))
                    h5py_img = h5[prop.format("image")][:].reshape((-1, 3))
                self.check_trajectory(
                    pos_unfolded_ref, image_box=h5py_img, pos_folded=h5py_pos)

    def test_folding(self):
        scenarios = [
            (0.400, "scenario: local cell resort only at first time step"),
            (0.040, "scenario: local cell resort in 1 cell at every time step"),
            (0.004, "scenario: local cell resort in 5 cells at every time step"),
            (0.000, "scenario: local cell resort in 10 cells at every time step"),
        ]
        with tempfile.TemporaryDirectory() as temp_dir_name:
            path = pathlib.Path(temp_dir_name).resolve()
            for le_flag in [False, True]:
                for skin, description in scenarios:
                    le_qualifier = "without" if le_flag else "with"
                    scenario = f"{description}; {le_qualifier} Lees-Edwards"
                    self.check_folding(scenario, skin, le_flag, path)


if __name__ == "__main__":
    ut.main()
