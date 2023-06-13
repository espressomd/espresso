#
# Copyright (C) 2021-2022 The ESPResSo project
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

import unittest as ut
import unittest_decorators as utx
import espressomd
import espressomd.electrostatics
import espressomd.electrostatic_extensions
import numpy as np


@utx.skipIfMissingFeatures(["ELECTROSTATICS"])
class Test(ut.TestCase):
    system = espressomd.System(box_l=[20., 20., 20.])
    system.cell_system.skin = 0.4
    system.time_step = 0.01

    def tearDown(self):
        self.system.part.clear()
        self.system.actors.clear()
        self.system.thermostat.turn_off()
        self.system.integrator.set_vv()

    def add_icc_particles(self):
        pos = [[0., 0., 0.], [1., 0., 0.]]
        q = [0.0001, -0.0001]
        p_slice = self.system.part.add(pos=pos, q=q)
        areas = self.system.box_l[0] * self.system.box_l[1] / 2. * np.ones(2)
        normals = 2 * [(0., 0., 1.)]
        return p_slice, normals, areas

    def setup_icc_particles_and_solver(self, **icc_params):
        part_slice, normals, areas = self.add_icc_particles()
        icc = espressomd.electrostatic_extensions.ICC(
            n_icc=len(part_slice),
            normals=normals,
            areas=areas,
            epsilons=np.ones_like(areas),
            first_id=part_slice.id[0],
            **icc_params
        )
        return icc, part_slice

    def valid_p3m_parameters(self):
        return {"prefactor": 1., "mesh": 32, "cao": 7, "accuracy": 1e-5,
                "r_cut": 1.25625, "alpha": 1.50505, "tune": False,
                "check_neutrality": False}

    def test_getters_and_setters(self):
        part_slice, normals, areas = self.add_icc_particles()

        params = {"n_icc": len(part_slice),
                  "normals": normals,
                  "areas": areas,
                  "epsilons": np.ones_like(areas),
                  "first_id": part_slice.id[0]}

        icc = espressomd.electrostatic_extensions.ICC(**params)
        icc_params = icc.get_params()
        for key, value in params.items():
            np.testing.assert_allclose(value, np.copy(icc_params[key]))
            with self.assertRaisesRegex(RuntimeError, f"Parameter '{key}' is read-only"):
                setattr(icc, key, 5)

    def test_invalid_parameters(self):
        part_slice, normals, areas = self.add_icc_particles()

        valid_params = {"n_icc": len(part_slice),
                        "normals": normals,
                        "areas": areas,
                        "epsilons": np.ones_like(areas),
                        "first_id": part_slice.id[0]}

        invalid_params = [({"n_icc": -1}, "Parameter 'n_icc' must be >= 1"),
                          ({"n_icc": 0}, "Parameter 'n_icc' must be >= 1"),
                          ({"first_id": -1}, "Parameter 'first_id' must be >= 0"),
                          ({"max_iterations": -1},
                           "Parameter 'max_iterations' must be > 0"),
                          ({"convergence": -1.},
                           "Parameter 'convergence' must be > 0"),
                          ({"relaxation": -0.1},
                           "Parameter 'relaxation' must be >= 0 and <= 2"),
                          ({"relaxation": 2.1},
                           "Parameter 'relaxation' must be >= 0 and <= 2"),
                          ({"eps_out": -1.}, "Parameter 'eps_out' must be > 0"),
                          ({"ext_field": 0.}, 'A single value was given but 3 were expected'), ]

        for kwargs, error in invalid_params:
            params = valid_params.copy()
            params.update(kwargs)
            with self.assertRaisesRegex(ValueError, error):
                espressomd.electrostatic_extensions.ICC(**params)

    @utx.skipIfMissingFeatures(["P3M"])
    def test_exceptions_small_r_cut(self):
        icc, _ = self.setup_icc_particles_and_solver(max_iterations=1)
        p3m = espressomd.electrostatics.P3M(
            prefactor=1., mesh=32, cao=7, accuracy=1e-5, r_cut=0.01875,
            alpha=0.005, tune=False, check_neutrality=False)
        self.system.actors.add(p3m)
        self.system.actors.add(icc)

        with self.assertRaisesRegex(Exception, "ICC found zero electric field on a charge"):
            self.system.integrator.run(0)

    @utx.skipIfMissingFeatures(["P3M"])
    def test_exceptions_large_r_cut(self):
        icc, (_, p) = self.setup_icc_particles_and_solver(
            max_iterations=1, convergence=10.)
        p3m = espressomd.electrostatics.P3M(
            check_complex_residuals=False, **self.valid_p3m_parameters())

        self.system.actors.add(p3m)
        self.system.actors.add(icc)
        with self.assertRaisesRegex(Exception, f"Particle with id {p.id} has a charge .+ that is too large for the ICC algorithm"):
            p.q = 1e8
            self.system.integrator.run(0)

        self.system.actors.remove(icc)
        self.system.part.clear()
        icc, (_, p) = self.setup_icc_particles_and_solver(max_iterations=1)
        self.system.actors.add(icc)
        with self.assertRaisesRegex(Exception, "ICC failed to converge in the given number of maximal steps"):
            p.q = 1.
            self.system.integrator.run(0)
        with self.assertRaisesRegex(Exception, "ICC found zero electric charge on a particle"):
            p.q = 0.
            self.system.integrator.run(0)

    @utx.skipIfMissingGPU()
    @utx.skipIfMissingFeatures(["P3M"])
    def test_exceptions_gpu(self):
        icc, _ = self.setup_icc_particles_and_solver()
        p3m = espressomd.electrostatics.P3MGPU(**self.valid_p3m_parameters())

        self.system.actors.add(p3m)
        with self.assertRaisesRegex(RuntimeError, "ICC does not work with P3MGPU"):
            self.system.actors.add(icc)
        self.assertEqual(len(self.system.actors), 1)
        self.system.integrator.run(0)

        self.system.actors.clear()
        elc = espressomd.electrostatics.ELC(
            actor=p3m, gap_size=5., maxPWerror=1e-3)
        self.system.actors.add(elc)
        with self.assertRaisesRegex(RuntimeError, "ICC does not work with P3MGPU"):
            self.system.actors.add(icc)
        self.assertEqual(len(self.system.actors), 1)
        self.system.integrator.run(0)

    @utx.skipIfMissingFeatures(["P3M"])
    def test_exceptions_elc(self):
        icc, _ = self.setup_icc_particles_and_solver()
        p3m = espressomd.electrostatics.P3M(**self.valid_p3m_parameters())
        elc = espressomd.electrostatics.ELC(
            actor=p3m, gap_size=5., maxPWerror=1e-3, pot_diff=-3.,
            delta_mid_top=-1., delta_mid_bot=-1., const_pot=True)

        self.system.actors.add(elc)
        with self.assertRaisesRegex(RuntimeError, "ICC conflicts with ELC dielectric contrast"):
            self.system.actors.add(icc)
        self.assertEqual(len(self.system.actors), 1)
        self.system.integrator.run(0)
        self.system.actors.clear()

        # valid ELC actor should pass sanity checks
        elc = espressomd.electrostatics.ELC(
            actor=p3m, gap_size=5., maxPWerror=1e-3)
        self.system.actors.add(elc)
        self.system.actors.add(icc)
        self.system.part.clear()
        self.system.integrator.run(0)

    def test_exceptions_dh(self):
        icc, _ = self.setup_icc_particles_and_solver()
        solver = espressomd.electrostatics.DH(
            prefactor=2., kappa=0.8, r_cut=1.2)

        self.system.actors.add(solver)
        with self.assertRaisesRegex(RuntimeError, "ICC does not work with DebyeHueckel"):
            self.system.actors.add(icc)
        self.assertEqual(len(self.system.actors), 1)
        self.system.integrator.run(0)

    def test_exceptions_rf(self):
        icc, _ = self.setup_icc_particles_and_solver()
        solver = espressomd.electrostatics.ReactionField(
            prefactor=1., kappa=2., epsilon1=1., epsilon2=2., r_cut=2.)

        self.system.actors.add(solver)
        with self.assertRaisesRegex(RuntimeError, "ICC does not work with ReactionField"):
            self.system.actors.add(icc)
        self.assertEqual(len(self.system.actors), 1)
        self.system.integrator.run(0)

    @utx.skipIfMissingFeatures(["NPT", "P3M"])
    def test_exceptions_npt(self):
        icc, _ = self.setup_icc_particles_and_solver()
        p3m = espressomd.electrostatics.P3M(**self.valid_p3m_parameters())

        self.system.actors.add(p3m)
        self.system.thermostat.set_npt(kT=1., gamma0=2., gammav=0.004, seed=42)
        self.system.integrator.set_isotropic_npt(ext_pressure=2., piston=0.001)
        with self.assertRaisesRegex(RuntimeError, "ICC does not work in the NPT ensemble"):
            self.system.actors.add(icc)
        self.assertEqual(len(self.system.actors), 1)
        self.system.integrator.run(0)


if __name__ == "__main__":
    ut.main()
