#
# Copyright (C) 2026 The ESPResSo project
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
import numpy as np
import espressomd
import espressomd.electrostatics
import espressomd.electrostatic_extensions


class ForceReduction(ut.TestCase):
    """
    Repeated force calculations must give the same forces. The kernels
    scatter forces (and P3M charges) into per-thread replicas, which are
    summed and zeroed after every force calculation; stale replica entries
    would add the previous force calculation's contributions to the next one.
    Run with several OpenMP threads to have more than one replica (the
    CMake entry asks for 2; with ESPRESSO_TEST_NT=1 there is only one
    replica and only the single-replica path is tested). The
    particles are fixed, so that integration steps recompute the forces of
    an unchanged configuration.
    """
    box_l = 12.
    system = espressomd.System(box_l=3 * [box_l])
    system.time_step = 0.01
    system.cell_system.skin = 0.4

    def setUp(self):
        self.system.cell_system.set_regular_decomposition()

    def tearDown(self):
        self.system.electrostatics.clear()
        self.system.part.clear()
        self.system.non_bonded_inter[0, 0].lennard_jones.deactivate()

    def add_particles(self, n, charged=False):
        rng = np.random.default_rng(42)
        # jittered lattice: no close contacts
        n_side = int(np.ceil(n**(1. / 3.)))
        a = self.box_l / n_side
        g = (np.arange(n_side) + 0.5) * a
        pos = np.array(np.meshgrid(g, g, g, indexing="ij")).reshape(3, -1).T
        pos = pos[:n] + rng.uniform(-0.1 * a, 0.1 * a, (n, 3))
        fix = n * [3 * [True]]
        if charged:
            q = np.where(np.arange(n) % 2 == 0, 1., -1.)
            return self.system.part.add(pos=pos, q=q, fix=fix)
        return self.system.part.add(pos=pos, fix=fix)

    def lj_forces_brute_force(self, eps, sig, cut):
        partcls = self.system.part.all()
        pos = np.copy(partcls.pos)
        d = pos[:, np.newaxis, :] - pos[np.newaxis, :, :]
        d -= self.box_l * np.round(d / self.box_l)
        r2 = np.sum(d**2, axis=-1)
        np.fill_diagonal(r2, np.inf)
        mask = r2 < cut**2
        r2 = np.where(mask, r2, 1.)
        sr6 = (sig**2 / r2)**3
        fmag = np.where(mask, 24. * eps * (2. * sr6**2 - sr6) / r2, 0.)
        return np.sum(fmag[:, :, np.newaxis] * d, axis=1)

    def check_repeated_forces(self, check):
        system = self.system
        system.integrator.run(0, recalc_forces=True)
        check()
        # force calculations in integration steps, without and with
        # interleaved energy and pressure calculations (they also reset the
        # force buffers)
        system.integrator.run(3)
        check()
        system.analysis.energy()
        system.integrator.run(3)
        check()
        system.analysis.pressure()
        system.integrator.run(3)
        check()

    @utx.skipIfMissingFeatures(["EXTERNAL_FORCES"])
    def test_lj(self):
        system = self.system
        eps, sig, cut = 1., 1., 2.5
        system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=eps, sigma=sig, cutoff=cut, shift=0.)
        partcls = self.add_particles(800)

        def check():
            np.testing.assert_allclose(
                np.copy(partcls.f), self.lj_forces_brute_force(eps, sig, cut),
                atol=1e-10, rtol=0.)

        self.check_repeated_forces(check)
        # moved particles: Verlet list rebuild with unchanged buffer extents
        partcls.pos = partcls.pos + 0.1
        self.check_repeated_forces(check)
        # new particles: force buffers and replicas are reallocated
        system.part.add(pos=[[0.1, 0.1, 0.1], [0.9, 0.1, 0.2]],
                        fix=2 * [3 * [True]])
        partcls = system.part.all()
        self.check_repeated_forces(check)
        partcls.pos = partcls.pos - 0.05
        self.check_repeated_forces(check)

    @utx.skipIfMissingFeatures(["P3M", "EXTERNAL_FORCES"])
    def test_p3m(self):
        system = self.system
        partcls = self.add_particles(216, charged=True)
        system.electrostatics.solver = espressomd.electrostatics.P3M(
            prefactor=1., accuracy=1e-4, mesh=[16, 16, 16], cao=5,
            r_cut=2.5, alpha=1.0, tune=False)
        system.integrator.run(0, recalc_forces=True)
        f_ref = np.copy(partcls.f)
        self.assertGreater(np.max(np.abs(f_ref)), 0.1)

        def check():
            np.testing.assert_allclose(np.copy(partcls.f), f_ref,
                                       atol=1e-12 * np.max(np.abs(f_ref)),
                                       rtol=0.)

        self.check_repeated_forces(check)

    @utx.skipIfMissingFeatures(["P3M", "EXTERNAL_FORCES"])
    def test_icc(self):
        # ICC computes its own forces through the force buffers between
        # two force calculations of the integrator
        system = self.system
        n_side = 6
        g = (np.arange(n_side) + 0.5) * self.box_l / n_side
        x, y = np.meshgrid(g, g)
        wall = np.stack((x.flatten(), y.flatten(),
                         np.full(n_side**2, 1.)), axis=-1)
        # ICC needs non-zero initial charges
        q_init = np.where(np.arange(n_side**2) % 2 == 0, 1e-4, -1e-4)
        icc_partcls = system.part.add(pos=wall, q=q_init,
                                      fix=n_side**2 * [3 * [True]])
        ions = system.part.add(pos=[[5., 6., 4.], [7., 6., 8.]], q=[1., -1.],
                               fix=2 * [3 * [True]])
        system.electrostatics.solver = espressomd.electrostatics.P3M(
            prefactor=1., accuracy=1e-4, mesh=[16, 16, 16], cao=5,
            r_cut=2.5, alpha=1.0, tune=False, check_neutrality=False)
        normals = np.zeros((n_side**2, 3))
        normals[:, 2] = 1.
        system.electrostatics.extension = \
            espressomd.electrostatic_extensions.ICC(
                n_icc=n_side**2, first_id=icc_partcls.id[0],
                normals=normals, areas=np.full(n_side**2, 4.),
                epsilons=np.full(n_side**2, 10.),
                sigmas=np.zeros(n_side**2), convergence=1e-9,
                max_iterations=200, eps_out=1., relaxation=0.75,
                ext_field=[0., 0., 0.])
        system.integrator.run(0, recalc_forces=True)
        f_ref = np.copy(ions.f)
        self.assertGreater(np.max(np.abs(f_ref)), 1e-3)

        def check():
            np.testing.assert_allclose(np.copy(ions.f), f_ref,
                                       atol=1e-6 * np.max(np.abs(f_ref)),
                                       rtol=0.)

        self.check_repeated_forces(check)


if __name__ == "__main__":
    ut.main()
