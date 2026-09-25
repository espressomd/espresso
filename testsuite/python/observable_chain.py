#
# Copyright (C) 2019-2026 The ESPResSo project
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
import unittest as ut
import espressomd
import numpy as np
import espressomd.observables
import sys
import subprocess


def cos_persistence_angles(positions):
    """ Python implementation for PersistenceAngles observable.

    """
    no_of_bonds = positions.shape[0] - 1
    no_of_angles = no_of_bonds - 1
    bond_vecs = positions[1:] - positions[:-1]
    bond_vecs = np.divide(bond_vecs, np.linalg.norm(
        bond_vecs, axis=1)[:, np.newaxis])
    angles = np.zeros(no_of_angles)
    for i in range(no_of_angles):
        average = 0.0
        for j in range(no_of_angles - i):
            average += np.dot(bond_vecs[j], bond_vecs[j + i + 1])
        angles[i] = average / (no_of_angles - i)
    return angles


class ObservableTests(ut.TestCase):
    n_tries = 50
    n_parts = 5
    box_l = 5.
    system = espressomd.System(box_l=3 * [box_l])
    system.periodicity = [True, True, True]
    system.time_step = 0.01
    system.cell_system.skin = 0.2 * box_l

    def setUp(self):
        for i in range(self.n_parts):
            self.system.part.add(pos=[1 + i, 1 + i, 1 + i], id=i)
        self.partcls = self.system.part.all()

    def tearDown(self):
        self.system.part.clear()

    def test_ParticleDistances(self):
        """
        Check ParticleDistances, for a particle pair and for a chain.
        """
        pids = list(range(self.n_parts))
        obs_single = espressomd.observables.ParticleDistances(ids=[0, 1])
        obs_chain = espressomd.observables.ParticleDistances(ids=pids)
        # take periodic boundaries into account: bond length cannot exceed
        # half the box size along the smallest axis
        min_dim = np.min(self.system.box_l)
        max_bond_length = min_dim / 2.01

        for _ in range(self.n_tries):
            # build polymer
            pos = np.zeros((self.n_parts, 3), dtype=float)
            pos[0] = np.random.uniform(low=0, high=min_dim, size=3)
            for i in range(1, self.n_parts):
                pos[i] = pos[i - 1] + np.random.uniform(
                    low=0, high=max_bond_length, size=3)
            self.partcls.pos = pos
            # expected values
            distances = np.linalg.norm(pos[1:] - pos[:-1], axis=1)
            # observed values
            self.system.integrator.run(0)
            res_obs_single = obs_single.calculate()
            res_obs_chain = obs_chain.calculate()
            # checks
            self.assertEqual(np.prod(res_obs_single.shape), 1)
            self.assertEqual(np.prod(res_obs_chain.shape), self.n_parts - 1)
            self.assertAlmostEqual(res_obs_single[0], distances[0], places=9)
            np.testing.assert_array_almost_equal(
                res_obs_chain, distances, decimal=9,
                err_msg="Data did not agree for observable ParticleDistances")

        # check exceptions
        for i in range(2):
            with self.assertRaises(RuntimeError):
                espressomd.observables.ParticleDistances(ids=np.arange(i))

    def test_BondAngles(self):
        """
        Check BondAngles, for a particle triple and for a chain.
        """
        pids = list(range(self.n_parts))
        obs_single = espressomd.observables.BondAngles(ids=[0, 1, 2])
        obs_chain = espressomd.observables.BondAngles(ids=pids)
        # take periodic boundaries into account: bond length cannot exceed
        # half the box size along the smallest axis
        min_dim = np.min(self.system.box_l)
        max_bond_length = min_dim / 2.01

        for _ in range(self.n_tries):
            # build polymer
            pos = np.zeros((self.n_parts, 3), dtype=float)
            pos[0] = np.random.uniform(low=0, high=min_dim, size=3)
            for i in range(1, self.n_parts):
                pos[i] = pos[i - 1] + np.random.uniform(
                    low=0, high=max_bond_length, size=3)
            self.partcls.pos = pos
            # expected values
            v1 = pos[:-2] - pos[1:-1]
            v2 = pos[2:] - pos[1:-1]
            l1 = np.linalg.norm(v1, axis=1)
            l2 = np.linalg.norm(v2, axis=1)
            angles = np.arccos((v1 * v2).sum(1) / l1 / l2)
            # observed values
            self.system.integrator.run(0)
            res_obs_single = obs_single.calculate()
            res_obs_chain = obs_chain.calculate()
            # checks
            self.assertEqual(np.prod(res_obs_single.shape), 1)
            self.assertEqual(np.prod(res_obs_chain.shape), self.n_parts - 2)
            self.assertAlmostEqual(res_obs_single[0], angles[0], places=9)
            np.testing.assert_array_almost_equal(
                res_obs_chain, angles, decimal=9,
                err_msg="Data did not agree for observable BondAngles")

        # check exceptions
        for i in range(3):
            with self.assertRaises(RuntimeError):
                espressomd.observables.BondAngles(ids=np.arange(i))

    def test_BondDihedrals(self):
        """
        Check BondDihedrals, for a particle quadruple and for a chain.
        """
        def rotate_vector(v, k, phi):
            """Rotates vector v around unit vector k by angle phi.
            Uses Rodrigues' rotation formula."""
            vrot = v * np.cos(phi) + np.cross(k, v) * \
                np.sin(phi) + k * np.dot(k, v) * (1.0 - np.cos(phi))
            return vrot

        def rotate_particle(p2, p3, p4, phi):
            """Rotates particle p4 around the axis formed by the bond
            between p2 and p3."""
            k = p3 - p2
            k /= np.linalg.norm(k)
            return p3 + rotate_vector(p4 - p3, k, phi)

        def calculate_dihedral(a, b, c, d):
            v1 = b - a
            v2 = c - b
            v3 = d - c
            b1 = np.cross(v1, v2)
            b2 = np.cross(v2, v3)
            u2 = v2 / np.linalg.norm(v2)
            return np.arctan2(np.dot(np.cross(b1, b2), u2), np.dot(b1, b2))

        def place_particles(bl, offset):
            """Place 5 particles in the XY plane with bond length `bl` and
            bond angle = 120 degrees. The chain is then shifted by `offset`."""
            phi = 2 * np.pi / 3
            pos = np.zeros((self.n_parts, 3), dtype=float)
            pos[0] = [bl * np.cos(phi), bl * np.sin(phi), 0.]
            pos[1] = [0., 0., 0.]
            pos[2] = [bl, 0., 0.]
            pos[3] = pos[2] + [bl * np.cos(np.pi - phi), bl * np.sin(phi), 0.]
            pos[4] = pos[3] + [bl, 0., 0.]
            pos += offset
            self.partcls.pos = pos
            return pos

        pids = list(range(self.n_parts))
        obs_single = espressomd.observables.BondDihedrals(ids=pids[:4])
        obs_chain = espressomd.observables.BondDihedrals(ids=pids)

        # test multiple angles, take periodic boundaries into account
        p0, p4 = self.system.part.by_ids([0, 4])
        for bond_length in [0.1, self.box_l / 2.0]:
            for offset in [1.0, self.box_l / 2.0]:
                for phi in np.arange(0, np.pi, np.pi / 6):
                    # place particles and keep list of unfolded positions
                    pos = place_particles(bond_length, 3 * [offset])
                    # rotate the 1st particle
                    p0.pos = pos[0] = rotate_particle(*pos[1:4, :][::-1],
                                                      phi=phi)
                    # rotate the 5th particle
                    p4.pos = pos[4] = rotate_particle(*pos[2:5, :], phi=phi)
                    # expected values
                    dih1 = calculate_dihedral(*pos[0:4, :][::-1])
                    dih2 = calculate_dihedral(*pos[1:5, :])
                    # observed values
                    self.system.integrator.run(0)
                    res_obs_single = obs_single.calculate()
                    res_obs_chain = obs_chain.calculate()
                    # checks
                    self.assertEqual(np.prod(res_obs_single.shape), 1)
                    self.assertEqual(
                        np.prod(res_obs_chain.shape),
                        self.n_parts - 3)
                    self.assertAlmostEqual(res_obs_single[0], dih1, places=9)
                    np.testing.assert_array_almost_equal(
                        res_obs_chain, [dih1, dih2], decimal=9,
                        err_msg="Data did not agree for observable BondDihedrals")

        # check exceptions
        for i in range(4):
            with self.assertRaises(RuntimeError):
                espressomd.observables.BondDihedrals(ids=np.arange(i))

    def test_CosPersistenceAngles(self):
        # First test: compare with python implementation
        self.system.part.clear()
        partcls = self.system.part.add(pos=np.array(
            [np.linspace(0, self.system.box_l[0], 20)] * 3).T + np.random.random((20, 3)))
        obs = espressomd.observables.CosPersistenceAngles(
            ids=partcls.id)
        np.testing.assert_allclose(
            obs.calculate(), cos_persistence_angles(partcls.pos))
        self.system.part.clear()
        # Second test: place particles with fixed angles and check that the
        # result of PersistenceAngle.calculate()[i] is i*phi
        delta_phi = np.radians(4)
        for i in range(10):
            pos = [np.cos(i * delta_phi), np.sin(i * delta_phi), 0.0]
            self.system.part.add(pos=pos)
        new_partcls = self.system.part.all()
        obs = espressomd.observables.CosPersistenceAngles(
            ids=new_partcls.id)
        expected = np.arange(1, 9) * delta_phi
        np.testing.assert_allclose(obs.calculate(), np.cos(expected))

        # check exceptions
        for i in range(3):
            with self.assertRaises(RuntimeError):
                espressomd.observables.CosPersistenceAngles(ids=np.arange(i))

    def test_nonexistent_id_is_rejected(self):
        """
        Bug #23: the documented "(existing) particle" precondition for
        observable ids is unenforced. Passing a nonexistent id to a
        get_all_particle_positions-family chain observable (here
        ParticleDistances) makes fetch_particles silently drop the id, so
        detail::get_argsort produces an out-of-range index that
        detail::get_all_particle_positions consumes as an out-of-bounds read
        (assert-abort in debug builds, garbage/SIGSEGV in release).

        Correct behavior is a clear RuntimeError at calculate() time.
        """
        # Only particles 0..n_parts-1 exist; 99 does not.
        nonexistent = 99
        self.assertNotIn(nonexistent, list(range(self.n_parts)))
        obs = espressomd.observables.ParticleDistances(
            ids=[0, 1, 2, nonexistent])
        self.system.integrator.run(0)
        with self.assertRaises(RuntimeError):
            obs.calculate()
        # Removing the offending particle is not possible here (it never
        # existed), but valid ids must still work afterwards.
        obs_ok = espressomd.observables.ParticleDistances(ids=[0, 1, 2])
        np.testing.assert_array_equal(obs_ok.calculate().shape, (2,))

    @ut.skipIf(
        system.cell_system.get_state()["n_nodes"] > 1,
        "subprocess re-initializes espressomd; only safe outside an MPI run")
    def test_nonexistent_id_does_not_crash_subprocess(self):
        """
        Bug #23 (crash-safety): the out-of-bounds read on the unfixed code can
        abort or segfault the whole process. Run the trigger out-of-process so
        a crash surfaces as a non-zero/negative return code rather than as a
        clean Python exception. After the fix the child must exit cleanly and
        report that the nonexistent id was REJECTED with an exception.

        This runs only single-rank; the multi-rank collective existence check
        is exercised by ``test_nonexistent_id_is_rejected``.
        """
        child_script = r"""
import espressomd
import espressomd.observables

system = espressomd.System(box_l=[10., 10., 10.])
system.time_step = 0.01
system.cell_system.skin = 0.4
for i in range(3):
    system.part.add(pos=[1. + i, 1. + i, 1. + i], id=i)
system.integrator.run(0)
try:
    obs = espressomd.observables.ParticleDistances(ids=[0, 1, 2, 99])
    res = obs.calculate()
except Exception as err:
    print("REJECTED:" + type(err).__name__)
else:
    print("NO_REJECTION:" + str(res))
"""
        # Run the child with the same interpreter/environment as the parent
        # (pypresso sets PYTHONPATH in os.environ, inherited by the child).
        proc = subprocess.run(
            [sys.executable, "-c", child_script],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines=True, timeout=300)
        out = proc.stdout + proc.stderr
        self.assertGreaterEqual(
            proc.returncode, 0,
            msg=f"child killed by signal {-proc.returncode} (OOB/crash); "
            f"output was:\n{out}")
        self.assertEqual(
            proc.returncode, 0,
            msg=f"child exited with code {proc.returncode}; output:\n{out}")
        self.assertIn(
            "REJECTED:", proc.stdout,
            msg=f"nonexistent observable id was not rejected with an "
            f"exception (OOB read / silently-wrong result); "
            f"output:\n{out}")


if __name__ == "__main__":
    ut.main()
