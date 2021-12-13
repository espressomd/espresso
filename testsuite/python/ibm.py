#
# Copyright (C) 2013-2019 The ESPResSo project
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
import itertools
import unittest as ut
import unittest_decorators as utx

import espressomd
import espressomd.interactions
import espressomd.virtual_sites


@utx.skipIfMissingFeatures(['VIRTUAL_SITES_INERTIALESS_TRACERS'])
class IBM(ut.TestCase):
    '''Test IBM implementation with a Langevin thermostat.'''
    system = espressomd.System(box_l=3 * [8.])
    system.time_step = 0.06
    system.cell_system.skin = 0.1

    def tearDown(self):
        self.system.part.clear()
        self.system.actors.clear()
        self.system.thermostat.turn_off()

    def compute_dihedral_angle(self, pos0, pos1, pos2, pos3):
        # first normal vector
        n1 = np.cross((pos1 - pos0), (pos2 - pos0))
        n2 = np.cross((pos2 - pos0), (pos3 - pos0))

        norm1 = np.linalg.norm(n1)
        norm2 = np.linalg.norm(n2)
        n1 = n1 / norm1
        n2 = n2 / norm2

        cos_alpha = min(1, np.dot(n1, n2))
        alpha = np.arccos(cos_alpha)
        return alpha

    def test_tribend(self):
        # two triangles with bending interaction
        # move nodes, should relax back

        system = self.system
        system.virtual_sites = espressomd.virtual_sites.VirtualSitesInertialessTracers()
        system.thermostat.set_langevin(kT=0, gamma=10, seed=1)

        # Add four particles
        p0 = system.part.add(pos=[4, 4, 4])
        p1 = system.part.add(pos=[4, 4, 5])
        p2 = system.part.add(pos=[4, 5, 5])
        p3 = system.part.add(pos=[4, 5, 4])

        # Add first triel, weak modulus
        tri1 = espressomd.interactions.IBM_Triel(
            ind1=p0.id, ind2=p1.id, ind3=p2.id, elasticLaw="Skalak", k1=0.1, k2=0, maxDist=2.4)
        system.bonded_inter.add(tri1)
        p0.add_bond((tri1, p1, p2))

        # Add second triel, strong modulus
        tri2 = espressomd.interactions.IBM_Triel(
            ind1=p0.id, ind2=p2.id, ind3=p3.id, elasticLaw="Skalak", k1=10, k2=0, maxDist=2.4)
        system.bonded_inter.add(tri2)
        p0.add_bond((tri2, p2, p3))

        # Add bending
        tribend = espressomd.interactions.IBM_Tribend(
            ind1=p0.id, ind2=p1.id, ind3=p2.id, ind4=p3.id, kb=1, refShape="Initial")
        system.bonded_inter.add(tribend)
        p0.add_bond((tribend, p1, p2, p3))

        # twist
        system.part.all().pos = system.part.all().pos + np.random.random((4, 3))

        # Perform integration
        system.integrator.run(200)
        angle = self.compute_dihedral_angle(p0.pos, p1.pos, p2.pos, p3.pos)
        self.assertLess(angle, 2E-2)

        # IBM doesn't implement energy and pressure kernels.
        energy = self.system.analysis.energy()
        pressure = self.system.analysis.pressure()
        self.assertAlmostEqual(energy['bonded'], 0., delta=1e-10)
        self.assertAlmostEqual(pressure['bonded'], 0., delta=1e-10)

    def test_triel(self):
        system = self.system
        system.virtual_sites = espressomd.virtual_sites.VirtualSitesInertialessTracers()
        system.thermostat.set_langevin(kT=0, gamma=1, seed=1)

        # Add particles: 0-2 are not bonded, 3-5 are bonded
        non_bound = system.part.add(pos=[[4, 4, 4], [4, 4, 5], [4, 5, 5]])

        p3 = system.part.add(pos=[1, 4, 4])
        p4 = system.part.add(pos=[1, 4, 5])
        p5 = system.part.add(pos=[1, 5, 5])
        all_partcls = system.part.all()

        # Add triel for 3-5
        tri = espressomd.interactions.IBM_Triel(
            ind1=p3.id, ind2=p4.id, ind3=p5.id, elasticLaw="Skalak", k1=15,
            k2=0, maxDist=2.4)
        system.bonded_inter.add(tri)
        p3.add_bond((tri, p4, p5))

        all_partcls.pos = all_partcls.pos + np.array((
            (0, 0, 0), (1, -.2, .3), (1, 1, 1),
            (0, 0, 0), (1, -.2, .3), (1, 1, 1)))

        distorted_pos = np.copy(non_bound.pos)

        system.integrator.run(110)
        dist1bound = system.distance(p3, p4)
        dist2bound = system.distance(p3, p5)

        # check bonded particles. Distance should restore to initial config
        self.assertAlmostEqual(dist1bound, 1, delta=0.05)
        self.assertAlmostEqual(dist2bound, np.sqrt(2), delta=0.05)

        # check not bonded particles. Positions should still be distorted
        np.testing.assert_allclose(np.copy(non_bound.pos), distorted_pos)

    def test_volcons(self):
        '''Check volume conservation forces on a simple mesh (cube).'''
        system = self.system
        system.virtual_sites = espressomd.virtual_sites.VirtualSitesOff()
        system.thermostat.set_langevin(kT=0, gamma=1, seed=1)

        # Place particles on a cube.
        positions = list(itertools.product((0, 1), repeat=3))
        positions = positions[:4] + positions[6:] + positions[4:6]
        positions = np.array(positions) - 0.5
        mesh_center_ref = np.copy(system.box_l) / 2.
        partcls = system.part.add(pos=positions + mesh_center_ref)

        # Divide the cube. All triangle normals must point inside the mesh.
        # Use the right hand rule to determine the order of the indices.
        triangles = [(0, 1, 2), (1, 3, 2),
                     (2, 3, 4), (3, 5, 4),
                     (4, 5, 6), (5, 7, 6),
                     (6, 7, 0), (7, 1, 0),
                     (0, 2, 4), (0, 4, 6),
                     (1, 5, 3), (1, 7, 5)]
        # Add triangle bonds that don't contribute to the force (infinite
        # elasticity). These bonds are needed to calculate the mesh volume.
        for id1, id2, id3 in triangles:
            bond = espressomd.interactions.IBM_Triel(
                ind1=id3, ind2=id2, ind3=id1, elasticLaw="Skalak", k1=0., k2=0., maxDist=3)
            system.bonded_inter.add(bond)
            system.part.by_id(id1).add_bond((bond, id2, id3))

        # Add volume conservation force.
        KAPPA_V = 0.01
        volCons = espressomd.interactions.IBM_VolCons(
            softID=15, kappaV=KAPPA_V)
        system.bonded_inter.add(volCons)
        for p in system.part:
            p.add_bond((volCons,))

        # Run the integrator to initialize the mesh reference volume.
        system.integrator.run(0, recalc_forces=True)
        self.assertAlmostEqual(volCons.current_volume(), 1., delta=1e-10)

        # The restorative force is zero at the moment.
        np.testing.assert_almost_equal(np.copy(partcls.f), 0.)

        # Double the cube dimensions. The volume increases by a factor of 8.
        partcls.pos = 2. * positions + mesh_center_ref
        system.integrator.run(0, recalc_forces=True)
        self.assertAlmostEqual(volCons.current_volume(), 8., delta=1e-10)

        # Reference forces for that particular mesh geometry.
        ref_forces = 1.75 * KAPPA_V * np.array(
            [(1, 2, 2), (2, 1, -2), (2, -1, 1), (1, -2, -1),
             (-1, -2, 2), (-2, -1, -2), (-2, 1, 1), (-1, 2, -1)])
        np.testing.assert_almost_equal(
            np.copy(partcls.f), ref_forces)

        # IBM doesn't implement energy and pressure kernels.
        energy = self.system.analysis.energy()
        pressure = self.system.analysis.pressure()
        self.assertAlmostEqual(energy['bonded'], 0., delta=1e-10)
        self.assertAlmostEqual(pressure['bonded'], 0., delta=1e-10)

        # Check the cube is shrinking. The geometrical center shouldn't move.
        volume_diff_ref = 0.1  # arbitrary, but should work for the given setup
        # warmup
        system.integrator.run(10)
        # sampling
        previous_volume = volCons.current_volume()
        for _ in range(10):
            system.integrator.run(5)
            current_volume = volCons.current_volume()
            volume_diff = previous_volume - current_volume
            self.assertLess(current_volume, previous_volume)
            self.assertGreater(volume_diff, volume_diff_ref)
            previous_volume = current_volume
            mesh_center = np.mean(partcls.pos, axis=0)
            np.testing.assert_allclose(mesh_center, mesh_center_ref, rtol=1e-3)
        # Halve the cube dimensions. The volume decreases by a factor of 8.
        partcls.pos = 0.5 * positions + mesh_center_ref
        system.integrator.run(0, recalc_forces=True)
        self.assertAlmostEqual(volCons.current_volume(), 1. / 8., delta=1e-10)

        # Check the cube is expanding. The geometrical center shouldn't move.
        volume_diff_ref = 0.005  # arbitrary, but should work for the given setup
        # warmup
        system.integrator.run(40)
        # sampling
        previous_volume = volCons.current_volume()
        for _ in range(10):
            system.integrator.run(5)
            current_volume = volCons.current_volume()
            volume_diff = current_volume - previous_volume
            self.assertGreater(current_volume, previous_volume)
            self.assertGreater(volume_diff, volume_diff_ref)
            previous_volume = current_volume
            mesh_center = np.mean(partcls.pos, axis=0)
            np.testing.assert_allclose(mesh_center, mesh_center_ref, rtol=1e-3)


if __name__ == "__main__":
    ut.main()
