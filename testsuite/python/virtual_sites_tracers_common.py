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
import espressomd
from espressomd import lb, shapes, lbboundaries
import numpy as np
try:
    from espressomd.virtual_sites import VirtualSitesInertialessTracers, VirtualSitesOff
except ImportError:
    pass
from espressomd.utils import handle_errors


class VirtualSitesTracersCommon:
    box_height = 10.
    box_lw = 8.
    system = espressomd.System(box_l=(box_lw, box_lw, box_height))
    system.time_step = 0.05
    system.cell_system.skin = 0.1

    def reset_lb(self, ext_force_density=(0, 0, 0)):
        self.system.actors.clear()
        self.system.lbboundaries.clear()
        self.lbf = self.LBClass(
            kT=0.0, agrid=1, dens=1, visc=1.8,
            tau=self.system.time_step, ext_force_density=ext_force_density)
        self.system.actors.add(self.lbf)
        self.system.thermostat.set_lb(
            LB_fluid=self.lbf,
            act_on_virtual=False,
            gamma=1)

        # Setup boundaries
        walls = [lbboundaries.LBBoundary() for k in range(2)]
        walls[0].set_params(shape=shapes.Wall(normal=[0, 0, 1], dist=0.5))
        walls[1].set_params(shape=shapes.Wall(normal=[0, 0, -1],
                                              dist=-self.box_height - 0.5))

        for wall in walls:
            self.system.lbboundaries.add(wall)

        handle_errors("setup")

    def test_aa_method_switching(self):
        # Virtual sites should be disabled by default
        self.assertIsInstance(self.system.virtual_sites, VirtualSitesOff)

        # Switch implementation
        self.system.virtual_sites = VirtualSitesInertialessTracers()
        self.assertIsInstance(
            self.system.virtual_sites, VirtualSitesInertialessTracers)

    def test_advection(self):
        self.reset_lb(ext_force_density=[0.1, 0, 0])
        # System setup
        system = self.system

        system.virtual_sites = VirtualSitesInertialessTracers()

        # Establish steady state flow field
        system.part.add(id=0, pos=(0, 5.5, 5.5), virtual=True)
        system.integrator.run(400)

        system.part[0].pos = (0, 5.5, 5.5)
        system.time = 0

        # Perform integration
        for _ in range(3):
            system.integrator.run(100)
            # compute expected position
            X = self.lbf.get_interpolated_velocity(
                system.part[0].pos)[0] * system.time
            self.assertAlmostEqual(
                system.part[0].pos[0] / X - 1, 0, delta=0.005)

    def compute_angle(self):
        system = self.system
        pos0 = system.part[0].pos
        pos1 = system.part[1].pos
        pos2 = system.part[2].pos
        pos3 = system.part[3].pos

        # first normal vector
        n1 = np.cross((pos1 - pos0), (pos2 - pos0))
        n2 = np.cross((pos2 - pos0), (pos3 - pos0))

        norm1 = np.linalg.norm(n1)
        norm2 = np.linalg.norm(n2)
        n1 = n1 / norm1
        n2 = n2 / norm2

        alpha = np.arccos(np.dot(n1, n2))
        return alpha

    def test_tribend(self):
        self.system.actors.clear()
        # two triangles with bending interaction
        # move nodes, should relax back

        system = self.system
        system.virtual_sites = VirtualSitesInertialessTracers()

        system.part.clear()

        # Add four particles
        system.part.add(id=0, pos=[5, 5, 5], virtual=True)
        system.part.add(id=1, pos=[5, 5, 6], virtual=True)
        system.part.add(id=2, pos=[5, 6, 6], virtual=True)
        system.part.add(id=3, pos=[5, 6, 5], virtual=True)

        # Add first triel, weak modulus
        from espressomd.interactions import IBM_Triel
        tri1 = IBM_Triel(
            ind1=0, ind2=1, ind3=2, elasticLaw="Skalak", k1=0.1, k2=0, maxDist=2.4)
        system.bonded_inter.add(tri1)
        system.part[0].add_bond((tri1, 1, 2))

        # Add second triel
        tri2 = IBM_Triel(
            ind1=0, ind2=2, ind3=3, elasticLaw="Skalak", k1=10, k2=0, maxDist=2.4)
        system.bonded_inter.add(tri2)
        system.part[0].add_bond((tri2, 2, 3))

        # Add bending
        from espressomd.interactions import IBM_Tribend
        tribend = IBM_Tribend(
            ind1=0, ind2=1, ind3=2, ind4=3, kb=1, refShape="Initial")
        system.bonded_inter.add(tribend)
        system.part[0].add_bond((tribend, 1, 2, 3))

        # twist
        system.part[1].pos = [5.2, 5, 6]

        self.reset_lb()

        # Perform integration
        last_angle = self.compute_angle()
        for _ in range(6):
            system.integrator.run(430)
            angle = self.compute_angle()
            self.assertLess(angle, last_angle)
            last_angle = angle
        self.assertLess(angle, 0.03)

    def test_triel(self):
        self.system.actors.clear()
        system = self.system
        system.virtual_sites = VirtualSitesInertialessTracers()
        system.virtual_sites = VirtualSitesInertialessTracers()

        system.part.clear()
        # Add particles: 0-2 are non-bonded, 3-5 are weakly bonded, 6-8 are
        # strongly bonded
        system.part.add(id=0, pos=[5, 5, 5], virtual=True)
        system.part.add(id=1, pos=[5, 5, 6], virtual=True)
        system.part.add(id=2, pos=[5, 6, 6], virtual=True)

        system.part.add(id=3, pos=[2, 5, 5], virtual=True)
        system.part.add(id=4, pos=[2, 5, 6], virtual=True)
        system.part.add(id=5, pos=[2, 6, 6], virtual=True)

        system.part.add(id=6, pos=[4, 7, 7], virtual=True)
        system.part.add(id=7, pos=[4, 7, 8], virtual=True)
        system.part.add(id=8, pos=[4, 8, 8], virtual=True)

        # Add triel, weak modulus for 3-5
        from espressomd.interactions import IBM_Triel
        triWeak = IBM_Triel(
            ind1=3, ind2=4, ind3=5, elasticLaw="Skalak", k1=5, k2=0, maxDist=2.4)
        system.bonded_inter.add(triWeak)
        system.part[3].add_bond((triWeak, 4, 5))

        # Add triel, strong modulus for 6-8
        triStrong = IBM_Triel(
            ind1=6, ind2=7, ind3=8, elasticLaw="Skalak", k1=15, k2=0, maxDist=2.4)
        system.bonded_inter.add(triStrong)
        system.part[6].add_bond((triStrong, 7, 8))

        self.reset_lb(ext_force_density=[0.1, 0, 0])
        # Perform integration
        system.integrator.run(4500)

        # For the cpu variant, check particle velocities
        if isinstance(self.lbf, lb.LBFluid):  # as opposed to LBFluidGPU
            for p in system.part:
                np.testing.assert_allclose(
                    np.copy(p.v), np.copy(
                        self.lbf.get_interpolated_velocity(p.pos)),
                    atol=2E-2)
        # get new shapes
        dist1non = np.linalg.norm(
            np.array(system.part[1].pos - system.part[0].pos))
        dist2non = np.linalg.norm(
            np.array(system.part[2].pos - system.part[0].pos))

        dist1weak = np.linalg.norm(
            np.array(system.part[3].pos - system.part[4].pos))
        dist2weak = np.linalg.norm(
            np.array(system.part[3].pos - system.part[5].pos))

        dist1strong = np.linalg.norm(
            np.array(system.part[6].pos - system.part[7].pos))
        dist2strong = np.linalg.norm(
            np.array(system.part[6].pos - system.part[8].pos))

        print("** Distances: non-bonded, weak, strong, expected")
        print(str(dist1non) + "    " + str(dist1weak)
              + "     " + str(dist1strong) + "    1")
        print(str(dist2non) + "    " + str(dist2weak)
              + "     " + str(dist2strong) + "    1.414")

        # test:
        # non-bonded should move apart by the flow (control group)
        # weakly-bonded should stretch somewhat
        # strongly-bonded should basically not stretch
        self.assertGreater(dist1non, 1.5)
        self.assertAlmostEqual(dist1weak, 1, delta=0.2)
        self.assertAlmostEqual(dist1strong, 1, delta=0.04)

        self.assertGreater(dist2non, 2)
        self.assertAlmostEqual(dist2weak, np.sqrt(2), delta=0.3)
        self.assertAlmostEqual(dist2strong, np.sqrt(2), delta=0.1)

    def test_zz_without_lb(self):
        """Check behaviour without lb. Ignore non-virtual particles, complain on
        virtual ones.

        """
        self.reset_lb()
        system = self.system
        system.virtual_sites = VirtualSitesInertialessTracers()
        system.actors.clear()
        system.part.clear()
        p = system.part.add(pos=(0, 0, 0))
        system.integrator.run(1)
        p.virtual = True
        with self.assertRaises(Exception):
            system.integrator.run(1)
