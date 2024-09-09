#
# Copyright (C) 2013-2022 The ESPResSo project
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
import espressomd.interactions
import numpy as np


@utx.skipIfMissingFeatures("COLLISION_DETECTION")
class CollisionDetection(ut.TestCase):

    """Tests collision detection interface and exceptions."""

    system = espressomd.System(box_l=[10.0, 10.0, 10.0])
    bond_angle_resolution = 3
    for i in range(bond_angle_resolution):
        system.bonded_inter.add(espressomd.interactions.AngleHarmonic(
            bend=1., phi0=float(i) / float(bond_angle_resolution - 1) * np.pi))
    bond_harmonic = espressomd.interactions.HarmonicBond(k=5., r_0=0.1)
    bond_angle = espressomd.interactions.AngleHarmonic(bend=1., phi0=np.pi)
    bond_dihe = espressomd.interactions.Dihedral(bend=1., mult=1, phase=0.)
    system.bonded_inter.add(bond_harmonic)
    system.bonded_inter.add(bond_angle)
    system.bonded_inter.add(bond_dihe)
    system.time_step = 0.01
    valid_coldet_params = {
        "bind_centers": {"distance": 0.1, "bond_centers": bond_harmonic},
        "bind_at_point_of_collision": {
            "bond_vs": bond_harmonic, "bond_centers": bond_harmonic,
            "part_type_vs": 1, "distance": 0.1, "vs_placement": 0.1
        },
        "glue_to_surface": {
            "distance": 0.1, "distance_glued_particle_to_vs": 0.02,
            "bond_centers": bond_harmonic, "bond_vs": bond_harmonic,
            "part_type_vs": 1, "part_type_to_attach_vs_to": 0,
            "part_type_to_be_glued": 2, "part_type_after_glueing": 3
        },
    }
    classes = {"bind_centers": espressomd.collision_detection.BindCenters,
               "bind_at_point_of_collision": espressomd.collision_detection.BindAtPointOfCollision,
               "glue_to_surface": espressomd.collision_detection.GlueToSurface}

    def tearDown(self):
        self.system.collision_detection.protocol = None

    def test_00_interface_and_defaults(self):
        self.assertIsNone(self.system.collision_detection.protocol)
        self.assertIsNone(self.system.collision_detection.call_method("none"))

    def check_stored_parameters(self, mode, **kwargs):
        """
        Check if collision detection stored parameters match input values.
        """
        parameters = self.system.collision_detection.protocol.get_params()
        parameters_ref = self.valid_coldet_params[mode].copy()
        parameters_ref.update(kwargs)
        for key, value_ref in parameters_ref.items():
            if isinstance(value_ref, float):
                self.assertAlmostEqual(parameters[key], value_ref, delta=1e-10)
            else:
                self.assertEqual(parameters[key], value_ref)

    def set_coldet(self, mode, **invalid_params):
        """
        Instantiate collision detection with one or more incorrect parameters.
        """
        params = self.valid_coldet_params.get(mode).copy()
        params.update(invalid_params)
        self.system.collision_detection.protocol = self.classes[mode](**params)

    def test_off(self):
        self.system.collision_detection.protocol = None
        with self.assertRaisesRegex(ValueError, "Parameter 'distance' must be > 0"):
            self.set_coldet("bind_centers", distance=-2.)
        # check if original parameters have been preserved
        self.assertIsNone(self.system.collision_detection.protocol)

    def test_bind_centers(self):
        self.set_coldet("bind_centers", distance=0.5)
        with self.assertRaisesRegex(ValueError, "Parameter 'distance' must be > 0"):
            self.set_coldet("bind_centers", distance=-2.)
        with self.assertRaisesRegex(ValueError, "Parameter 'distance' must be > 0"):
            self.set_coldet("bind_centers", distance=0.)
        with self.assertRaisesRegex(ValueError, "Bond in parameter list was not added to the system"):
            bond = espressomd.interactions.HarmonicBond(k=1., r_0=0.1)
            self.set_coldet("bind_centers", bond_centers=bond)
        with self.assertRaisesRegex(RuntimeError, "The bond type to be used for binding particle centers needs to be a pair bond"):
            self.set_coldet("bind_centers", bond_centers=self.bond_angle)
        with self.assertRaisesRegex(RuntimeError, "Parameter 'distance' is missing"):
            self.system.collision_detection.protocol = espressomd.collision_detection.BindCenters(
                bond_centers=self.bond_harmonic)
        # check if original parameters have been preserved
        self.check_stored_parameters("bind_centers", distance=0.5)

    @utx.skipIfMissingFeatures("VIRTUAL_SITES_RELATIVE")
    def test_bind_at_point_of_collision(self):
        self.set_coldet("bind_at_point_of_collision", distance=0.5)
        with self.assertRaisesRegex(ValueError, "Parameter 'distance' must be > 0"):
            self.set_coldet("bind_at_point_of_collision", distance=-2.)
        with self.assertRaisesRegex(ValueError, "Parameter 'distance' must be > 0"):
            self.set_coldet("bind_at_point_of_collision", distance=0.)
        with self.assertRaisesRegex(ValueError, "Parameter 'vs_placement' must be between 0 and 1"):
            self.set_coldet("bind_at_point_of_collision", vs_placement=-0.01)
        with self.assertRaisesRegex(ValueError, "Parameter 'vs_placement' must be between 0 and 1"):
            self.set_coldet("bind_at_point_of_collision", vs_placement=1.01)
        with self.assertRaisesRegex(ValueError, "Bond in parameter list was not added to the system"):
            bond = espressomd.interactions.HarmonicBond(k=1., r_0=0.1)
            self.set_coldet("bind_at_point_of_collision", bond_vs=bond)
        with self.assertRaisesRegex(RuntimeError, "bond type to be used for binding virtual sites needs to be a pair bond"):
            self.set_coldet(
                "bind_at_point_of_collision", bond_vs=self.bond_dihe)
        with self.assertRaisesRegex(ValueError, "type for virtual sites needs to be >=0"):
            self.set_coldet("bind_at_point_of_collision", part_type_vs=-1)
        # check if original parameters have been preserved
        self.check_stored_parameters(
            "bind_at_point_of_collision", distance=0.5)

    @utx.skipIfMissingFeatures("VIRTUAL_SITES")
    def test_bind_at_point_of_collision_norotation(self):
        if not espressomd.has_features("VIRTUAL_SITES_RELATIVE"):
            with self.assertRaisesRegex(RuntimeError, "Missing features VIRTUAL_SITES_RELATIVE"):
                self.set_coldet("bind_at_point_of_collision")

    @utx.skipIfMissingFeatures("VIRTUAL_SITES_RELATIVE")
    def test_glue_to_surface(self):
        self.set_coldet("glue_to_surface", distance=0.5)
        with self.assertRaisesRegex(ValueError, "Parameter 'distance' must be > 0"):
            self.set_coldet("glue_to_surface", distance=-2.)
        with self.assertRaisesRegex(ValueError, "Parameter 'distance' must be > 0"):
            self.set_coldet("glue_to_surface", distance=0.)
        with self.assertRaisesRegex(ValueError, "type for virtual sites needs to be >=0"):
            self.set_coldet("glue_to_surface", part_type_vs=-1)
        with self.assertRaisesRegex(ValueError, "type to be glued needs to be >=0"):
            self.set_coldet("glue_to_surface", part_type_to_be_glued=-1)
        with self.assertRaisesRegex(ValueError, "type to attach the virtual site to needs to be >=0"):
            self.set_coldet("glue_to_surface", part_type_to_attach_vs_to=-1)
        with self.assertRaisesRegex(ValueError, "type after gluing needs to be >=0"):
            self.set_coldet("glue_to_surface", part_type_after_glueing=-1)
        # check if original parameters have been preserved
        self.check_stored_parameters("glue_to_surface", distance=0.5)


if __name__ == "__main__":
    ut.main()
