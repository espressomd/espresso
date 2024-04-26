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
import espressomd.particle_data
import espressomd.propagation
import numpy as np
import collections
import itertools


class ParticleProperties(ut.TestCase):

    # Particle id to work on
    pid = 17

    # Error tolerance when comparing arrays/tuples...
    tol = 1E-9

    # Handle for espresso system
    system = espressomd.System(box_l=[100.0, 100.0, 100.0])
    system.cell_system.skin = 0
    system.time_step = 0.01

    f1 = espressomd.interactions.FeneBond(k=1, d_r_max=2)
    system.bonded_inter.add(f1)
    f2 = espressomd.interactions.FeneBond(k=1, d_r_max=2)
    system.bonded_inter.add(f2)

    def setUp(self):
        self.partcl = self.system.part.add(id=self.pid, pos=(0, 0, 0))

    def tearDown(self):
        self.system.part.clear()

    def generateTestForVectorProperty(_propName, *_values):
        """Generates test cases for vectorial particle properties such as
        position, velocity...
        1st arg: name of the property (e.g., "pos"),
        2nd array: value to be used for testing. Has to be numpy.array of floats
        """
        # This is executed, when generateTestForVectorProperty() is called
        propName = _propName
        values = _values

        def func(self):
            # This code is run at the execution of the generated function.
            # It will use the state of the variables in the outer function,
            # which was there, when the outer function was called
            for value in values:
                setattr(self.partcl, propName, value)
                np.testing.assert_allclose(
                    np.array(getattr(self.partcl, propName)), value,
                    err_msg=propName + ": value set and value gotten back differ.",
                    atol=self.tol)

        return func

    def generateTestForScalarProperty(_propName, *_values):
        """Generates test cases for scalar particle properties such as
        type, mass, charge...
        1st arg: name of the property (e.g., "type"),
        2nd array: value to be used for testing. int or float
        """
        # This is executed, when generateTestForVectorProperty() is called
        propName = _propName
        values = _values

        def func(self):
            # This code is run at the execution of the generated function.
            # It will use the state of the variables in the outer function,
            # which was there, when the outer function was called
            for value in values:
                setattr(self.partcl, propName, value)
                self.assertEqual(getattr(self.partcl, propName),
                                 value, propName + ": value set and value gotten back differ.")

        return func

    test_pos = generateTestForVectorProperty("pos", np.array([0.1, 0.2, 0.3]))
    test_v = generateTestForVectorProperty("v", np.array([0.2, 0.3, 0.4]))
    test_f = generateTestForVectorProperty("f", np.array([0.2, 0.3, 0.7]))
    test_type = generateTestForScalarProperty("type", int(3))
    test_mol_id = generateTestForScalarProperty("mol_id", int(3))
    test_propagation = generateTestForScalarProperty(
        "propagation", int(espressomd.propagation.Propagation.TRANS_LANGEVIN))

    def test_invalid_propagation(self):
        Propagation = espressomd.propagation.Propagation
        with self.assertRaisesRegex(ValueError, "propagation combination not accepted"):
            self.partcl.propagation = (
                Propagation.TRANS_LANGEVIN | Propagation.TRANS_LANGEVIN_NPT | Propagation.TRANS_LB_MOMENTUM_EXCHANGE)

    def test_propagation_enum(self):
        Propagation = espressomd.propagation.Propagation
        flags_core = self.system.call_method("get_propagation_modes_enum")
        flags_si = {e.name: e.value for e in Propagation}
        # Python 3.11+ skips empty bitfields during enum iteration
        flags_si["NONE"] = Propagation.NONE
        self.assertEqual(flags_si, flags_core)

    test_bonds_property = generateTestForScalarProperty(
        "bonds", ((f1, 1), (f2, 2)))

    if espressomd.has_features(["MASS"]):
        test_mass = generateTestForScalarProperty("mass", 1.3)

    if espressomd.has_features(["ROTATION"]):
        test_rotation = generateTestForVectorProperty(
            "rotation", *itertools.product((True, False), repeat=3))
        test_omega_lab = generateTestForVectorProperty(
            "omega_lab", np.array([4., 2., 1.]))
        test_omega_body = generateTestForVectorProperty(
            "omega_body", np.array([4., 72., 1.]))
        test_torque_lab = generateTestForVectorProperty(
            "torque_lab", np.array([4., 72., 3.7]))
        # The tested value has to be normalized!
        test_quat = generateTestForVectorProperty(
            "quat", np.array([0.5, 0.5, 0.5, 0.5]))

        def test_director(self):
            """
            Test `director`. When set, it should get normalized.

            """
            sample_vector = np.array([0.5, -0.4, 1.3])
            sample_vector_normalized = sample_vector / \
                np.linalg.norm(sample_vector)

            setattr(self.partcl, "director", sample_vector)
            np.testing.assert_allclose(
                np.array(getattr(self.partcl, "director")),
                sample_vector_normalized
            )

        if espressomd.has_features(["THERMOSTAT_PER_PARTICLE"]):
            if espressomd.has_features(["PARTICLE_ANISOTROPY"]):
                test_gamma = generateTestForVectorProperty(
                    "gamma", np.array([2., 9., 0.23]))

                def test_gamma_single(self):
                    self.partcl.gamma = 17.4
                    np.testing.assert_array_equal(
                        np.copy(self.partcl.gamma),
                        np.array([17.4, 17.4, 17.4]),
                        "gamma: value set and value gotten back differ.")
            else:
                test_gamma = generateTestForScalarProperty("gamma", 17.3)

            if espressomd.has_features(["PARTICLE_ANISOTROPY"]):
                test_gamma_rot = generateTestForVectorProperty(
                    "gamma_rot", np.array([5., 10., 0.33]))

                def test_gamma_rot_single(self):
                    self.partcl.gamma_rot = 15.4
                    np.testing.assert_array_equal(
                        np.copy(self.partcl.gamma_rot),
                        np.array([15.4, 15.4, 15.4]),
                        "gamma_rot: value set and value gotten back differ.")
            else:
                test_gamma_rot = generateTestForScalarProperty(
                    "gamma_rot", 14.23)

    if espressomd.has_features(["ELECTROSTATICS"]):
        test_charge = generateTestForScalarProperty("q", -19.7)

    if espressomd.has_features(["EXTERNAL_FORCES"]):
        test_ext_force = generateTestForVectorProperty(
            "ext_force", [0.1, 0.2, 0.3])
        test_fix = generateTestForVectorProperty("fix", [True, False, True])

    if espressomd.has_features(["EXTERNAL_FORCES", "ROTATION"]):
        test_ext_torque = generateTestForVectorProperty(
            "ext_torque", [.4, .5, .6])

    if espressomd.has_features(["DIPOLES"]):
        test_dip = generateTestForVectorProperty(
            "dip", np.array([0.5, -0.5, 3]))
        test_dipm = generateTestForScalarProperty("dipm", -9.7)

    @utx.skipIfMissingFeatures(["VIRTUAL_SITES_RELATIVE"])
    def test_vs_relative(self):
        self.system.part.add(id=0, pos=(0, 0, 0))
        p1 = self.system.part.add(id=1, pos=(0, 0, 0))
        self.assertFalse(p1.is_virtual())
        p1.vs_relative = (0, 5.0, (0.5, -0.5, -0.5, -0.5))
        p1.vs_quat = [1, 2, 3, 4]
        self.assertTrue(p1.is_virtual())
        np.testing.assert_array_equal(p1.vs_quat, [1, 2, 3, 4])
        res = p1.vs_relative
        self.assertEqual(res[0], 0, f"vs_relative: {res}")
        self.assertEqual(res[1], 5.0, f"vs_relative: {res}")
        np.testing.assert_allclose(
            res[2], np.array((0.5, -0.5, -0.5, -0.5)),
            err_msg=f"vs_relative: {res}", atol=self.tol)
        # check exceptions
        error_msg = r"attribute 'vs_relative' of 'ParticleHandle' must take the form \[id, distance, quaternion\]"
        with self.assertRaisesRegex(ValueError, error_msg):
            p1.vs_relative = (0, 5.0)
        with self.assertRaisesRegex(ValueError, error_msg):
            p1.vs_relative = ('0', 5.0, (1, 0, 0, 0))
        with self.assertRaisesRegex(ValueError, error_msg):
            p1.vs_relative = (0, '5', (1, 0, 0, 0))
        with self.assertRaisesRegex(ValueError, error_msg):
            p1.vs_relative = (0, 5.0, (1, 0, 0))
        with self.assertRaisesRegex(ValueError, error_msg):
            p1.vs_relative = (0, 5.0, (0, 0, 0, 0))
        with self.assertRaisesRegex(ValueError, "attribute 'vs_quat' of 'ParticleHandle' must be non-zero"):
            p1.vs_quat = [0, 0, 0, 0]

    @utx.skipIfMissingFeatures(["EXCLUSIONS"])
    def test_invalid_exclusions(self):
        pid1 = self.pid
        pid2 = self.pid + 1

        p1 = self.partcl
        with self.assertRaisesRegex(RuntimeError, rf"Particles cannot exclude themselves \(id {self.pid}\)"):
            p1.add_exclusion(pid1)
        with self.assertRaisesRegex(RuntimeError, f"Particle with id {pid2} not found"):
            p1.add_exclusion(pid2)

        self.system.part.add(id=pid2, pos=(0, 0, 0))
        with self.assertRaisesRegex(RuntimeError, f"Particle with id {pid2} is not in exclusion list of particle with id {pid1}"):
            p1.delete_exclusion(pid2)
        with self.assertRaisesRegex(RuntimeError, f"Particle with id {pid2} is already in exclusion list of particle with id {pid1}"):
            p1.add_exclusion(pid2)
            p1.add_exclusion(pid2)

    @utx.skipIfMissingFeatures(["ROTATION"])
    def test_contradicting_properties_quat(self):
        invalid_combinations = [
            {'quat': [1., 1., 1., 1.], 'director': [1., 1., 1.]},
        ]
        if espressomd.has_features(["DIPOLES"]):
            invalid_combinations += [
                {'dip': [1., 1., 1.], 'dipm': 1.},
                {'dip': [1., 1., 1.], 'quat': [1., 1., 1., 1.]},
                {'dip': [1., 1., 1.], 'director': [1., 1., 1.]},
            ]
        # check all methods that can instantiate particles
        for kwargs in invalid_combinations:
            for make_new_particle in (
                    self.system.part.add, espressomd.particle_data.ParticleHandle):
                with self.assertRaises(ValueError):
                    make_new_particle(pos=[0., 0., 0.], **kwargs)

    @utx.skipIfMissingFeatures(["ROTATION"])
    def test_invalid_quat(self):
        system = self.system
        with self.assertRaisesRegex(ValueError, "attribute 'quat' of 'ParticleHandle' must be non-zero"):
            system.part.add(pos=[0., 0., 0.], quat=[0., 0., 0., 0.])

    @utx.skipIfMissingFeatures("ELECTROSTATICS")
    def test_particle_selection(self):
        self.system.part.clear()
        positions = ((0.2, 0.3, 0.4), (0.4, 0.2, 0.3), (0.7, 0.7, 0.7))
        charges = [0, 1E-6, -1, 1]

        # Place particles
        i = 0
        for pos in positions:
            for q in charges:
                self.system.part.add(pos=pos, q=q, id=i)
                i += 2

        # Scalar property
        res = self.system.part.select(q=0)
        self.assertEqual(len(res.id), len(positions))
        for p in res:
            self.assertAlmostEqual(p.q, 0, places=13)

        # Vectorial property
        res = self.system.part.select(pos=(0.2, 0.3, 0.4))
        self.assertEqual(len(res.id), len(charges))
        for p in res:
            np.testing.assert_allclose(
                (0.2, 0.3, 0.4), np.copy(p.pos), atol=1E-12)

        # Two criteria
        res = self.system.part.select(pos=(0.2, 0.3, 0.4), q=0)
        self.assertEqual(tuple(res.id), (0,))

        # Empty result
        res = self.system.part.select(q=17)
        self.assertEqual(tuple(res.id), ())
        # User-specified criterion
        res = self.system.part.select(lambda p: p.pos[0] < 0.5)
        np.testing.assert_equal(sorted(res.id), np.arange(0, 16, 2, dtype=int))

    def test_image_box(self):
        self.system.part.clear()
        pos = 1.5 * self.system.box_l
        p = self.system.part.add(pos=pos)
        np.testing.assert_equal(np.copy(p.image_box), [1, 1, 1])

    def test_particle_numbering(self):
        """
        Instantiate particles on different nodes and check they are
        numbered sequentially.
        """
        system = self.system
        offset = system.box_l / 100.
        # clear the cached value for largest particle id
        system.part.clear()
        # update the cached value to 20
        p_id_start = self.pid + 3
        system.part.add(pos=offset, id=p_id_start)
        # check the cached value is updated regardless of the MPI node
        # where particles are inserted
        for i in range(3):
            pos = [0., 0., 0.]
            pos[i] = -offset[i]  # guaranteed to hit different MPI domains
            p = self.system.part.add(pos=pos)
            self.assertEqual(p.id, p_id_start + i + 1)
        # removing the particle with highest id should update the cached value
        p.remove()
        self.assertEqual(system.part.highest_particle_id, p_id_start + 2)

        # update the cached value to 32
        p_id_start += 10
        p0 = system.part.add(pos=offset, id=p_id_start + 0)
        p1 = system.part.add(pos=offset, id=p_id_start + 1)
        p2 = system.part.add(pos=offset, id=p_id_start + 2)
        # invalidate the cache by introducing a gap at 31 and then removing 32
        p1.remove()
        p2.remove()
        # the cache should now be 30
        self.assertEqual(system.part.highest_particle_id, p0.id)
        p3 = system.part.add(pos=offset)
        self.assertEqual(p3.id, p0.id + 1)

    def test_invalid_particle_ids(self):
        with self.assertRaisesRegex(ValueError, f"Particle {self.partcl.id} already exists"):
            self.system.part.add(id=self.partcl.id, pos=self.partcl.pos)
        self.system.part.clear()
        handle_to_non_existing_particle = self.system.part.by_id(42)
        with self.assertRaisesRegex(RuntimeError, "Particle node for id 42 not found"):
            handle_to_non_existing_particle.type
        for i in range(1, 10):
            p = espressomd.particle_data.ParticleHandle(id=-i)
            with self.assertRaisesRegex(ValueError, f"Invalid particle id: {-i}"):
                p.node
            with self.assertRaisesRegex(ValueError, f"Invalid particle id: {-i}"):
                p.remove()
            with self.assertRaisesRegex(ValueError, f"Invalid particle id: {-i}"):
                self.system.part.add(pos=[0., 0., 0.], id=-i)

    def test_invalid_particle_creation(self):
        err_msg = r"add\(\) takes either a dictionary or a bunch of keyword args"
        with self.assertRaisesRegex(ValueError, err_msg):
            self.system.part.add(1)
        with self.assertRaisesRegex(ValueError, err_msg):
            self.system.part.add(1, 2)
        with self.assertRaisesRegex(ValueError, err_msg):
            self.system.part.add(1, id=2)
        with self.assertRaisesRegex(ValueError, err_msg):
            self.system.part.add([1, 2])
        for coord in (0, 1, 2):
            for value in (float('nan'), float('inf'), -float('inf')):
                with self.assertRaisesRegex(ValueError, "Particle position must be finite"):
                    pos = [0., 0., 0.]
                    pos[coord] = value
                    self.system.part.add(pos=pos)

    def test_invalid_particle_attributes(self):
        p = self.partcl
        err_msg = "attribute '{}' of 'ParticleHandle' {}"
        with self.assertRaisesRegex(ValueError, err_msg.format("type", "must be an integer >= 0")):
            p.type = -1
        with self.assertRaisesRegex(ValueError, err_msg.format("mol_id", "must be an integer >= 0")):
            p.mol_id = -1
        if espressomd.has_features("MASS"):
            for mass in [0., -1., -2.]:
                with self.assertRaisesRegex(ValueError, err_msg.format("mass", "must be a float > 0")):
                    p.mass = mass

    def test_missing_features(self):
        def check(feature, prop, throwing_values, valid_value=None):
            if not espressomd.has_features(feature):
                if valid_value is not None:
                    # this should not throw
                    setattr(self.partcl, prop, valid_value)
                for throwing_value in throwing_values:
                    with self.assertRaisesRegex(RuntimeError, f"Feature {feature} not compiled in"):
                        setattr(self.partcl, prop, throwing_value)

        check("MASS", "mass", [1.1, 0., -1.], 1.)
        check("ELECTROSTATICS", "q", [1., -1.], 0.)

    def test_parallel_property_setters(self):
        system = self.system
        system.part.clear()
        partcls = system.part.add(
            pos=system.box_l * np.random.random((100, 3)))

        # Copy individual properties of particle 0
        print(
            "If this test hangs, there is an mpi deadlock in a particle property setter.")
        for attr in espressomd.particle_data.particle_attributes:
            # Uncomment to identify guilty property
            # print(attr)

            assert hasattr(system.part.by_id(0), attr), \
                "Inconsistency between ParticleHandle and particle_data.particle_attributes"
            try:
                setattr(partcls, attr, getattr(system.part.by_id(0), attr))
            except RuntimeError as err:
                self.assertEqual(str(err), f"Parameter '{attr}' is read-only.")
            # Cause a different mpi callback to uncover deadlock immediately
            _ = getattr(partcls, attr)

    def test_remove_particle(self):
        """Tests that if a particle is removed,
        it no longer exists and bonds to the removed particle are
        also removed."""

        p1 = self.system.part.by_id(self.pid)
        p2 = self.system.part.add(pos=p1.pos, bonds=[(self.f1, p1.id)])

        p1.remove()
        self.assertFalse(self.system.part.exists(self.pid))
        self.assertEqual(len(p2.bonds), 0)

    def test_bonds(self):
        """Tests bond addition and removal."""

        p1 = self.system.part.by_id(self.pid)
        p2 = self.system.part.add(pos=p1.pos)
        inactive_bond = espressomd.interactions.FeneBond(k=1, d_r_max=2)
        p2.add_bond([self.f1, p1])
        with self.assertRaisesRegex(RuntimeError, "already exists on particle"):
            p2.add_bond([self.f1, p1.id])
        with self.assertRaisesRegex(RuntimeError, "already exists on particle"):
            p2.add_bond((self.f1, p1))
        with self.assertRaisesRegex(Exception, "1st element of Bond has to be of type BondedInteraction or int"):
            p2.add_bond(('self.f1', p1))
        with self.assertRaisesRegex(ValueError, "Bond partners have to be of type integer or ParticleHandle"):
            p2.add_bond((self.f1, '1'))
        with self.assertRaisesRegex(ValueError, r"Bond FeneBond\(.+?\) needs 1 partner"):
            p2.add_bond((self.f1, p1, p2))
        with self.assertRaisesRegex(Exception, "The bonded interaction has not yet been added to the list of active bonds in ESPResSo"):
            p2.add_bond((inactive_bond, p1))
        p2.delete_bond([self.f1, p1])
        with self.assertRaisesRegex(RuntimeError, "doesn't exist on particle"):
            p2.delete_bond([self.f1, p1])
        with self.assertRaisesRegex(ValueError, "Bond partners have to be of type integer or ParticleHandle"):
            p2.delete_bond((self.f1, 'p1'))

    def test_zz_remove_all(self):
        for p in self.system.part.all():
            p.remove()
        self.system.part.add(
            pos=np.random.random((100, 3)) * self.system.box_l,
            id=np.arange(100, dtype=int))
        ids = self.system.part.all().id
        np.random.shuffle(ids)
        for pid in ids:
            self.system.part.by_id(pid).remove()
        p = self.system.part.by_id(17)
        with self.assertRaises(RuntimeError):
            p.remove()

    def test_coord_fold_corner_cases(self):
        system = self.system
        system.time_step = .5
        system.cell_system.set_regular_decomposition(use_verlet_lists=False)
        system.cell_system.skin = 0
        system.min_global_cut = 3
        system.part.clear()
        p1 = system.part.add(
            pos=3 * [np.nextafter(0., -1.)], v=system.box_l / 3)
        print(p1.pos)
        p2 = system.part.add(
            pos=np.nextafter(system.box_l, 2 * system.box_l), v=system.box_l / 3)
        print(p2.pos)
        p3 = system.part.add(
            pos=np.nextafter(system.box_l, (0, 0, 0)), v=system.box_l / 3)
        print(p3.pos)
        p4 = system.part.add(
            pos=3 * [np.nextafter(0., 1.)], v=system.box_l / 3)
        print(p4.pos)
        system.integrator.run(3)
        for p in system.part:
            for i in range(3):
                self.assertGreaterEqual(p.pos_folded[i], 0)
                self.assertLess(p.pos_folded[i], system.box_l[i])

        # Force resort
        system.part.add(pos=(0, 0, 0))
        system.integrator.run(9)
        for p in system.part:
            for i in range(3):
                self.assertGreaterEqual(p.pos_folded[i], 0)
                self.assertLess(p.pos_folded[i], system.box_l[i])

    def test_particle_list(self):
        self.assertEqual(str(self.system.part), "ParticleList([17])")

    def test_particle_slice(self):
        """Tests operations on slices of particles"""

        system = self.system

        # Empty slice
        system.part.clear()
        all_partcls_empty = system.part.all()
        self.assertEqual(len(all_partcls_empty), 0)
        self.assertEqual(len(all_partcls_empty), 0)
        self.assertEqual(len(all_partcls_empty), 0)
        with self.assertRaises(AttributeError):
            all_partcls_empty.pos = ((1, 2, 3,),)

        # Slice containing particles
        ids = [1, 4, 6, 3, 8, 9, 5]
        pos = np.random.random((len(ids), 3))
        system.part.add(id=ids, pos=pos)

        # All particles
        all_partcls = system.part.all()
        self.assertEqual(len(all_partcls), len(ids))
        np.testing.assert_equal(all_partcls.id, sorted(ids))
        np.testing.assert_equal(all_partcls.pos, pos[np.argsort(ids)])

        # Access via slicing
        np.testing.assert_equal(system.part.by_ids(range(3, 6)).id,
                                [i for i in sorted(ids) if i >= 3 and i < 6])
        np.testing.assert_equal(system.part.by_ids(range(6, 3, -1)).id,
                                [i for i in sorted(ids, key=lambda i: -i) if i > 3 and i <= 6])

        # Setting particle properties on a slice
        system.part.by_ids(range(9, 10)).pos = (0, 0, 0)
        np.testing.assert_equal(system.part.all().pos,
                                [pos[i] if ids[i] != 9 else [0, 0, 0] for i in np.argsort(ids)])

        # Slice access via explicit list of ids
        np.testing.assert_equal(system.part.by_ids(ids[1:4]).id, ids[1:4])
        # Check that ids passed in an explicit list must exist
        with self.assertRaisesRegex(IndexError, "Particle does not exist: 99"):
            system.part.by_ids([99, 3])
        # Check that wrong types are not accepted
        with self.assertRaisesRegex(RuntimeError, "it contains a value that is not convertible to 'int'"):
            system.part.by_ids([ids[0], 1.2])

    def test_to_dict(self):
        self.system.part.clear()
        p = self.system.part.add(
            pos=np.random.uniform(size=(10, 3)) * self.system.box_l)
        pp = str(p)
        pdict = p.to_dict()
        p.remove()
        p = self.system.part.add(pdict)
        self.assertEqual(str(self.system.part.select()), pp)
        p.remove()
        p = self.system.part.add(collections.OrderedDict(pdict))
        self.assertEqual(str(self.system.part.select()), pp)
        p.remove()

    def test_update(self):
        self.system.part.clear()
        p = self.system.part.add(pos=0.5 * self.system.box_l)
        # cannot change id (to avoid corrupting caches in the core)
        with self.assertRaisesRegex(RuntimeError, "Cannot change particle id"):
            p.update({'id': 1})
        with self.assertRaisesRegex(RuntimeError, "Cannot change particle id"):
            self.system.part.all().update({'id': 1})
        # check value change
        new_pos = [1., 2., 3.]
        p.update({'pos': new_pos})
        np.testing.assert_almost_equal(p.pos, new_pos)
        # updating self should not change anything
        pdict = p.to_dict()
        del pdict['id']
        p.update(pdict)
        new_pdict = p.to_dict()
        del new_pdict['id']
        self.assertEqual(str(new_pdict), str(pdict))
        with self.assertRaisesRegex(RuntimeError, "Parameter 'test' is missing"):
            p.call_method("set_param_parallel", name="test")


if __name__ == "__main__":
    ut.main()
