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

import unittest as ut
import unittest_decorators as utx
import unittest.mock
import espressomd
import numpy as np
import contextlib

Calculator = unittest.mock.MagicMock()
with contextlib.suppress(ImportError):
    import espressomd.plugins.ase
    from ase.calculators.calculator import Calculator

# Global system instance shared between test classes
system = espressomd.System(box_l=[20., 20., 20.])


@utx.skipIfMissingFeatures(["EXTERNAL_FORCES", "MASS", "ELECTROSTATICS"])
@utx.skipIfMissingModules("ase.calculators")
class ASEInterfaceTest(ut.TestCase):
    """test suite for the ASE interface focusing on update_ase() method."""

    def setUp(self):
        """Set up system with particles having various properties."""
        system.time_step = 0.01
        system.box_l = [20., 20., 20.]  # Ensure correct box size

        # Add particles with positions, charges, masses, velocities
        system.part.add(
            pos=[1., 2., 3.], q=1.0, mass=2.0, v=[0.1, 0.2, 0.3], type=0)
        system.part.add(
            pos=[4., 5., 6.], q=-0.5, mass=1.5, v=[0.4, 0.5, 0.6], type=1)
        system.part.add(
            pos=[7., 8., 9.], q=2.0, mass=3.0, v=[0.7, 0.8, 0.9], type=0)

    def tearDown(self):
        """Clean up after each test."""
        system.part.clear()

    def test_ase_interface_instantiation_basic(self):
        """Test basic ASE interface creation without optional exports."""
        ase_interface = espressomd.plugins.ase.ASEInterface(
            system=system,
            particle_slice=system.part.all()
        )

        # Check that interface is properly initialized
        self.assertIsNotNone(ase_interface.atoms)
        self.assertEqual(ase_interface.export_charges, False)
        self.assertEqual(ase_interface.export_masses, False)
        self.assertEqual(ase_interface.export_momenta, False)
        self.assertEqual(ase_interface.export_forces, False)

        # Verify atoms object has correct basic properties
        atoms = ase_interface.atoms
        self.assertEqual(len(atoms), 3)
        np.testing.assert_array_equal(atoms.numbers, [0, 1, 0])
        np.testing.assert_allclose(
            atoms.positions, [[1., 2., 3.], [4., 5., 6.], [7., 8., 9.]])

    def test_ase_interface_instantiation_with_charges(self):
        """Test ASE interface creation with charge export enabled."""
        ase_interface = espressomd.plugins.ase.ASEInterface(
            system=system,
            particle_slice=system.part.all(),
            export_charges=True
        )

        self.assertEqual(ase_interface.export_charges, True)
        # Verify charges are set in ASE atoms
        atoms = ase_interface.atoms
        np.testing.assert_allclose(
            atoms.get_initial_charges(), [1.0, -0.5, 2.0])

    def test_ase_interface_instantiation_with_masses(self):
        """Test ASE interface creation with mass export enabled."""
        ase_interface = espressomd.plugins.ase.ASEInterface(
            system=system,
            particle_slice=system.part.all(),
            export_masses=True
        )

        self.assertEqual(ase_interface.export_masses, True)
        # Verify masses are set in ASE atoms
        atoms = ase_interface.atoms
        np.testing.assert_allclose(atoms.get_masses(), [2.0, 1.5, 3.0])

    def test_ase_interface_instantiation_with_momenta(self):
        """Test ASE interface creation with momentum export enabled."""
        ase_interface = espressomd.plugins.ase.ASEInterface(
            system=system,
            particle_slice=system.part.all(),
            export_momenta=True
        )

        self.assertEqual(ase_interface.export_momenta, True)
        # Verify momenta are set in ASE atoms (momentum = mass * velocity)
        atoms = ase_interface.atoms
        expected_momenta = np.array([[0.1 * 2.0, 0.2 * 2.0, 0.3 * 2.0],
                                    [0.4 * 1.5, 0.5 * 1.5, 0.6 * 1.5],
                                    [0.7 * 3.0, 0.8 * 3.0, 0.9 * 3.0]])
        np.testing.assert_allclose(atoms.get_momenta(), expected_momenta)

    def test_ase_interface_instantiation_with_forces(self):
        """Test ASE interface creation with force export enabled."""
        # Set some forces on ESPResSo particles
        system.part.by_id(0).f = [0.1, 0.2, 0.3]
        system.part.by_id(1).f = [0.4, 0.5, 0.6]
        system.part.by_id(2).f = [0.7, 0.8, 0.9]

        ase_interface = espressomd.plugins.ase.ASEInterface(
            system=system,
            particle_slice=system.part.all(),
            export_forces=True
        )

        self.assertEqual(ase_interface.export_forces, True)
        # Verify forces are set in ASE atoms
        atoms = ase_interface.atoms
        expected_forces = np.array([[0.1, 0.2, 0.3],
                                    [0.4, 0.5, 0.6],
                                    [0.7, 0.8, 0.9]])
        np.testing.assert_allclose(atoms.arrays["forces"], expected_forces)

    def test_ase_interface_instantiation_all_exports(self):
        """Test ASE interface creation with all exports enabled."""
        ase_interface = espressomd.plugins.ase.ASEInterface(
            system=system,
            particle_slice=system.part.all(),
            export_charges=True,
            export_masses=True,
            export_momenta=True,
            export_forces=True
        )

        # Check all export flags
        self.assertEqual(ase_interface.export_charges, True)
        self.assertEqual(ase_interface.export_masses, True)
        self.assertEqual(ase_interface.export_momenta, True)
        self.assertEqual(ase_interface.export_forces, True)

        # Verify all properties are set
        atoms = ase_interface.atoms
        np.testing.assert_allclose(
            atoms.get_initial_charges(), [1.0, -0.5, 2.0])
        np.testing.assert_allclose(atoms.get_masses(), [2.0, 1.5, 3.0])
        expected_momenta = np.array(
            [[0.2, 0.4, 0.6], [0.6, 0.75, 0.9], [2.1, 2.4, 2.7]])
        np.testing.assert_allclose(atoms.get_momenta(), expected_momenta)
        expected_forces = np.array([[0., 0., 0.],
                                    [0., 0., 0.],
                                    [0., 0., 0.]])
        np.testing.assert_allclose(atoms.arrays["forces"], expected_forces)

    def test_update_ase_basic_no_exports(self):
        """Test update_ase() with no exports enabled - positions and types should update."""
        ase_interface = espressomd.plugins.ase.ASEInterface(
            system=system,
            particle_slice=system.part.all()
        )

        # Change particle positions
        system.part.by_id(0).pos = [10., 11., 12.]
        system.part.by_id(1).pos = [13., 14., 15.]
        system.part.by_id(2).pos = [16., 17., 18.]

        # Change particle types
        system.part.by_id(0).type = 1
        system.part.by_id(1).type = 0
        system.part.by_id(2).type = 1

        # Update ASE
        ase_interface.update_ase()

        # Check that positions are updated
        atoms = ase_interface.atoms
        expected_positions = system.part.all().pos
        np.testing.assert_allclose(atoms.positions, expected_positions)
        # Check that types are updated
        expected_types = system.part.all().type
        np.testing.assert_array_equal(atoms.numbers, expected_types)

    def test_update_ase_charges_not_constant(self):
        """Test update_ase() with charges enabled and not assumed constant."""
        ase_interface = espressomd.plugins.ase.ASEInterface(
            system=system,
            particle_slice=system.part.all(),
            export_charges=True,
            assume_constant_charges=False
        )

        # Change particle charges
        system.part.all().q = [5., -2., 3.5]

        # Update ASE (charges should be updated)
        ase_interface.update_ase()

        # Check that charges are updated
        atoms = ase_interface.atoms
        np.testing.assert_allclose(
            atoms.get_initial_charges(), np.copy(system.part.all().q))

    def test_update_ase_charges_assumed_constant(self):
        """Test update_ase() with charges enabled but assumed constant."""
        ase_interface = espressomd.plugins.ase.ASEInterface(
            system=system,
            particle_slice=system.part.all(),
            export_charges=True,
            assume_constant_charges=True
        )

        # Store original charges
        original_charges = ase_interface.atoms.get_initial_charges().copy()

        # Change particle charges
        system.part.all().q = [5., -2., 3.5]
        # Update ASE (charges should NOT be updated due to
        # assume_constant_charges)
        ase_interface.update_ase()

        # Check that charges are NOT updated (should still be original)
        atoms = ase_interface.atoms
        np.testing.assert_allclose(
            atoms.get_initial_charges(),
            original_charges)

    def test_update_ase_masses_not_constant(self):
        """Test update_ase() with masses enabled and not assumed constant."""
        ase_interface = espressomd.plugins.ase.ASEInterface(
            system=system,
            particle_slice=system.part.all(),
            export_masses=True,
            assume_constant_masses=False
        )

        # Change particle masses
        system.part.all().mass = [10.0, 5.5, 7.2]

        # Update ASE (masses should be updated)
        ase_interface.update_ase()

        # Check that masses are updated
        atoms = ase_interface.atoms
        np.testing.assert_allclose(
            atoms.get_masses(), np.copy(system.part.all().mass))

    def test_update_ase_masses_assumed_constant(self):
        """Test update_ase() with masses enabled but assumed constant."""
        ase_interface = espressomd.plugins.ase.ASEInterface(
            system=system,
            particle_slice=system.part.all(),
            export_masses=True,
            assume_constant_masses=True
        )

        # Store original masses
        original_masses = ase_interface.atoms.get_masses().copy()

        # Change particle masses
        system.part.by_id(2).mass = 7.2

        # Update ASE (masses should NOT be updated due to
        # assume_constant_masses)
        ase_interface.update_ase()

        # Check that masses are NOT updated (should still be original)
        atoms = ase_interface.atoms
        np.testing.assert_allclose(atoms.get_masses(), original_masses)

    def test_update_ase_momenta_always_updates(self):
        """Test update_ase() with momenta enabled - should always update (no skip option)."""
        ase_interface = espressomd.plugins.ase.ASEInterface(
            system=system,
            particle_slice=system.part.all(),
            export_momenta=True
        )

        # Change particle velocities
        system.part.by_id(0).v = [1.0, 1.1, 1.2]
        system.part.by_id(1).v = [2.0, 2.1, 2.2]
        system.part.by_id(2).v = [3.0, 3.1, 3.2]

        # Update ASE
        ase_interface.update_ase()

        # Check that momenta are updated (momentum = mass * velocity)
        atoms = ase_interface.atoms
        expected_momenta = np.array([[1.0 * 2.0, 1.1 * 2.0, 1.2 * 2.0],
                                    [2.0 * 1.5, 2.1 * 1.5, 2.2 * 1.5],
                                    [3.0 * 3.0, 3.1 * 3.0, 3.2 * 3.0]])
        np.testing.assert_allclose(atoms.get_momenta(), expected_momenta)

    def test_update_ase_types_assumed_constant(self):
        """Test update_ase() with types assumed constant."""
        ase_interface = espressomd.plugins.ase.ASEInterface(
            system=system,
            particle_slice=system.part.all(),
            assume_constant_types=True
        )

        # Store original types
        original_types = ase_interface.atoms.numbers.copy()

        # Change particle types
        system.part.by_id(0).type = 1
        system.part.by_id(1).type = 0
        system.part.by_id(2).type = 1

        # Update ASE (types should NOT be updated due to
        # assume_constant_types)
        ase_interface.update_ase()

        # Check that types are NOT updated (should still be original)
        atoms = ase_interface.atoms
        np.testing.assert_array_equal(atoms.numbers, original_types)

    def test_update_ase_forces(self):
        """Test update_ase() with forces enabled."""
        # Set some forces on ESPResSo particles
        ase_interface = espressomd.plugins.ase.ASEInterface(
            system=system,
            particle_slice=system.part.all(),
            export_forces=True
        )

        # Change particle forces
        system.part.by_id(0).f = [1.1, 1.2, 1.3]
        system.part.by_id(1).f = [2.1, 2.2, 2.3]
        system.part.by_id(2).f = [3.1, 3.2, 3.3]

        # Update ASE
        ase_interface.update_ase()

        # Check that forces are updated in ASE atoms
        atoms = ase_interface.atoms
        expected_forces = np.array([[1.1, 1.2, 1.3],
                                    [2.1, 2.2, 2.3],
                                    [3.1, 3.2, 3.3]])
        np.testing.assert_allclose(atoms.arrays["forces"], expected_forces)

    def test_error_handling_no_atoms(self):
        """Test that update_ase() raises error when atoms is None."""
        ase_interface = espressomd.plugins.ase.ASEInterface(
            system=system,
            particle_slice=system.part.all()
        )

        # Set atoms to None to simulate uninitialized state
        ase_interface.atoms = None

        # Should raise RuntimeError
        with self.assertRaisesRegex(RuntimeError, "atoms object not initialized"):
            ase_interface.update_ase()

    def test_folded_positions_default(self):
        """Test that folded positions are used by default (when not specified)."""
        # Place a particle outside the box
        system.part.by_id(0).pos = [25., 25., 25.]  # Box is [20, 20, 20]

        # Create interface WITHOUT specifying use_folded_positions (should
        # default to True)
        ase_interface = espressomd.plugins.ase.ASEInterface(
            system=system,
            particle_slice=system.part.all()
        )

        # Verify the default value is True
        self.assertEqual(ase_interface.use_folded_positions, True)

        # Get folded position from ESPResSo
        folded_pos = np.copy(system.part.by_id(0).pos_folded)

        # Check that ASE atoms use folded position
        atoms = ase_interface.atoms
        np.testing.assert_allclose(atoms.positions[0], folded_pos)

        # Folded position should be inside box
        self.assertTrue(np.all(folded_pos >= 0))
        self.assertTrue(np.all(folded_pos < system.box_l))

    def test_folded_positions_true(self):
        """Test that folded positions are used when use_folded_positions=True."""
        # Place a particle outside the box
        system.part.by_id(0).pos = [25., 25., 25.]  # Box is [20, 20, 20]

        # Create interface with folded positions explicitly set
        ase_interface = espressomd.plugins.ase.ASEInterface(
            system=system,
            particle_slice=system.part.all(),
            use_folded_positions=True
        )

        # Get folded position from ESPResSo
        folded_pos = np.copy(system.part.by_id(0).pos_folded)

        # Check that ASE atoms use folded position
        atoms = ase_interface.atoms
        np.testing.assert_allclose(atoms.positions[0], folded_pos)

        # Folded position should be inside box
        self.assertTrue(np.all(folded_pos >= 0))
        self.assertTrue(np.all(folded_pos < system.box_l))

    def test_folded_positions_false(self):
        """Test that unfolded positions are used when use_folded_positions=False."""
        # Place a particle outside the box
        system.part.by_id(0).pos = [25., 25., 25.]  # Box is [20, 20, 20]

        # Create interface with unfolded positions
        ase_interface = espressomd.plugins.ase.ASEInterface(
            system=system,
            particle_slice=system.part.all(),
            use_folded_positions=False
        )

        # Get unfolded position from ESPResSo
        unfolded_pos = np.copy(system.part.by_id(0).pos)

        # Check that ASE atoms use unfolded position
        atoms = ase_interface.atoms
        np.testing.assert_allclose(atoms.positions[0], unfolded_pos)

        # Unfolded position should match what we set (outside box)
        np.testing.assert_allclose(unfolded_pos, [25., 25., 25.])

    def test_folded_positions_update(self):
        """Test that folded positions are updated correctly during update_ase()."""
        # Create interface with folded positions
        ase_interface = espressomd.plugins.ase.ASEInterface(
            system=system,
            particle_slice=system.part.all(),
            use_folded_positions=True
        )

        # Move particle outside box
        system.part.by_id(0).pos = [22., 23., 24.]

        # Update ASE
        ase_interface.update_ase()

        # Check that folded position is used
        folded_pos = np.copy(system.part.by_id(0).pos_folded)
        atoms = ase_interface.atoms
        np.testing.assert_allclose(atoms.positions[0], folded_pos)

        # Folded position should be inside box
        self.assertTrue(np.all(atoms.positions[0] >= 0))
        self.assertTrue(np.all(atoms.positions[0] < system.box_l))


class FixedForceCalculator(Calculator):
    """ASE calculator that returns fixed forces for testing."""

    implemented_properties = ["energy", "forces"]

    def __init__(self, forces_array, energy=0.0):
        """
        Initialize with fixed forces array and optional energy.

        Parameters
        ----------
        forces_array : np.ndarray
            Array of forces with shape (n_atoms, 3)
        energy : float, optional
            Total energy value to return
        """
        Calculator.__init__(self)
        self.forces_array = np.array(forces_array)
        self.energy_value = energy

    def calculate(self, atoms=None, properties=("forces",),
                  system_changes=("positions",)):
        """Calculate forces and energy."""
        Calculator.calculate(self, atoms, properties, system_changes)

        if "forces" in properties:
            self.results["forces"] = self.forces_array.copy()
        if "energy" in properties:
            self.results["energy"] = self.energy_value

    def get_forces(self, atoms):  # pylint: disable=unused-argument
        """Get forces directly."""
        return self.forces_array.copy()


@utx.skipIfMissingFeatures("EXTERNAL_FORCES")
@utx.skipIfMissingModules("ase.calculators")
class ASEIntegrationTest(ut.TestCase):
    """Test suite for ASE interface integration functionality."""

    def setUp(self):
        """Set up system with particles for integration tests."""
        system.time_step = 0.01
        system.box_l = [20., 20., 20.]  # Ensure correct box size

        # Add particles with different masses and initial positions
        system.part.add(pos=[5., 5., 5.], mass=1.0, v=[0., 0., 0.], type=0)
        system.part.add(pos=[10., 10., 10.], mass=2.0, v=[0., 0., 0.], type=1)
        system.part.add(pos=[15., 15., 15.], mass=0.5, v=[0., 0., 0.], type=0)

        # Set up velocity verlet integrator
        system.integrator.set_vv()

    def tearDown(self):
        """Clean up after each test."""
        system.part.clear()
        system.integrator.set_vv()

    def test_newton_second_law_integration(self):
        """Test 5-step integration to verify Newton's second law."""
        # Create ASE interface
        system.cell_system.skin = 0.2
        system.integrator.set_vv()
        ase_interface = espressomd.plugins.ase.ASEInterface(
            system=system,
            particle_slice=system.part.all(),
            export_masses=True
        )

        # Define constant forces for testing Newton's second law
        # F = ma, so we expect acceleration = F/m
        test_forces = np.array([
            # Force on particle 0 (mass=1.0), expected a = [1.0, 0.0, 0.0]
            [1.0, 0.0, 0.0],
            # Force on particle 1 (mass=2.0), expected a = [2.0, 0.0, 0.0]
            [4.0, 0.0, 0.0],
            # Force on particle 2 (mass=0.5), expected a = [2.0, 0.0, 0.0]
            [1.0, 0.0, 0.0]
        ])

        # Create calculator with fixed forces
        calculator = FixedForceCalculator(test_forces)

        # Store initial positions and velocities
        initial_pos = np.copy(system.part.all().pos)
        initial_vel = np.copy(system.part.all().v)

        # Integrate for 5 steps
        steps_performed = ase_interface.integrate(5, calculator)

        # Verify 5 steps were performed
        self.assertEqual(steps_performed, 5)

        # Get final positions and velocities
        final_pos = np.copy(system.part.all().pos)
        final_vel = np.copy(system.part.all().v)

        # Expected accelerations based on F = ma
        expected_accel = np.array([
            [1.0, 0.0, 0.0],   # F=1.0, m=1.0 -> a=1.0
            [2.0, 0.0, 0.0],   # F=4.0, m=2.0 -> a=2.0
            [2.0, 0.0, 0.0]    # F=1.0, m=0.5 -> a=2.0
        ])

        # For constant acceleration, v_final = v_initial + a*t
        dt = system.time_step
        total_time = 5 * dt
        expected_final_vel = initial_vel + expected_accel * total_time

        # For constant acceleration, x_final = x_initial + v_initial*t +
        # 0.5*a*t^2
        expected_final_pos = (initial_pos + initial_vel * total_time +
                              0.5 * expected_accel * total_time**2)

        # Test velocities (allow some numerical tolerance)
        np.testing.assert_allclose(final_vel, expected_final_vel, rtol=1e-10,
                                   err_msg="Final velocities don't match Newton's second law prediction")

        # Test positions (allow some numerical tolerance)
        np.testing.assert_allclose(final_pos, expected_final_pos, rtol=1e-10,
                                   err_msg="Final positions don't match kinematic prediction")

    def test_newton_second_law_varying_forces(self):
        """Test 5-step integration with different forces in each time step."""
        # Create ASE interface
        system.cell_system.skin = 0.2
        system.integrator.set_vv()
        ase_interface = espressomd.plugins.ase.ASEInterface(
            system=system,
            particle_slice=system.part.all(),
            export_masses=True
        )

        # Define different forces for each of the 5 time steps
        # Each entry is forces for [particle0, particle1, particle2] at that
        # time step
        forces_per_step = np.array([
            [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.5, 0.0, 0.0]],  # Step 1
            [[0.0, 1.0, 0.0], [0.0, 4.0, 0.0], [0.0, 1.0, 0.0]],  # Step 2
            [[0.0, 0.0, 1.0], [0.0, 0.0, 2.0], [0.0, 0.0, 2.0]],  # Step 3
            [[-1.0, 0.0, 0.0], [-2.0, 0.0, 0.0], [-0.5, 0.0, 0.0]],  # Step 4
            [[2.0, 1.0, 0.5], [4.0, 2.0, 1.0], [1.0, 0.5, 1.0]]   # Step 5
        ])

        # Store initial conditions
        initial_pos = np.copy(system.part.all().pos)
        initial_vel = np.copy(system.part.all().v)
        masses = np.copy(system.part.all().mass)

        # Track positions and velocities at each step for verification
        positions = [initial_pos.copy()]
        velocities = [initial_vel.copy()]

        dt = system.time_step

        # Integrate step by step with different forces
        for step in range(5):
            # Create calculator with forces for this step
            calculator = FixedForceCalculator(forces_per_step[step])

            # Perform one integration step
            steps_performed = ase_interface.integrate(1, calculator)
            self.assertEqual(steps_performed, 1)

            # Record state after this step
            positions.append(np.copy(system.part.all().pos))
            velocities.append(np.copy(system.part.all().v))

        # Verify Newton's second law for each step
        # For each step: F = ma, so a = F/m
        # Using velocity Verlet: v(t+dt) = v(t) + a*dt, x(t+dt) = x(t) +
        # v(t)*dt + 0.5*a*dt^2

        for step in range(5):
            current_forces = forces_per_step[step]
            # F/m for each particle
            expected_accel = current_forces / masses[:, np.newaxis]

            # Expected velocity after this step
            expected_vel = velocities[step] + expected_accel * dt

            # Expected position after this step
            expected_pos = (positions[step] + velocities[step] * dt +
                            0.5 * expected_accel * dt**2)

            # Compare with actual results
            np.testing.assert_allclose(
                velocities[step + 1], expected_vel, rtol=1e-10,
                err_msg=f"Step {step + 1}: Velocities don't match Newton's second law prediction"  # nopep8
            )

            np.testing.assert_allclose(
                positions[step + 1], expected_pos, rtol=1e-10,
                err_msg=f"Step {step + 1}: Positions don't match kinematic prediction"  # nopep8
            )

        # Additional verification: check that different forces produced
        # different accelerations
        for step in range(4):  # Compare adjacent steps
            forces_1 = forces_per_step[step]
            forces_2 = forces_per_step[step + 1]

            # Forces should be different between steps
            self.assertFalse(np.allclose(forces_1, forces_2),
                             f"Forces in steps {step + 1} and {step + 2} should be different")

            # This should result in different velocity changes
            accel_1 = forces_1 / masses[:, np.newaxis]
            accel_2 = forces_2 / masses[:, np.newaxis]

            vel_change_1 = velocities[step + 1] - velocities[step]
            vel_change_2 = velocities[step + 2] - velocities[step + 1]

            expected_vel_change_1 = accel_1 * dt
            expected_vel_change_2 = accel_2 * dt

            np.testing.assert_allclose(
                vel_change_1, expected_vel_change_1, rtol=1e-10)
            np.testing.assert_allclose(
                vel_change_2, expected_vel_change_2, rtol=1e-10)

    def test_steepest_descent_integration(self):
        """Test steepest descent integration with max_displacement constraint."""
        system.cell_system.skin = 0.2
        # Create ASE interface
        ase_interface = espressomd.plugins.ase.ASEInterface(
            system=system,
            particle_slice=system.part.all()
        )

        # Define forces pointing in positive x direction with different magnitudes
        # These will be used to test max_displacement limit
        large_forces = np.array([
            [100.0, 0.0, 0.0],  # Large force on particle 0
            [50.0, 0.0, 0.0],   # Medium force on particle 1
            [200.0, 0.0, 0.0]   # Very large force on particle 2
        ])

        calculator = FixedForceCalculator(large_forces)

        # Set up steepest descent integrator with small max_displacement
        max_displacement = 0.1
        system.integrator.set_steepest_descent(
            f_max=0,  # No force threshold, run for specified steps
            gamma=1.0,
            max_displacement=max_displacement
        )

        # Store initial positions
        initial_pos = np.copy(system.part.all().pos)

        # Run one integration step
        steps_performed = ase_interface.integrate(1, calculator)

        # Verify one step was performed
        self.assertEqual(steps_performed, 1)

        # Get final positions
        final_pos = np.copy(system.part.all().pos)

        # Calculate actual displacements
        displacements = final_pos - initial_pos
        displacement_magnitudes = np.linalg.norm(displacements, axis=1)

        # All displacements should be limited by max_displacement
        for i, disp_mag in enumerate(displacement_magnitudes):
            self.assertLessEqual(disp_mag, max_displacement + 1e-10,
                                 f"Particle {i} displacement {disp_mag:.6f} exceeds max_displacement {max_displacement}")

        # All displacements should be in the positive x direction (force
        # direction)
        for i, disp in enumerate(displacements):
            self.assertGreater(
                disp[0], 0, f"Particle {i} should move in positive x direction")
            self.assertAlmostEqual(
                disp[1], 0, places=10, msg=f"Particle {i} should not move in y direction")
            self.assertAlmostEqual(
                disp[2], 0, places=10, msg=f"Particle {i} should not move in z direction")

        # For particles with large forces, displacement should equal max_displacement
        # (since gamma * F  >> max_displacement, the displacement is limited)
        np.testing.assert_allclose(displacement_magnitudes, max_displacement, rtol=1e-10,
                                   err_msg="All particles should move exactly max_displacement distance")

    def test_steepest_descent_small_forces(self):
        """Test steepest descent with small forces where max_displacement is not limiting."""
        system.cell_system.skin = 0.2
        # Create ASE interface
        ase_interface = espressomd.plugins.ase.ASEInterface(
            system=system,
            particle_slice=system.part.all()
        )

        # Define small forces
        small_forces = np.array([
            [0.01, 0.0, 0.0],   # Small force on particle 0
            [0.02, 0.0, 0.0],   # Small force on particle 1
            [0.005, 0.0, 0.0]   # Very small force on particle 2
        ])

        calculator = FixedForceCalculator(small_forces)

        # Set up steepest descent with large max_displacement but small gamma
        max_displacement = 1.0  # Large max displacement
        gamma = 0.1  # Small gamma
        system.integrator.set_steepest_descent(
            f_max=0,
            gamma=gamma,
            max_displacement=max_displacement
        )

        # Store initial positions
        initial_pos = np.copy(system.part.all().pos)

        # Run one integration step
        steps_performed = ase_interface.integrate(1, calculator)
        self.assertEqual(steps_performed, 1)

        # Get final positions
        final_pos = np.copy(system.part.all().pos)
        displacements = final_pos - initial_pos
        displacement_magnitudes = np.linalg.norm(displacements, axis=1)

        # Calculate expected displacements based on steepest descent formula:
        # displacement = min(|gamma * F, max_displacement) * F/|F|
        force_magnitudes = np.linalg.norm(small_forces, axis=1)
        expected_displacements = gamma * force_magnitudes

        # All expected displacements should be less than max_displacement
        for expected_disp in expected_displacements:
            self.assertLess(expected_disp, max_displacement)

        # Actual displacements should match expected (not limited by
        # max_displacement)
        np.testing.assert_allclose(displacement_magnitudes, expected_displacements, rtol=1e-10,
                                   err_msg="Displacements should match gamma* when not limited by max_displacement")


if __name__ == "__main__":
    ut.main()
