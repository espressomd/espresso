#
# Copyright (C) 2025 The ESPResSo project
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

"""
Tests the Pairwise Distance observable and the ContactTimes accumulator.
"""

import unittest as ut
import unittest_decorators as utx
import numpy as np
import espressomd
import espressomd.accumulators
import espressomd.observables


class ContactTimeTest(ut.TestCase):

    @classmethod
    def calculate_minimum_image_distance(cls, coords1, coords2, size_box):
        """
        Calculates the minimum distance between two particles in the simulation box,
        taking into account the minimum image convention and periodic boundary conditions.

        Args:
            coords1(`list`): array with the cartesian coordinates of the first particle.
            coords2(`list`): array with the cartesian coordinates of the second particle.
            size_box(`list`): xyz size of the simulation box.

        Returns:
            dist(`float`): minimum distance between the particles
        """
        # Sanity check, check that all input variables have the same size
        assert len(coords1) == len(size_box)
        assert len(coords2) == len(size_box)

        n_axes = len(size_box)
        coords = {"particle1": coords1,
                  "particle2": coords2}

        # Fold particles
        for particle in coords.keys():
            for axis in range(n_axes):
                coords[particle][axis] = coords[particle][axis] % size_box[axis]

        # Calculate minimum distance
        dist2 = 0.
        for axis in range(n_axes):
            dist = abs(coords["particle1"][axis] -
                       coords["particle2"][axis])
            if dist > size_box[axis] * 0.5:
                dist -= size_box[axis]
            dist2 += dist**2
        return np.sqrt(dist2)

    @classmethod
    def calculate_contact_times(
            cls, system, accumulator_pos, contact_threshold):
        """
        Calculate the contact times from the time series of the coordinates.

        Args:
            accumulator_pos(`espressomd.accumulators.TimeSeries`): accumulator with the time series of the particles.
            contact_threshold(`float`): threshold distance under which particles are considered to be in contact.
            system(`espressomd.System`): instance of the espresso simulation system.
        """
        # Get the coordinates
        coords = accumulator_pos.time_series()
        # Check for the initial contacts (if any)
        n_ids = len(coords[0])
        contacts = np.zeros((n_ids, n_ids), dtype=bool)
        time_first_contact = np.zeros((n_ids, n_ids), dtype=float)
        # Calculate contact times:
        frame_time = 0
        contact_times = []
        distances = []
        for coords_frame in coords:
            frame_dists = []
            for index1 in range(n_ids):
                for index2 in range(index1 + 1, n_ids):
                    dist = cls.calculate_minimum_image_distance(
                        coords_frame[index1], coords_frame[index2], system.box_l)
                    frame_dists.append(dist)
                    if dist < contact_threshold:
                        # They are in contact now
                        if not contacts[index1][index2]:
                            # but they were not contact before!
                            contacts[index1][index2] = True
                            time_first_contact[index1][index2] = frame_time
                    else:
                        # They are not in contact
                        if contacts[index1][index2]:
                            # but they were in contact before!
                            # Calculate the total contact time
                            contact_time = frame_time - \
                                time_first_contact[index1][index2] - \
                                system.time_step
                            if contact_time < system.time_step:
                                # Avoid rounding errors
                                contact_time = 0
                            contacts[index1][index2] = False
                            contact_times.append(float(contact_time))
            frame_time += system.time_step
            distances.append(frame_dists)
        return (np.array(contact_times), np.array(distances))

    @utx.skipIfMissingFeatures(["LENNARD_JONES"])
    def test_contact_time(self):
        # System: LJ fluid
        contact_threshold = 1.
        seed = 23
        N = 40

        # Simulation parameters
        N_steps = 1400
        box_l = 20
        system = espressomd.System(box_l=[box_l] * 3)
        system.time_step = 0.01
        system.cell_system.skin = 1
        rng = np.random.default_rng(seed)
        system.part.add(pos=rng.random((N, 3)) * np.copy(system.box_l))
        ids = np.copy(system.part.all().id)

        # relax the system
        system.integrator.set_steepest_descent(f_max=1e-3,
                                               gamma=1.,
                                               max_displacement=0.1)
        system.integrator.run(1000)

        # setup LD integration
        system.thermostat.set_langevin(kT=1., gamma=1., seed=seed)
        system.integrator.set_vv()

        # Setup LJ interaction
        system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=1., sigma=1., shift="auto", cutoff=3.)

        # Setup observables to track pairwise distances and particle positions
        pairwise_dist_obs = espressomd.observables.PairwiseDistances(ids=ids,
                                                                     target_ids=ids)
        particle_pos_obs = espressomd.observables.ParticlePositions(ids=ids)

        # Setup the accumulators to track the contact times and the time series
        contact_time_accumulator = espressomd.accumulators.ContactTimes(
            obs=pairwise_dist_obs, delta_N=1, contact_threshold=contact_threshold)

        accumulator_pos = espressomd.accumulators.TimeSeries(obs=particle_pos_obs,
                                                             delta_N=1)
        accumulator_dist = espressomd.accumulators.TimeSeries(obs=pairwise_dist_obs,
                                                              delta_N=1)
        system.auto_update_accumulators.add(contact_time_accumulator)
        system.auto_update_accumulators.add(accumulator_dist)
        system.auto_update_accumulators.add(accumulator_pos)

        # Run the simulation
        system.integrator.run(N_steps)

        # Get the pairwise distances calculated on the fly
        distances_test = accumulator_dist.time_series()

        # Get the contact time calculated on the fly
        contact_times_test = contact_time_accumulator.contact_times()
        # Calculate distances and contact times directly from the time series
        contact_times_ref, distances_benchmark = self.calculate_contact_times(
            system, accumulator_pos, contact_threshold)

        np.testing.assert_allclose(actual=distances_test,
                                   desired=distances_benchmark,
                                   err_msg="The pairwise distances calculated on the fly are not consistent with the ones calculated from the trajectory",
                                   verbose=True, atol=1e-12, rtol=1e-7)

        np.testing.assert_allclose(actual=contact_times_test,
                                   desired=contact_times_ref,
                                   err_msg="The contact times calculated on the fly are not consistent with the ones calculated from the trajectory",
                                   verbose=True, atol=1e-12, rtol=1e-7)

        np.testing.assert_array_equal(contact_time_accumulator.shape(), [9])
        np.testing.assert_array_equal(pairwise_dist_obs.shape(),
                                      [(N * (N - 1)) // 2])

        # clear accumulator
        contact_time_accumulator.clear()
        self.assertEqual(len(contact_time_accumulator.contact_times()), 0)
        np.testing.assert_array_equal(contact_time_accumulator.shape(), [0])
        contact_time_accumulator.update()
        self.assertEqual(len(contact_time_accumulator.contact_times()), 0)
        np.testing.assert_array_equal(contact_time_accumulator.shape(), [0])

    def test_interface(self):
        ref_ids1 = [1, 5, 2]
        ref_ids2 = [5, 7, 9]
        obs = espressomd.observables.PairwiseDistances(
            ids=ref_ids1, target_ids=ref_ids2)
        acc = espressomd.accumulators.ContactTimes(
            obs=obs, delta_N=1, contact_threshold=0.2)
        np.testing.assert_array_equal(obs.ids, ref_ids1)
        np.testing.assert_array_equal(obs.target_ids, ref_ids2)
        self.assertAlmostEqual(acc.contact_threshold, 0.2, delta=1e-7)
        self.assertIsNone(obs.call_method("unknown"))
        self.assertIsNone(acc.call_method("unknown"))
        # empty set case
        obs = espressomd.observables.PairwiseDistances(ids=[], target_ids=[])
        acc = espressomd.accumulators.ContactTimes(
            obs=obs, delta_N=1, contact_threshold=0.2)
        res = acc.contact_times()
        self.assertEqual(len(res), 0)

    def test_exceptions(self):
        with self.assertRaisesRegex(ValueError, "Attribute 'contact_threshold' must be >= 0"):
            # check input data in accumulators.ContactTimes
            obs = espressomd.observables.PairwiseDistances(
                ids=[1], target_ids=[2])
            espressomd.accumulators.ContactTimes(
                obs=obs, delta_N=1, contact_threshold=-1.)


if __name__ == "__main__":
    ut.main()
