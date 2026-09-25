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
import numpy as np
import espressomd


class RecalcForcesTest(ut.TestCase):
    """
    Regression test for the ``recalc_forces`` flag of ``integrator.run``.

    With a Langevin thermostat active, a prior ``run(N>=1)`` performs a force
    calculation and leaves ``propagation.recalc_forces == false``. A
    subsequent ``run(0, recalc_forces=True)`` must, per the documented
    contract, recompute the forces unconditionally. Because the Langevin
    thermostat draws a stochastic force from a Philox counter that advanced
    during the prior ``run(1)``, a genuine recompute re-draws the random
    force at the new counter value, so the forces *change*.

    If the ``recalc_forces`` flag is dropped, the recompute is skipped because
    ``propagation.recalc_forces`` is false after the prior run, leaving the
    forces byte-identical (stale).
    """

    def test_recalc_forces_recomputes_langevin_force(self):
        system = espressomd.System(box_l=[10., 10., 10.])
        system.time_step = 0.01
        system.cell_system.skin = 0.4

        rng = np.random.default_rng(42)
        for _ in range(20):
            system.part.add(pos=rng.random(3) * system.box_l)

        system.thermostat.set_langevin(kT=1., gamma=1., seed=123)

        # Prior run(1): advances the Philox counter and leaves
        # propagation.recalc_forces == false.
        system.integrator.run(1)
        forces_before = np.copy(system.part.all().f)

        # The trigger: request an explicit recompute of forces. A correct
        # implementation reads recalc_forces=True and re-draws the Langevin
        # force at the advanced counter. Crucially, no particle / constraint /
        # interaction property is touched between the two runs, since that
        # would set propagation.recalc_forces=true and mask the bug.
        system.integrator.run(0, recalc_forces=True)
        forces_after = np.copy(system.part.all().f)

        # A genuine recompute redraws the stochastic force, so the forces must
        # differ. If recalc_forces is dropped, no recompute happens and the
        # two arrays are byte-identical.
        self.assertFalse(
            np.array_equal(forces_before, forces_after),
            msg="recalc_forces=True did not recompute forces: the Langevin "
                "stochastic force is byte-identical before and after "
                "run(0, recalc_forces=True), so forces were not recomputed.")


if __name__ == "__main__":
    ut.main()
