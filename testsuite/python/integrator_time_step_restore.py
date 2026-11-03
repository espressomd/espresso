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
import pickle

import espressomd


class IntegratorTimeStepRestore(ut.TestCase):
    """
    Regression test for the ``on_bind_system`` sentinel-skip guard in
    ``IntegratorHandle``. A :class:`System` whose ``time_step`` was never set
    keeps the ``-1.`` sentinel and the default Velocity Verlet integrator
    (``INTEG_METHOD_NVT``). Serializing and deserializing such a system must
    succeed: the guard is supposed to *skip* re-applying ``time_step == -1.``
    during ``on_bind_system`` so that ``System::set_time_step(-1.)`` is never
    called (it throws ``std::domain_error`` for non-positive time steps).
    """

    def test_restore_unset_time_step(self):
        system = espressomd.System(box_l=[1., 1., 1.])
        # time_step must still be the -1. sentinel (never set by the user)
        self.assertEqual(system.integrator.time_step, -1.)

        # round-trip via pickle: this triggers
        #   System._restore_object -> do_construct (else-branch) ->
        #   do_set_parameter("integrator", deserialized handle) ->
        #   IntegratorHandle::on_bind_system -> sentinel-skip guard
        state = pickle.dumps(system)
        restored = pickle.loads(state)

        # restore must succeed and preserve the unset sentinel
        self.assertEqual(restored.integrator.time_step, -1.)


if __name__ == "__main__":
    ut.main()
