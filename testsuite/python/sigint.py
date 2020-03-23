#
# Copyright (C) 2019 The ESPResSo project
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
import signal
import subprocess
import time


class SigintTest(ut.TestCase):

    def setUp(self):
        self.process = subprocess.Popen(
            ['@CMAKE_BINARY_DIR@/pypresso',
             '@CMAKE_CURRENT_BINARY_DIR@/sigint_child.py'])

    def test_signal_handling(self):
        self.process.send_signal(signal.SIGINT)
        # Wait for the signal to arrive and one integration step to be finished
        time.sleep(1)
        self.assertIsNotNone(self.process.poll())


if __name__ == '__main__':
    ut.main()
