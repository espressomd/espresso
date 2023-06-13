#
# Copyright (C) 2019-2022 The ESPResSo project
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
import sys
import pathlib
import os


EXPECTED_TRACEBACK_ENDING = """ in handle_sigint
    signal.raise_signal(signal.Signals.SIGINT)
KeyboardInterrupt
"""


class SigintTest(ut.TestCase):

    script = str(pathlib.Path(__file__).parent / 'sigint_child.py')

    def check_signal_handling(self, process, sig):
        # send signal
        process.send_signal(sig)
        # capture stderr and return code (negative of signum)
        stdout, stderr = process.communicate(input=None, timeout=6.)
        assert stdout is None
        traceback = stderr.decode()
        return_code = process.poll()
        signum = -return_code
        self.assertEqual(signum, sig.value)
        if sig == signal.Signals.SIGTERM:
            self.assertEqual(traceback, "")
        elif sig == signal.Signals.SIGINT:
            self.assertIn(" self.integrator.run(", traceback)
            self.assertTrue(traceback.endswith(EXPECTED_TRACEBACK_ENDING),
                            msg=f"Traceback failed string match:\n{traceback}")

    def test_signal_handling(self):
        signals = [signal.Signals.SIGINT, signal.Signals.SIGTERM]
        processes = []
        # open asynchronous processes with non-blocking read access on stderr
        for _ in range(len(signals)):
            process = subprocess.Popen([sys.executable, self.script],
                                       stderr=subprocess.PIPE)
            os.set_blocking(process.stderr.fileno(), False)
            processes.append(process)

        # wait for the script to reach the integration loop
        time.sleep(0.5)
        for process, sig in zip(processes, signals):
            tick = time.time()
            while True:
                message = process.stderr.readline().decode()
                if message == "start of integration loop\n":
                    # wait for the script to enter the integrator run method
                    time.sleep(0.1)
                    # send signal and check process behavior
                    self.check_signal_handling(process, sig)
                    break
                tock = time.time()
                assert tock - tick < 8., "subprocess timed out"
                time.sleep(0.1)


if __name__ == '__main__':
    ut.main()
