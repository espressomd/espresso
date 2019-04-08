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
        time.sleep(1)
                   # Wait for the signal to arrive and one integration step to
                   # be finished
        self.assertFalse(self.process.poll() is None)

if __name__ == '__main__':
    ut.main()
