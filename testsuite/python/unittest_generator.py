#
# Copyright (C) 2022 The ESPResSo project
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

import sys
import inspect
import pathlib


class TestGenerator:
    """
    Class factory to generate test classes on-the-fly. Parse the test name
    from the command line arguments to create the corresponding test case.

    Use this class as follows::

        import unittest as ut
        import unittest_generator as utg

        config = utg.TestGenerator()
        modes = config.get_modes()

        class Test(ut.TestCase):
            def test_modes(self):
                self.assertIn("P3M.CPU", modes)

        if __name__ == "__main__":
            config.bind_test_class(Test)
            ut.main()

    Invoke the test file with a test name of the form "Test_{suffix}__{modes}"
    to dynamically bind to the ``Test`` class and populate the local variable
    ``modes``. The suffix is optional. Example:

    .. code-block:: none

        $ python3 test.py Test__lb_cpu_ascii__p3m_cpu
        .
        ----------------------------------------------------------------------
        Ran 1 test in 0.000s

        OK

    To make the test fail, do:

    .. code-block:: none

        $ python3 test.py Test__lb_cpu_ascii__p3m_gpu
        F
        ======================================================================
        FAIL: test_modes (__main__.Test)
        ----------------------------------------------------------------------
        AssertionError: 'P3M.CPU' not found in {'LB', 'LB.CPU', 'LB.ASCII', 'P3M', 'P3M.GPU'}

        ----------------------------------------------------------------------
        Ran 1 test in 0.000s

        FAILED (failures=1)

    Since the test name contains information that is required to run the test,
    the test file can no longer be invoked without a test name as argument:

    .. code-block:: none

        $ python3 test.py
        AssertionError: please provide a test name as argument,
        like 'Test_lb_cpu__p3m_cpu' (got ['test.py'])

    At the CMake level, configure CTest to invoke the test file with a test
    name as argument. An optional suffix can be added to disambiguate two
    tests that take the same modes but run with e.g. a different number of
    MPI ranks. This framework provides a solution to the Python ticket
    `39283 <https://bugs.python.org/issue39283>`__.

    Although CMake provides a function ``configure_file()`` that substitutes
    CMake variables into test files that are copied to the build folder, this
    is only useful for variables that are known at build time. It also forces
    us to copy the original test file multiple times with a different filepath
    in the build folder, making it difficult for Python code coverage tools to
    identify which file in the root folder matches the recorded coverage
    information.

    """

    def __init__(self):
        self.main_module = inspect.getmodule(inspect.stack()[1][0])
        self.test_name = None
        self.test_feat = None
        prefix = 'Test_'
        for arg in sys.argv:
            if arg.startswith(prefix):
                self.test_name = arg
                self.test_feat = arg.split('__', 1)[1]
                self.test_idx = arg.split('_', 1)[1].lstrip('_')
                break
        err_msg = f"please provide a test name as argument, like '{prefix}lb_cpu__p3m_cpu' (got {sys.argv})"
        assert self.test_name is not None, err_msg

    def bind_test_class(self, base_class):
        """
        Dynamically re-bind an existing test class in the main module
        using the test name passed as command line argument. When running
        the test file, the original test class name is displayed instead
        of the new name.
        """
        setattr(self.main_module, self.test_name, base_class)

    def get_modes(self):
        """
        Generate the list of modes for the test name found in the command line.
        A mode consists of a feature with zero or more options; separate
        features with 2 underscores and options with 1 underscore (options
        can appear in any order). For example, "p3m_cpu__lb_cpu_ascii"
        generates modes P3M, P3M.CPU, LB, LB.CPU, LB.ASCII.
        """
        modes = set()
        for item in self.test_feat.upper().split('__'):
            feature, *options = item.split('_')
            for option in options:
                modes.add(f"{feature}.{option}")
            modes.add(feature)
        return modes

    @classmethod
    def recursive_unlink(cls, root):
        """
        Delete files in a folder recursively but preserve the tree structure.
        """
        if root.exists():
            for filepath in root.iterdir():
                if filepath.is_file():
                    filepath.unlink()

    def cleanup_old_checkpoint(self):
        """
        Remove the contents of the checkpoint directory if it exists.
        The directory itself and its subfolder structure are preserved
        since they will typically be created soon afterwards (risk of
        race condition on file systems with latency).
        """
        args = self.get_checkpoint_params()
        root = pathlib.Path(args["checkpoint_path"]) / args["checkpoint_id"]
        self.recursive_unlink(root)

    def get_checkpoint_params(self):
        """
        Generate parameters to instantiate an ESPResSo checkpoint file.
        """
        return {"checkpoint_id": f"checkpoint_{self.test_idx}",
                "checkpoint_path": str(pathlib.Path(__file__).parent)}
