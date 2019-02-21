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

import unittest as ut
import importlib_wrapper as iw
import sys


class importlib_wrapper(ut.TestCase):

    def test_substitute_variable_values(self):
        str_inp = "n_steps=5000\nnsteps == 5\n"
        str_exp = "n_steps = 10; _n_steps__original=5000\nnsteps == 5\n"
        str_out = iw.substitute_variable_values(str_inp, n_steps=10)
        self.assertEqual(str_out, str_exp)
        # test exceptions
        str_inp = "n_steps=5000\nnsteps == 5\n"
        self.assertRaises(AssertionError, iw.substitute_variable_values,
                          str_inp, other_var=10)
        str_inp = "other_var == 5\n"
        self.assertRaises(AssertionError, iw.substitute_variable_values,
                          str_inp, other_var=10)

    def test_set_cmd(self):
        original_sys_argv = list(sys.argv)
        sys.argv = [0, "test"]
        # test substitutions
        str_inp = "import sys\nimport argparse"
        str_exp = "import sys\nsys.argv = ['a.py', '1', '2']\nimport argparse"
        str_out, sys_argv = iw.set_cmd(str_inp, "a.py", (1, 2))
        self.assertEqual(str_out, str_exp)
        self.assertEqual(sys_argv, [0, "test"])
        str_inp = "import argparse"
        str_exp = "import argparse\nimport sys\nsys.argv = ['a.py', '1', '2']"
        str_out, sys_argv = iw.set_cmd(str_inp, "a.py", ["1", 2])
        self.assertEqual(str_out, str_exp)
        self.assertEqual(sys_argv, [0, "test"])
        # test exceptions
        str_inp = "import re"
        self.assertRaises(AssertionError, iw.set_cmd, str_inp, "a.py", (1, 2))
        # restore sys.argv
        sys.argv = original_sys_argv

    def test_disable_matplotlib_gui(self):
        str_inp = "\nimport matplotlib as mp\nmp.use('PS')\n"
        str_exp = ("\nimport matplotlib as _mpl;_mpl.use('Agg');"
                   "import matplotlib as mp\n\n")
        str_out = iw.disable_matplotlib_gui(str_inp)
        self.assertEqual(str_out, str_exp)

    def test_set_random_seeds(self):
        # ESPResSo seed
        str_es_sys = "system = espressomd.System(box_l=[box_l] * 3)\n"
        str_inp = str_es_sys + "system.set_random_state_PRNG()"
        str_exp = str_es_sys + "system.set_random_state_PRNG()"
        str_out = iw.set_random_seeds(str_inp)
        self.assertEqual(str_out, str_exp)
        str_inp = str_es_sys + "system.random_number_generator_state = 7 * [0]"
        str_exp = str_es_sys + "system.set_random_state_PRNG();" + \
            " _random_seed_es__original = 7 * [0]"
        str_out = iw.set_random_seeds(str_inp)
        self.assertEqual(str_out, str_exp)
        str_inp = str_es_sys + "system.seed = 42"
        str_exp = str_es_sys + "system.set_random_state_PRNG();" + \
            " _random_seed_es__original = 42"
        str_out = iw.set_random_seeds(str_inp)
        self.assertEqual(str_out, str_exp)
        # NumPy seed
        str_lambda = "(lambda *args, **kwargs: None)"
        str_inp = "\nnp.random.seed(seed=system.seed)"
        str_exp = "\n_random_seed_np = " + str_lambda + "(seed=system.seed)"
        str_out = iw.set_random_seeds(str_inp)
        self.assertEqual(str_out, str_exp)
        str_inp = "\nnumpy.random.seed(42)"
        str_exp = "\n_random_seed_np = " + str_lambda + "(42)"
        str_out = iw.set_random_seeds(str_inp)
        self.assertEqual(str_out, str_exp)

    def test_mock_es_visualization(self):
        mock_module = "unittest.mock" if sys.version_info >= (3, 3) else "mock"
        statement = "import espressomd.visualization"
        expected = """
try:
    import espressomd.visualization
    if hasattr(espressomd.visualization.mayaviLive, 'deferred_ImportError') or \\
       hasattr(espressomd.visualization.openGLLive, 'deferred_ImportError'):
        raise ImportError()
except ImportError:
    from {} import MagicMock
    import espressomd
    espressomd.visualization = MagicMock()
""".format(mock_module)
        self.assertEqual(iw.mock_es_visualization(statement), expected[1:])

        statement = "import espressomd.visualization as test"
        expected = """
try:
    import espressomd.visualization as test
    if hasattr(test.mayaviLive, 'deferred_ImportError') or \\
       hasattr(test.openGLLive, 'deferred_ImportError'):
        raise ImportError()
except ImportError:
    from {} import MagicMock
    import espressomd
    test = MagicMock()
""".format(mock_module)
        self.assertEqual(iw.mock_es_visualization(statement), expected[1:])

        statement = "from espressomd import visualization"
        expected = """
try:
    from espressomd import visualization
    if hasattr(visualization.mayaviLive, 'deferred_ImportError') or \\
       hasattr(visualization.openGLLive, 'deferred_ImportError'):
        raise ImportError()
except ImportError:
    from {} import MagicMock
    import espressomd
    visualization = MagicMock()
""".format(mock_module)
        self.assertEqual(iw.mock_es_visualization(statement), expected[1:])

        statement = "from espressomd import visualization as test"
        expected = """
try:
    from espressomd import visualization as test
    if hasattr(test.mayaviLive, 'deferred_ImportError') or \\
       hasattr(test.openGLLive, 'deferred_ImportError'):
        raise ImportError()
except ImportError:
    from {} import MagicMock
    import espressomd
    test = MagicMock()
""".format(mock_module)
        self.assertEqual(iw.mock_es_visualization(statement), expected[1:])

        statement = "from espressomd import visualization_mayavi"
        expected = """
try:
    from espressomd import visualization_mayavi
except ImportError:
    from {} import MagicMock
    import espressomd
    visualization_mayavi = MagicMock()
""".format(mock_module)
        self.assertEqual(iw.mock_es_visualization(statement), expected[1:])

        statement = "from espressomd import visualization_mayavi as test"
        expected = """
try:
    from espressomd import visualization_mayavi as test
except ImportError:
    from {} import MagicMock
    import espressomd
    test = MagicMock()
""".format(mock_module)
        self.assertEqual(iw.mock_es_visualization(statement), expected[1:])

        statement = "from espressomd.visualization_mayavi import mayaviLive"
        expected = """
try:
    from espressomd.visualization_mayavi import mayaviLive
except ImportError:
    from {} import MagicMock
    import espressomd
    mayaviLive = MagicMock()
""".format(mock_module)
        self.assertEqual(iw.mock_es_visualization(statement), expected[1:])

        statement = "from espressomd.visualization_mayavi import mayaviLive as test"
        expected = """
try:
    from espressomd.visualization_mayavi import mayaviLive as test
except ImportError:
    from {} import MagicMock
    import espressomd
    test = MagicMock()
""".format(mock_module)
        self.assertEqual(iw.mock_es_visualization(statement), expected[1:])

        statement = "from espressomd.visualization import openGLLive"
        expected = """
try:
    from espressomd.visualization import openGLLive
    if hasattr(openGLLive, 'deferred_ImportError'):
        raise openGLLive.deferred_ImportError
except ImportError:
    from {} import MagicMock
    import espressomd
    openGLLive = MagicMock()
""".format(mock_module)
        self.assertEqual(iw.mock_es_visualization(statement), expected[1:])

        statement = "from espressomd.visualization import openGLLive as test"
        expected = """
try:
    from espressomd.visualization import openGLLive as test
    if hasattr(test, 'deferred_ImportError'):
        raise test.deferred_ImportError
except ImportError:
    from {} import MagicMock
    import espressomd
    test = MagicMock()
""".format(mock_module)
        self.assertEqual(iw.mock_es_visualization(statement), expected[1:])

        # test exceptions
        statements_without_namespace = [
            "from espressomd.visualization import *",
            "from espressomd.visualization_opengl import *",
            "from espressomd.visualization_mayavi import *"
        ]
        for s in statements_without_namespace:
            self.assertRaises(AssertionError, iw.mock_es_visualization, s)


if __name__ == "__main__":
    ut.main()
