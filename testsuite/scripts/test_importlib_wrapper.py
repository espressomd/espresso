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
import ast


class importlib_wrapper(ut.TestCase):

    def test_substitute_variable_values(self):
        str_inp = "n_steps=5000\nnsteps == 5\n"
        str_exp = "n_steps = 10; _n_steps__original=5000\nnsteps == 5\n"
        str_out = iw.substitute_variable_values(str_inp, n_steps=10)
        self.assertEqual(str_out, str_exp)
        str_out = iw.substitute_variable_values(str_inp, n_steps='10',
                                                strings_as_is=True)
        self.assertEqual(str_out, str_exp)
        str_inp = "N=5000\nnsteps == 5\n"
        str_exp = "N = 10\nnsteps == 5\n"
        str_out = iw.substitute_variable_values(str_inp, N=10, keep_original=0)
        self.assertEqual(str_out, str_exp)
        # test exceptions
        str_inp = "n_steps=5000\nnsteps == 5\n"
        self.assertRaises(AssertionError, iw.substitute_variable_values,
                          str_inp, other_var=10)
        str_inp = "other_var == 5\n"
        self.assertRaises(AssertionError, iw.substitute_variable_values,
                          str_inp, other_var=10)
        str_inp = "var, other_var = 5, 6\n"
        self.assertRaises(AssertionError, iw.substitute_variable_values,
                          str_inp, var=10)

    def test_set_cmd(self):
        original_sys_argv = list(sys.argv)
        sys.argv = [0, "test"]
        # test substitutions
        str_inp = "import sys\nimport argparse"
        str_exp = "import sys;sys.argv = ['a.py', '1', '2'];" + str_inp
        str_out, sys_argv = iw.set_cmd(str_inp, "a.py", (1, 2))
        self.assertEqual(str_out, str_exp)
        self.assertEqual(sys_argv, [0, "test"])
        str_inp = "import argparse"
        str_exp = "import sys;sys.argv = ['a.py', '1', '2'];" + str_inp
        str_out, sys_argv = iw.set_cmd(str_inp, "a.py", ["1", 2])
        self.assertEqual(str_out, str_exp)
        self.assertEqual(sys_argv, [0, "test"])
        # test exceptions
        str_inp = "import re"
        self.assertRaises(AssertionError, iw.set_cmd, str_inp, "a.py", (1, 2))
        # restore sys.argv
        sys.argv = original_sys_argv

    def test_disable_matplotlib_gui(self):
        str_inp = "if 1:\n\timport matplotlib as mp\nmp.use('PS')\n"
        str_exp = ("if 1:\n\timport matplotlib as _mpl;_mpl.use('Agg');"
                   "import matplotlib as mp\n#mp.use('PS')\n")
        str_out = iw.disable_matplotlib_gui(str_inp)
        self.assertEqual(str_out, str_exp)
        str_inp = "if 1:\n    import matplotlib.pyplot as plt\nplt.ion()\n"
        str_exp = ("if 1:\n    import matplotlib as _mpl;_mpl.use('Agg');"
                   "import matplotlib.pyplot as plt\n#plt.ion()\n")
        str_out = iw.disable_matplotlib_gui(str_inp)
        self.assertEqual(str_out, str_exp)
        str_inp = "if 1:\n\tget_ipython(\n).run_line_magic('matplotlib', 'x')\n"
        str_exp = "if 1:\n#\tget_ipython(\n#).run_line_magic('matplotlib', 'x')\n"
        str_out = iw.disable_matplotlib_gui(str_inp)
        self.assertEqual(str_out, str_exp)

    def test_set_random_seeds(self):
        # NumPy seed
        str_np = "import numpy as np\n"
        str_lambda = "(lambda *args, **kwargs: None)"
        str_inp = str_np + "np.random.seed(seed=42)"
        str_exp = str_np + "np.random.seed;_random_seed_np = " + \
            str_lambda + "(seed=42)"
        str_out = iw.set_random_seeds(str_inp)
        self.assertEqual(str_out, str_exp)
        str_np = "import numpy.random as npr\n"
        str_inp = str_np + "npr.seed(42)"
        str_exp = str_np + "npr.seed;_random_seed_np = " + str_lambda + "(42)"
        str_out = iw.set_random_seeds(str_inp)
        self.assertEqual(str_out, str_exp)

    def test_mock_es_visualization(self):
        statement = "import espressomd.visualization"
        expected = """
try:
    import espressomd.visualization
    if hasattr(espressomd.visualization.mayaviLive, 'deferred_ImportError') or \\
       hasattr(espressomd.visualization.openGLLive, 'deferred_ImportError'):
        raise ImportError()
except ImportError:
    from unittest.mock import MagicMock
    import espressomd
    espressomd.visualization = MagicMock()
"""
        self.assertEqual(iw.mock_es_visualization(statement), expected[1:])

        statement = "import espressomd.visualization as test"
        expected = """
try:
    import espressomd.visualization as test
    if hasattr(test.mayaviLive, 'deferred_ImportError') or \\
       hasattr(test.openGLLive, 'deferred_ImportError'):
        raise ImportError()
except ImportError:
    from unittest.mock import MagicMock
    import espressomd
    test = MagicMock()
"""
        self.assertEqual(iw.mock_es_visualization(statement), expected[1:])

        statement = "import espressomd.visualization, espressomd.visualization as test"
        expected = """
try:
    import espressomd.visualization
    if hasattr(espressomd.visualization.mayaviLive, 'deferred_ImportError') or \\
       hasattr(espressomd.visualization.openGLLive, 'deferred_ImportError'):
        raise ImportError()
except ImportError:
    from unittest.mock import MagicMock
    import espressomd
    espressomd.visualization = MagicMock()
try:
    import espressomd.visualization as test
    if hasattr(test.mayaviLive, 'deferred_ImportError') or \\
       hasattr(test.openGLLive, 'deferred_ImportError'):
        raise ImportError()
except ImportError:
    from unittest.mock import MagicMock
    import espressomd
    test = MagicMock()
"""
        self.assertEqual(iw.mock_es_visualization(statement), expected[1:])

        statement = "from espressomd import visualization"
        expected = """
try:
    from espressomd import visualization
    if hasattr(visualization.mayaviLive, 'deferred_ImportError') or \\
       hasattr(visualization.openGLLive, 'deferred_ImportError'):
        raise ImportError()
except ImportError:
    from unittest.mock import MagicMock
    import espressomd
    visualization = MagicMock()
"""
        self.assertEqual(iw.mock_es_visualization(statement), expected[1:])

        statement = "from espressomd import visualization as test"
        expected = """
try:
    from espressomd import visualization as test
    if hasattr(test.mayaviLive, 'deferred_ImportError') or \\
       hasattr(test.openGLLive, 'deferred_ImportError'):
        raise ImportError()
except ImportError:
    from unittest.mock import MagicMock
    import espressomd
    test = MagicMock()
"""
        self.assertEqual(iw.mock_es_visualization(statement), expected[1:])

        statement = "from espressomd import visualization_mayavi"
        expected = """
try:
    from espressomd import visualization_mayavi
except ImportError:
    from unittest.mock import MagicMock
    import espressomd
    visualization_mayavi = MagicMock()
"""
        self.assertEqual(iw.mock_es_visualization(statement), expected[1:])

        statement = "from espressomd import visualization_mayavi as test"
        expected = """
try:
    from espressomd import visualization_mayavi as test
except ImportError:
    from unittest.mock import MagicMock
    import espressomd
    test = MagicMock()
"""
        self.assertEqual(iw.mock_es_visualization(statement), expected[1:])

        statement = "from espressomd.visualization_mayavi import mayaviLive"
        expected = """
try:
    from espressomd.visualization_mayavi import mayaviLive
except ImportError:
    from unittest.mock import MagicMock
    import espressomd
    mayaviLive = MagicMock()
"""
        self.assertEqual(iw.mock_es_visualization(statement), expected[1:])

        statement = "from espressomd.visualization_mayavi import mayaviLive as test"
        expected = """
try:
    from espressomd.visualization_mayavi import mayaviLive as test
except ImportError:
    from unittest.mock import MagicMock
    import espressomd
    test = MagicMock()
"""
        self.assertEqual(iw.mock_es_visualization(statement), expected[1:])

        statement = "from espressomd.visualization_mayavi import a as b, mayaviLive"
        expected = """
try:
    from espressomd.visualization_mayavi import a as b
except ImportError:
    from unittest.mock import MagicMock
    import espressomd
    b = MagicMock()
try:
    from espressomd.visualization_mayavi import mayaviLive
except ImportError:
    from unittest.mock import MagicMock
    import espressomd
    mayaviLive = MagicMock()
"""
        self.assertEqual(iw.mock_es_visualization(statement), expected[1:])

        statement = "from espressomd.visualization import openGLLive"
        expected = """
try:
    from espressomd.visualization import openGLLive
    if hasattr(openGLLive, 'deferred_ImportError'):
        raise openGLLive.deferred_ImportError
except ImportError:
    from unittest.mock import MagicMock
    import espressomd
    openGLLive = MagicMock()
"""
        self.assertEqual(iw.mock_es_visualization(statement), expected[1:])

        statement = "from espressomd.visualization import openGLLive as test"
        expected = """
try:
    from espressomd.visualization import openGLLive as test
    if hasattr(test, 'deferred_ImportError'):
        raise test.deferred_ImportError
except ImportError:
    from unittest.mock import MagicMock
    import espressomd
    test = MagicMock()
"""
        self.assertEqual(iw.mock_es_visualization(statement), expected[1:])

        # test exceptions
        statements_without_namespace = [
            "from espressomd.visualization import *",
            "from espressomd.visualization_opengl import *",
            "from espressomd.visualization_mayavi import *"
        ]
        for s in statements_without_namespace:
            self.assertRaises(ValueError, iw.mock_es_visualization, s)

    def test_matplotlib_pyplot_visitor(self):
        import_stmt = [
            'import matplotlib',
            'import matplotlib as mpl',
            'import matplotlib.pyplot',
            'import matplotlib.pyplot as plt',
            'import matplotlib.pyplot.figure as fig',
            'import ast, matplotlib',
            'import ast, matplotlib as mpl',
            'import ast, matplotlib.pyplot',
            'import ast, matplotlib.pyplot as plt',
            'import ast, matplotlib.pyplot.figure as fig',
            'from matplotlib import pyplot',
            'from matplotlib import pyplot as plt',
            'from matplotlib.pyplot import figure',
            'matplotlib.pyplot.ion()',
            'matplotlib.pyplot.ioff()',
            'plt.ion()',
            'mpl.use("PS")',
            'matplotlib.use("Agg")',
            'get_ipython().run_line_magic("matplotlib", "notebook")',
        ]
        tree = ast.parse('\n'.join(import_stmt))
        v = iw.GetMatplotlibPyplot()
        v.visit(tree)
        # find first line where matplotlib is imported
        self.assertEqual(v.matplotlib_first, 1)
        # find all aliases for matplotlib
        expected_mpl_aliases = ['matplotlib', 'mpl', 'matplotlib', 'mpl']
        self.assertEqual(v.matplotlib_aliases, expected_mpl_aliases)
        # find all aliases for matplotlib.pyplot
        expected_plt_aliases = [
            'matplotlib.pyplot', 'mpl.pyplot', 'matplotlib.pyplot', 'plt',
            'matplotlib.pyplot', 'mpl.pyplot', 'matplotlib.pyplot', 'plt',
            'pyplot', 'plt',
        ]
        self.assertEqual(v.pyplot_aliases, expected_plt_aliases)
        expected_plt_paths = {('matplotlib', 'pyplot'), ('mpl', 'pyplot'),
                              ('plt',), ('pyplot',)}
        self.assertEqual(v.pyplot_paths, expected_plt_paths)
        # find lines interactive mode, backend setup and magic functions
        self.assertEqual(v.pyplot_interactive_linenos, [14, 16])
        self.assertEqual(v.matplotlib_backend_linenos, [17, 18])
        self.assertEqual(v.ipython_magic_linenos, [19])

    def test_prng_seed_espressomd_system_visitor(self):
        import_stmt = [
            'sys0 = espressomd.System()   # nothing: espressomd not imported',
            'import espressomd as es1',
            'import espressomd.system as es2',
            'import espressomd.System as s1, espressomd.system.System as s2',
            'from espressomd import System as s3, electrostatics',
            'from espressomd.system import System as s4',
            'sys1 = es1.System()',
            'sys2 = es1.system.System()',
            'sys3 = es2.System()',
            'sys4 = s1()',
            'sys5 = s2()',
            'sys6 = s3()',
            'sys7 = s4()',
            'import numpy as np',
            'import numpy.random as npr1',
            'from numpy import random as npr2',
            'np.random.seed(1)',
            'npr1.seed(1)',
            'npr2.seed(1)',
        ]
        tree = ast.parse('\n'.join(import_stmt))
        v = iw.GetPrngSeedEspressomdSystem()
        v.visit(tree)
        # find all aliases for espressomd.system.System
        expected_es_sys_aliases = {'es1.System', 'es1.system.System',
                                   'es2.System', 's1', 's2', 's3', 's4'}
        self.assertEqual(v.es_system_aliases, expected_es_sys_aliases)
        # find all variables of type espressomd.system.System
        expected_es_sys_objs = set('sys' + str(i) for i in range(1, 8))
        self.assertEqual(v.variable_system_aliases, expected_es_sys_objs)
        # find all seeds setup
        self.assertEqual(v.numpy_seeds, [17, 18, 19])
        # test exceptions
        str_es_sys_list = [
            'import espressomd.System',
            'import espressomd.system.System',
            'from espressomd import System',
            'from espressomd.system import System',
        ]
        exception_stmt = [
            's, var = System(), 5',
            'class A:\n\ts = System()',
            'def A():\n\ts = System()',
        ]
        for str_es_sys in str_es_sys_list:
            for str_stmt in exception_stmt:
                for alias in ['', ' as EsSystem']:
                    str_import = str_es_sys + alias + '\n'
                    alias = str_import.split()[-1]
                    code = str_import + str_stmt.replace('System', alias)
                    v = iw.GetPrngSeedEspressomdSystem()
                    tree = ast.parse(code)
                    err_msg = v.__class__.__name__ + \
                        ' should fail on ' + repr(code)
                    with self.assertRaises(AssertionError, msg=err_msg):
                        v.visit(tree)

    def test_delimit_statements(self):
        lines = [
            'a = 1 # NEWLINE becomes NL after a comment',
            'print("""',
            '',
            '""")',
            '',
            'b = 1 +\\',
            '3 + (',
            '4)',
            'if True:',
            '   c = 1',
        ]
        source_code = '\n'.join(lines)
        linenos_exp = {1: 1, 2: 4, 6: 8, 9: 9, 10: 10}
        linenos_out = iw.delimit_statements(source_code)
        self.assertEqual(linenos_out, linenos_exp)


if __name__ == "__main__":
    ut.main()
