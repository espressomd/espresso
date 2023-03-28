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

import unittest as ut
import importlib_wrapper as iw
import sys
import ast
import pathlib
import tempfile
import espressomd


class importlib_wrapper(ut.TestCase):

    def test_substitute_variable_values(self):
        str_inp = "n_steps=(\n5000)\nnsteps == 5\n"
        str_exp = "n_steps = 10; _n_steps__original=(\n5000)\nnsteps == 5\n"
        str_out = iw.substitute_variable_values(str_inp, n_steps=10)
        self.assertEqual(str_out, str_exp)
        str_out = iw.substitute_variable_values(str_inp, n_steps='10',
                                                strings_as_is=True)
        self.assertEqual(str_out, str_exp)
        str_inp = "N=(\n5000)\nnsteps == 5\n"
        str_exp = "N = 10\n\nnsteps == 5\n"
        str_out = iw.substitute_variable_values(str_inp, N=10,
                                                keep_original=False)
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
        path = pathlib.Path("a.py")
        # test substitutions
        str_inp = "import sys\nimport argparse"
        str_exp = f"import sys;sys.argv = ['{path.name}', '1', '2'];" + str_inp
        str_out, sys_argv = iw.set_cmd(str_inp, path, (1, 2))
        self.assertEqual(str_out, str_exp)
        self.assertEqual(sys_argv, [0, "test"])
        str_inp = "import argparse"
        str_exp = f"import sys;sys.argv = ['{path.name}', '1', '2'];" + str_inp
        str_out, sys_argv = iw.set_cmd(str_inp, path, ["1", 2])
        self.assertEqual(str_out, str_exp)
        self.assertEqual(sys_argv, [0, "test"])
        # test exceptions
        str_inp = "import re"
        self.assertRaises(AssertionError, iw.set_cmd, str_inp, path, (1, 2))
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

    def test_mock_es_visualization(self):
        statement = "import espressomd.visualization"
        expected = f"""
try:
    {statement}
except ImportError:
    import unittest.mock
    import espressomd
    espressomd.visualization = unittest.mock.MagicMock()
"""
        self.assertEqual(iw.mock_es_visualization(statement), expected[1:])

        statement = "import espressomd.visualization as test"
        expected = f"""
try:
    {statement}
except ImportError:
    import unittest.mock
    import espressomd
    test = unittest.mock.MagicMock()
"""
        self.assertEqual(iw.mock_es_visualization(statement), expected[1:])

        statement = "import espressomd.visualization, espressomd.visualization as test"
        expected = """
try:
    import espressomd.visualization
except ImportError:
    import unittest.mock
    import espressomd
    espressomd.visualization = unittest.mock.MagicMock()
try:
    import espressomd.visualization as test
except ImportError:
    import unittest.mock
    import espressomd
    test = unittest.mock.MagicMock()
"""
        self.assertEqual(iw.mock_es_visualization(statement), expected[1:])

        statement = "from espressomd import visualization"
        expected = f"""
try:
    {statement}
except ImportError:
    import unittest.mock
    import espressomd
    visualization = unittest.mock.MagicMock()
"""
        self.assertEqual(iw.mock_es_visualization(statement), expected[1:])

        statement = "from espressomd import visualization as test"
        expected = f"""
try:
    {statement}
except ImportError:
    import unittest.mock
    import espressomd
    test = unittest.mock.MagicMock()
"""
        self.assertEqual(iw.mock_es_visualization(statement), expected[1:])

        statement = "from espressomd.visualization import openGLLive"
        expected = f"""
try:
    {statement}
except ImportError:
    import unittest.mock
    import espressomd
    openGLLive = unittest.mock.MagicMock()
"""
        self.assertEqual(iw.mock_es_visualization(statement), expected[1:])

        statement = "from espressomd.visualization import openGLLive as test"
        expected = f"""
try:
    {statement}
except ImportError:
    import unittest.mock
    import espressomd
    test = unittest.mock.MagicMock()
"""
        self.assertEqual(iw.mock_es_visualization(statement), expected[1:])

        statement = "from espressomd.visualization import a as b, openGLLive"
        expected = """
try:
    from espressomd.visualization import a as b
except ImportError:
    import unittest.mock
    import espressomd
    b = unittest.mock.MagicMock()
try:
    from espressomd.visualization import openGLLive
except ImportError:
    import unittest.mock
    import espressomd
    openGLLive = unittest.mock.MagicMock()
"""
        self.assertEqual(iw.mock_es_visualization(statement), expected[1:])

        # test exceptions
        self.assertRaises(ValueError, iw.mock_es_visualization,
                          "from espressomd.visualization import *")

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

    def test_ipython_magics(self):
        lines = [
            '%matplotlib inline',
            '%%matplotlib notebook',
            'import matplotlib.pyplot as plt',
            '__doc__="%matplotlib inline"',
        ]
        code = '\n'.join(lines)
        code_protected = iw.protect_ipython_magics(code)
        code_protected_ref = f'\n{code}'.replace(
            '\n%', '\n#_IPYTHON_MAGIC_%')[1:]
        self.assertEqual(code_protected, code_protected_ref)
        code_deprotected = iw.deprotect_ipython_magics(code_protected)
        self.assertEqual(code_deprotected, code)

    def test_configure_and_import(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            path_in = pathlib.Path(temp_dir) / "sample.py"
            path_out = pathlib.Path(temp_dir) / "sample_test_processed.py"
            path_features = pathlib.Path(temp_dir) / "sample_impossible.py"

            # test importing a simple sample
            sys.argv.append("42")
            sys_argv_ref = list(sys.argv)
            path_in.write_text("""
import sys
from argparse import ArgumentParser
value = 42
argv = list(sys.argv)
import espressomd.visualization
""")
            sample, _ = iw.configure_and_import(
                path_in, move_to_script_dir=False, cmd_arguments=["TestCase"],
                gpu=False, script_suffix="test", value=43)
            self.assertEqual(sys.argv, sys_argv_ref)
            self.assertEqual(sample.argv, [path_in.name, "TestCase"])
            self.assertTrue(path_out.exists(), f"File {path_out} not found")
            self.assertEqual(sample.value, 43)
            self.assertIn(
                "espressomd.visualization = unittest.mock.MagicMock()",
                path_out.read_text())

            # test importing a sample that relies on features not compiled in
            inactive_features = set(
                espressomd.all_features()) - set(espressomd.features())
            if inactive_features:
                path_features.write_text(f"""
import espressomd
espressomd.assert_features({list(inactive_features)})
""")
                module, _ = iw.configure_and_import(path_features)
                self.assertIsInstance(module, ut.mock.MagicMock)

                # importing a valid sample is impossible after a failed import
                sample, _ = iw.configure_and_import(
                    path_in, move_to_script_dir=True, script_suffix="test",
                    cmd_arguments=["TestCase"], gpu=False, value=43)
                self.assertIsInstance(sample, ut.mock.MagicMock)


if __name__ == "__main__":
    ut.main()
