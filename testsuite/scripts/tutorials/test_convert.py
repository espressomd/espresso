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
import sys
import nbformat
import traceback
import pathlib

sys.path.insert(0, '@CMAKE_BINARY_DIR@/doc/tutorials')
import convert


def skipIfMissingModules(x): return x


try:
    import yaml  # pylint: disable=unused-import
    import autopep8  # pylint: disable=unused-import
except ImportError:
    skipIfMissingModules = ut.skip(
        "Python modules 'yaml' or 'autopep8' not available, skipping test!")
else:
    def skipIfMissingModules(x): return x


class HtmlRunner(ut.TestCase):
    """
    Test the :file:`doc/tutorials/convert.py` script. A new Jupyter
    notebook and a new python script are created, and both are supplied to
    convert.py, which will include the python script in a new code cell,
    substitute global variables, run the code and then save the result in
    a new notebook. The input notebook contains IPython magic commands and
    imports matplotlib and an ESPResSo visualizer, all of which require
    special treatment.
    """

    cell_md_src = '''
Cannot plot in the same cell where `matplotlib` is imported for the first
time, so CI/CD needs to split the code cell after the first matplotlib
import statement. IPython magic commands `%matplotlib` must be set to inline.
'''.strip()

    cell_py_src = '''
import numpy as np
%matplotlib notebook
import matplotlib as mpl  # split here
import matplotlib.pyplot as plt  # don't split
try:
    from espressomd.visualization import openglLive
except ImportError:
    mpl.use('Agg')  # running in CI without graphical output
plt.ion()
global_var = 5
plt.plot([1, 2], [3, global_var])
plt.show()
'''.strip()

    nb_metadata = {
        "kernelspec": {
            "display_name": "Python 3",
            "language": "python",
            "name": "python3"},
        "language_info": {
            "codemirror_mode": {"name": "ipython", "version": 3},
            "file_extension": ".py",
            "mimetype": "text/x-python",
            "name": "python",
            "nbconvert_exporter": "python",
            "pygments_lexer": "ipython3",
            "version": ".".join(map(str, sys.version_info[:3]))}
    }

    def run_command(self, cmd, output=None):
        error_msg = ("Could not run @CMAKE_BINARY_DIR@/pypresso "
                     "@CMAKE_BINARY_DIR@/doc/tutorials/convert.py "
                     f"{' '.join(cmd)}")
        try:
            args = convert.parser.parse_args(cmd)
            args.callback(args)
        except BaseException:
            traceback.print_exc()
            self.fail(error_msg)
        if output is not None:
            self.assertTrue(output.exists(), f"File {output} not created")

    def test_html_wrapper(self):
        root = pathlib.Path("@CMAKE_CURRENT_BINARY_DIR@")
        f_input = root / "test_convert_notebook.ipynb"
        f_output = root / "test_convert_notebook.run.ipynb"
        f_script = root / "test_convert_script.py"
        # setup
        f_output.unlink(missing_ok=True)
        f_script.write_text('global_var = 5')
        with open(f_input, 'w', encoding='utf-8') as f:
            nb = nbformat.v4.new_notebook(metadata=self.nb_metadata)
            cell_md = nbformat.v4.new_markdown_cell(source=self.cell_md_src)
            nb['cells'].append(cell_md)
            cell_code = nbformat.v4.new_code_cell(source=self.cell_py_src)
            nb['cells'].append(cell_code)
            nbformat.write(nb, f)
        # run command and check for errors
        cmd = ['ci',
               '--input', str(f_input),
               '--output', str(f_output),
               '--scripts', str(f_script),
               '--substitutions', 'global_var=20',
               '--execute']
        self.run_command(cmd, f_output)
        # read processed notebook
        with open(f_output, encoding='utf-8') as f:
            nb_output = nbformat.read(f, as_version=4)
        # the first Markdown cell must be identical
        self.assertEqual(nb_output['cells'][0]['cell_type'], 'markdown')
        self.assertEqual(nb_output['cells'][0]['source'], self.cell_md_src)
        # the first Python cell must have been split
        self.assertEqual(nb_output['cells'][1]['cell_type'], 'code')
        lines = (self.cell_py_src
                 .replace('%matplotlib notebook', '%matplotlib inline')
                 .replace('global_var = 5', 'global_var = 20')
                 ).split('\n')
        self.assertEqual(nb_output['cells'][1]['source'], '\n'.join(lines[:3]))
        self.assertEqual(nb_output['cells'][2]['source'], '\n'.join(lines[3:]))
        # the cell should have produced a plot
        graphical_plots = True
        try:
            from espressomd.visualization import openglLive  # pylint: disable=unused-import
        except ImportError:
            graphical_plots = False  # running in CI without graphical output
        if graphical_plots:
            outputs = nb_output['cells'][2]['outputs']
            self.assertTrue(outputs, 'cell has no output')
            self.assertIn('image/png', outputs[0]['data'])
            self.assertGreater(len(outputs[0]['data']['image/png']), 6000)
        # check the external script was correctly inserted
        self.assertEqual(nb_output['cells'][3]['cell_type'], 'markdown')
        self.assertEqual(nb_output['cells'][3]['source'],
                         'Solution from test_convert_script.py')
        self.assertEqual(nb_output['cells'][4]['cell_type'], 'code')
        self.assertEqual(nb_output['cells'][4]['source'], 'global_var = 20')

    def test_exercise2_plugin(self):
        root = pathlib.Path("@CMAKE_CURRENT_BINARY_DIR@")
        f_input = root / "test_convert_exercise2.ipynb"
        f_output = root / "test_convert_exercise2.run.ipynb"
        # setup
        f_output.unlink(missing_ok=True)
        with open(f_input, 'w', encoding='utf-8') as f:
            nb = nbformat.v4.new_notebook(metadata=self.nb_metadata)
            # question with 2 answers and an empty cell
            cell_md = nbformat.v4.new_markdown_cell(source='Question 1')
            cell_md['metadata']['solution2_first'] = True
            cell_md['metadata']['solution2'] = 'shown'
            nb['cells'].append(cell_md)
            code = '```python\n1\n```'
            cell_md = nbformat.v4.new_markdown_cell(source=code)
            cell_md['metadata']['solution2'] = 'shown'
            cell_md['metadata']['key'] = 'value'
            nb['cells'].append(cell_md)
            cell_md = nbformat.v4.new_markdown_cell(source='1b')
            cell_md['metadata']['solution2'] = 'shown'
            nb['cells'].append(cell_md)
            cell_code = nbformat.v4.new_code_cell(source='')
            nb['cells'].append(cell_code)
            # question with 1 answer and a non-empty cell
            cell_md = nbformat.v4.new_markdown_cell(source='Question 2')
            cell_md['metadata']['solution2_first'] = True
            cell_md['metadata']['solution2'] = 'hidden'
            nb['cells'].append(cell_md)
            code = '```python\n2\nimport matplotlib.pyplot\nglobal_var = 5\n```'
            cell_md = nbformat.v4.new_markdown_cell(source=code)
            cell_md['metadata']['solution2'] = 'hidden'
            nb['cells'].append(cell_md)
            cell_code = nbformat.v4.new_code_cell(source='3')
            nb['cells'].append(cell_code)
            nbformat.write(nb, f)
        # run command and check for errors
        cmd = ['ci',
               '--input', str(f_input),
               '--output', str(f_output),
               '--substitutions', 'global_var=20',
               '--exercise2', '--remove-empty-cells']
        self.run_command(cmd, f_output)
        # read processed notebook
        with open(f_output, encoding='utf-8') as f:
            nb_output = nbformat.read(f, as_version=4)
        # check cells
        cells = iter(nb_output['cells'])
        cell = next(cells)
        self.assertEqual(cell['cell_type'], 'markdown')
        self.assertEqual(cell['source'], 'Question 1')
        cell = next(cells)
        self.assertEqual(cell['cell_type'], 'code')
        self.assertEqual(cell['source'], '1')
        self.assertEqual(cell['metadata']['key'], 'value')
        cell = next(cells)
        self.assertEqual(cell['cell_type'], 'markdown')
        self.assertEqual(cell['source'], '1b')
        cell = next(cells)
        self.assertEqual(cell['cell_type'], 'markdown')
        self.assertEqual(cell['source'], 'Question 2')
        cell = next(cells)
        self.assertEqual(cell['cell_type'], 'code')
        self.assertEqual(cell['source'], '2\nimport matplotlib.pyplot')
        cell = next(cells)
        self.assertEqual(cell['cell_type'], 'code')
        self.assertEqual(cell['source'], 'global_var = 20')
        cell = next(cells)
        self.assertEqual(cell['cell_type'], 'code')
        self.assertEqual(cell['source'], '3')
        self.assertEqual(next(cells, 'EOF'), 'EOF')

    def test_exercise2_conversion(self):
        root = pathlib.Path("@CMAKE_CURRENT_BINARY_DIR@")
        f_input = root / "test_convert_exercise2_conversion.ipynb"
        # setup
        with open(f_input, 'w', encoding='utf-8') as f:
            nb = nbformat.v4.new_notebook(metadata=self.nb_metadata)
            # question and code answer
            cell_md = nbformat.v4.new_markdown_cell(source='Question 1')
            cell_md['metadata']['solution2_first'] = True
            cell_md['metadata']['solution2'] = 'hidden'
            nb['cells'].append(cell_md)
            code = '```python\n1\n```'
            cell_md = nbformat.v4.new_markdown_cell(source=code)
            cell_md['metadata']['solution2'] = 'hidden'
            cell_md['metadata']['key'] = 'value'
            nb['cells'].append(cell_md)
            nbformat.write(nb, f)
        # run command and check for errors
        cmd = ['exercise2', '--to-py', str(f_input)]
        self.run_command(cmd, f_input)
        # read processed notebook
        with open(f_input, encoding='utf-8') as f:
            nb_output = nbformat.read(f, as_version=4)
        # check cells
        cells = iter(nb_output['cells'])
        cell = next(cells)
        self.assertEqual(cell['cell_type'], 'markdown')
        self.assertEqual(cell['source'], 'Question 1')
        cell = next(cells)
        self.assertEqual(cell['cell_type'], 'code')
        self.assertEqual(cell['source'], '1')
        self.assertEqual(cell['metadata']['solution2'], 'shown')
        self.assertEqual(cell['metadata']['key'], 'value')
        self.assertEqual(next(cells, 'EOF'), 'EOF')
        # run command and check for errors
        cmd = ['exercise2', '--to-md', str(f_input)]
        self.run_command(cmd, f_input)
        # read processed notebook
        with open(f_input, encoding='utf-8') as f:
            nb_output = nbformat.read(f, as_version=4)
        # check cells
        cells = iter(nb_output['cells'])
        cell = next(cells)
        self.assertEqual(cell['cell_type'], 'markdown')
        self.assertEqual(cell['source'], 'Question 1')
        cell = next(cells)
        self.assertEqual(cell['cell_type'], 'markdown')
        self.assertEqual(cell['source'], '```python\n1\n```')
        self.assertEqual(cell['metadata']['solution2'], 'hidden')
        self.assertEqual(cell['metadata']['key'], 'value')
        self.assertEqual(next(cells, 'EOF'), 'EOF')

    @skipIfMissingModules
    def test_exercise2_autopep8(self):
        root = pathlib.Path("@CMAKE_CURRENT_BINARY_DIR@")
        f_input = root / "test_convert_exercise2_autopep8.ipynb"
        # setup
        with open(f_input, 'w', encoding='utf-8') as f:
            nb = nbformat.v4.new_notebook(metadata=self.nb_metadata)
            # question and code answer
            cell_md = nbformat.v4.new_markdown_cell(source='Question 1')
            cell_md['metadata']['solution2_first'] = True
            cell_md['metadata']['solution2'] = 'hidden'
            nb['cells'].append(cell_md)
            code = '```python\n\nif 1: #comment\n  print( [5+1,4])\n\n```'
            cell_md = nbformat.v4.new_markdown_cell(source=code)
            cell_md['metadata']['solution2'] = 'hidden'
            nb['cells'].append(cell_md)
            nbformat.write(nb, f)
        # run command and check for errors
        cmd = ['exercise2', '--pep8', str(f_input)]
        self.run_command(cmd)
        # read processed notebook
        with open(f_input, encoding='utf-8') as f:
            nb_output = nbformat.read(f, as_version=4)
        # check cells
        cells = iter(nb_output['cells'])
        cell = next(cells)
        self.assertEqual(cell['cell_type'], 'markdown')
        self.assertEqual(cell['source'], 'Question 1')
        cell = next(cells)
        self.assertEqual(cell['cell_type'], 'markdown')
        self.assertEqual(
            cell['source'],
            '```python\nif 1:  # comment\n    print([5 + 1, 4])\n```')
        self.assertEqual(cell['metadata']['solution2'], 'hidden')
        self.assertEqual(next(cells, 'EOF'), 'EOF')


if __name__ == "__main__":
    ut.main()
