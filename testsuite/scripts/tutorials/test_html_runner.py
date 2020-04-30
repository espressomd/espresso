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
import os
import sys
import subprocess
import nbformat


class HtmlRunner(ut.TestCase):
    """
    Test the :file:`doc/tutorials/html_runner.py` script. A new Jupyter
    notebook and a new python script are created, and both are supplied to
    html_runner.py, which will include the python script in a new code cell,
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
    from espressomd.visualization_opengl import openglLive
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

    def test_html_wrapper(self):
        f_input = '@CMAKE_CURRENT_BINARY_DIR@/test_html_runner_notebook.ipynb'
        f_output = '@CMAKE_CURRENT_BINARY_DIR@/test_html_runner_notebook.run.ipynb'
        f_script = '@CMAKE_CURRENT_BINARY_DIR@/test_html_runner_script.py'
        # setup
        if os.path.isfile(f_output):
            os.remove(f_output)
        with open(f_script, 'w') as f:
            f.write('global_var = 5')
        with open(f_input, 'w', encoding='utf-8') as f:
            nb = nbformat.v4.new_notebook(metadata=self.nb_metadata)
            cell_md = nbformat.v4.new_markdown_cell(source=self.cell_md_src)
            nb['cells'].append(cell_md)
            cell_code = nbformat.v4.new_code_cell(source=self.cell_py_src)
            nb['cells'].append(cell_code)
            nbformat.write(nb, f)
        # run command
        cmd = ['@CMAKE_BINARY_DIR@/pypresso',
               '@CMAKE_BINARY_DIR@/doc/tutorials/html_runner.py',
               '--input', f_input,
               '--output', f_output,
               '--scripts', f_script,
               '--substitutions', 'global_var=20']
        completedProc = subprocess.run(cmd)
        # check the command ran without any error
        self.assertEqual(completedProc.returncode, 0, 'non-zero return code')
        self.assertTrue(os.path.isfile(f_output), f_output + ' not created')
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
            from espressomd.visualization_opengl import openglLive  # pylint: disable=unused-import
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
                         'Solution from test_html_runner_script.py')
        self.assertEqual(nb_output['cells'][4]['cell_type'], 'code')
        self.assertEqual(nb_output['cells'][4]['source'], 'global_var = 20')


if __name__ == "__main__":
    ut.main()
