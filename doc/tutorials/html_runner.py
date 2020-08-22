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
"""
This script processes Jupyter notebooks. External Python scripts
can be inserted as new code cells (e.g. solutions to exercises).
The notebook may also be executed, if necessary with modified
global variables to reduce runtime. The processed notebook can
then be converted to HTML externally.
"""

import argparse

parser = argparse.ArgumentParser(description='Process Jupyter notebooks.',
                                 epilog=__doc__)
parser.add_argument('--input', type=str, nargs=1, required=True,
                    help='Path to the original Jupyter notebook')
parser.add_argument('--output', type=str, nargs=1,
                    help='Path to the processed Jupyter notebook')
parser.add_argument('--substitutions', nargs='*',
                    help='Variables to substitute')
parser.add_argument('--scripts', nargs='*',
                    help='Scripts to insert in new cells')
parser.add_argument('--execute', action='store_true',
                    help='Run the script')
args = parser.parse_args()

import nbformat
import re
import os
import ast
import sys
import uuid
sys.path.append('@CMAKE_SOURCE_DIR@/testsuite/scripts')
import importlib_wrapper as iw


def get_code_cells(nb):
    return [c['source'] for c in nb['cells'] if c['cell_type'] == 'code']


def set_code_cells(nb, new_cells):
    i = 0
    for c in nb['cells']:
        if c['cell_type'] == 'code':
            c['source'] = new_cells[i]
            i += 1


def add_cell_from_script(nb, filepath):
    """
    Create new code cell at the end of a notebook and populate it with
    the content of a script.
    """
    with open(filepath, encoding='utf-8') as f:
        code = f.read()
    # remove ESPResSo copyright header
    m = re.search('# Copyright \(C\) [\d\-,]+ The ESPResSo project\n.+?'
                  'If not, see <http://www\.gnu\.org/licenses/>\.\n', code, re.DOTALL)
    if m and all(x.startswith('#') for x in m.group(0).strip().split('\n')):
        code = re.sub('^(#\n)+', '', code.replace(m.group(0), ''), re.M)
    # strip first component in relative paths
    code = re.sub('(?<=[\'\"])\.\./', './', code)
    # create new cells
    filename = os.path.relpath(os.path.realpath(filepath))
    if len(filename) > len(filepath):
        filename = filepath
    cell_md = nbformat.v4.new_markdown_cell(source='Solution from ' + filename)
    nb['cells'].append(cell_md)
    cell_code = nbformat.v4.new_code_cell(source=code.strip())
    nb['cells'].append(cell_code)


def disable_plot_interactivity(nb):
    """
    Replace all occurrences of the magic command ``%matplotlib notebook``
    by ``%matplotlib inline``.
    """
    for cell in nb['cells']:
        if cell['cell_type'] == 'code' and 'matplotlib' in cell['source']:
            cell['source'] = re.sub('^%matplotlib +notebook',
                                    '%matplotlib inline',
                                    cell['source'], flags=re.M)


def split_matplotlib_cells(nb):
    """
    If a cell imports matplotlib, split the cell to keep the
    import statement separate from the code that uses matplotlib.
    This prevents a known bug in the Jupyter backend which causes
    the plot object to be represented as a string instead of a canvas
    when created in the cell where matplotlib is imported for the
    first time (https://github.com/jupyter/notebook/issues/3523).
    """
    for i in range(len(nb['cells']) - 1, -1, -1):
        cell = nb['cells'][i]
        if cell['cell_type'] == 'code' and 'matplotlib' in cell['source']:
            code = iw.protect_ipython_magics(cell['source'])
            # split cells after matplotlib imports
            mapping = iw.delimit_statements(code)
            tree = ast.parse(code)
            visitor = iw.GetMatplotlibPyplot()
            visitor.visit(tree)
            if visitor.matplotlib_first:
                code = iw.deprotect_ipython_magics(code)
                lines = code.split('\n')
                lineno_end = mapping[visitor.matplotlib_first]
                split_code = '\n'.join(lines[lineno_end:]).lstrip('\n')
                if split_code:
                    new_cell = nbformat.v4.new_code_cell(source=split_code)
                    nb['cells'].insert(i + 1, new_cell)
                lines = lines[:lineno_end]
                nb['cells'][i]['source'] = '\n'.join(lines).rstrip('\n')


def execute_notebook(nb, src, cell_separator):
    """
    Run the notebook in a python3 kernel. The ESPResSo visualizers are
    disabled to prevent the kernel from crashing and to allow running
    the notebook in a CI environment.
    """
    import nbconvert.preprocessors
    notebook_dirname = os.path.dirname(notebook_filepath)
    # disable OpenGL/Mayavi GUI
    src_no_gui = iw.mock_es_visualization(src)
    # update notebook with new code
    set_code_cells(nb, src_no_gui.split(cell_separator))
    # execute notebook
    ep = nbconvert.preprocessors.ExecutePreprocessor(
        timeout=20 * 60, kernel_name='python3')
    ep.preprocess(nb, {'metadata': {'path': notebook_dirname}})
    # restore notebook with code before the GUI removal step
    set_code_cells(nb, src.split(cell_separator))


notebook_filepath = args.input[0]
if args.output:
    notebook_filepath_edited = args.output[0]
else:
    notebook_filepath_edited = notebook_filepath + '~'

# parse original notebook
with open(notebook_filepath, encoding='utf-8') as f:
    nb = nbformat.read(f, as_version=4)

# add new cells containing the solutions
if args.scripts:
    for filepath in args.scripts:
        add_cell_from_script(nb, filepath)

# disable plot interactivity
disable_plot_interactivity(nb)

# guard against a jupyter bug involving matplotlib
split_matplotlib_cells(nb)

if args.substitutions or args.execute:
    # substitute global variables
    cell_separator = '\n##{}\n'.format(uuid.uuid4().hex)
    src = cell_separator.join(get_code_cells(nb))
    new_values = args.substitutions or []
    parameters = dict(x.split('=', 1) for x in new_values)
    src = iw.substitute_variable_values(src, strings_as_is=True,
                                        keep_original=False, **parameters)
    set_code_cells(nb, src.split(cell_separator))
    if args.execute:
        execute_notebook(nb, src, cell_separator)

# write edited notebook
with open(notebook_filepath_edited, 'w', encoding='utf-8') as f:
    nbformat.write(nb, f)
