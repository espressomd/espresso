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
This script runs Jupyter notebooks and writes the output to new notebooks.
Global variables can be edited to reduce runtime. External Python scripts
can be inserted as new code cells (e.g. solutions to exercises).
The output notebooks can then be converted to HTML externally.
"""

import argparse

parser = argparse.ArgumentParser(description='Process Jupyter notebooks.',
                                 epilog=__doc__)
parser.add_argument('--input', type=str,
                    help='Path to the original Jupyter notebook')
parser.add_argument('--output', type=str, nargs='?',
                    help='Path to the processed Jupyter notebook')
parser.add_argument('--substitutions', nargs='*',
                    help='Variables to substitute')
parser.add_argument('--scripts', nargs='*',
                    help='Scripts to insert in new cells')
args = parser.parse_args()

import nbformat
from nbconvert.preprocessors import ExecutePreprocessor
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


notebook_filepath = args.input
notebook_filepath_edited = args.output or args.input + '~'
notebook_dirname = os.path.dirname(notebook_filepath)
new_values = args.substitutions or []
new_cells = args.scripts or []

# parse original notebook
with open(notebook_filepath, encoding='utf-8') as f:
    nb = nbformat.read(f, as_version=4)

# add new cells containing the solutions
for filepath in new_cells:
    with open(filepath, encoding='utf-8') as f:
        code = f.read()
    # remove ESPResSo copyright header
    m = re.search('# Copyright \(C\) \d+(?:-\d+)? The ESPResSo project\n.+?'
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


# disable plot interactivity
for i in range(len(nb['cells'])):
    cell = nb['cells'][i]
    if cell['cell_type'] == 'code' and 'matplotlib' in cell['source']:
        cell['source'] = re.sub('^%matplotlib +notebook', '%matplotlib inline',
                                cell['source'], flags=re.M)


# if matplotlib is used in this script, split cell to keep the import
# statement separate and avoid a known bug in the Jupyter backend which
# causes the plot object to be represented as a string instead of a
# canvas when created in the cell where matplotlib is imported for the
# first time (https://github.com/jupyter/notebook/issues/3523)
for i in range(len(nb['cells'])):
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
            break

# substitute global variables and disable OpenGL/Mayavi GUI
cell_separator = '\n##{}\n'.format(uuid.uuid4().hex)
src = cell_separator.join(get_code_cells(nb))
parameters = dict(x.split('=', 1) for x in new_values)
src = iw.substitute_variable_values(src, strings_as_is=True,
                                    keep_original=False, **parameters)
src_no_gui = iw.mock_es_visualization(src)

# update notebook with new code
set_code_cells(nb, src_no_gui.split(cell_separator))

# execute notebook
ep = ExecutePreprocessor(timeout=20 * 60, kernel_name='python3')
ep.preprocess(nb, {'metadata': {'path': notebook_dirname}})

# restore notebook with code before the GUI removal step
set_code_cells(nb, src.split(cell_separator))

# write edited notebook
with open(notebook_filepath_edited, 'w', encoding='utf-8') as f:
    nbformat.write(nb, f)
