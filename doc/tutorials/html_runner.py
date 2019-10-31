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
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor
import re
import os
import sys
import uuid
import argparse
sys.path.append('@CMAKE_SOURCE_DIR@/testsuite/scripts')
from importlib_wrapper import substitute_variable_values, mock_es_visualization

parser = argparse.ArgumentParser(description='Process IPython notebooks.')
parser.add_argument('--input', type=str,
                    help='Path to the original IPython notebook')
parser.add_argument('--output', type=str, nargs='?',
                    help='Path to the processed IPython notebook')
parser.add_argument('--substitutions', nargs='*',
                    help='Variables to substitute')
parser.add_argument('--scripts', nargs='*',
                    help='Scripts to insert in new cells')
args = parser.parse_args()


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
    # if matplotlib is used in this script, split cell to keep the import
    # statement separate and avoid a know bug in the Jupyter backend which
    # causes the plot object to be represented as a string instead of a
    # canvas when created in the cell where matplotlib is imported
    # (https://github.com/jupyter/notebook/issues/3523)
    if 'import matplotlib' in code:
        cells_code = re.split('^((?:|.*\n)import matplotlib.*?)\n', code,
                              maxsplit=1, flags=re.DOTALL)[1:]
    else:
        cells_code = [code]
    # create new cells
    cell_md = nbformat.v4.new_markdown_cell(source='Solution from ' + filepath)
    nb['cells'].append(cell_md)
    for cell_code in cells_code:
        cell_code = nbformat.v4.new_code_cell(source=cell_code.strip())
        nb['cells'].append(cell_code)

# substitute global variables and disable OpenGL/Mayavi GUI
cell_separator = '\n##{}\n'.format(uuid.uuid4().hex)
src = cell_separator.join(get_code_cells(nb))
parameters = dict(x.split('=', 1) for x in new_values)
src = substitute_variable_values(src, strings_as_is=True, keep_original=False,
                                 **parameters)
src_no_gui = mock_es_visualization(src)

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
