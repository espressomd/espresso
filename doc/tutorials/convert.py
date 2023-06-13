#
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
#
"""
This script processes Jupyter notebooks. External Python scripts
can be inserted as new code cells (e.g. solutions to exercises).
Hidden solutions from the ``exercise2`` plugin can be converted
to code cells. The notebook may also be executed, if necessary
with modified global variables to reduce runtime. The processed
notebook can then be converted to HTML externally.
"""

import argparse
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
    m = re.search('# Copyright \\(C\\) [\\d\\-,]+ The ESPResSo project\n.+?'
                  'If not, see <http://www\\.gnu\\.org/licenses/>\\.\n', code, re.DOTALL)
    if m and all(x.startswith('#') for x in m.group(0).strip().split('\n')):
        code = re.sub('^(#\n)+', '', code.replace(m.group(0), ''), re.M)
    # strip first component in relative paths
    code = re.sub(r'(?<=[\'\"])\.\./', './', code)
    # create new cells
    filename = os.path.relpath(os.path.realpath(filepath))
    if len(filename) > len(filepath):
        filename = filepath
    cell_md = nbformat.v4.new_markdown_cell(source='Solution from ' + filename)
    nb['cells'].append(cell_md)
    cell_code = nbformat.v4.new_code_cell(source=code.strip())
    nb['cells'].append(cell_code)


def remove_empty_cells(nb):
    for i in range(len(nb['cells']) - 1, 0, -1):
        cell = nb['cells'][i]
        if cell['source'].strip() == '':
            nb['cells'].pop(i)


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
                    if 'id' not in cell and 'id' in new_cell:
                        del new_cell['id']
                    nb['cells'].insert(i + 1, new_cell)
                lines = lines[:lineno_end]
                nb['cells'][i]['source'] = '\n'.join(lines).rstrip('\n')


def convert_exercise2_to_code(nb):
    """
    Walk through the notebook cells and convert exercise2 Markdown cells
    containing fenced python code to exercise2 code cells.
    """
    for i, cell in enumerate(nb['cells']):
        if 'solution2' in cell['metadata']:
            cell['metadata']['solution2'] = 'shown'
        # convert solution markdown cells into code cells
        if cell['cell_type'] == 'markdown' and 'solution2' in cell['metadata'] \
                and 'solution2_first' not in cell['metadata']:
            lines = cell['source'].strip().split('\n')
            if lines[0].strip() == '```python' and lines[-1].strip() == '```':
                source = '\n'.join(lines[1:-1]).strip()
                nb['cells'][i] = nbformat.v4.new_code_cell(source=source)
                nb['cells'][i]['metadata'] = cell['metadata']
                nb['cells'][i]['metadata']['solution2'] = 'shown'
                if 'id' in nb['cells'][i]:
                    del nb['cells'][i]['id']


def convert_exercise2_to_markdown(nb):
    """
    Walk through the notebook cells and convert exercise2 Python cells
    to exercise2 Markdown cells using a fenced code block.
    """
    for i, cell in enumerate(nb['cells']):
        if 'solution2' in cell['metadata']:
            cell['metadata']['solution2'] = 'hidden'
        # convert solution code cells into markdown cells
        if cell['cell_type'] == 'code' and 'solution2' in cell['metadata']:
            content = '```python\n' + cell['source'] + '\n```'
            nb['cells'][i] = nbformat.v4.new_markdown_cell(source=content)
            nb['cells'][i]['metadata'] = cell['metadata']
            nb['cells'][i]['metadata']['solution2'] = 'hidden'
            if 'id' in nb['cells'][i]:
                del nb['cells'][i]['id']


def convert_exercise2_to_jupyterlab(nb):
    """
    Walk through the notebook cells and convert exercise2 Markdown cells
    containing fenced python code to a JupyterLab-compatible format.
    As of 2022, there is no equivalent of exercise2 for JupyterLab
    ([chart](https://jupyterlab-contrib.github.io/migrate_from_classical.html)),
    but a similar effect can be obtained with basic HTML.

    This also converts a notebook to Notebook Format 4.5. ESPResSo notebooks
    cannot be saved in 4.5 format since both Jupyter Notebook and JupyterLab
    overwrite the cell ids with random strings after each save, which is a
    problem for version control. The notebooks need to be converted to the
    4.5 format to silence JSON parser errors in JupyterLab.
    """
    jupyterlab_tpl = """\
<details {0} style="margin: 0.8em 4em;">\
<summary style="cursor: pointer; margin-left: -3em;">Show solution</summary>
<div style="margin-bottom: 2em;"></div>

{1}
<div style="margin-top: 2em;"></div>
</details>\
"""
    for i, cell in enumerate(nb['cells']):
        # convert solution markdown cells into code cells
        if cell['cell_type'] == 'markdown' and 'solution2' in cell['metadata'] \
                and 'solution2_first' not in cell['metadata']:
            lines = cell['source'].strip().split('\n')
            shown = 'open=""' if cell['metadata']['solution2'] == 'shown' else ''
            if lines[0].strip() == '```python' and lines[-1].strip() == '```':
                source = jupyterlab_tpl.format(shown, '\n'.join(lines).strip())
                nb['cells'][i] = nbformat.v4.new_markdown_cell(source=source)
        # convert cell to notebook format 4.5
        if 'id' not in cell:
            cell = uuid.uuid4().hex[:8]

    # change to notebook format 4.5
    current_version = (nb['nbformat'], nb['nbformat_minor'])
    assert current_version >= (4, 0)
    if current_version < (4, 5):
        nb['nbformat_minor'] = 5


def convert_exercise2_to_vscode_jupyter(nb):
    """
    Walk through the notebook cells and convert exercise2 Markdown cells
    containing fenced python code to a VS Code Jupyter-compatible format.
    As of 2022, there is no equivalent of exercise2 for VS Code Jupyter.
    """
    convert_exercise2_to_jupyterlab(nb)


def apply_autopep8(nb):
    import yaml
    import autopep8

    def get_autopep8_options():
        options = {'aggressive': 0, 'ignore': [], 'max_line_length': 120}
        with open('@CMAKE_SOURCE_DIR@/.pre-commit-config.yaml') as f:
            pre_config = yaml.safe_load(f)
        for repo in pre_config['repos']:
            for hook in repo['hooks']:
                if hook['id'] == 'autopep8':
                    for arg in hook['args']:
                        if arg == '--aggressive':
                            options['aggressive'] += 1
                        elif arg.startswith('--ignore='):
                            options['ignore'] = arg.split('=', 1)[0].split(',')
                    return options
        return options

    pep8_opts = get_autopep8_options()
    for cell in nb['cells']:
        source = None
        header = ''
        footer = ''
        if cell['cell_type'] == 'code':
            source = cell['source']
        elif cell['cell_type'] == 'markdown' and 'solution2' in cell['metadata']:
            lines = cell['source'].strip().split('\n')
            if lines[0].strip() == '```python' and lines[-1].strip() == '```':
                source = '\n'.join(lines[1:-1])
                header = lines[0] + '\n'
                footer = '\n' + lines[-1]
        if source is not None:
            source = autopep8.fix_code(source, options=pep8_opts).strip()
            cell['source'] = header + source + footer


def execute_notebook(nb, src, cell_separator, notebook_filepath):
    """
    Run the notebook in a python3 kernel. The ESPResSo visualizers are
    disabled to prevent the kernel from crashing and to allow running
    the notebook in a CI environment.
    """
    import nbconvert.preprocessors
    notebook_dirname = os.path.dirname(notebook_filepath)
    # disable OpenGL GUI
    src_no_gui = iw.mock_es_visualization(src)
    # update notebook with new code
    set_code_cells(nb, src_no_gui.split(cell_separator))
    # execute notebook
    ep = nbconvert.preprocessors.ExecutePreprocessor(
        timeout=20 * 60, kernel_name='python3')
    ep.preprocess(nb, {'metadata': {'path': notebook_dirname}})
    # restore notebook with code before the GUI removal step
    set_code_cells(nb, src.split(cell_separator))


def handle_ci_case(args):
    notebook_filepath = args.input
    if args.output:
        notebook_filepath_edited = args.output
    else:
        notebook_filepath_edited = notebook_filepath + '~'

    # parse original notebook
    with open(notebook_filepath, encoding='utf-8') as f:
        nb = nbformat.read(f, as_version=4)

    # add new cells containing the solutions
    if args.scripts:
        for filepath in args.scripts:
            add_cell_from_script(nb, filepath)

    # convert solution cells to code cells
    if args.exercise2:
        convert_exercise2_to_code(nb)

    # remove empty cells (e.g. those below exercise2 cells)
    if args.remove_empty_cells:
        remove_empty_cells(nb)

    # disable plot interactivity
    disable_plot_interactivity(nb)

    # guard against a jupyter bug involving matplotlib
    split_matplotlib_cells(nb)

    if args.substitutions or args.execute:
        # substitute global variables
        cell_separator = f'\n##{uuid.uuid4().hex}\n'
        src = cell_separator.join(get_code_cells(nb))
        new_values = args.substitutions or []
        parameters = dict(x.split('=', 1) for x in new_values)
        src = iw.substitute_variable_values(src, strings_as_is=True,
                                            keep_original=False, **parameters)
        set_code_cells(nb, src.split(cell_separator))

    if args.execute:
        execute_notebook(nb, src, cell_separator, args.input)

    # write edited notebook
    with open(notebook_filepath_edited, 'w', encoding='utf-8') as f:
        nbformat.write(nb, f, version=nbformat.NO_CONVERT)


def handle_exercise2_case(args):
    # parse original notebook
    with open(args.input, encoding='utf-8') as f:
        nb = nbformat.read(f, as_version=4)

    if args.to_md:
        convert_exercise2_to_markdown(nb)
    elif args.to_jupyterlab:
        convert_exercise2_to_jupyterlab(nb)
    elif args.to_vscode_jupyter:
        convert_exercise2_to_vscode_jupyter(nb)
    elif args.to_py:
        convert_exercise2_to_code(nb)
    elif args.pep8:
        convert_exercise2_to_code(nb)
        apply_autopep8(nb)
        convert_exercise2_to_markdown(nb)
    elif args.remove_empty_cells:
        remove_empty_cells(nb)

    # write edited notebook
    with open(args.input, 'w', encoding='utf-8') as f:
        nbformat.write(nb, f, version=nbformat.NO_CONVERT)


parser = argparse.ArgumentParser(description='Process Jupyter notebooks.',
                                 epilog=__doc__)
subparsers = parser.add_subparsers(help='Submodules')
# CI module
parser_ci = subparsers.add_parser(
    'ci', help='module for CI (variable substitution, code execution, etc.)')
parser_ci.add_argument('--input', type=str, required=True,
                       help='path to the original Jupyter notebook')
parser_ci.add_argument('--output', type=str,
                       help='path to the processed Jupyter notebook')
parser_ci.add_argument('--substitutions', nargs='*',
                       help='variables to substitute')
parser_ci.add_argument('--scripts', nargs='*',
                       help='scripts to insert in new cells')
parser_ci.add_argument('--exercise2', action='store_true',
                       help='convert exercise2 solutions into code cells')
parser_ci.add_argument('--remove-empty-cells', action='store_true',
                       help='remove empty cells')
parser_ci.add_argument('--execute', action='store_true',
                       help='run the notebook')
parser_ci.set_defaults(callback=handle_ci_case)
# exercise2 module
parser_exercise2 = subparsers.add_parser(
    'exercise2', help='module for exercise2 conversion (Markdown <-> Python)')
parser_exercise2.add_argument('input', type=str, help='path to the Jupyter '
                              'notebook (in-place conversion)')
group_exercise2 = parser_exercise2.add_mutually_exclusive_group(required=True)
group_exercise2.add_argument('--to-md', action='store_true',
                             help='convert solution cells to Markdown')
group_exercise2.add_argument('--to-jupyterlab', action='store_true',
                             help='convert solution cells to JupyterLab')
group_exercise2.add_argument('--to-vscode-jupyter', action='store_true',
                             help='convert solution cells to VS Code Jupyter')
group_exercise2.add_argument('--to-py', action='store_true',
                             help='convert solution cells to Python')
group_exercise2.add_argument('--pep8', action='store_true',
                             help='apply autopep8 formatting')
group_exercise2.add_argument('--remove-empty-cells', action='store_true',
                             help='remove empty cells')
parser_exercise2.set_defaults(callback=handle_exercise2_case)


if __name__ == "__main__":
    args = parser.parse_args()
    args.callback(args)
