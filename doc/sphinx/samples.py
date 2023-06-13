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

import os
import ast
import textwrap

samples_dir = '/home/thilo/code/espresso2/espresso/samples/'


def get_docstring(filenames):
    docstrings = {}
    for filename in filenames:
        with open(samples_dir + filename) as f:
            src = f.read()
        tree = ast.parse(src)
        module = next(ast.walk(tree))
        docstrings[filename] = ast.get_docstring(module)
    return docstrings


# extract docstrings
samples = [x for x in os.listdir(samples_dir) if x.endswith('.py')]
samples += ['immersed_boundary/sampleImmersedBoundary.py',
            'object_in_fluid/motivation.py']
docstrings = get_docstring(samples)


# write documentation
sphinx_tpl = '* :file:`{}`\n{}\n'
with open('/home/thilo/code/espresso2/espresso/doc/sphinx/samples.rst', 'w') as f:
    for filename in sorted(docstrings, key=lambda x: x.lower()):
        docstring = (docstrings[filename] or '').replace('ESPResSo', '|es|')
        paragraph = textwrap.indent(docstring, prefix='    ')
        f.write(sphinx_tpl.format(filename, paragraph))
