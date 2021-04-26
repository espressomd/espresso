# Copyright (C) 2020 The ESPResSo project
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

import sys
import nbformat
import importlib
import unittest as ut

sys.path.insert(0, '@CMAKE_SOURCE_DIR@/maintainer/CI')
module = importlib.import_module('jupyter_warnings')


class Test(ut.TestCase):

    cell_md_src = '''
ignored: http://espressomd.org/wordpress/documentation/
ignored: http://espressomd.org/wordpress/unknown_folder/
valid:   https://espressomd.github.io/doc/index.html
valid:   https://espressomd.github.io/doc/index.html#python-module-documentation
valid:   https://espressomd.github.io/doc/index.html?highlight=highlander#python-module-documentation
valid:   https://espressomd.github.io/doc/index.html?highlight=highlander
invalid: https://espressomd.github.io/doc/index.html#unknown_anchor
invalid: https://espressomd.github.io/doc/unknown_file.html
invalid: [footnote 1](#unknown-footnote-1)
'''

    def test_detect_invalid_urls(self):
        nb = nbformat.v4.new_notebook()
        cell_md = nbformat.v4.new_markdown_cell(source=self.cell_md_src)
        nb['cells'].append(cell_md)
        ref_issues = [
            'https://espressomd.github.io/doc/index.html has no anchor "unknown_anchor"',
            'https://espressomd.github.io/doc/unknown_file.html does not exist',
            'notebook has no anchor "unknown-footnote-1"'
        ]
        issues = module.detect_invalid_urls(nb, '@CMAKE_BINARY_DIR@')
        self.assertEqual(issues, ref_issues)


if __name__ == "__main__":
    ut.main()
