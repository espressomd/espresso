#
# Copyright (C) 2020-2022 The ESPResSo project
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

import sys
import nbformat
import nbconvert
import importlib
import unittest as ut

sys.path.insert(0, '@CMAKE_SOURCE_DIR@/maintainer/CI')
module = importlib.import_module('jupyter_warnings')


class Test(ut.TestCase):

    cell_md_src = '''
ignored: https://espressomd.org/wordpress/documentation/
ignored: https://espressomd.org/wordpress/unknown_folder/
valid:   https://espressomd.github.io/doc/index.html
valid:   https://espressomd.github.io/doc/index.html#python-modules
valid:   https://espressomd.github.io/doc/index.html?highlight=highlander#python-modules
valid:   https://espressomd.github.io/doc/index.html?highlight=highlander
invalid: https://espressomd.github.io/doc/index.html#unknown_anchor
invalid: https://espressomd.github.io/doc/unknown_file.html
invalid: [footnote 1](#unknown-footnote-1)
invalid: [resource](file:///home/espresso/image.png)
'''

    def test_detect_invalid_urls(self):
        nbconvert.HTMLExporter.mathjax_url = "file:///usr/share/javascript/mathjax/MathJax.js?config=Safe"
        nbconvert.HTMLExporter.require_js_url = "file:///usr/share/javascript/requirejs/require.min.js"
        html_exporter = nbconvert.HTMLExporter()
        nb = nbformat.v4.new_notebook()
        cell_md = nbformat.v4.new_markdown_cell(source=self.cell_md_src)
        nb['cells'].append(cell_md)
        ref_issues = [
            '"https://espressomd.github.io/doc/index.html" has no anchor "unknown_anchor"',
            '"https://espressomd.github.io/doc/unknown_file.html" does not exist',
            'notebook has no anchor "unknown-footnote-1"',
            '"file:///home/espresso/image.png" is an absolute path to a local file',
            '"file:///usr/share/javascript/requirejs/require.min.js" is an absolute path to a local file',
            '"file:///usr/share/javascript/mathjax/MathJax.js?config=Safe" is an absolute path to a local file',
        ]
        issues = module.detect_invalid_urls(
            nb, build_root='@CMAKE_BINARY_DIR@', html_exporter=html_exporter)
        self.assertEqual(issues, ref_issues)


if __name__ == "__main__":
    ut.main()
