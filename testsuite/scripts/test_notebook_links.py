#
# Copyright (C) 2020-2026 The ESPResSo project
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
import time
import pathlib
import tempfile
import importlib
import contextlib
import unittest as ut

missing_modules = True
with contextlib.suppress(ImportError):
    import nbformat
    import nbconvert
    sys.path.insert(0, '@CMAKE_SOURCE_DIR@/maintainer/parsing')
    module = importlib.import_module('notebook_links')
    missing_modules = False


@ut.skipIf(missing_modules, "Jupyter-related dependencies are unavailable")
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
    html_src = '''\
<!DOCTYPE html>
<html lang="en" data-content_root="./">
  <head>
    <meta charset="utf-8" />
  </head>
  <body>
    <section id="python-modules">
    <h2>Python modules</h2>
  </body>
</html>
'''

    def test_detect_invalid_urls(self):
        temp_dir = tempfile.TemporaryDirectory()
        temp_root = pathlib.Path(temp_dir.name).resolve()
        temp_file = temp_root / 'doc' / 'sphinx' / 'html' / 'index.html'
        temp_file.parent.mkdir(parents=True, exist_ok=False)
        temp_file.write_text(self.html_src, encoding="utf-8")
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
        # absorb file system latency
        for _ in range(10):
            time.sleep(0.01)
            if temp_file.exists():
                break
        assert temp_file.exists(), f"excessive write latency to '{temp_file}'"
        # check validation
        validator = module.NotebookLinksValidator(
            build_root=str(temp_root), html_exporter=html_exporter)
        issues = validator.validate(nb)
        temp_dir.cleanup()
        self.assertEqual(issues, ref_issues)


if __name__ == "__main__":
    ut.main()
