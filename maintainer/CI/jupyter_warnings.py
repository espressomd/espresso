#!/usr/bin/env python3
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
"""
Generate a summary of invalid URLs found in Jupyter notebooks that are
pointing to sections of the notebooks or to the online Sphinx documentation.
"""

import sys
import pathlib
import nbformat

module_path = pathlib.Path(__file__).resolve().parent.parent / 'parsing'
sys.path.insert(0, str(module_path))
import notebook_links


if __name__ == '__main__':
    error_code = 0
    nb_filepaths = sorted(pathlib.Path().glob('doc/tutorials/*/*.ipynb'))
    assert len(nb_filepaths) != 0, 'no Jupyter notebooks could be found!'
    validator = notebook_links.NotebookLinksValidator()
    for nb_filepath in nb_filepaths:
        with open(nb_filepath, encoding='utf-8') as f:
            nb = nbformat.read(f, as_version=4)
        issues = validator.validate(nb)
        if issues:
            error_code = 1
            print(f'In notebook {nb_filepath.name}:', file=sys.stderr)
            for issue in issues:
                print(f'* {issue}', file=sys.stderr)
    if not error_code:
        print('Found no URLs in Jupyter notebooks requiring fixing.')
    exit(error_code)
