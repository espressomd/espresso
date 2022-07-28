#!/usr/bin/env python3
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

"""List all tutorial files to deploy (PDF, HTML and figures)"""

import pathlib
import lxml.html

deploy_list = list(pathlib.Path().glob('*/*.pdf'))
for filepath in pathlib.Path().glob('*/*.html'):
    deploy_list.append(filepath)
    # extract all figures
    with open(filepath, encoding='utf-8') as f:
        html = lxml.html.parse(f)
    figures = filter(lambda src: not src.startswith('data:image'),
                     html.xpath('//img/@src'))
    deploy_list += list(map(lambda src: filepath.parent / src, figures))

with open('deploy_list.txt', 'w') as f:
    f.write('\n'.join(map(str, deploy_list)))
