#!/usr/bin/env python3
#
# Copyright (C) 2019-2023 The ESPResSo project
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
List all tutorial files to deploy (PDF, HTML and figures) and write
their filepaths to a text file.
Update metadata of Jupyter notebooks written in the HTML format.
"""

import pathlib
import lxml.etree
import lxml.html
import re

re_title_head = re.compile("<title[ >].+?</title>")

deploy_list = list(pathlib.Path().glob("*/*.pdf"))
for filepath in pathlib.Path().glob("*/*.html"):
    deploy_list.append(filepath)
    root = lxml.html.parse(filepath)
    # extract all figures
    figures = filter(lambda src: not src.startswith("data:image"),
                     root.xpath("//img/@src"))
    deploy_list += list(map(lambda src: filepath.parent / src, figures))
    # update metadata
    try:
        first_title = root.xpath("/html/body//h1")[0]
        metadata = root.xpath("/html/head")[0]
        meta_title_old = root.xpath("/html/head/title")[0]
        meta_title_new = lxml.etree.SubElement(
            metadata, "title", attrib=meta_title_old.attrib)
        meta_title_new.text = first_title.text
        new_title = lxml.html.tostring(meta_title_new, encoding=str)
        with open(filepath, "r+") as f:
            content = f.read()
            assert len(re_title_head.findall(content)) == 1
            content = re_title_head.sub(lambda m: new_title, content, 1)
            f.seek(0)
            f.truncate()
            f.write(content)
    except Exception as err:
        print(f"could not process '{str(filepath)}':")
        print(f"{type(err).__name__}: {err}")

with open("deploy_list.txt", "w") as f:
    f.write("\n".join(map(str, deploy_list)))
