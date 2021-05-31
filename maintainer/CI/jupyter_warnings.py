#!/usr/bin/env python3
#
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
#
"""
Generate a summary of invalid URLs found in Jupyter notebooks that are
pointing to sections of the notebooks or to the online Sphinx documentation.
"""

import os
import sys
import glob

import lxml.etree
import nbformat
import nbconvert


sphinx_docs = {}


def detect_invalid_urls(nb, sphinx_root='.'):
    '''
    Find all links. Check that links to the Sphinx documentation are valid
    (the target HTML files exist and contain the anchors). These links are
    checked against the Sphinx documentation in the current build directory.
    Links to sections of the notebook (table of contents, bibliography, etc.)
    are also checked.

    Parameters
    ----------
    nb: :obj:`nbformat.notebooknode.NotebookNode`
        Jupyter notebook to process

    Returns
    -------
    :obj:`list`
        List of warnings formatted as strings.
    '''
    # convert notebooks to HTML
    html_exporter = nbconvert.HTMLExporter()
    html_exporter.template_name = 'classic'
    html_string = html_exporter.from_notebook_node(nb)[0]
    # parse HTML
    html_parser = lxml.etree.HTMLParser()
    root = lxml.etree.fromstring(html_string, parser=html_parser)
    # process all links
    espressomd_website_root = 'https://espressomd.github.io/doc/'
    sphinx_html_root = os.path.join(sphinx_root, 'doc', 'sphinx', 'html')
    broken_links = []
    for link in root.xpath('//a'):
        url = link.attrib.get('href', '')
        if url.startswith(espressomd_website_root):
            # extract anchor
            anchor = None
            if '#' in url:
                url, anchor = url.split('#', 1)
            # strip query
            if '?' in url:
                url = url.split('?', 1)[0]
            # check file exists
            basename = url.split(espressomd_website_root, 1)[1]
            filepath = os.path.join(sphinx_html_root, basename)
            if not os.path.isfile(filepath):
                broken_links.append(f'{url} does not exist')
                continue
            # check anchor exists
            if anchor is not None:
                if filepath not in sphinx_docs:
                    sphinx_docs[filepath] = lxml.etree.parse(
                        filepath, parser=html_parser)
                doc = sphinx_docs[filepath]
                nodes = doc.xpath(f'//*[@id="{anchor}"]')
                if not nodes:
                    broken_links.append(f'{url} has no anchor "{anchor}"')
        elif url.startswith('#'):
            # check anchor exists
            anchor = url[1:]
            nodes = root.xpath(f'//*[@id="{anchor}"]')
            if not nodes:
                broken_links.append(f'notebook has no anchor "{anchor}"')
    return broken_links


if __name__ == '__main__':
    error_code = 0
    for nb_filepath in sorted(glob.glob('doc/tutorials/*/*.ipynb')):
        with open(nb_filepath, encoding='utf-8') as f:
            nb = nbformat.read(f, as_version=4)
        issues = detect_invalid_urls(nb)
        if issues:
            error_code = 1
            basename = os.path.basename(nb_filepath)
            print(f'In notebook {basename}:', file=sys.stderr)
            for issue in issues:
                print(f'* {issue}', file=sys.stderr)
    if not error_code:
        print('Found no URLs in Jupyter notebooks requiring fixing.')
    exit(error_code)
