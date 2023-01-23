#!/usr/bin/env python3
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
"""
Generate a summary of invalid URLs found in Jupyter notebooks that are
pointing to sections of the notebooks or to the online Sphinx documentation.
"""

import sys
import pathlib

import lxml.etree
import nbformat
import nbconvert


sphinx_docs = {}


def detect_invalid_urls(nb, build_root='.', html_exporter=None):
    '''
    Find all links. Check that links to the Sphinx documentation are valid
    (the target HTML files exist and contain the anchors). These links are
    checked against the Sphinx documentation in the current build directory.
    Links to sections of the notebook (table of contents, bibliography, etc.)
    are also checked.

    Parameters
    ----------
    nb: :obj:`nbformat.notebooknode.NotebookNode`
        Jupyter notebook to process.
    build_root: :obj:`str`
        Path to the ESPResSo build directory. The Sphinx files will be
        searched in :file:`doc/sphinx/html`.
    html_exporter: :obj:`nbformat.HTMLExporter`
        Custom NB convert HTML exporter.

    Returns
    -------
    :obj:`list`
        List of warnings formatted as strings.
    '''
    # convert notebooks to HTML
    if html_exporter is None:
        html_exporter = nbconvert.HTMLExporter()
    html_exporter.template_name = 'classic'
    html_string = html_exporter.from_notebook_node(nb)[0]
    # parse HTML
    html_parser = lxml.etree.HTMLParser()
    root = lxml.etree.fromstring(html_string, parser=html_parser)
    # process all links
    espressomd_website_root = 'https://espressomd.github.io/doc/'
    sphinx_html_root = pathlib.Path(build_root) / 'doc' / 'sphinx' / 'html'
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
            filepath = sphinx_html_root / basename
            if not filepath.is_file():
                broken_links.append(f'"{url}" does not exist')
                continue
            # check anchor exists
            if anchor is not None:
                if filepath not in sphinx_docs:
                    sphinx_docs[filepath] = lxml.etree.parse(
                        str(filepath), parser=html_parser)
                doc = sphinx_docs[filepath]
                nodes = doc.xpath(f'//*[@id="{anchor}"]')
                if not nodes:
                    broken_links.append(f'"{url}" has no anchor "{anchor}"')
        elif url.startswith('#'):
            # check anchor exists
            anchor = url[1:]
            nodes = root.xpath(f'//*[@id="{anchor}"]')
            if not nodes:
                broken_links.append(f'notebook has no anchor "{anchor}"')
        elif url.startswith('file:///'):
            broken_links.append(f'"{url}" is an absolute path to a local file')
    for link in root.xpath('//script'):
        url = link.attrib.get('src', '')
        if url.startswith('file:///'):
            broken_links.append(f'"{url}" is an absolute path to a local file')
    return broken_links


if __name__ == '__main__':
    error_code = 0
    nb_filepaths = sorted(pathlib.Path().glob('doc/tutorials/*/*.ipynb'))
    assert len(nb_filepaths) != 0, 'no Jupyter notebooks could be found!'
    for nb_filepath in nb_filepaths:
        with open(nb_filepath, encoding='utf-8') as f:
            nb = nbformat.read(f, as_version=4)
        issues = detect_invalid_urls(nb)
        if issues:
            error_code = 1
            print(f'In notebook {nb_filepath.name}:', file=sys.stderr)
            for issue in issues:
                print(f'* {issue}', file=sys.stderr)
    if not error_code:
        print('Found no URLs in Jupyter notebooks requiring fixing.')
    exit(error_code)
