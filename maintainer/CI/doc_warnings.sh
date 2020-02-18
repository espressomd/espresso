#!/usr/bin/env bash
# Copyright (C) 2019 The ESPResSo project
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


# This script searches for :class:`foo` links to non-existent 'foo'
# classes/functions in the output of the Sphinx documentation. When
# broken links are found in .rst files, the putative line numbers are
# printed, otherwise the .html files are printed without line number.

[ -z "${srcdir}" ] && srcdir=$(realpath ..)

# Pattern for :class:`foo` commands to non-existent 'foo' classes/functions.
# They are formatted by Sphinx just like regular links to functions, but are
# not enclosed within <a href="..."></a> tags. Sphinx doesn't use line
# wrapping, so these broken links can be found via text search. The first
# negative lookahead filters out common Python types (for performance reasons).
regex_sphinx_broken_link='<code class=\"xref py py-[a-z]+ docutils literal notranslate\"><span class=\"pre\">(?!(int|float|complex|bool|str|bytes|array|bytearray|memoryview|object|list|tuple|range|slice|dict|set|frozenset|(?:numpy\.|np\.)?(?:nd)?array)<)[^<>]+?</span></code>(?!</a>)'

if [ ! -f doc/sphinx/html/index.html ]; then
    echo "Please run Sphinx first."
    exit 1
fi

n_warnings=0
grep -qrP --include='*.html' --exclude-dir=_modules "${regex_sphinx_broken_link}" doc/sphinx/html/
if [ "${?}" = "0" ]; then
    rm -f doc_warnings.log~
    touch doc_warnings.log~
    found="false"
    grep -rPno --include='*.html' --exclude-dir=_modules "${regex_sphinx_broken_link}" doc/sphinx/html/ | sort | uniq | while read -r line; do
        # extract link target
        reference=$(echo "${line}" | sed -r 's|^.+<span class="pre">(.+)</span></code>$|\1|' | sed 's/()$//')
        lineno=$(echo "${line}" | cut -d ':' -f 2)
        # skip if broken link refers to a standard Python type or to a
        # class/function from an imported module other than espressomd
        is_standard_type_or_module="false"
        grep -Pq '^([a-zA-Z0-9_]+Error|[a-zA-Z0-9_]*Exception|(?!espressomd\.)[a-zA-Z0-9_]+\.[a-zA-Z0-9_\.]+)$' <<< "${reference}"
        [ "${?}" = "0" ] && is_standard_type_or_module="true"
        # private objects are not documented and cannot be linked
        is_private="false"
        grep -Pq "(^_|\._)" <<< "${reference}"
        [ "${?}" = "0" ] && is_private="true"
        # filter out false positives
        if [ "${is_standard_type_or_module}" = "true" ] || [ "${is_private}" = "true" ]; then
            continue
        fi
        if [ "${found}" = "false" ]; then
            echo "The Sphinx documentation contains broken links:"
        fi
        found="true"
        # locate the .rst file containing the broken link
        filepath_html=$(echo "${line}" | cut -d ':' -f 1)
        filepath_rst=$(echo "${filepath_html}" | sed 's|/html/|/|')
        filepath_rst="${filepath_rst%.html}.rst"
        if [ -f "${filepath_rst}" ]; then
            # look for the reference
            grep -q -F "\`${reference}\`" "${filepath_rst}"
            if [ "${?}" = "0" ]; then
                grep --color -FnHo -m 1 "\`${reference}\`" "${filepath_rst}" | tee -a doc_warnings.log~
                continue
            fi
            # if not found, check if reference was shortened, for example
            # :class:`~espressomd.system.System` outputs a link named `System`
            grep -q -P "^[a-zA-Z0-9_]+$" <<< "${reference}"
            if [ "${?}" = "0" ]; then
                grep -q -P "\`~.+${reference}\`" "${filepath_rst}"
                if [ "${?}" = "0" ]; then
                    grep --color -PnHo -m 1 "\`~.+${reference}\`" "${filepath_rst}" | tee -a doc_warnings.log~
                    continue
                fi
            fi
        fi
        # if not in a .rst file, show the .html file
        echo "${filepath_html}:${lineno}:\`${reference}\`" | tee -a doc_warnings.log~
    done
    # generate log file
    n_warnings=$(wc -l < doc_warnings.log~)
    echo "The Sphinx documentation contains ${n_warnings} broken links:" > doc_warnings.log
    cat doc_warnings.log~ >> doc_warnings.log
    rm doc_warnings.log~
fi

# Find malformed reSt roles, appearing as raw text in the HTML output:
#    * unparsed roles, e.g. ":cite:`UnknownKey`"
#    * broken math formula, e.g. ":math:` leading whitespace`"
#    * incorrect syntax, e.g. "obj:float:"
#    * incorrect numpydoc syntax, e.g. ":rtype:`float`"
# They are difficult to predict, so we leave them to the user's discretion
grep -qrP --include='*.html' --exclude-dir=_modules '(:py)?:[a-z]+:' doc/sphinx/html/
if [ "${?}" = "0" ]; then
    echo "Possibly errors:"
    grep -rP --color --include='*.html' --exclude-dir=_modules '(:py)?:[a-z]+:' doc/sphinx/html/
fi

if [ "${CI}" != "" ]; then
    "${srcdir}/maintainer/gh_post_docs_warnings.py" sphinx ${n_warnings} doc_warnings.log || exit 1
fi

if [ "${n_warnings}" = "0" ]; then
    echo "Found no broken link requiring fixing."
    exit 0
else
    exit 1
fi
