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

[ -z "${srcdir}" ] && srcdir=`realpath ..`

# Pattern for :class:`foo` commands to non-existent 'foo' classes/functions.
# They are formatted by Sphinx just like regular links to functions, but are
# not enclosed within <a href="..."></a> tags. Sphinx doesn't use line
# wrapping, so these broken links can be found via text search. The first
# negative lookahead filters out common Python types (for performance reasons).
regex_sphinx_broken_link='<code class=\"xref py py-[a-z]+ docutils literal notranslate\"><span class=\"pre\">(?!(int|float|bool|str|object|list|tuple|dict)<)[^<>]+?</span></code>(?!</a>)'

# list of espresso modules not compiled in CI (visualization, scafacos)
regex_ignored_es_features_ci='(visualization|[a-z]+\.[sS]cafacos)'

if [ ! -f doc/sphinx/html/index.html ]; then
    echo "Please run Sphinx first."
    exit 1
fi

n_warnings=0
grep -qrP --include \*.html --exclude-dir=_modules "${regex_sphinx_broken_link}" doc/sphinx/html/
if [ $? = "0" ]; then
    rm -f doc_warnings.log~
    touch doc_warnings.log~
    found="false"
    grep -rPno --include \*.html --exclude-dir=_modules "${regex_sphinx_broken_link}" doc/sphinx/html/ | sort | uniq | while read -r line
    do
        # extract link target
        reference=$(echo "${line}" | sed -r 's|^.+<span class="pre">(.+)</span></code>$|\1|' | sed  's/()$//')
        # skip if broken link refers to a standard Python type or to a
        # class/function from an imported module other than espressomd
        is_standard_type_or_module="false"
        grep -Pq '^([a-zA-Z0-9_]+Error|[a-zA-Z0-9_]*Exception|(?!espressomd\.)[a-zA-Z0-9_]+\.[a-zA-Z0-9_\.]+)$' <<< "${reference}"
        [ "$?" = "0" ] && is_standard_type_or_module="true"
        # skip espresso modules not compiled in CI (visualization, scafacos)
        is_es_feature_skipped="false"
        if [ "${CI}" != "" ]; then
            grep -Pq "^espressomd\.${regex_ignored_es_features_ci}" <<< "${reference}"
            [ "$?" = "0" ] && is_es_feature_skipped="true"
        fi
        # private objects are not documented and cannot be linked
        is_private="false"
        grep -Pq "(^_|\._)" <<< "${reference}"
        [ "$?" = "0" ] && is_private="true"
        # filter out false positives
        if [ ${is_standard_type_or_module} = "true" ] || [ ${is_es_feature_skipped} = "true" ] || [ ${is_private} = "true" ]; then
            continue
        fi
        if [ ${found} = "false" ]; then
            echo "The Sphinx documentation contains broken links:"
        fi
        found="true"
        # locate the .rst file containing the broken link
        filepath_html=$(echo "${line}" | sed -r 's/^(.+?):[0-9]+:<code.+$/\1/')
        filepath_rst=$(echo "${filepath_html}" | sed 's|/html/|/|')
        filepath_rst="${filepath_rst%.html}.rst"
        if [ -f "${filepath_rst}" ]; then
            grep -q -F "\`${reference}\`" "${filepath_rst}"
            if [ $? = "0" ]; then
                grep --color -FnHo -m 1 "\`${reference}\`" "${filepath_rst}" | tee -a doc_warnings.log~
                continue
            fi
        fi
        # if not in a .rst file, show the .html file without line number
        echo "${filepath_html}:\`${reference}\`" | tee -a doc_warnings.log~
    done
    # generate log file
    n_warnings=$(cat doc_warnings.log~ | wc -l)
    echo "The Sphinx documentation contains ${n_warnings} broken links:" > doc_warnings.log
    cat doc_warnings.log~ >> doc_warnings.log
    rm doc_warnings.log~
    # warn user about ignored features in CI
    grep -Pq "espressomd\.${regex_ignored_es_features_ci}" doc_warnings.log
    if [ "$?" = "0" ] && [ "${CI}" = "" ]; then
        echo "(Note that features visualization and Scafacos are ignored in CI)"
    fi
fi

# Find malformed reSt roles, appearing as raw text in the HTML output:
#    * unparsed roles, e.g. ":cite:`UnknownKey`"
#    * broken math formula, e.g. ":math:` leading whitespace`"
#    * incorrect syntax, e.g. "obj:float:"
#    * incorrect numpydoc syntax, e.g. ":rtype:`float`"
# They are difficult to predict, so we leave them to the user's discretion
grep -qrP --include \*.html --exclude-dir=_modules '(:py)?:[a-z]+:' doc/sphinx/html/
if [ $? = "0" ]; then
    echo "Possibly errors:"
    grep -rP --color --include \*.html --exclude-dir=_modules '(:py)?:[a-z]+:' doc/sphinx/html/
fi

if [ "${CI}" != "" ]; then
    "${srcdir}/maintainer/gh_post_docs_warnings.py" sphinx ${n_warnings} doc_warnings.log
fi

if [ ${n_warnings} = "0" ]; then
    echo "Found no broken link requiring fixing."
    exit 0
else
    exit 1
fi
