#!/usr/bin/env sh
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


# This script generates a summary of all Doxygen warnings that require fixing.

DOXYGEN=$(cmake -LA -N | grep "DOXYGEN_EXECUTABLE:FILEPATH" | cut -d'=' -f2-)
if ! hash "${DOXYGEN}" 2>/dev/null; then
    echo "No doxygen found."
    exit 2
fi
dox_version=$("${DOXYGEN}" --version)
dox_version_supported=false
for supported_version in 1.8.11 1.8.13 1.8.17; do
    if [ "${dox_version}" = "${supported_version}" ]; then
        dox_version_supported=true
        break
    fi
done
if [ "${dox_version_supported}" = false ]; then
    echo "Doxygen version ${dox_version} not fully supported" >&2
    if [ ! -z "${CI}" ]; then
        exit 1
    fi
    echo "Proceeding anyway"
fi

[ -z "${srcdir}" ] && srcdir=$(realpath ..)

# prepare Doxyfile for the XML version
cp doc/doxygen/Doxyfile doc/doxygen/Doxyfile.bak
cat >> doc/doxygen/Doxyfile <<EOF
QUIET                  = YES
WARNINGS               = YES
WARN_FORMAT            = "doxygen:\$file:\$line: \$text"
WARN_LOGFILE           = warnings.log

# Generate warnings for potential errors in the documentation, such as not
# documenting some parameters in a documented function, or documenting
# parameters that don't exist or using markup commands wrongly.

WARN_IF_DOC_ERROR      = YES

# Warn about wrong or incomplete parameter documentation, but not about
# the absence of documentation.

WARN_NO_PARAMDOC       = NO

# Generate XML output instead of HTML to ignore latex- and dot-related issues
# and save time
GENERATE_XML           = YES
GENERATE_HTML          = NO
EOF

# compile the XML version to get the list of warnings
make doxygen

# restore Doxyfile
mv doc/doxygen/Doxyfile.bak doc/doxygen/Doxyfile

# print enabled features
if [ "${CI}" != "" ]; then
    cat doc/doxygen/doxy-features
fi

# find @param without description
grep -Prn '^[\ \t]*(?:\*?[\ \t]+)?[@\\]t?param(?:\[[in, out]+\])?[\ \t]+[a-zA-Z0-9_\*]+[\ \t]*$' "${srcdir}/src" > doc/doxygen/empty-params.log

# remove old summary
rm -f dox_warnings_summary.log

# process warnings
"${srcdir}/maintainer/CI/dox_warnings.py" || exit 2

# print summary
cat dox_warnings_summary.log
n_warnings=$(cat dox_warnings_summary.log | wc -l)

if [ "${CI}" != "" ]; then
    "${srcdir}/maintainer/gh_post_docs_warnings.py" doxygen "${n_warnings}" dox_warnings_summary.log || exit 2
fi

if [ "${n_warnings}" = "0" ]; then
    echo "Found no warning requiring fixing."
    exit 0
else
    if [ "${CI}" = "" ] && [ "${dox_version_supported}" = false ]; then
        echo "Doxygen version ${dox_version} not fully supported." >&2
        echo "This list of warnings may contain false positives." >&2
    fi
    exit 2
fi
