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


# This script generates a log of all Doxygen warnings that require fixing.

dox_version=`doxygen --version`
if [ "${dox_version}" != "1.8.13" ] && [ "${dox_version}" != "1.8.11" ]; then
    echo "Doxygen version ${dox_version} not supported"
    exit 1
fi

[ -z "${srcdir}" ] && srcdir=`realpath ..`

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

rm -f dox_warnings.log

# process warnings
${srcdir}/maintainer/CI/dox_warnings.py
n_warnings=$?

if [ ! -f dox_warnings.log ]; then
    exit 1
fi

# print logfile
cat dox_warnings.log

if [ "${CI}" != "" ]; then
    "${srcdir}/maintainer/gh_post_docs_warnings.py" doxygen ${n_warnings} dox_warnings.log
fi

if [ ${n_warnings} = "0" ]; then
    echo "Found no warning requiring fixing."
    exit 0
else
    exit 1
fi
