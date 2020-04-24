#!/usr/bin/env python3
#
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
#
import re
import os

# collect list of Doxygen warnings
with open('doc/doxygen/warnings.log') as f:
    content = f.read()

raw_warnings = re.findall(
    r'(?:^|\n)doxygen:(.+?):(\d+): warning: (.+?)(?=\n\S|\n*$)',
    content, re.DOTALL)

# collect list of empty @param and @tparam blocks
with open('doc/doxygen/empty-params.log') as f:
    content = f.read().strip()

if content:
    source_code_ext = set(['.hpp', '.cpp', '.hh', '.cc', '.h', '.c', '.cuh',
                           '.cu', '.dox'])
    for line in content.strip().split('\n'):
        m = re.search(r'^(.+):(\d+):[\s\*]*([@\\]t?param).*\s(\S+)\s*$', line)
        filepath, lineno, paramtype, varname = m.groups()
        ext = os.path.splitext(filepath)[1]
        if ext.lower() not in source_code_ext:
            continue
        warning = ('argument \'{0}\' of {1} has no description, either add one'
                   ' or remove {1}'.format(varname, paramtype))
        raw_warnings.append((filepath, lineno, warning))

# remove duplicated warnings
n_all = len(raw_warnings)
raw_warnings = {(filepath, int(lineno), warning.split('\n')[0]):
                set(warning.split('\n')[1:])
                for filepath, lineno, warning in raw_warnings}
n_unique_raw = len(raw_warnings)

# filter out non-critical warnings
warnings = {}
for (filepath, lineno, warning), warning_list in raw_warnings.items():
    if re.search(r'^member \S+ belongs to two different groups\. '
                 r'The second one found here will be ignored\.$', warning):
        # happens when a function is declared in a group in the .hpp file but
        # defined in another group in the .cpp file; this is usually caused by
        # the "Private functions" and "Exported functions" groups in .hpp files
        continue
    if re.search(r'^documented symbol `\S+\' was not declared or defined\.$',
                 warning):
        # known bug, fixed in 1.8.16
        continue
    if re.search('^no uniquely matching class member found for $', warning):
        # known bug, not fixed yet
        continue
    if re.search(
            '^The following parameters? of .+ (is|are) not documented:$', warning):
        # non-critical warning, enforcing it would encourage bad behavior, i.e.
        # inserting "@param argument" without a description to silence the
        # warning, when in reality the warning is silenced because the text on
        # the following line is captured and becomes the argument description
        continue
    filepath = re.sub(r'^.*(?=src/)', '', filepath)
    if filepath not in warnings:
        warnings[filepath] = {}
    warnings[filepath][(lineno, warning)] = warning_list

n_unique = sum(map(len, warnings.values()))
if n_unique == 0:
    with open('dox_warnings_summary.log', 'w') as f:
        pass
    exit()

# generate a log file
with open('dox_warnings_summary.log', 'w') as f:
    f.write('The Doxygen documentation generated {} unique warnings (total: {},'
            ' ignored: {}):\n'.format(n_unique, n_all, n_unique_raw - n_unique))
    for filepath in sorted(warnings.keys()):
        f.write(filepath + ':\n')
        for (lineno, warning) in sorted(warnings[filepath].keys()):
            warning_list = warnings[filepath][(lineno, warning)]
            s = re.sub(r'\(.*\)', '()', warning)
            if warning_list:
                s += ': ' + ', '.join(x.strip() for x in warning_list)
            f.write('  line {}: {}\n'.format(lineno, s))
