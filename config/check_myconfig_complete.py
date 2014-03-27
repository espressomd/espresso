# Copyright (C) 2012,2013 Olaf Lenz
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
# Check whether all features used in the code are defined
#
from __future__ import print_function
import sys, os, re, fileinput
sys.path.append(os.path.join(sys.path[0], '..', '..', 'config'))

import featuredefs

if len(sys.argv) < 3:
    print("Usage: %s DEFFILE [FILE...]" % sys.argv[0])
    exit(2)

print("Checking for completeness of features in test configurations...")

fdefs = featuredefs.defs(sys.argv[1])

featurefound = set()
featurere = re.compile('^#define (\w+)')
for line in fileinput.input(sys.argv[2:]):
    res = featurere.match(line)
    if res is not None:
        feature = res.group(1)
        featurefound.add(feature)

unused = fdefs.features.difference(featurefound)
unused = unused.difference(fdefs.notestfeatures)
if len(unused) > 0:
    for feature in unused:
        print("check_myconfig_complete: %s is not used" % feature)
else:
    print("check_myconfig_complete: All features are used!")
        

