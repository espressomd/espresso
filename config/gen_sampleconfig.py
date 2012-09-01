# Copyright (C) 2012 Olaf Lenz
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
# This script appends the sample list of features to the file 
#   myconfig-sample.h.
#
import sys, featuredefs, time, string, fileinput

if len(sys.argv) != 2:
    print >> sys.stderr, "Usage: %s DEFFILE" % sys.argv[0]
    exit(2)

deffilename = sys.argv[1]

#print "Reading definitions from " + deffilename + "..."
defs = featuredefs.defs(deffilename)
#print "Done."

#print "Writing " + hfilename + "..."
featuresdone = set()

for line in fileinput.input(deffilename):
    line = line.strip()

    # Handle empty and comment lines
    if len(line) == 0:
        print 
        continue
    elif line.startswith('#'):
        continue
    elif line.startswith('//') or line.startswith('/*'):
        print line
        continue

    # Tokenify line
    feature = line.split(None, 1)[0]

    if feature in defs.features and feature not in featuresdone:
        print '//#define %s' % feature
        featuresdone.add(feature)
