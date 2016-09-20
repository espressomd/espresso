#
# Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
#   Max-Planck-Institute for Polymer Research, Theory Group
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
import espressomd
import numpy
import cPickle as pickle
import sys

# Import system properties
system = espressomd.System()
pickle.load(open("system_config", "r"))

rdf_bins = 100
rdf_list_a = [0]
rdf_list_b = [1]

# Initialize the total RDF with zero
cnt = 0
avg_rdf = numpy.zeros(rdf_bins)

for filename in sys.argv[1:]:
    # Import of particle properties
    try:
        pfile = open(filename, "r")
        pickle.load(pfile)
        pfile.close()
    except Exception as e:
        print("Could not import particles from '{}'.".format(filename))
        print(e)
        continue

    # Calculate rdf
    r, rdf = system.analysis.rdf(rdf_type='rdf', type_list_a=rdf_list_a,
                                 type_list_b=rdf_list_b, r_min=0.9,
                                 r_max=system.box_l[0] / 2., r_bins=rdf_bins)
    avg_rdf += rdf

    system.part[:].remove()

    cnt += 1

avg_rdf *= 1. / cnt

# Write averaged data to plot file
plot = open("rdf.data", "w")
plot.write("# r\trdf(r)\n")
for i in xrange(len(rdf)):
    plot.write("{0}\t{1}\n".format(r[i], avg_rdf[i]))
plot.close()
