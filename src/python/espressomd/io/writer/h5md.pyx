#  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
#  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
#    Max-Planck-Institute for Polymer Research, Theory Group
#
#  This file is part of ESPResSo.
#
#  ESPResSo is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ESPResSo is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.


import sys

cdef class H5md:
    cdef File* c_h5md

    def __cinit__(self, filename):
        self.c_h5md = new File(filename, sys.argv[0])
    def __dealloc__(self):
        del self.c_h5md
    def write(self, arg):
        cmap = {'t': self.c_h5md.W_TYPE,
                'v': self.c_h5md.W_V,
                'f': self.c_h5md.W_F,
                'p': self.c_h5md.W_POS,
                'm': self.c_h5md.W_MASS}
        if type(arg) is str:
            i = reduce(lambda a, b: a | b, map(lambda c: cmap[c], arg))
            self.c_h5md.Write(i)
        elif type(arg) is int:
            self.c_h5md.Write(arg)
    def close(self):
        self.c_h5md.Close()
