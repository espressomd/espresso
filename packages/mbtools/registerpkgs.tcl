# 
# tcl script to create pkgIndex for all packages inside the packages directory
#
# If you add a new package you should ensure that the package gets its
# own directory and that you add the name of that directory to the
# list "pkgs" below. 
#
# Inside the directory you should always have one main tcl file (with
# same name as directory) that contains the package provide command
# and gives the version.  The minor version number should represent
# bug fixes.  The middle version should represent new features and the
# major version should be changed if there are compatibility issues
# between versions. 
#
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

set pkgs { system_generation analysis utils }

foreach package $pkgs {
    pkg_mkIndex -verbose $package 
}


