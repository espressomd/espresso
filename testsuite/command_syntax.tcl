# Copyright (C) 2010,2011,2012,2013 The ESPResSo project
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

source "tests_common.tcl"

setmd box_l 10.0 10.0 10.0

if { [has_feature "LENNARD_JONES"]} {
    # test variants of using lennard-jones
    inter 0 0 lennard-jones 1.0 1.0 1.12246
    inter 0 1 lennard-jones 1.0 1.0 1.12246 auto 0.0
    inter 1 0 lennard-jones 1.0 1.0 1.12246 auto 0.0 0.5
}

if { [has_feature "LENNARD_JONES"] && [has_feature "LJCOS"]} {
    # test using several interactions in one command
   inter 0 0 lennard-jones 1.0 1.0 1.12246 auto 0 0 0 lj-cos 1.0 1.0 2.0 0.0
}

if { [has_feature "DPD"] } {
    thermostat dpd 1.0 1.0 1.0
}

exit 0
