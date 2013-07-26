# Copyright (C) 2010,2011,2012,2013 The ESPResSo project
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
#############################################################
#                                                             #
# init.tcl                                                    #
# ========                                                    #
#                                                             #
# Used for initializing tcl-scripts powered by the 'Espresso' #
# MPI parallel molecular dynamics program package.            #
#                                                             #
# It is a very good idea to include all custom tcl-scripts    #
# here which may enhance other programs running 'Espresso'    #
# by providing additional tcl-commands/-functions.            #
#                                                             #
###############################################################
#
# here do everything you want to do upon initialization.
# e. g.

# Read user defined settings
if { [file exists "~/.espressorc" ] } {
    source ~/.espressorc
}

# see whether the user wants the startup
# message. Can also be permanently switched off
# in "~/.espressorc"
set quiet 0
foreach arg $argv { if {$arg == "-quiet"} { set quiet 1 } }

if {!$quiet} {
    puts stderr "*******************************************************"
    puts stderr "*                                                     *"
    puts stderr "*                    - ESPResSo -                     *"
    puts stderr "*                    ============                     *"
    puts stderr "*        A Parallel Molecular Dynamics Program        *"
    puts stderr "*                                                     *"
    puts stderr "* (c) 2010,2011,2012,2013                             *"
    puts stderr "*     The ESPResSo project                            *"
    puts stderr "*                                                     *"
    puts stderr "* (c) 2002,2003,2004,2005,2006,2007,2008,2009,2010    *"
    puts stderr "*     Max-Planck-Institute for Polymer Research       *"
    puts stderr "*     Mainz, Germany                                  *"
    puts stderr "*                                                     *"
    puts stderr "*******************************************************"
    puts stderr ""
    puts stderr "This is [code_info version]."
    puts stderr ""
    puts stderr "ESPResSo is free software: you can redistribute it and/or modify"
    puts stderr "it under the terms of the GNU General Public License as published by"
    puts stderr "the Free Software Foundation, either version 3 of the License, or"
    puts stderr "(at your option) any later version."
    puts stderr ""
    puts stderr "ESPResSo is distributed in the hope that it will be useful,"
    puts stderr "but WITHOUT ANY WARRANTY; without even the implied warranty of"
    puts stderr "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the"
    puts stderr "GNU General Public License for more details."
    puts stderr "" 
    puts stderr "You should have received a copy of the GNU General Public License"
    puts stderr "along with this program.  If not, see <http://www.gnu.org/licenses/>."
    puts stderr ""
}

# Check to see if the user specified a packages directory
if { [ catch { set dummy $env(ESPRESSO_PACKAGES)} ] } {
    # If not then guess it
    cd ..
    lappend auto_path "[pwd]/packages/mbtools/"
    lappend auto_path "[pwd]/packages/mmsg/"
    lappend auto_path "[pwd]/packages/"
    catch { cd "scripts" }
} else {
    # set the packages directory to the user specified value
    lappend auto_path "$env(ESPRESSO_PACKAGES)/"
    lappend auto_path "$env(ESPRESSO_PACKAGES)/mmsg/"
    lappend auto_path "$env(ESPRESSO_PACKAGES)/mbtools/"
}


# Include useful tcl-scripts providing new functions etc.

source countBonds.tcl
source auxiliary.tcl
source blockfile_support.tcl
source pdb.tcl
source polymer.tcl
source statistics.tcl
source ABHmath.tcl
source vtf.tcl
source vtk.tcl
source dielectrics.tcl
source object_in_fluid.tcl
