# Copyright (C) 2010,2012,2013 The ESPResSo project
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
# blockfile_convert.tcl [-out <filename.vtf> [-vsf <vsfargs>] [-vcf <vcfargs>]]
#          <file 0> <file 1> ....
#
# convert one or more blockfiles to vtf. If output is not given, the
# configurations are sent to VMD. Each blockfile can contain one or
# more particle blocks, which will be treated as individual frames. If
# output goes to VTF, arguments can be given to writevsf/writevcf, for
# example -vsf "radius '0 0.5'" -vcf "folded"
#
#####################################################

# to where should we output?
set outfilename ""
set vsfargs ""
set vcfargs ""

if {[lindex $argv 0] == "-out"} {
    set outfilename [lindex $argv 1] 
    set argv [lrange $argv 2 end]

    while {[lindex $argv 0] != ""} {
	switch [lindex $argv 0] {
	    "-vsf" {
		set vsfargs [lindex $argv 1] 
		set argv [lrange $argv 2 end]
	    }
	    "-vcf" {
		set vcfargs [lindex $argv 1] 
		set argv [lrange $argv 2 end]
	    }
	    default {
		break
	    }
	}
    }
}

proc process_particles {} {
    # called every time we have a new set of particles in order to
    # write out the configuration, either to VMD or vtf.  On first
    # call, this opens a connection to VMD or the vtf file to write
    # to.

    global out

    # open on first call
    if {$out == ""} {
	global outfilename vsfargs
	if {$outfilename == ""} {
	    prepare_vmd_connection "blockfile_convert_tmp" 10000 1
	    set out "<vmd>"
	} {
	    # open vtf file and prepare with topology
	    set out [open $outfilename "w"]
	    eval writevsf $out $vsfargs
	}
    }

    # output configuration
    if { $out == "<vmd>" } {
	imd positions -unfolded
    } {
	global vcfargs
	eval writevcf $out $vcfargs
    }
}

# don't read any Tcl variables, we can't process them anyways
set blockfile_tclvariable_whitelist ""

# loop all blockfiles
set out ""
foreach infilename $argv {
    set in [open "$infilename" "r"]

    while {![eof $in]} {
	set block [blockfile $in read auto]
	switch $block {
	    "eof"       { break }
	    "particles" { process_particles }
	}
    }
}
