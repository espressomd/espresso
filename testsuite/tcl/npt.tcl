# Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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
source "tests_common.tcl"

require_feature "NPT"
require_feature "LENNARD_JONES"
# incompatible with NPT
require_feature "ROTATION_PER_PARTICLE" off

puts "----------------------------------------"
puts "- Testcase npt.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "----------------------------------------"

proc average { data { from 0 } { to 0 } } {
    set sum 0
    if { $to == 0 || $to >= [llength $data] } { 
	set to [expr [llength $data] - 1 ] 
    }
    if { $from < 0 || $from >= [llength $data] } {
	return 0
    }    
    	for { set i $from } { $i < [expr $to + 1] } { incr i } {
	    set sum [expr $sum + [lindex $data $i] ]
	}
    
    
    return [expr $sum/(1.0*[expr $to - $from + 1 ])]
}

proc stdev { data { from 0 } { to 0 } } {
    set sum 0
    if { $to == 0 || $to >= [llength $data] } { 
	set to [expr [llength $data] - 1 ] 
    }
    if { $from < 0 || $from >= [llength $data] } {
	return 0
    }
    
    set av [average $data $from $to]
    
    for { set i $from } { $i < [expr $to + 1] } { incr i } {
	set sum [expr $sum + ([lindex $data $i] - $av)* ([lindex $data $i] - $av)]
    }
    
    set var [expr $sum/(1.0*[expr $to - $from  ])]
    return [expr sqrt($var)]
}

set epsilon 1e-4
setmd time_step 1
setmd skin 0

proc read_data {file} {
    set f [open $file "r"]
    while {![eof $f]} { blockfile $f read auto}
    close $f
}

proc write_data {file} {
    global energy pressure
    set f [open $file "w"]
    set energy [analyze energy total]
    set pressure [analyze pressure total]
    blockfile $f write tclvariable {energy pressure}
    blockfile $f write variable box_l
    blockfile $f write particles {id type pos f q}
    close $f
}

read_data "npt_lj_system.data"

    ############## lj-specific part



# ----------- Integration Parameters before warmup -----------#
setmd time_step 0.01
setmd skin      0.4
thermostat set langevin 1.0 1.0 

set p_ext 2.0
set piston_mass 0.0001
set gamma_0 0.5
set gamma_v 0.001
set blen [lindex [setmd box_l] 0]

inter 0 0 lennard-jones 1.0 1.0 1.12246 0.25 0.0

integrate set npt_isotropic $p_ext $piston_mass 
#-cubic_box
thermostat set npt_isotropic 1.0 $gamma_0 $gamma_v 

set avp 0
set avpi 0
set avbl 0
set avbli 0

puts -nonewline "Running NPT Test for uncharged lj system "
for { set t 0 } { $t < 1000 } { incr t } {
    integrate 10
    set totprs [analyze pressure total]
    set avp [expr $avp + $totprs]
    incr avpi

    set Linst [lindex [setmd box_l] 0]

    set avbl [expr $avbl + $Linst ]
    incr avbli

    if { $t%100==0 } {
	puts -nonewline ". "
	flush stdout
    }

    if { $t > 500 } {
	lappend storedV [expr pow($Linst,3)]
#	set vkappa [analyze Vkappa]
    }
}

set Vvar [expr pow([stdev $storedV],2)]
set compressibility [expr $Vvar/([average $storedV]*1.0)]
#set compressibility [expr $vkappa/([average $storedV]*1.0)]

if { [expr abs($avp/(1.0*$avpi) - 2.0) ] > 0.2 } {
    error "ERROR: Average pressure <P> = [expr $avp/(1.0*$avpi)] deviates from imposed pressure P = 2.0"
} else { puts "Pressure deviations: [expr abs($avp/(1.0*$avpi) - 2.0) ] acceptable" }  

set Vvar [expr pow([stdev $storedV],2)]
set compressibility [expr $Vvar/([average $storedV]*1.0)]
if { [expr abs($compressibility - 0.209) ] > 0.05 } {
    error "ERROR: Compressibility <K> = $compressibility deviates from known value "
} else { puts "Compressibility deviations: [expr abs($compressibility-0.2093)] acceptable" } 

# --------------------------------------------------------- #
# -- Code to calculate theoretical value of Compressibility #
# --------------------------------------------------------- #

# First set the system back to nvt with a slightly smaller box length

#setmd box_l 5.65 5.65 5.65 
#integrate set nvt
#thermostat set langevin 1.0 0.5
#set f_p [open "TvsP_low" w] 
#for { set t 0 } { $t < 1000 } { incr t } {
#    integrate 10
#    set totprs [analyze pressure total]
#    if { $t > 500 } {
#	lappend storedPlow $totprs
#    }
#}

# Then set the system back to  a slightly larger box length

#setmd box_l 5.8 5.8 5.8 
#for { set t 0 } { $t < 1000 } { incr t } {
#    integrate 10
#    set totprs [analyze pressure total]
#    if { $t > 500 } {
#	lappend storedPhigh $totprs
#    }
#}

#set Pdiff [expr [average $storedPhigh] - [average $storedPlow]]
#puts $Pdiff
#puts [average $storedV]
#set Tcompressibility [expr -(pow(5.8,3)-pow(5.65,3))/($Pdiff*[average $storedV])]
#set Tcompressibility [expr -(pow(5.8,3)-pow(5.65,3))/($Pdiff*187.6405)]
#puts $Tcompressibility


# here you can create the necessary snapshot
#    write_data "npt_lj_system.data"

# ---------------------------------------------------------------------------#
# ---------------------------------------------------------------------------#
# Test NPT with ELECTROSTATICS  #

# Somebody who understands this please feel free to write it #

#} res ] } {
#    error_exit $res
#}

exit 0



