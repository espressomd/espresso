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

# check the electristatics algorithms
source "tests_common.tcl"

require_feature "ELECTROSTATICS"

puts "---------------------------------------------------------------"
puts "- Testcase coulomb_cloud_wall.tcl running on [format %02d [setmd n_nodes]] nodes:"
puts "---------------------------------------------------------------"

thermostat off
setmd time_step 0.01
setmd skin 0.05

set file "coulomb_cloud_wall_system.data"
set f [open $file "r"]

while {![eof $f]} { blockfile $f read auto}
    close $f

for { set i 0 } { $i <= [setmd max_part] } { incr i } {
  set F($i) [part $i print f]
}

# write vtf file
# set vtffile [open "es_cloud_wall.vtf" "w"]
# writevsf $vtffile
# for { set i 0 } { $i < 299 } { incr i } {
#     if { $i % 2 == 0 } {
# 	puts $vtffile "atom $i name S radius 0.1"
#     } else {
# 	puts $vtffile "atom $i name O radius 0.1"
#     }
# }
# writevcf $vtffile
# close $vtffile

# To retune, use
#puts [inter coulomb 1.0 p3m tune accuracy 1.e-3]

# If you want to use specific P3M params, use
# inter coulomb 1.0 p3m 2.000000 32 6 1.322773
# inter coulomb n_interpol 0

set methods {}
set setups {}
set accuracies {}

if { [ has_feature "SCAFACOS" ] && [lsearch [scafacos_methods] "ewald"] } then {
    proc setup_scafacos_ewald {} {
        inter coulomb 1.0 scafacos ewald ewald_r_cut 1.001 tolerance_field 1e-4
    }
    lappend methods "SCAFACOS_ewald"
    lappend setups "setup_scafacos_ewald"
    lappend accuracies 1e-3
}

if { [ has_feature "SCAFACOS" ] && [lsearch [scafacos_methods] "p3m"] } then {
    proc setup_scafacos_p3m {} {
        inter coulomb 1.0 scafacos p3m p3m_r_cut 1.001 p3m_grid 64 p3m_cao 7 p3m_alpha 2.70746
    }
    lappend methods "SCAFACOS_p3m"
    lappend setups "setup_scafacos_p3m"
    lappend accuracies 1e-3
}

if { [ has_feature "SCAFACOS" ] && [lsearch [scafacos_methods] "p2nfft"] } then {
    proc setup_scafacos_p2nfft {} {
        inter coulomb 1.0 scafacos p2nfft p2nfft_r_cut 1.001 tolerance_field 1e-4
    }
    lappend methods "SCAFACOS_p2nfft"
    lappend setups "setup_scafacos_p2nfft"
    lappend accuracies 1e-3
}

if { [ has_feature "P3M" ] } then {
    proc setup_p3m {} {
        inter coulomb 1.0 p3m 1.001 64 7 2.70746
        inter coulomb epsilon metallic n_interpol 32768 mesh_off 0.5 0.5 0.5 
    }
    lappend methods "P3M" 
    lappend setups "setup_p3m"
    lappend accuracies 1e-3
}

if { [ has_feature "CUDA" ] && [ has_feature "P3M" ] } then {
    proc setup_p3m_gpu {} {
        inter coulomb 1.0 p3m gpu 1.001 64 7 2.70746
        inter coulomb epsilon metallic n_interpol 32768 mesh_off 0.5 0.5 0.5 
    }
    lappend methods "P3M-GPU"
    lappend setups "setup_p3m_gpu"
    lappend accuracies 1e-3
}

if { [ has_feature "CUDA" ] && [ has_feature "EWALD_GPU" ] } then {
    proc setup_ewald_gpu {} {
	inter coulomb 1.0 ewaldgpu 2.70746 64 1.001
    }
    lappend methods "EWALD-GPU"
    lappend setups "setup_ewald_gpu"
    lappend accuracies 1e-3
}

if { ! [ llength $methods ] } {
    puts "No Coulomb methods to test."
    ignore_exit
}
puts "Methods to test: $methods"

foreach method $methods setup $setups accuracy $accuracies {
    puts "Testing $method..."
    eval $setup
    puts [inter coulomb]

    if { [catch {
        for { set pid 0 } { $pid <= [setmd max_part] } { incr pid } {
            part $pid f 0.0 0.0 0.0
        }        
        integrate 0 recalc_forces

        set rmsf 0
        for { set pid 0 } { $pid <= [setmd max_part] } { incr pid } {
            set resF [part $pid print f]
            set tgtF $F($pid)
            set dx [expr ([lindex $resF 0] - [lindex $tgtF 0])]
            set dy [expr ([lindex $resF 1] - [lindex $tgtF 1])]
            set dz [expr ([lindex $resF 2] - [lindex $tgtF 2])]
            set rmsf [expr $rmsf + $dx*$dx + $dy*$dy + $dz*$dz]
        }

        set rmsf [expr sqrt($rmsf / [setmd n_part])]
        puts [format "  rms_force_error=%e" $rmsf]
        if { $rmsf > $accuracy } then {
            error [format \
                       "coulomb_cloud_wall: rms_force_error=%e larger than accuracy=%e" \
                       $rmsf $accuracy]
        }
    } res ] } { error_exit $res }

}
exit 0
