# Copyright (C) 2010,2011,2012 The ESPResSo project
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
#                                                           #
# Basic tests of the Lattice Boltzmann implementation       #
#                                                           #
# 1) check conservation of fluid mass                       #
# 2) check conservation of total momentum                   #
# 3) measure temperature of colloid (no auto-check so far)  #
#                                                           #
#############################################################
set tcl_precision 6 ; # we can't expect more from the gpu...
source "tests_common.tcl"
require_feature "LB_GPU"
require_feature "ADRESS" off

puts "----------------------------------------"
puts "- Testcase $argv0 running on [format %02d [setmd n_nodes]] nodes  -"
puts "----------------------------------------"

#############################################################
# Procedures                                                #
#############################################################

proc read_data {file} {
    set f [open $file "r"]
    while {![eof $f]} {blockfile $f read auto}
    close $f
}


proc write_data {file} {
    set f [open $file "w"]
    blockfile $f write interactions
    blockfile $f write particles {id type q pos v f}
    blockfile $f write bonds
    close $f
}

proc setnodes { u } {
   global box_l
   for { set i 0 } { $i < $box_l } { incr i } {
    for { set j 0 } { $j < $box_l } { incr j } {
     for { set k 0 } { $k < $box_l } { incr k } {
        lbnode $i $j $k set u $u 0. 0. 
     }
    }
   }
}

proc saveparts { file } { 
   set f [ open $file "w" ] 
   puts $f [part 0 print pos v f]
   close $f
} 

proc saveobs { file } {
   set f [ open $file "w" ] 
   puts $f "[analyze fluid momentum]  [analyze fluid mass] [analyze fluid temperature ]"
   close $f
}

proc computediff_obs { file } {
   set f [open $file "r"] 
   set df 0
   gets $f line 
   set data  "[analyze fluid momentum]  [analyze fluid mass] [analyze fluid temperature ]"
   for { set i 0 } { $i < 5 } { incr i } { 
       set val1 [ lindex $data $i ]
       set val2 [ lindex $line $i ]
       set val [expr ($val1-$val2)*($val1-$val2) ]
       set df [expr $df + sqrt($val) ]
   }
   close $f
   return $df
}




proc savenodes { file } {
   global box_l
   set f [open $file "w"] 
   for { set i 0 } { $i < $box_l } { incr i } {
    for { set j 0 } { $j < $box_l } { incr j } {
     for { set k 0 } { $k < $box_l } { incr k } {
   	puts $f "$i $j $k  [lbnode $i $j $k print u ]"
     }
    }
   }
   close $f
}

proc computediff_part { file } {
   set f [open $file "r"] 
   set df 0
   gets $f line 
   set data  [ part 0 print pos v f]
   for { set i 0 } { $i < 9 } { incr i } { 
       set val1 [ lindex $data $i ]
       set val2 [ lindex $line $i ]
       set val [expr ($val1-$val2)*($val1-$val2) ]
       set df [expr $df + sqrt($val) ]
   }
   close $f
   return $df
}

proc printall { } {
   global box_l
   for { set i 0 } { $i < $box_l } { incr i } {
    for { set j 0 } { $j < $box_l } { incr j } {
     for { set k 0 } { $k < $box_l } { incr k } {
	puts "$i $j $k : [lbnode $i $j $k print u]";
     }
    }
   }
}


proc computediff { file } {
   global box_l
   set f [open $file "r"] 
   set df 0
   for { set i 0 } { $i < $box_l } { incr i } {
    for { set j 0 } { $j < $box_l } { incr j } {
     for { set k 0 } { $k < $box_l } { incr k } {
   	gets $f line 
	set ii [lindex $line 0 ] ; set jj [lindex $line 1 ] ; set kk [lindex $line 2 ] ;
	set fx [lindex $line 3 ] ; set fy [lindex $line 4 ] ; set fz [lindex $line 5 ] ;
        if { $ii != $i || $jj != $j || $kk != $k } { 
		puts "Error reading sc.dat, file format not compatible"
        }
	set vel [lbnode $i $j $k print u] ; set velx [ lindex $vel 0 ] ; set vely [ lindex $vel 1 ] ; set velz [ lindex $vel 2 ]
        set df [ expr $df + sqrt (( $fx - $velx) * ($fx - $velx) + ($fy - $vely) * ($fy - $vely) + ($fz - $velz) * ($fz - $velz) )]
     }
    }
   }
   set df [ expr $df / pow($box_l,3)  ]
   close $f
   return $df
}


proc testoutput { df } { 
   if { $df < 1e-10 } { 
	puts "OK -- error = $df"
   } else { 
      if { $df < 1e-5 } { 
	puts "ALMOST OK -- error = $df"
      } else { 
	puts "FAILED -- error = $df"
      }
   }
}



# Integration parameters
#############################################################
set int_steps     100
set int_times     10

set time_step     0.005
set tau           0.005

set agrid         1.0

set box_l         2.0

set dens          0.85
set viscosity     30.0
set friction      2.0

set temp          1.0

set skin          0.5

set mom_prec      1.e-2
set mass_prec     1.e-8
set temp_confidence 10

# Other parameters
#############################################################
set generate 0
if { $argc > 0 } { 
      puts "in"
      set generate [ lindex $argv 0 ]
      puts  "generate  = $generate"
}


#############################################################
# System Setup                                              #
#############################################################

# make real random draw
#set cmd "t_random seed"
#for {set i 0} {$i < [setmd n_nodes]} { incr i } {
#    lappend cmd [expr [pid] + $i] }
#eval $cmd

setmd time_step $time_step
setmd skin $skin

# Simulation box
#############################################################
setmd box_l $box_l $box_l $box_l
setmd periodic 1 1 1

thermostat off

#############################################################
lbfluid gpu dens $dens visc $viscosity agrid $agrid tau $tau
thermostat off

puts -nonewline "Test 1 : "
set filen "sc.1.dat"
setnodes  0.1
integrate 100

if { $generate == 1 } { 
   savenodes $filen
   puts "Data stored in $filen"
  # printall

} else { 
   testoutput [computediff $filen]
  # printall
}

puts -nonewline "Test 2 : "
set filen "sc.2.dat"
setnodes 0.1
lbnode 0 0 0 set u -.1 0. 0.
integrate 100
setnodes 0.1132
if { $generate == 1 } { 
   savenodes $filen
   puts "Data stored in $filen"

} else { 
   testoutput [computediff $filen]
}

thermostat off
set filen "sc.3.dat"
set filep "sc.3p.dat"
setnodes 0.1
thermostat  langevin 1.0 .1
integrate 0
part 0 pos .1 .1 .1 v -0.1 0. 0. f 0 0 0 
integrate 100
setnodes 0.1123

if { $generate == 1 } { 
   savenodes $filen
   saveparts $filep
   puts "Data stored in $filen and $filep"
} else { 
puts -nonewline "Test 3 : "
   testoutput [computediff $filen]
puts -nonewline "Test 3p: "
   testoutput [computediff_part $filep]
}

thermostat off
lbfluid gpu dens $dens visc $viscosity agrid $agrid tau $tau
set filen "sc.4.dat"
set filep "sc.4p.dat"
setnodes 0.1
integrate 0
lbfluid friction 0.1
thermostat lb 1.0
part 0 pos 1 1 1 v -0.3 0. 0. f 0 0 0
integrate 100
setnodes 0.1123

if { $generate == 1 } { 
   savenodes $filen
   saveparts $filep
   puts "Data stored in $filen and $filep"
} else { 
puts -nonewline "Test 4 : "
   testoutput [computediff $filen]
puts -nonewline "Test 4p: "
   testoutput [computediff_part $filep]
#   puts "[part 0 print pos v f] "
}

set filen "sc.5.dat"
if { $generate == 1 } { 
   saveobs $filen
   puts "Data stored in $filen"
} else { 
puts -nonewline "Test 5 : "
   testoutput [computediff_obs $filen]
#   puts "[part 0 print pos v f] "
#   puts "[analyze fluid momentum]  [analyze fluid mass] [analyze fluid temperature ]"
}


exit 0

#############################################################
