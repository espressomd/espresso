#############################################################
#                                                           #
#  Analyze for:                                             #
#  Bundle simulation: Ion distribution                      #
#                                                           #
#############################################################
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

puts " "
puts "======================================================="
puts "=                  pe_analyze.tcl                     ="
puts "=     Analysis script for pe_solution.tcl             ="
puts "======================================================="
puts " "


# Simulation to analyze 
# Assumes Configuration files in blockfileformat: $name$ident.####
set name  "pe_solution"
set ident "_s3"

# If time_step and skin are not set in the blockfile, we set it here.
setmd time_step 0.0125
setmd skin      0.4

# Specify which configurations to analyze
set anaident ""
set start 0
set end   100
set step  1

# Ion distributions
set sections 5
set r_min    1.0
set r_max    300.0
set bins     100
set int_flag 1

set tcl_precision 8

# preparations
set config_cnt 0
set sec_cnt 1

set sec_size [expr ($end-$start+1)/$sections]
puts "$sec_size"

# read in configuration and add up the analyzation values
for { set config $start } { $config <= $end } { incr config $step} {

    # Read in configuration
    set file [open "$name$ident.[format %04d $config]" r]
    while { [blockfile $file read auto] != "eof" } {}
    close $file

    # Do analyzations
    set dist [analyze distribution {2} {0 1} $r_min $r_max $bins 1 $int_flag]

    # Add up analyzation values in result
    if { $config == $start } {
	set result "[lindex $dist 1]"
    } else {
	if { ($config-$start)%$sec_size == 0 } { 
	    incr sec_cnt 
	    for { set i 0 } {$i < [llength $result]} {incr i } {
		set temp2 [lindex $result $i]
		lappend temp2 0
		lset result $i $temp2
	    }
	}
 	set temp "[lindex $dist 1]"
	for { set i 0 } {$i < [llength $result]} {incr i } {
	    set ind_t [list $i 1]
	    set ind_r [list $i $sec_cnt]
	    lset result $ind_r [expr [lindex $result $ind_r] + [lindex $temp $ind_t] ]
	}
    }
    puts -nonewline "Analyze ($start-$end) now: $name$ident.[format %04d $config], sec_cnt=$sec_cnt\r"
    flush stdout
    incr config_cnt
}
# average the analyzation values
for { set i 0 } {$i < [llength $result]} {incr i } {
    for { set j 1 } {$j <= $sections} {incr j } {
	set ind [list $i $j]
	lset result $ind [expr [lindex $result $ind]/($config_cnt/$sections)]
    }
}

# write results to file
set file [open "$name$ident.iondis$anaident.dat" w]
puts $file "\# r\tP_1(r)\tP_2(r)\t... "
for { set i 0 } {$i < [llength $result]} {incr i } {
    puts $file "[lindex $result $i]"
}
close $file

# create iondistributions
set file [open "gnu_script.gnu" w]
    if { $int_flag == 1 } {
	puts $file "set title \"Integrated Iondistribution ($name$ident)\""
	puts $file "set yrange \[0:1\]"
	puts $file "set ylabel \"P(r)\""
	puts $file "set log x"
    } else {
	puts $file "set title \"Ion density ($name$ident)\""
	puts $file "set ylabel \"rho\""
	puts $file "set log xy"
    }
puts $file "set xlabel \"r\""
puts $file "set xrange \[$r_min:$r_max\]"
puts $file "set key left top"
set cs $start
set ce [expr $start+$sec_size-1]
    if { $int_flag == 1 } {
	puts $file "p '$name$ident.iondis$anaident.dat' u 1:2 tit \"section 0: $cs-$ce\"w lp"
    } else {
	puts $file "p '$name$ident.iondis$anaident.dat' u 1:(\$2/(3.14*\$1**2)) tit \"section 0: $cs-$ce\"w lp"
    }
for { set j 1 } {$j < $sections} {incr j } {
    set cs [expr $cs+$sec_size]
    set ce [expr $ce+$sec_size]
    puts $file "rep '$name$ident.iondis$anaident.dat' u 1:[expr $j+2] tit \"section $j: $cs-$ce\"w lp"
}
puts $file "set term postscript"
puts $file "set output '$name$ident.iondis$anaident.ps'"
puts $file "rep"
close $file

exec gnuplot gnu_script.gnu &
exec sleep 1
exec gv $name$ident.iondis$anaident.ps -scale -2 &

#create differences 
if { $sections > 1 } {
    set file [open "gnu_script.gnu" w]
    if { $int_flag == 1 } {
	puts $file "set title \"Comparison: Integrated Iondistributions of different sections ($name$ident)\""
    } else {
	puts $file "set title \"Comparison: Iondistributions of different sections ($name$ident)\""
    }
    puts $file "set xlabel \"r\""
    puts $file "set ylabel \"P_n(r)/P_m(r)\""
    puts $file "set log x"
    puts $file "set xrange \[0.7:70\]"
    puts $file "set yrange \[0.8:1.2\]"
    puts $file "p '$name$ident.iondis$anaident.dat' u 1:(\$2/\$3) tit \"sec. 0 / sec. 1\"w lp"
    for { set j 2 } {$j < $sections} {incr j } {
	puts $file "rep '$name$ident.iondis$anaident.dat' u 1:(\$[expr $j+1]/\$[expr $j+2]) tit \"sec. [expr $j-1] / sec. $j\"w lp"
}   
    puts $file "set term postscript"
    puts $file "set output '$name$ident.iondis_diff$anaident.ps'"
    puts $file "rep"
    close $file
    
    exec gnuplot gnu_script.gnu &
    exec sleep 1
    exec gv $name$ident.iondis_diff$anaident.ps -scale -2 &
}   

puts "\nFinished."
exit
