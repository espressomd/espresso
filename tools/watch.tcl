#!/bin/sh
# tricking... the line after a these comments are interpreted as standard shell script \
    exec $ESPRESSO_SOURCE/Espresso $0 $*
#############################################################
#                                                           #
#  Watch Configurations as animation in vmd                 #
#                                                           #
#############################################################
#
# Copyright (C) 2010 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany
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
puts "=                     watch.tcl                       ="
puts "======================================================="
puts " "
puts " "

# Default values
set dir   "./"
set start 0
set end   -1
set step  1
# flag wether you want to delete the created psf/pdb files
set format "%04d"
set delete 1

# Parse command line options
if {$argc > 0} { set dir    [lindex $argv 0] }
if {$argc > 1} { set start  [lindex $argv 1] }
if {$argc > 2} { set end    [lindex $argv 2] }
if {$argc > 3} { set step   [lindex $argv 3] }
if {$argc > 4} { set format   [lindex $argv 4] }
if {$argc > 5} { set delete [lindex $argv 5] }

if { $dir == "?" ||  $dir == "h" || $dir == "-h" || $dir == "help" } {
    puts "Usage: watch.tcl <directory> <first> <last> <step> <delete>"
    puts "Description: Watch Espresso configurations with VMD"
    puts "<directory>  Complete path, where the configurations are."
    puts "<first>      Number of first configuration (default 0)."
    puts "<last>       Number of last configuration (default auto=-1)."
    puts "<step>       Step size (default 1)."
    puts "<format>     configuration name format (default \"%04d\")."
    puts "<delete>     Wether to delete the created pdb files (1) or not (0) (default 1)."
    exit
}
# try to get configuration names 

set tmp_name [exec ls -C1 $dir | grep "[format $format $start]"]

set choice  [exec ls -C1 $dir | grep "[format $format $start]" | wc -l ]
if { $choice > 1 } {
    puts "Found more than one matching file:"
    puts "\n$tmp_name\n"
    puts -nonewline "Specify additional pattern: "; flush stdout
    gets stdin pattern
    set tmp_name [exec ls -C1 $dir | grep "[format $format $start]" | grep "$pattern"]
    if { [string length $tmp_name] == 0 } {
	puts "No match found"
	exit
    }
}

set length  [string length $tmp_name]
set name [string replace $tmp_name [expr $length-5] $length]
# check directory ending
if { [string index $dir end] != "/" } { set dir "$dir/" }
# Check existing configurations
if { $end == -1 } {
    set flag 1
    set i $start
    while { $flag == 1 } {
	if { ![file exists "$dir$name.[format $format $i]"] } { set flag 0 }
	set i [expr $i+$step]
    }
    set end [expr $i-2*$step]
}
# Tell what will be done
puts "Watch configs in $dir with VMD"
puts "$name.$start-$end step $step (del=$delete)"
# Read in configurations and create psf/pdb files
set name "$dir$name"
set pdb_cnt 0
for { set config $start } { $config <= $end } { incr config $step} {
    puts -nonewline "Read $name.[format $format $config] - create pdb\r" 
    flush stdout
    set file [open "$name.[format $format $config]" r]
    while { [blockfile $file read auto] != "eof" } {}
    close $file

    if { $config == $start } { writepsf "$name.vmd.psf" }
    writepdb "$name.vmd[format $format $pdb_cnt].pdb"
    incr pdb_cnt
}
# Create vmd script
set vmd_file [open "vmd_movie.script" "w"]
puts $vmd_file "loadseries $name.vmd"
#puts $vmd_file "mol modstyle 0 0 CPK 1.000000 0.300000 8.000000 6.000000"
#puts $vmd_file "mol modcolor 0 0 SegName"
#puts $vmd_file "nearclip set 0.01"
#puts $vmd_file "logfile vmd.log"
#puts $vmd_file "scale by 1.7"
puts $vmd_file "animate forward"
puts $vmd_file "logfile off"
close $vmd_file

puts "\nStart VMD"
exec vmd -e vmd_movie.script

# wait until vmd has read in the pdb files
set vmd_finish 0
while { $vmd_finish == 0 } {
    if { [file exists "vmd.log" ] } { set vmd_finish 1 }
}
exec sleep 2
# Clean up
if { $delete != 0 } {
    puts "Delete psf/pdb files"
    set pdb_cnt 0
    for { set config $start } { $config <= $end } { incr config $step} {
	if { $config == $start } { exec rm "$name.vmd.psf" }
	exec rm "$name.vmd[format $format $pdb_cnt].pdb"
	incr pdb_cnt
    }
}
puts "Clean up working directory"
exec rm vmd.log
exec rm vmd_movie.script

puts "\nFinished\n"
exit
