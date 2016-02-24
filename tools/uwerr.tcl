# Copyright (C) 2012,2013,2014,2015,2016 The ESPResSo project
# Copyright (C) 2006,2007,2008,2009,2010 Olaf Lenz
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
# Estimate errors

proc usage {} {
    puts "usage:"
    puts "  uwerr.tcl \[OPTIONS...\] FILE COL \[COL...\]"
    puts ""
    puts "  OPTIONS:"
    puts "    -first \$first"
    puts "    -last \$last"
    puts "    -plot"
    puts "    -stepcol \$stepcol"
    exit 1
}

set argc [llength $argv]
if { $argc < 2 } then { usage }

set first 0
set last -1
set stepcol 0
set plot 0

set cols {}
for { set argnum 0 } { $argnum < [llength $argv] } { incr argnum } {
    set arg [ lindex $argv $argnum ]
    set val [ lindex $argv [expr {$argnum + 1}]]
    switch -glob -- $arg {
	"-first" { set first $val; incr argnum }
	"-last" { set last $val; incr argnum }
	"-stepcol" { set stepcol $val; incr argnum }
	"-plot" { set plot 1 }
	"-*" { puts "unknown option: $arg"; usage }
	default {
	    set filename [lindex $argv $argnum]
	    set cols [lrange $argv [expr $argnum + 1] end]
	    break
	}
    }
}

set ocols $cols

puts  "Using uwerr on file $filename, columns $cols."
if { [string length $filename] < 1 || [llength $cols] < 1 } then { usage }

set file [open $filename r]
set colnames {}

# read comment lines at the beginning of the file
# the last comment line before the first data line is interpreted as
# name line

# the number of the column
while { ! [eof $file] } {
    set line [string trim [gets $file]]
    if { [string index $line 0] != "\#" } then {
	break
    } else {
	# found a comment line before the first data line
	# interpret it as label line
	set colnames [split [string range $line 1 end]]

	# identify $cols
	set i 0
	foreach col $cols {
	    # test whether $col is a name
	    if { ! [string is integer $col] } then {
		set j 0
		foreach colname $colnames {
		    if { $colname == $col } then {
			lset cols $i $j
			break
		    }
		    incr j
		}
	    } 
	    incr i
	}
	# identify $stepcol
	if { ! [string is integer $stepcol] } then {
	    set j 0
	    foreach colname $colnames {
		if { $colname == $stepcol } then {
		    lset stepcol $j
		    break
		}
		incr j
	    }
	}
    }
}

# test whether cols is useful
foreach col $cols {
    if { ! [string is integer $col] } then {
	puts "Error: Column $col does not exist!"
	if { [llength $colnames] > 0 } then {
	    puts "Column names: $colnames"
	}
	exit 2
    }
}

if { ! [string is integer $stepcol] } then {
    puts "Error: Step column $stepcol does not exist!"
    if { [llength $colnames > 0] } then {
	puts "Column names: $colnames"
    }
    exit 2
}

# read the data (measurement file)
puts "Reading columns $cols between steps $first and $last from $filename..."
set data {}
set N 0

set nfields 0
set lineno 1
while { ! [eof $file] } {
    if { [string index $line 0] != "\#" } then {
	set record [ split $line ]
	# check the number of fields
	if { $nfields == 0 } then {
	    set nfields [llength $record]
	} elseif { $nfields != [llength $record] } then {
	    puts "WARNING: Wrong number of fields in line $lineno!"
	}
	# check the fields
	foreach field $record {
	    if { [string length $field] == 0 } then {
		puts "WARNING: Empty field in line $lineno!"
	    }
	}

	set step [lindex $record $stepcol]
	puts -nonewline "$step\r"
	flush stdout
	if { $step >= $first && \
		 ( $step <= $last || $last == -1 ) } \
	    then { lappend data $record; incr N }
    }
    set line [string trim [gets $file]]
    incr lineno
}

close $file

puts "Got $N data samples."

set res ""
foreach col $cols ocol $ocols {
    # execute uwerr on the data
    if { $plot } then {
	set results [uwerr $data $N [expr $col + 1] plot]
    } else {
	set results [uwerr $data $N [expr $col + 1]]
    }
    set mean [lindex $results 0]
    set error [lindex $results 1]
    set errerr [lindex $results 2]
    set actime [lindex $results 3]
    set acerr [lindex $results 4]
    
    puts "column=$ocol ($col)"
    puts "mean=$mean"
    puts "error=$error"
    puts "error of the error=$errerr"
    puts "autocorrelation time=$actime"
    puts "autocorrelation time error=$acerr"
    puts ""

    set res "$res\t$mean\t$error"
}

puts "--------------------------------------------------"
puts "DATA:$res"
