#!/bin/sh
# tricking...\
    PLATFORM=`uname -s`;
# OSF1 \
    if test $PLATFORM = OSF1; then  exec dmpirun -np 8 $PLATFORM/tcl_md $0 $*
# AIX \
    elif test $PLATFORM = AIX; then exec poe $PLATFORM/tcl_md $0 $* -procs 8
# Linux \
    else lamboot; exec /usr/lib/lam/bin/mpirun -np 8 -nsigs $PLATFORM/tcl_md $0 $*;

# \
    fi;

##################################################
# settings
##################################################
set npart 100
set tcl_precision 5
set writemethod ascii

##################################################
# subroutines
##################################################

# writing ascii output
##################################################
# index list in output of part *:
# pos   1  2  3
# type  5
# q     7
# v     9 10 11
# f    13 14 15

proc awritemd {file} {
    set npart [setmd npart]
    puts $file "\# part | type | p_x p_y p_z"
    for {set p 0} {$p < $npart} {incr p} { if {! [catch {part $p} l]} { puts $file "$p | [lindex $l 5] | [lindex $l 1] [lindex $l 2] [lindex $l 3]" } }
}

proc areadmd {file} {
    # consistency check
    if {[gets $file] != "\# part | type | p_x p_y p_z"} {
	error "$file not in MDASCII format"
	exit
    }
    while {! [eof $file]} {
	gets $file l
	if {$l != ""} { part [lindex $l 0] pos [lindex $l 4] [lindex $l 5] [lindex $l 6] type [lindex $l 2] }
    }
}

# data initialization
##################################################
puts "nproc = [setmd nproc]"
setmd box_l 10.0 8.0 16.0
puts "box =\{[setmd box]\}"
if {[setmd nproc] == 8} {
    setmd procgrid 2 2 2
}
puts "grid = \{[setmd proc]\}"
setmd niatypes 1
puts "niatypes = [setmd niatypes]"

# setup interactions
##################################################
inter 0 0 lennard-jones 1 1 1.2 0 0
inter 1 1 lennard-jones 1 1 1.2 0 0
inter 2 2 lennard-jones 1 1 1.2 0 0
inter 0 1 lennard-jones 1 1 1.2 0 0
inter 0 2 lennard-jones 3 1 1.2 0 0
inter 1 2 lennard-jones 2 1 1.2 0 0
puts "nptypes = [setmd nptypes]"

# Set up particle positions
##################################################
set read 0
if {0 && [file exists "config"]} {
    set compressed ""
    set f [open "|gzip -cd config" r]
    if {[gets $f] != ""} { set compressed .gz } { set f [open "config" r] }
    set l [gets $f]
    switch -glob $l {
	"MDASCII" {
	    set method ascii1
	    set etime [time {
		areadmd $f
	    }]
	}
	"MD01*" {
	    set method binary
	    close $f
	    set f [open "config" r]
	    set etime [time {
		readmd $f
	    }]
	}
	default {
	    error "format not recognized in file $f"
	}
    }
    close $f

    set npart [setmd npart]
    puts "read $npart particles, $method$compressed time $etime"
    set read 1
}

if {$read == 0} {
    set etime [time {
	part 0 pos 1.0 2.0 1.0 type 0 q 1 v 20 10 5 f 0 0 0
	puts "part 0 = [part 0]"
	part 1 pos 1.0 1.0 1.0 type 2 q 1 v 0 0 0 f 10 20 30
	puts "part 1 = [part 1]"
	part 2 pos 10.0 1 1 type 1 q -2   v 0 0 0 f 0 0 0
	puts "part 2 = [part 2]"
	part 3 pos 6.0 1 1 type 3 q -1    v 0 0 0 f 0 0 0
	puts "part 3 = [part 3]"
	part 4 pos 3.0 1 1 type 0 q 1     v 0 0 0 f 0 0 0
	puts "part 4 = [part 4]"
	
	# pump up
	##################################################
#	for {set i 0} { $i < 300 } { incr i} {
#	    puts "test $i [expr srand($i,$i)]"
#	}
	for {set i 5} { $i < $npart } { incr i} {
	    if {[expr $i % 100 == 0]} {
		puts "adding part $i"
	    }
	    part $i pos [expr 10*rand()] [expr 10*rand()] [expr 10*rand()] \
		q [ expr ($i % 2 == 1) ? -1 : 1 ] \
		type 0
	}
    }]
    set npart [setmd npart]
    puts "setup $npart particles, time $etime"
}

# write test config
##################################################
set etime ""
if {$writemethod == "ascii"} {
    set etime [time {
	set f [open "config" w]
	puts $f "MDASCII"
	awritemd $f
	close $f
    }]
}

if {$writemethod == "ascii.gz"} {
    set etime [time {
	set f [open "|gzip -c - >config" w]
	puts $f "MDASCII"
	awritemd $f
	close $f
    }]
}

if {$writemethod == "binary"} {
    set etime [time {
	set f [open "config" w]
	writemd $f type posx posy posz
	close $f
    }]
}

if {$etime != ""} {
    puts "$writemethod time $etime"
}

# integration
##################################################
integrate init
integrate 2
integrate exit

# exit
##################################################
puts "finished"
exit
