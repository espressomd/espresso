#!/usr/bin/tclsh8.4

set f [open [lindex $argv 0] "r"]

while {![eof $f]} {
    set line [gets $f]
    set line [string map {\" <} $line]

    if {[lindex $line 1] == "alloc"} {
	set size [lindex $line 2]
	set new [lindex $line 4]
	if {$new != "(nil)"} {
	    set allocated($new) "$size at [lindex $line 6]"
	}
    }
    if {[lindex $line 1] == "realloc"} {
	set old [lindex $line 2]
	if {$old != "(nil)"} {
	    unset allocated($old)
	}
	set size [lindex $line 6]
	set new [lindex $line 4]
	if {$new != "(nil)"} {
	    set allocated($new) "$size at [lindex $line 8]"
	}
    }
    if {[lindex $line 1] == "free"} {
	set old [lindex $line 2]
	if {$old != "(nil)"} {
	    unset allocated($old)
	}
    }
}

set pl ""
foreach {m f} [array get allocated] {
    lappend pl "$m $f"
}

set pl [lsort -index 0 -integer $pl]

foreach p $pl {
    set m [lindex $p 0]
    set f [lrange $p 1 end]
    puts "$m: $f"
}