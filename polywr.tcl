# expects the particles to be numbered as 0...maxpart consecutively !!!!!
proc polywrite {name} {
    set out [open $name "w"]

    puts $out "# Praefix-String.......: keiner"
    puts $out "# Anzahl der Polymere..: 1"
    puts $out "# Ladungen pro Polymer.: 1"
    puts $out "# Ladungsabstand.......: 1"
    puts $out "# Valenz der Polymer-Gegenionen.: 1"
    puts $out [format "# Teilchenzahl.........: %d" [expr [setmd maxpart] + 1]]
    puts $out [format "# Boxlaenge............: %3f" [lindex [setmd box_l] 0]]

    for {set p 0} { $p <= [setmd maxpart]} {incr p} {
	set pd [part $p]
	# pos v q
	puts $out "[lrange $pd 1 3] [lrange $pd 9 11] [lindex $pd 7]"
    }

    close $out
}