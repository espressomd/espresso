#!/bin/sh
# tricking... the line after a these comments are interpreted as standard shell script \
    PLATFORM=`uname -s`; if [ "$1" != "" ]; then NP=$1; else NP=2; fi
# OSF1 \
    if test $PLATFORM = OSF1; then  exec dmpirun -np $NP $ESPRESSO_SOURCE/$PLATFORM/Espresso $0 $*
# AIX \
    elif test $PLATFORM = AIX; then exec poe $ESPRESSO_SOURCE/$PLATFORM/Espresso $0 $* -procs $NP
# Linux \
    else export EF_ALLOW_MALLOC_0=1; exec mpirun -np $NP -nsigs $ESPRESSO_SOURCE/$PLATFORM/Espresso $0 $*;
# \
    fi;

#############################################################
#                                                           #
#  Demo of a buckyball                                      #
#                                                           #
#                                                           #
#  Created:       26.08.2003 by AxA & BAM                   #
#                                                           #
#############################################################

set case2timeout 120

#############################################################
#  Setup GUI                                                #
#############################################################

# popup window if desired
proc popup {text} {
    toplevel .tag
    button .tag.btn -text "text" -command "destroy .tag"
    pack .tag.btn
    wm geometry .tag +60+0
    tkwait window .tag
    update
}

wm geometry . -0+0

##### case 1
frame .case1 -relief raised -border 2
pack .case1 -expand 1 -fill both -in .

# title
frame .case1.title
label .case1.title.l -justify left -text "Schrumpfen"
button .case1.title.b -text "Start" -command Case1Start
pack .case1.title.l -fill both -side left -in .case1.title
pack .case1.title.b -side right -in .case1.title
pack .case1.title -pady {0 10} -fill x -side top -in .case1

# sliders
label .case1.label -justify left -text "Stärke der Elektrostatik"
scale .case1.slider -orient h -from 0 -to 15 \
    -resolution [expr 10/100.] -command Case1BjerrumChange
pack .case1.slider .case1.label -fill x -in .case1

label .case1.label2 -justify left -text "Temperatur"
scale .case1.slider2 -orient h -from 0 -to 2 \
    -resolution [expr 10/100.] -command Case1TempChange
pack .case1.slider2 .case1.label2 -fill x -in .case1

# radiobuttons
frame .case1.butFH
label .case1.labelF -justify left -text "FENE"
label .case1.labelH -justify left -text "Harmonisch"
radiobutton .case1.butF -variable butFH -value F -command Case1PotChange
radiobutton .case1.butH -variable butFH -value H -command Case1PotChange
pack .case1.butF .case1.labelF -fill both -side left  -in .case1.butFH
pack .case1.labelH .case1.butH -fill both -side right -in .case1.butFH
pack .case1.butFH  -pady {0 10} -fill x -side bottom  -in .case1

set disabledfg [.case1.label cget -disabledforeground]
set normalfg  [.case1.label cget -foreground]

##### case 2
frame .case2 -relief raised -border 2
pack .case2 -expand 1 -fill both -in .

# title
frame .case2.title
label .case2.title.l -justify left -text "Irrgarten"
button .case2.title.b -text "Start" -command Case2Start
pack .case2.title.l -fill both -side left -in .case2.title
pack .case2.title.b -side right -in .case2.title
pack .case2.title -pady {0 10} -fill x -side top -in .case2

# sliders
label .case2.label1 -justify left -text "Ladung Wand 1"
scale .case2.slider1 -orient h -from -1 -to 1 \
    -resolution [expr 2/100.] -command Case2Wall1ChargeChange
pack .case2.slider1 .case2.label1 -fill x -in .case2

label .case2.label2 -justify left -text "Ladung Wand 2"
scale .case2.slider2 -orient h -from -1 -to 1 \
    -resolution [expr 2/100.] -command Case2Wall2ChargeChange
pack .case2.slider2 .case2.label2 -fill x -in .case2

label .case2.label3 -justify left -text "Ladung Wand 3"
scale .case2.slider3 -orient h -from -1 -to 1 \
    -resolution [expr 2/100.] -command Case2Wall3ChargeChange
pack .case2.slider3 .case2.label3 -fill x -in .case2

label .case2.status -justify left -text ""
pack .case2.status -fill x -in .case2

proc disableCase1 {} {
    global disabledfg butFH
    .case1.slider set 8
    .case1.label configure -state disabled
    .case1.slider configure -state disabled -foreground $disabledfg 
    .case1.slider2 set 0
    .case1.label2 configure -state disabled
    .case1.slider2 configure -state disabled -foreground $disabledfg 
    set butFH F
    .case1.labelF configure -state disabled
    .case1.butF configure -state disabled -foreground $disabledfg 
    .case1.labelH configure -state disabled
    .case1.butH configure -state disabled -foreground $disabledfg 
}

proc enableCase1 {} {
    global normalfg butFH
    .case1.label configure -state normal
    .case1.slider configure -state normal -foreground $normalfg 
    .case1.slider set 8
    .case1.label2 configure -state normal
    .case1.slider2 configure -state normal -foreground $normalfg 
    .case1.slider2 set 0
    .case1.labelF configure -state normal
    .case1.butF configure -state normal -foreground $normalfg 
    .case1.labelH configure -state normal
    .case1.butH configure -state normal -foreground $normalfg 
    set butFH F

    .case2.status configure -text ""
}

proc disableCase2 {} {
    global disabledfg displayclock
    .case2.slider1 set 0
    .case2.label1 configure -state disabled
    .case2.slider1 configure -state disabled -foreground $disabledfg
    .case2.slider2 set 0
    .case2.label2 configure -state disabled
    .case2.slider2 configure -state disabled -foreground $disabledfg
    .case2.slider3 set 0
    .case2.label3 configure -state disabled
    .case2.slider3 configure -state disabled -foreground $disabledfg

    .case2.status configure -text ""
    set displayclock 0
}

proc enableCase2 {} {
    global normalfg
    .case2.label1 configure -state normal
    .case2.slider1 configure -state normal -foreground $normalfg
    .case2.slider1 set 0
    .case2.label2 configure -state normal
    .case2.slider2 configure -state normal -foreground $normalfg
    .case2.slider2 set 0
    .case2.label3 configure -state normal
    .case2.slider3 configure -state normal -foreground $normalfg
    .case2.slider3 set 0
}

proc disableStarts {} {
    .case1.title.b configure -state disabled
    .case2.title.b configure -state disabled
}

proc enableStarts {} {
    .case1.title.b configure -state normal
    .case2.title.b configure -state normal
}

#############################################################
#      Helpers                                              #
#############################################################

proc read_highscore { } {
    set f [open highscore.txt r]
    while {[eof $f]==0} { if { [expr [gets $f inp] > 0] } { lappend highscore [list [lindex $inp 0] [concat [lrange $inp 1 end]]] } }
    close $f
    return $highscore
}

proc write_highscore { highscore } {
    set f [open highscore.txt w]
    for {set i 0} {$i < [llength $highscore]} {incr i} {
	set tmp [lindex $highscore $i]
	puts $f "[lindex $tmp 0]  [lindex $tmp 1]"
	puts "[lindex $tmp 0]  [lindex $tmp 1]"
    }
    close $f
}

proc imd_reconnect {case} {
    imd disconnect

    puts "before [glob -nocomplain .lock*]"
    eval "exec rm -f [glob -nocomplain .lock*]"
    puts "after [glob -nocomplain .lock*]"

    while { [catch {imd connect 10000} res] } {
	.case2.status configure -text "Waiting for port to unbind...\nerror: $res"
	update
	after 500
    }
    exec vmd -pos 100 0 -e $case.script &
    while { [glob -nocomplain ".lock*" ] == "" } {
	puts "waiting for vmd to come up"
	after 500
    }
    while { [catch {imd positions} res] } {
	puts "connect failed: $res, retrying."
	after 500
    }
}

#############################################################
#      Events                                               #
#############################################################

proc Case1Start {} {
    global run_sim

    .case2.status configure -text "Starte..."    
    if { $run_sim == 1 } { set restart 1 } else { set restart 0 }
    set run_sim 0
    update

    disableStarts
    disableCase1
    disableCase2

    part deleteall
    set inp [open "case1.blk" r]
    while { [blockfile $inp read auto] != "eof" } {}
    close $inp
    puts "Read [setmd n_part] particles..."

    if { $restart == 0 } { imd_reconnect case1 }

    after 3000

    enableCase1
    enableStarts
    set run_sim 1

    .case2.status configure -text "Fertig"; update
}

set displayclock 0
proc Case2Start {} {
    global run_sim case2stime displayclock N0 N1 N2 N3 N4 Rx Ry Rz Rt

    .case2.status configure -text "Starte..."    
    if { [expr abs($run_sim)] == 2 } { set restart 1 } else { set restart 0 }
    set run_sim 0
    update

    disableStarts
    disableCase1
    disableCase2

    part deleteall
    set inp [open "case2.blk" r]
    while { [blockfile $inp read auto] != "eof" } {}
    close $inp
    for {set i $N0} {$i < [expr $N0+$N1+$N2+$N3+$N4]} {incr i} { part $i fix }
    puts "Read & fixed [setmd n_part] particles..."

    if { $restart == 0 } { imd_reconnect case2  }

    set run_sim 2
    Case2Wall1ChargeChange 0
    Case2Wall2ChargeChange 0
    Case2Wall3ChargeChange 0

    after 3000
    .case2.status configure -text "...3..."
    update; after 1000
    .case2.status configure -text "...2..."
    update; after 1000
    .case2.status configure -text "...1..."
    update; after 1000

    enableCase2
    enableStarts

    .case2.status configure -text "Los!"
    set case2stime [clock seconds]
    set displayclock 1
}

proc Case1BjerrumChange {l} {
    global run_sim
    if { $run_sim == 1} { 
	inter coulomb $l dh 0 [expr 2*$l]
	puts "Changed Bjerrum length to $l."
    }
}

proc Case1TempChange {T} {
    global run_sim
    if { $run_sim == 1} { 
	setmd temp $T
	puts "Changed temperature to $T."
    }
}

proc Case1PotChange { } {
    global run_sim butFH
    if { $run_sim == 1} {
	if { $butFH == "F" } {
	    set tsik [setmd temp]; setmd temp 0.0
	    inter 0 harmonic [expr 4*[lindex [inter 0] 2]] [expr 0.25*[lindex [inter 0] 3]]
	    integrate 100; imd positions; setmd temp $tsik; integrate 100; imd positions
	    inter 0 FENE [expr 0.125*[lindex [inter 0] 2]] [expr 8*[lindex [inter 0] 3]] 
	}
	if { $butFH == "H" } { inter 0 harmonic [expr 2*[lindex [inter 0] 2]] [expr 0.5*[lindex [inter 0] 3]] }
    }
}

proc Case2Wall1ChargeChange {c} {
    global run_sim N0 N1
    if { $run_sim == 2} { 
	for {set i [expr $N0]} {$i < [expr $N0+$N1]} {incr i} { part $i q $c }
	puts "Charged Wall 1 with $c (with Bjerrum = [lindex [lindex [inter coulomb] 0] 1])."
    }
}

proc Case2Wall2ChargeChange {c} {
    global run_sim N0 N1 N2
    if { $run_sim == 2 } { 
	for {set i [expr $N0+$N1]} {$i < [expr $N0+$N1+$N2]} {incr i} { part $i q $c } 
	puts "Charged Wall 2 with $c (with Bjerrum = [lindex [lindex [inter coulomb] 0] 1])."
    }
}

proc Case2Wall3ChargeChange {c} {
    global run_sim N0 N1 N2 N3
    if { $run_sim == 2} {
	set ind [expr $N0+$N1+$N2]; set le [expr round(sqrt($N3))]; set maxdist [expr $le/sqrt(2)]
	for {set i 0} {$i < $N3} {incr i} { 
	    set x [expr $i % $le]; set y [expr $i / $le]; set dist [expr sqrt(pow(0.5*$N3-$x,2)+pow(0.5*$N3-$y,2))]
	    part $ind q [expr $c*sqrt($dist/$maxdist)]; incr ind
	}
	puts "Charged Wall 3 with 0...$c (with Bjerrum = [lindex [lindex [inter coulomb] 0] 1])."
    }
}

#############################################################
#      Integration                                          #
#############################################################

disableCase1
disableCase2

set displayclock 0
set run_sim 0
set firework 250

while { 1 } {
    if { $run_sim > 0 } {
	integrate 10
	while { [catch {imd positions} res] } { puts "positions failed: $res, retrying."; after 500 }
    } {
	after 100
	while { [catch {imd listen 1} res] } { puts "listen failed: $res, retrying."; after 500 }
    }
    if { $displayclock } {
	set nbh [analyze nbhood $Rx $Ry $Rz [expr 1.25*$Rt]]
	if { [llength $nbh] > $N4 } {
	    set etime [expr [clock seconds] - $case2stime]
	    disableCase2; set run_sim -2
	    .case2.status configure -text "Gewonnen nach [clock format $etime -format "%M:%S"]!"; update
	    for {set i [expr $N0]} {$i < [expr $N0+$N1+$N2+$N3 +$N4]} {incr i} { part $i unfix }
	    setmd temp 1.0
	    for {set i 0} {$i < $firework} {incr i} {
		integrate 10; while { [catch {imd positions} res] } { puts "positions failed: $res, retrying."; after 500 }
	    }
	    set highscore [read_highscore]
	    if { $etime < [lindex [lindex $highscore end] 0] } {
		for {set i 0} {$i < [llength $highscore]} {incr i} {
		    if { $etime < [lindex [lindex $highscore $i] 0] } {
			set victor [clock format [clock seconds]]
			set victory [list $etime $victor]
			break
		    }
		}
		for {set j $i} {$j < [llength $highscore]} {incr j} {
		    if { [expr $j+1] < [llength $highscore] } { set old [lindex $highscore $j] }
		    set highscore [lreplace $highscore $j $j $victory]
		    if { [expr $j+1] < [llength $highscore] } { set victory $old }
		}
		puts "New Highscore (Rank \#$i)!!!"
		write_highscore $highscore
	    } 
	} else {
	    set etime [expr [clock seconds] - $case2stime]
	    .case2.status configure -text [clock format $etime -format "%M:%S"]
	    if { $etime > $case2timeout } {
		.case2.status configure -text "Die Zeit ist um!"
		set displayclock 0
		set run_sim -2
	    }
	}
    }
    update
}