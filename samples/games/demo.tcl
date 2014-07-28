#!/usr/bin/env Espresso $*
#  This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
#  It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
#  and by which you are legally bound while utilizing this file in any form or way.
#  There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  You should have received a copy of that license along with this program;
#  if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
#  write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
#  Copyright (c) 2002-2013; all rights reserved unless otherwise stated.
# 
#############################################################
#                                                           #
#  Demo of a buckyball                                      #
#                                                           #
#                                                           #
#  Created:       26.08.2003 by AxA & BAM                   #
#                                                           #
#############################################################

set case2timeout 120

set startTemp 1.0
set startBjerrum 8.0

#############################################################
#  Initialization                                           #
#############################################################

set run_sim 0
thermostat langevin 1 .1

#############################################################
#  Setup GUI                                                #
#############################################################

bind . <Destroy> onExitHandler

# popup window if desired
proc popup {title text} {
    toplevel .tag
    label .tag.title -text "$title" -justify center -font TkCaptionFont
    text .tag.txt
    .tag.txt insert 1.0 "$text"
    .tag.txt configure -state disabled
    button .tag.btn -text "Ok" -command "destroy .tag"
    pack .tag.title .tag.txt .tag.btn
    wm geometry .tag +60+0
    wm title .tag $title
    tkwait window .tag
    update
}

# popup-input-window
proc ret_input {win} {
    global input
    set input [.tag.inp get]
    destroy $win
}
proc get_input { } {
    toplevel .tag
    label .tag.title -justify center -text "Herzlichen Glückwunsch!!!

Sie haben sich für einen Eintrag in der Highscore-Tabelle qualifiziert!
Bitte geben Sie nun Ihren Namen ein, so, wie er dort erscheinen soll:
"
    entry .tag.inp
    button .tag.btn -text "OK" -justify left -command "ret_input .tag"
    pack .tag.title -fill both -side top -in .tag
    pack .tag.inp -fill both -side top  -in .tag
    pack .tag.btn -side top
    wm geometry .tag +60+0
    tkwait window .tag
    update
}

proc set_status {stat} {
    .case2.status configure -text $stat
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
scale .case1.slider -orient h -from 0 -to 10 \
    -resolution [expr 10/100.] -command Case1BjerrumChange
pack .case1.slider .case1.label -fill x -in .case1

label .case1.label2 -justify left -text "Temperatur"
scale .case1.slider2 -orient h -from 0 -to 1.5 \
    -resolution [expr 2/100.] -command Case1TempChange
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
label .case2.label1 -justify left -text "Ladung Wand mitte"
scale .case2.slider1 -orient h -from -1.5 -to 1.5 \
    -resolution [expr 3/5.] -command Case2Wall1ChargeChange
pack .case2.slider1 .case2.label1 -fill x -in .case2

label .case2.label2 -justify left -text "Ladung Wand links"
scale .case2.slider2 -orient h -from -1.5 -to 1.5 \
    -resolution [expr 3/5.] -command Case2Wall2ChargeChange
pack .case2.slider2 .case2.label2 -fill x -in .case2

label .case2.label3 -justify left -text "Ladung Wand oben"
scale .case2.slider3 -orient h -from -1.5 -to 1.5 \
    -resolution [expr 3/5.] -command Case2Wall3ChargeChange
pack .case2.slider3 .case2.label3 -fill x -in .case2

label .case2.status -justify left -text ""
pack .case2.status -fill x -in .case2

proc disableCase1 {} {
    global disabledfg butFH startBjerrum startTemp
    .case1.slider set $startBjerrum
    .case1.label configure -state disabled
    .case1.slider configure -state disabled -foreground $disabledfg 
    .case1.slider2 set $startTemp
    .case1.label2 configure -state disabled
    .case1.slider2 configure -state disabled -foreground $disabledfg 
    set butFH F
    .case1.labelF configure -state disabled
    .case1.butF configure -state disabled -foreground $disabledfg 
    .case1.labelH configure -state disabled
    .case1.butH configure -state disabled -foreground $disabledfg 
}

proc enableCase1 {} {
    global normalfg butFH startBjerrum startTemp
    .case1.label configure -state normal
    .case1.slider configure -state normal -foreground $normalfg 
    .case1.label2 configure -state normal
    .case1.slider2 configure -state normal -foreground $normalfg 
    .case1.labelF configure -state normal
    .case1.butF configure -state normal -foreground $normalfg 
    .case1.labelH configure -state normal
    .case1.butH configure -state normal -foreground $normalfg 
    set butFH F

    set_status ""
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

    set_status ""
    set displayclock 0
}

proc enableCase2 {} {
    global normalfg
    .case2.label1 configure -state normal
    .case2.slider1 configure -state normal -foreground $normalfg
    .case2.label2 configure -state normal
    .case2.slider2 configure -state normal -foreground $normalfg
    .case2.label3 configure -state normal
    .case2.slider3 configure -state normal -foreground $normalfg
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
  while {[eof $f]==0} {
    if { [gets $f inp] > 0 } {
      lappend highscore [list [lindex $inp 0] [concat [lrange $inp 1 end]]]
    }
  }
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

proc disp_highscore { highscore } {
    set ret " \# | Zeit | Name \n=======================================================\n\n"
    for {set i 0} {$i < [llength $highscore]} {incr i} {
	set tmp [lindex $highscore $i]
	set ret "$ret[format %2d [expr $i+1]] | [format %3d [lindex $tmp 0]]s | [lindex $tmp 1]\n"
    }
    return $ret
}

proc onExitHandler {} {
    # terminate vmd
    eval "exec rm -f [glob -nocomplain .lock*]"    
}

proc imd_reconnect {case} {
    imd disconnect

    eval "exec rm -f [glob -nocomplain .lock*]"

    puts "creating $case.vtf"
    set f [open "$case.vtf" "w"]
    writevsf $f radius {0 2 1 2 2 3 3 2 4 2}
    writevcf $f
    close $f

    while { [catch {imd connect 10000} res] } {
	set_status  "Waiting for port to unbind...\nerror: $res"
	update
	after 500
    }
    exec vmd -size 1000 750 -e $case-vmd.script -startup ./vmdrc &
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
    global run_sim startTemp startBjerrum

    set_status "Starte..."    
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

    set run_sim 1

    enableCase1
    enableStarts

    Case1BjerrumChange $startBjerrum
    Case1TempChange $startTemp

    set_status "Fertig"
}

set displayclock 0
proc Case2Start {} {
    global run_sim case2stime displayclock N0 N1 N2 N3 N4 Rx Ry Rz Rt

    set_status "Starte..."    
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

    thermostat langevin 0 0.5

    for {set i $N0} {$i < [expr $N0+$N1+$N2+$N3+$N4]} {incr i} { part $i fix }
    puts "Read & fixed [setmd n_part] particles..."

    if { $restart == 0 } { imd_reconnect case2  }

    set run_sim 2
    Case2Wall1ChargeChange 0
    Case2Wall2ChargeChange 0
    Case2Wall3ChargeChange 0

    set_status "...3..."
    update; after 1000
    set_status "...2..."
    update; after 1000
    set_status "...1..."
    update; after 1000

    enableCase2
    enableStarts

    set_status "Los!"
    set case2stime [clock seconds]
    set displayclock 1
}

proc checkLambda {} {
    global Bjerrum temp
    if { $Bjerrum > 10 || $temp > 1 } {
	puts "Increasing friction and slowing down"
	thermostat langevin $temp 5
	integrate 50
    }
}

proc Case1BjerrumChange {l} {
    global run_sim Bjerrum
    set Bjerrum $l
    if { $run_sim == 1} { 
	inter coulomb $l dh 0 [expr 2*$l]
	puts "Changed Bjerrum length to $l."
	checkLambda
    }
}

proc Case1TempChange {T} {
    global run_sim temp
    set temp $T
    if { $run_sim == 1} { 
	thermostat langevin $T 1
	puts "Changed temperature to $T."
	checkLambda
    }
}

proc Case1PotChange { } {
    global run_sim butFH temp
    .case1.butH configure -state disabled
    .case1.butF configure -state disabled

    if { $run_sim == 1} {
	if { $butFH == "F" } {
	    set_status "Bindungen straffen"

	    set param1 [lindex [inter 0] 2]
	    set param2 [lindex [inter 0] 3]
	    set tg [setmd gamma]
	    foreach factor "2 4 6" {
		puts "fac $factor"
		inter 0 harmonic [expr $factor*$param1] [expr 1./$factor*$param2]
		integrate 50; imd positions
		integrate 50; imd positions
		integrate 50; imd positions
		integrate 50; imd positions
		# by this
		thermostat langevin $temp [expr 1./$factor]
	    }

	    set_status "Viskositaet erhöhen..."

	    foreach gamma "10 20 30 50" {
		puts "gamma $gamma"
		thermostat langevin $temp $gamma
		integrate 50; imd positions
		integrate 50; imd positions
		integrate 50; imd positions
		integrate 50; imd positions
	    }

	    set_status "Bindungstyp umschalten..."

	    inter 0 FENE [expr 0.5*$param1] [expr 2*$param2]

	    set_status "Viskositaet senken"

	    foreach gamma "40 35 30 20 10 $tg" {
		puts "gamma $gamma"
		thermostat langevin $temp $gamma
		integrate 50; imd positions
		integrate 50; imd positions
		integrate 50; imd positions
		integrate 50; imd positions
	    }
	}
	if { $butFH == "H" } { inter 0 harmonic [expr 2*[lindex [inter 0] 2]] [expr 0.5*[lindex [inter 0] 3]] }
    }

    .case1.butH configure -state active
    .case1.butF configure -state active
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
	set ind [expr $N0+$N1+$N2]; set le [expr round(sqrt($N3))]; set maxdist [expr $le/sqrt(2.)]
	for {set i 0} {$i < $N3} {incr i} { 
	    set x [expr $i % $le]; set y [expr $i / $le]; set dist [expr sqrt(pow(0.5*$le-$x,2)+pow(0.5*$le-$y,2))]
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
set firework 250

set etime 91; set i 4
set highscore [read_highscore]

while { 1 } {
    if { $run_sim > 0 } {
 	integrate 10
	while { [catch {imd positions} res] } { puts "positions failed: $res, retrying."; after 500 }
    } {
	after 100
	while { [catch {imd listen 1} res] } { puts "listen failed: $res, retrying."; after 500 }
    }
    if { $displayclock } {
	set hit 0
	for {set part 0} { $part < $N0 } { incr part } {
	    set pos [part $part pr pos]
	    set dy [expr [lindex $pos 1] - $Ry]
	    set dz [expr [lindex $pos 2] - $Rz]
	    if {[expr sqrt($dy*$dy + $dz*$dz)*0.9] < $Rt } {
		set hit 1
	    }
	}
	if { $hit } {
	    set etime [expr [clock seconds] - $case2stime]
	    disableCase2; set run_sim -2
	    set_status "Gewonnen nach [clock format $etime -format "%M:%S"]!"
	    for {set i [expr $N0]} {$i < [expr $N0+$N1+$N2+$N3+$N4]} {incr i} { part $i unfix }
	    thermostat langevin 1 1
	    for {set i 0} {$i < $firework} {incr i} {
		integrate 10; while { [catch {imd positions} res] } { puts "positions failed: $res, retrying."; after 500 }
	    }
	    set highscore [read_highscore]
	    if { $etime < [lindex [lindex $highscore end] 0] } {
		for {set i 0} {$i < [llength $highscore]} {incr i} {
		    if { $etime < [lindex [lindex $highscore $i] 0] } { break }
		}
		get_input; set victor "$input ([clock format [clock seconds] -format "%d.%m.%y @ %Rh"])"; # "%a, %d.%m.%Y @ %T"
		set victory [list $etime $victor]
		for {set j $i} {$j < [llength $highscore]} {incr j} {
		    if { [expr $j+1] < [llength $highscore] } { set old [lindex $highscore $j] }
		    set highscore [lreplace $highscore $j $j $victory]
		    if { [expr $j+1] < [llength $highscore] } { set victory $old }
		}
		incr i

		popup "HIGHSCORE" "\
Sie haben es tatsächlich geschafft, das geladene Objekt\ 
in einer Rekordzeit ins Ziel zu steuern - phantastisch!\ 
                                                       \ 
Mit Ihrem beachtlichen Wert von [format %3d $etime]s haben Sie sich den\ 
----->  [format %2d $i]. Platz in der Highscore-Liste   <-----   \ 
erobert und Ihren Eintrag redlich verdient!\ 
                                                       \ 
Herzlichen Glückwunsch!!!                              \ 
                                                       \ 
[disp_highscore $highscore]\
"

	       puts "New Highscore (Rank \#$i)!!!"
               write_highscore $highscore

               catch {exec /usr/local/bin/mplayer Logo-Espresso-Bat-HQ.avi >& /dev/null}
	    } else {
               popup "Gewonnen!" "\ 
Sie haben es tatsächlich geschafft, das geladene Objekt\ 
in der gegebenen Zeit ins Ziel zu steuern - super!     \ 
                                                       \ 
Beinahe hätte Ihre Zeit von [format %3d $etime]s auch für einen Eintrag\ 
in die Highscore-Liste gereicht -- vielleicht probieren\ 
Sie es daher gleich noch einmal?!                      \ 
                                                       \ 
[disp_highscore $highscore]\
"

               catch {exec /usr/local/bin/mplayer Logo-Espresso-Bat-HQ.avi >& /dev/null}
	    }
	} else {
	    set etime [expr [clock seconds] - $case2stime]
	    set_status [clock format $etime -format "%M:%S"]
	    if { $etime > $case2timeout } {
		set_status "Die Zeit ist um!"
		set displayclock 0
		set run_sim -2
                popup "Schade!" "\ 
Leider hat es diesmal nicht geklappt, das Ladungsobjekt\ 
in der gegebenen Zeit ins Ziel zu steuern - schade! :-(\ 
                                                       \ 
Aber mit ein wenig Übung schaffen Sie  es bestimmt auch\ 
sich mit einem Eintrag in die Highscore-Liste verewigen\ 
zu können -- vielleicht probieren Sie es daher sogleich\ 
noch einmal?!                                          \ 
                                                       \ 
[disp_highscore $highscore]\
"
	    }
	}
    }
    update
}
