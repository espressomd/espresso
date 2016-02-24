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
#
# procedures for loading and rendering molecular dynamics
# utils

global numps povraycmd connectcom pscom imgname
global scsrender scsnumber scsname scsstyle scsglob

# the name of the povray files (exactly one %d must be included 
# for the filecounter). Used for writing out and reading in for making a film
set scsname "povs/poly%04d.pov"
# if set to 1, scshot really renders
set scsrender 0
# number of next snapshot to render
set scsnumber 1
# style when rendering
set scsstyle {CPK 1.000000 0.300000 8.000000 6.000000}
# number of povray processes to spawn
set numps 2
# command to check if a process is alive
set pscom "ps -p %d"
# povray command. The first %s gives the name of the povray,
# a %*d the image counter
set povraycmd "x-povray +W768 +H576 -I%s -Otgas/poly%05d.tga >&/dev/null"
# command to make a film from the tgas. First %d gives number of images,
# %s the filmname
set connectcom "dmconvert -v -fqt -pvideo,comp=qt_anim,inrate=25 \
		-ntgas/poly\#\#\#\#\#.tga,start=1,end=%d \
                tgas/poly\#\#\#\#\#.tga %s.mov"

################### commands for loading files from poly2pdb

proc loadseries {name {step 1} {start 0} {end 10000}} {
    if [regexp "(.*)%.*d" $name dummy pref] {
	# load with an additional counter
	set fmt "$name%04d.pdb"
	set psffile "$pref.psf"
	if { ![file exists $psffile] } {
		set psffile [format "$name.psf" 0]
		if { ![file exists $psffile] } {
			error "no psf file found"
		}
	}
	mol load psf $psffile
	set ser $start
	while { ($ser < $end) &&
		[file exists [format $fmt $ser 0]]} {
	    set index 0
	    while { ($index < 10000) &&
		    [file exists [set file [format $fmt $ser $index]]] } {
		animate read pdb $file
		incr index $step
	    }
	    incr ser
	}
	return
    }
    if {! [file exists "$name.psf"] } { error "no psf file found" }
    mol load psf "$name.psf"
    set num $start
    while { ($num < $end) && \
		[file exists [set file [format "$name%04d.pdb" $num]]] } {
        animate read pdb $file
	incr num $step
    }
}

################### commands for rendering

global scsprevst
set scsprevst {Lines 1.0}

proc applystyle {style} {
    eval "mol modstyle 0 top $style"
}

proc scshot {{job render} {arg {}}} {
    global scsrender scsnumber scsname scsstyle scsprevst
    switch $job {
	start {
	    set scsrender 1
	    set scsnumber 1
	    if { $arg  != {} } {
		set scsname $arg
	    }
	    set scsprevst [molinfo top get "{rep 0}"]
	    applystyle $scsstyle
	}
	restart {
	    set scsrender 1
	    applystyle $scsstyle
	}
	stop {
	    set scsrender 0
	    applystyle $scsprevst
	}
	wind {
	    if { $scsrender != 1 } {
		puts "Warning: rewinding while not rendering"
	    }
	    switch $arg {
		{} {
		    incr scsnumber -1
		}
		start {
		    set scsnumber 1
		}
		last {
		    set scsnumber 1
		    while [file exists [format $scsname $scsnumber]] {
			incr scsnumber
		    }
		}
		default {
		    set scsnumber [expr $scsnumber + $arg]
		}
	    }
	}
	render {
	    # rendering
	    if $scsrender {
		render POV3 [format $scsname $scsnumber]
		incr scsnumber
	    } {
		display update
	    }
	}
	default {
	    error "unknown command $job"
	}
    }
}

################## transformation

proc vecdot {v1 v2} {
    lassign $v1 x y z q
    lassign $v2 u v w p
    return [expr $x*$u + $y*$v + $z*$w + $p*$q]
}

proc vecnorm {v} {
    return [expr sqrt([vecdot $v $v])]
}

proc vecnormalize {v} {
    return [vecscale $v [expr 1/[vecnorm $v]]]
}

proc veccross {b1 b2} {
    lassign $b1 b11 b12 b13
    lassign $b2 b21 b22 b23

    set b31 [expr $b12*$b23 - $b13*$b22]
    set b32 [expr $b13*$b21 - $b11*$b23]
    set b33 [expr $b11*$b22 - $b12*$b21]
    return [list $b31 $b32 $b33 0]
}

proc matmult {m1 m2} {
    lassign $m1 z1 z2 z3 z4
    lassign [transtranspose $m2] y1 y2 y3 y4
    return [list [list [vecdot $z1 $y1] [vecdot $z1 $y2] [vecdot $z1 $y3] [vecdot $z1 $y4]]\
		 [list [vecdot $z2 $y1] [vecdot $z2 $y2] [vecdot $z2 $y3] [vecdot $z2 $y4]]\
		 [list [vecdot $z3 $y1] [vecdot $z3 $y2] [vecdot $z3 $y3] [vecdot $z3 $y4]]\
		 [list [vecdot $z4 $y1] [vecdot $z4 $y2] [vecdot $z4 $y3] [vecdot $z4 $y4]]]
}

proc matscale {m d} {
    lassign $m v1 v2 v3 v4
    return [list [vecscale $v1 $d] [vecscale $v2 $d] \
		[vecscale $v3 $d] [vecscale $v4 $d]]
}

proc matadd {m1 m2} {
    lassign $m1 v11 v12 v13 v14
    lassign $m2 v21 v22 v23 v24    
    return [list [vecadd $v11 $v21] [vecadd $v12 $v22] \
		[vecadd $v13 $v23] [vecadd $v14 $v24]]
}

proc matsub {m1 m2} {
    lassign $m1 v11 v12 v13 v14
    lassign $m2 v21 v22 v23 v24    
    return [list [vecsub $v11 $v21] [vecsub $v12 $v22] \
		[vecsub $v13 $v23] [vecsub $v14 $v24]]
}

proc check_vectors {v1 v2 v3} {
    puts "cos: [vecdot $v1 $v2] [vecdot $v2 $v3] [vecdot $v1 $v3]"
    puts "len: [vecdot $v1 $v1] [vecdot $v2 $v2] [vecdot $v3 $v3]"
}

proc rotx {ang} {
    set c [expr cos($ang)]
    set s [expr sin($ang)]
    return [list "$c $s 0 0" "[expr -$s] $c 0 0" {0 0 1 0} {0 0 0 1}]
}

proc roty {ang} {
    set c [expr cos($ang)]
    set s [expr sin($ang)]
    return [list "$c 0 $s 0" {0 1 0 0} "[expr -$s] 0 $c 0" {0 0 0 1}]
}

proc rotz {ang} {
    set c [expr cos($ang)]
    set s [expr sin($ang)]
    return [list {1 0 0 0} "0 $c $s 0" "0 [expr -$s] $c 0" {0 0 0 1}]
}

proc gettrafo {} {
    # translation
    set cm [molinfo top get center_matrix]
    lassign [lindex $cm 0] v1 v2 v3 v4
    set cmv [list [lindex $v1 3] [lindex $v2 3] [lindex $v3 3]]
 
    # rotation
    set rm [lindex [molinfo top get rotate_matrix] 0]

    # scaling
    set sm [molinfo top get scale_matrix]
    set sf [lindex [lindex [lindex $sm 0] 0] 0]

    # global
    set gm [lindex [molinfo top get global_matrix] 0]
    return [list $cmv $rm $sf $gm]
}

proc settrafo {mats} {
    lassign $mats cmv rm sf gm

    # translation
    set cm [list "1 0 0 [lindex $cmv 0]" "0 1 0 [lindex $cmv 1]"\
		"0 0 1 [lindex $cmv 2]" "0 0 0 1"]

    # scaling
    set sm [list "$sf 0 0 0" "0 $sf 0 0" "0 0 $sf 0" {0 0 0 1}]
 
    molinfo top set {center_matrix rotate_matrix scale_matrix global_matrix} \
	[list $cm $rm $sm $gm]
}

proc transto {steps mats rotv {animfac 0}} {
    set PI 3.14159265358979323846
    set animcount 0.0
    set state [gettrafo]
    lassign $state cmv rm sf gm
    lassign $mats tcmv trm tsf tgm

    # most tricky part: rotation
    # first determine basis of the rotation plane per vector
    lassign $trm tx ty tz
    lassign $rm  sx sy sz

    # determine coordinate basis along one rotation
    set ambigous 0
    set txortho [vecsub $tx [vecscale $sx [vecdot $sx $tx]]]
    if { [vecnorm $txortho] > 1e-4 } {
	set b1 $sx
	set b2 [vecnormalize $txortho]
	set b3 [veccross $b1 $b2]

	set ta [expr atan2([vecdot $tx $b2], [vecdot $tx $b1])]

	# rotation of the second vector gives rotation of b3
	# sy in <b2, b3>
	set sb [expr atan2([vecdot $b3 $sy], [vecdot $b2 $sy])]
	# ty in <R b2,  R b3> = <R b2 , b3>
	set Rb2 [vecadd [vecscale $b2 [expr cos($ta)]] [vecscale $b1 [expr -sin($ta)]]]
	set tb [expr atan2([vecdot $b3 $ty], [vecdot $Rb2 $ty]) - $sb]
    } {
	# check wether x is mirrored - complain
	if { [vecnorm [vecsub $tx $sx]] > 1e-4 } {
	    set ambigous 1
	    set ambvec x
	} {
	    set tyortho [vecsub $ty [vecscale $sy [vecdot $sy $ty]]]
	    if { [vecnorm $tyortho] > 1e-4 } {
		set b1 $sy
		set b2 [vecnormalize $tyortho]
		
		set ta [expr atan2([vecdot $ty $b2], [vecdot $ty $b1])]
		
		# as x doesn't move, this is a rotation around x || b3
		set tb 0
	    } {
		if { [vecnorm [vecsub $ty $sy]] > 1e-4 } {
		    set ambigous 1
		    set ambvec y
		} {
		    # x AND y do not move - so no rotation at all
		    set b1 $sx
		    set b2 $sy
		    set b3 $sz
		    set ta 0
		    set tb 0
		}
	    }
	}
    }
    
    if { $ambigous } {
	puts "Ambigous rotation in $ambvec axis, please add keyframe. Skipping"
	settrafo [list $cmv $rm $sf $gm]
	return
    }

    # puts "rot_xy $ta rot_yz $tb"

    # either 1 0 0 in case 1 or 0 0 1 in case x doesn´t move
    set bx1 [vecdot $sx $b1]
    set bx2 [vecdot $sx $b2]
    set bx3 [vecdot $sx $b3]

    set by1 [vecdot $sy $b1]
    set by2 [vecdot $sy $b2]
    set by3 [vecdot $sy $b3]

    set bz1 [vecdot $sz $b1]
    set bz2 [vecdot $sz $b2]
    set bz3 [vecdot $sz $b3]
	
    set sc [expr 1.0 / $steps]
    set dcmv [vecscale [vecsub $tcmv $cmv] $sc]
    set dsf [expr ($tsf - $sf) * $sc]
    set dgm [matscale [matsub $tgm $gm] $sc]
    set da [expr $ta*$sc]
    set db [expr $tb*$sc]

    for {set s 0} { $s < $steps} {incr s} {
	# calculation transition rotation
	set wa [expr $da*$s]
	set wb [expr $db*$s]
	set ca [expr cos($wa)]
	set sa [expr sin($wa)]
	set cb [expr cos($wb)]
	set sb [expr sin($wb)]

	set Rb1 [vecadd [vecscale $b1 $ca]            [vecscale $b2 [expr $sa]]]
	set Rb2 [vecadd [vecscale $b1 [expr -$sa*$cb]] [vecscale $b2 [expr $ca*$cb]] [vecscale $b3 $sb]]
	set Rb3 [vecadd [vecscale $b1 [expr $sa*$sb]] [vecscale $b2 [expr -$ca*$sb]] [vecscale $b3 $cb]]
	
	set Rsx [vecadd [vecscale $Rb1 $bx1] [vecscale $Rb2 $bx2] [vecscale $Rb3 $bx3]]
	set Rsy [vecadd [vecscale $Rb1 $by1] [vecscale $Rb2 $by2] [vecscale $Rb3 $by3]]
	set Rsz [vecadd [vecscale $Rb1 $bz1] [vecscale $Rb2 $bz2] [vecscale $Rb3 $bz3]]

	set rm [list $Rsx $Rsy $Rsz "0 0 0 1"]

	# add additional axial rotations
	set rm [matmult [rotx [expr $s*[lindex $rotv 0]]] $rm]
	set rm [matmult [roty [expr $s*[lindex $rotv 1]]] $rm]
	set rm [matmult [rotz [expr $s*[lindex $rotv 2]]] $rm]
	
	set cmv [vecadd $cmv $dcmv]
	set sf [expr $sf + $dsf]
	set gm [matadd $gm $dgm]
	settrafo [list $cmv $rm $sf $gm]
	set animcount [expr $animcount + $animfac]
	while {$animcount >= 1.0} {
	    set animcount [expr $animcount - 1.0]
	    animate next
	}
	while {$animcount < 0.0} {
	    set animcount [expr $animcount + 1.0]
	    animate prev
	}
	scshot
    }
}

################### recording of transformations
global playlist plfrm plpos pldur plmat

# protect playlist for debugging...
if { ! [info exists playlist] } {
    set playlist {}
}
set plfrm 0
set plpos 0
set pldur 1
set plrx  0
set plry  0
set plrz  0

proc playlist {{job play} {arg {}}} {
    global playlist plfrm plpos pldur plrx plry plrz M_PI
    switch $job {
	save {
	    if {$arg == {}} {
		error "playlist save needs a filename to save to"
	    }
	    set fd [open $arg w]
	    foreach entry $playlist {
		puts $fd "$entry"
	    }
	    close $fd
	}
	load {
	    if {$arg == {}} {
		error "playlist load needs a filename to load"
	    }
	    set playlist ""
	    set fd [open $arg]
	    while {! [eof $fd]} {
		gets $fd line
		if {$line != ""} {
		    lappend playlist $line
		}
	    }
	    close $fd
	    set plfrm 0
	    set plpos 0
	    set pldur 1
	    set plrx  0
	    set plry  0
	    set plrz  0
	}
	play {
	    set len [expr [llength $playlist] - 1]
	    if {$len == 0} return
	    lassign [lindex $playlist 0] olddur oldmat oldrotv oldfrm
	    set plpos 0
	    while {$plpos < $len} {
		incr plpos
		set entry [lindex $playlist $plpos]
		lassign $entry pldur plmat plrotv plfrm
		animate goto $oldfrm
		settrafo $oldmat
		set animfac [expr double($plfrm - $oldfrm)/$pldur]
		lassign $plrotv rx ry rz
		set plrotv [list [expr $rx*2*$M_PI/$pldur] \
				[expr $ry*2*$M_PI/$pldur] \
				[expr $rz*2*$M_PI/$pldur]]
		transto $pldur $plmat $plrotv $animfac
		set olddur $pldur
		set oldfrm $plfrm
		set oldmat $plmat
	    }
	}
	editor {
	    showplwin
	}
	default {
	    error "unknown command $job"
	}
    }
}

proc showplwin {} {
    global playlist plpos pldur plfrm plrx plry plrz
    if [winfo exists .pled] {
	destroy .pled
    }
    toplevel .pled
    wm title .pled "Playlist Editor"

    frame .pled.f
    entry .pled.ep -textvariable plpos
    button .pled.pp -text "prev" -command {incr plpos -1; adjpos}
    button .pled.np -text "next" -command {incr plpos; adjpos}
    pack .pled.pp .pled.ep .pled.np -fill x -side left -in .pled.f

    frame .pled.ff
    label .pled.lf -text "frame   "
    entry .pled.ef -textvariable plfrm
    pack .pled.ef .pled.lf -fill x -side right -in .pled.ff

    frame .pled.fd
    label .pled.ld -text "duration"
    entry .pled.ed -textvariable pldur
    pack .pled.ed .pled.ld -fill x -side right -in .pled.fd

    frame .pled.frx
    label .pled.lrx -text "rock x"
    entry .pled.erx -textvariable plrx
    pack .pled.erx .pled.lrx -fill x -side right -in .pled.frx

    frame .pled.fry
    label .pled.lry -text "rock y"
    entry .pled.ery -textvariable plry
    pack .pled.ery .pled.lry -fill x -side right -in .pled.fry

    frame .pled.frz
    label .pled.lrz -text "rock z"
    entry .pled.erz -textvariable plrz
    pack .pled.erz .pled.lrz -fill x -side right -in .pled.frz

    button .pled.set -text "set" -command setplf
    button .pled.ins -text "insert" -command insplf
    button .pled.del -text "delete" -command delplf
    button .pled.pl -text "play" -command {playlist play}

    pack .pled.f .pled.ff .pled.fd .pled.frx .pled.fry .pled.frz\
	.pled.pl .pled.set .pled.ins .pled.del -fill x -side top -in .pled
    bind .pled.ep <Return> adjpos
    bind .pled.ef <Return> adjfrm
}

proc adjfrm {} {
    global plfrm
    animate goto $plfrm
}

proc adjpos {} {
    global playlist pldur plfrm plpos plmat plrx plry plrz
    set b [llength $playlist]
    if {$plpos < 0} {
	set plpos 0
    }
    if {$plpos > $b} {
	set plpos $b
    }
    set entry [lindex $playlist $plpos]
    if {[llength $entry] != 4} {
	return
    }
    lassign $entry dur mat rotv frm
    set pldur $dur
    set plfrm $frm
    set plmat $mat
    lassign $rotv rx ry rz
    set plrx $rx
    set plry $ry
    set plrz $rz
    animate goto $plfrm
    settrafo $plmat
}

proc insplf {} {
    global playlist plpos pldur plfrm plrx plry plrz
    set pre [lrange $playlist 0 [expr $plpos]]
    set post [lrange $playlist [expr $plpos] end]
    set playlist [concat $pre $post]
    
}

proc delplf {} {
    global playlist plpos pldur plfrm plrx plry plrz
    set pre [lrange $playlist 0 [expr $plpos - 1]]
    set post [lrange $playlist [expr $plpos + 1] end]
    set playlist [concat $pre $post]
    
}

proc setplf {} {
    global playlist plpos pldur plfrm plrx plry plrz
    set pre [lrange $playlist 0 [expr $plpos - 1]]
    set post [lrange $playlist [expr $plpos + 1] end]

    set newentry [list $pldur [gettrafo] [list $plrx $plry $plrz] $plfrm]
    lappend pre $newentry
    set playlist [concat $pre $post]
}

global makeprocesses
set makeprocesses {}

proc waitforps {} {
    global makeprocesses numps pscom
    while {[llength $makeprocesses] >= $numps} {
	set tmpps {}
	foreach ps $makeprocesses {
	    if {[catch "exec [format $pscom $ps]"] == 0} {
		lappend tmpps $ps
	    }
	}
	set makeprocesses $tmpps
	puts "waiting for a process to finish..."
	after 1000
    }
}

proc waitforfinish {} {
    global makeprocesses pscom
    while {[llength $makeprocesses] != 0} {
	set tmpps {}
	foreach ps $makeprocesses {
	    if {[catch "exec [format $pscom $ps]"] == 0} {
		lappend tmpps $ps
	    }
	}
	set makeprocesses $tmpps
	puts "waiting for all processes to finish"
	after 1000
    }
}

proc makefilm {film {filelist {}}} {
    global povraycmd connectcom scsname imgname
    puts "$filelist"
    if {$filelist == {}} {
	# use all images that fit the scsname mask
	set num 1
	while [file exists [set fn [format $scsname $num]]] {
	    lappend filelist $fn
	    incr num
	}
    }
    set num 1
    foreach ifile $filelist {
	global makeprocesses
	waitforps
	puts "rendering $ifile ($num)"
	lappend makeprocesses [lindex [eval \
          "exec [format $povraycmd $ifile $num] &"] 0]
	incr num
    }
    waitforfinish
    eval "exec [format $connectcom $num $film]"
}
