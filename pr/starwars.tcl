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
#  Demo of a polyelectrolyte star                           #
#                                                           #
#                                                           #
#  Created: 15.09.03  Hans Joerg Limbach                    #
#                                                           #
#############################################################

#############################################################
#  Initialization                                           #
#############################################################

source ../scripts/bundle.tcl

set name ""

#############################################################
#  Setup GUI                                                #
#############################################################

bind . <Destroy> onExitHandler

wm geometry . -0+0

# star
#############################################################
frame .star -relief raised -border 2
pack .star -expand 1 -fill both -in .

# title
frame .star.title
label .star.title.l -justify left -text "Polyelectrolyte Star Simulation"

pack .star.title.l -fill both -side left -in .star.title
pack .star.title -pady {0 10} -fill x -side top -in .star

# menu
frame .star.menu -relief raised -border 1
label  .star.menu.t1 -text "Main Menu"
button .star.menu.b1 -justify left -text "Setup System"     -command SetupStar
button .star.menu.b2 -justify left -text "Read in System"   -command ReadStar
button .star.menu.b3 -justify left -text "Test System"   -command TestStar
button .star.menu.b4 -justify left -text "Polyelectrolyte"   -command PolyelectrolyteTest

pack .star.menu.t1 .star.menu.b1 .star.menu.b2 .star.menu.b3 .star.menu.b4 -in .star.menu
pack .star.menu -pady {0 10} -fill x -side top -in .star

#############################################################
#      Events                                               #
#############################################################

proc TestStar {} {
    global name
    global n_arms
    global l_arms
    global c_dist
    global density
    set name teststar
    set n_arms 5
    set l_arms 40
    set c_dist 2
    set density 0.001

    CreateSystem
}

proc PolyelectrolyteTest {} {
    global name
    global n_arms
    global l_arms
    global c_dist
    global density
    set name polyelectrolyte
    set n_arms 2
    set l_arms 50
    set c_dist 2
    set density 0.001

    CreateSystem
}

proc CreateSystem {} {
    global name n_arms l_arms c_dist density rad

    frame .star.create -relief raised -border 1
    label .star.create.title -text "Create System"
    pack .star.create.title -in .star.create
    pack .star.create -pady {0 10} -fill x -side top -in .star

    # Calculate number of particles
    set n_ci   [expr $n_arms*floor($l_arms/$c_dist)]
    set n_part [expr $n_arms*$l_arms + 1 + $n_ci]
    # Volume
    set vol [expr $n_part/$density]
    set rad [expr pow((3.0*$vol)/(4.0*[PI]),1.0/3.0)]
    set box [expr 4.0*$rad+20.0]
    set cen [expr $box/2.0]

    label .star.create.t1 -justify left -text "Star: $n_arms arms of length $l_arms
Put $n_part Particles, $n_ci Counterions
in Sphere with radius $rad"	
    pack .star.create.t1 -in .star.create

    # Simulation Box
    setmd box_l $box $box $box
    constraint sphere center $cen $cen $cen radius $rad type 4
    setmd periodic  0 0 0

    # Integrator
    set_time_step -2
    set_skin 0.5
    setmd gamma 1.0
    setmd temperature 1.0
    button .star.create.int_but -text "Change Integrator settings" -command change_integrator
    pack .star.create.int_but -in .star.create

    # Particle types and interactions:
    # uncharged monomer 0, charge monomer 2, counterion 1, star center 3, wall 4
    for {set i 0} { $i <= 4} {incr i} {
	for {set j $i} { $j <= 4 } {incr j} {
	    inter $i $j lennard-jones 1.0 1.0 1.12246 1.0 0.0 
	}
    }
    global lj_eps lj_cut
    set lj_eps 1.75
    set lj_cut 2.5
    set_lj_eps $lj_eps
    set_bjerrum 1.0
    inter 0 fene 7.0 2.0
    button .star.create.inter_but -text "Change Interactions" -command change_inter
    pack .star.create.inter_but -in .star.create

    # Create Particles
    set part_id 0
    # Arms
    for {set i 0} { $i < $n_arms} {incr i} {
	set theta [expr     [PI]*[t_random]]
	set phi   [expr 2.0*[PI]*[t_random]]
	set posx  [expr $cen + sin($theta)*cos($phi)]
	set posy  [expr $cen + sin($theta)*sin($phi)]
	set posz  [expr $cen + cos($theta)]
	polymer 1 $l_arms 1.0 start $part_id pos $posx $posy $posz types 0 FENE 0
	for {set j 1} {$j <= $l_arms} {incr j} {
	    if {[expr $j%$c_dist] == 0 } { part $part_id type 2 q 1.0 }
	    incr part_id
	}
    }
    # Center
    part $part_id pos $cen $cen $cen type 3
    for {set i 0} { $i < $n_arms} {incr i} { part $part_id bond 0 [expr $i*$l_arms] }
    incr part_id

    # Counterions
    set center [list $cen $cen $cen]
    bundle_counterion_setup $n_ci $rad -1.0 $center 4 $part_id
    set part_id [expr $part_id+$n_ci]
    # Particle Control
    #        for {set i 0} { $i < $part_id} {incr i} {
    #   	puts [part $i]
    #      }
    polyBlockWrite "$name.start" 

    # Start Integration
    button .star.create.start_but -text "Start Simulation" -command simulation
    pack .star.create.start_but -in .star.create
}

proc get_Setup {win} {
    global name
    global n_arms
    global l_arms
    global c_dist
    global density
    set name [$win.inp1.i get]
    if {$name == ""} { 
	popup_message "System needs a name!"
	return
    }
    set n_arms [$win.inp2.i get]
    if { $n_arms < 1 || $n_arms > 12} { 
	popup_message "Number of Arms has to be between 1 and 12" 
	return
    }
    set l_arms [$win.inp3.i get]
    if {$l_arms < 1} {
	popup_message "Arm length has to be positive"
	return
    }
    set c_dist [$win.inp4.i get]
    if {$c_dist < 1} {
	popup_message "Charge distance has to be positive"
	return
    }
    set density [$win.inp5.i get]
    if {$density < 1e-8 || $density > 0.7} {
	popup_message "Density is restricted between 1e-8 and 0.7"
	return
    }

    CreateSystem

    destroy $win 
}

proc SetupStar  {} {
    frame .star.setup -relief raised -border 1
    label .star.setup.title -text "Setup System"

    frame .star.setup.inp1
    label .star.setup.inp1.t -justify left -text "Name:            "
    entry .star.setup.inp1.i
    pack  .star.setup.inp1.t -fill both -side left -in .star.setup.inp1
    pack  .star.setup.inp1.i -in .star.setup.inp1

    frame .star.setup.inp2
    label .star.setup.inp2.t -justify left -text "Arm Number:      "
    entry .star.setup.inp2.i
    pack  .star.setup.inp2.t -fill both -side left -in .star.setup.inp2
    pack  .star.setup.inp2.i -in .star.setup.inp2

    frame .star.setup.inp3
    label .star.setup.inp3.t -justify left -text "Arm Length:      "
    entry .star.setup.inp3.i
    pack  .star.setup.inp3.t -fill both -side left -in .star.setup.inp3
    pack  .star.setup.inp3.i -in .star.setup.inp3

    frame .star.setup.inp4
    label .star.setup.inp4.t -justify left -text "Charge Distance: "
    entry .star.setup.inp4.i
    pack  .star.setup.inp4.t -fill both -side left -in .star.setup.inp4
    pack  .star.setup.inp4.i -in .star.setup.inp4

    frame .star.setup.inp5
    label .star.setup.inp5.t -justify left -text "Density:         "
    entry .star.setup.inp5.i
    pack  .star.setup.inp5.t -fill both -side left -in .star.setup.inp5
    pack  .star.setup.inp5.i -in .star.setup.inp5

    button .star.setup.b1 -text "OK" -command "get_Setup .star.setup"

    pack .star.setup.title .star.setup.inp1 .star.setup.inp2 .star.setup.inp3 .star.setup.inp4 .star.setup.inp5 .star.setup.b1 -in .star.setup
    pack .star.setup -pady {0 10} -fill x -side top -in .star
    tkwait window .star.setup
    update
}

proc change_integrator {} {
    global time_step skin gamma temp
    toplevel .integrator
    label .integrator.title -text "Change Integrator parameters"
    pack .integrator.title -in .integrator

    frame .integrator.slider1 -relief raised -border 1
    label .integrator.slider1.t -text "log(time_step)"
    scale .integrator.slider1.s -orient h -from -3.5 -to -1.5 -resolution [expr 1/100.] -command set_time_step
    .integrator.slider1.s set $time_step
    pack .integrator.slider1.t .integrator.slider1.s -fill both -in .integrator.slider1
    pack .integrator.slider1 -in .integrator

    frame .integrator.slider2 -relief raised -border 1
    label .integrator.slider2.t -text "Skin"
    scale .integrator.slider2.s -orient h -from 0 -to 4 -resolution [expr 1/10.] -command set_skin
    .integrator.slider2.s set $skin
    pack .integrator.slider2.t .integrator.slider2.s -fill both -in .integrator.slider2
    pack .integrator.slider2 -in .integrator
 
    button .integrator.but1 -text "OK" -command "destroy .integrator"
    pack .integrator.but1 -in .integrator

    tkwait window .integrator
    update

}

proc set_time_step {ts} {
    global time_step
    setmd time_step [expr pow(10, $ts)]
    set time_step $ts
}

proc set_skin {sk} {
    global skin
    setmd skin $sk
    set skin $sk
}

proc change_inter {} {
    global lj_eps bjerrum lj_cut
    toplevel .inter
    label .inter.title -text "Change Interaction parameters"
    pack .inter.title -in .inter

    frame .inter.slider1 -relief raised -border 1
    label .inter.slider1.t -text "Solvent Quality (lj_eps)"
    scale .inter.slider1.s -orient h -from 0 -to 2 -resolution [expr 1/100.] -command set_lj_eps
    .inter.slider1.s set $lj_eps
    pack .inter.slider1.t .inter.slider1.s -fill both -in .inter.slider1
    pack .inter.slider1 -in .inter

    frame .inter.slider2 -relief raised -border 1
    label .inter.slider2.t -text "Poor Solvent range (lj_cut)"
    scale .inter.slider2.s -orient h -from 1.12246 -to 3.0 -resolution [expr 1/100.] -command set_lj_cut
    .inter.slider2.s set $lj_cut
    pack .inter.slider2.t .inter.slider2.s -fill both -in .inter.slider2
    pack .inter.slider2 -in .inter

    frame .inter.slider3 -relief raised -border 1
    label .inter.slider3.t -text "Electrostatics (Bjerrum length)"
    scale .inter.slider3.s -orient h -from 0 -to 10 -resolution [expr 1/10.] -command set_bjerrum
    .inter.slider3.s set $bjerrum
    pack .inter.slider3.t .inter.slider3.s -fill both -in .inter.slider3
    pack .inter.slider3 -in .inter
 
    button .inter.but1 -text "OK" -command "destroy .inter"
    pack .inter.but1 -in .inter

    tkwait window .inter
    update

}

proc set_lj_eps {val} {
    global lj_eps lj_cut
    set lj_eps $val
    set lj_shift [expr -4.0*$lj_eps*(pow($lj_cut,-12)-pow($lj_cut,-6))]
    inter 0 0 lennard-jones $lj_eps 1.0 $lj_cut $lj_shift 0
    inter 0 2 lennard-jones $lj_eps 1.0 $lj_cut $lj_shift 0
    inter 2 2 lennard-jones $lj_eps 1.0 $lj_cut $lj_shift 0
}

proc set_lj_cut {val} {
    global lj_eps lj_cut
    set lj_cut $val
    set lj_shift [expr -4.0*$lj_eps*(pow($lj_cut,-12)-pow($lj_cut,-6))]
    inter 0 0 lennard-jones $lj_eps 1.0 $lj_cut $lj_shift 0
    inter 0 2 lennard-jones $lj_eps 1.0 $lj_cut $lj_shift 0
    inter 2 2 lennard-jones $lj_eps 1.0 $lj_cut $lj_shift 0
}


proc set_bjerrum {val} {
    global bjerrum rad
    set bjerrum $val
    inter coulomb $bjerrum dh 0.0 [expr 2.0*$rad]

}

proc simulation {} {
    global name bjerrum lj_cut lj_eps n_arms l_arms

    frame .star.sim -relief raised -border 1
    label .star.sim.title -text "Simulation"
    pack .star.sim.title -in .star.sim
    pack .star.sim -in .star

    prepare_vmd_connection "$name" 0 1

    # Check if Warmup is necessary
    set warmup_dist 0.9
    set cap 20  
    set act_min_dist [analyze mindist]
    if { $act_min_dist < $warmup_dist } {
	frame .star.sim.warmup 
	label .star.sim.warmup.title -text "Warmup:"
	label .star.sim.warmup.md   -text "mindist: *"
	label .star.sim.warmup.time -text "time:    *"
	pack .star.sim.warmup.title .star.sim.warmup.md .star.sim.warmup.time -in .star.sim.warmup
	pack .star.sim.warmup -in .star.sim
	pack .star.sim -in .star
	
	set tmp_bjerrum $bjerrum
	set tmp_lj_cut $lj_cut
	set tmp_lj_eps $lj_eps
	set_bjerrum 0
	set_lj_eps 1.0
	set_lj_cut 1.12246
	inter ljforcecap $cap
    }
    while { $act_min_dist < $warmup_dist } {
	integrate 20
	imd positions
	set act_min_dist [analyze mindist]
	set time [setmd time]
	.star.sim.warmup.md   conf -text "mindist: $act_min_dist"
	.star.sim.warmup.time conf -text "time:  : $time"
	set cap [expr $cap+5]
	inter ljforcecap $cap  
    }
    if { $cap > 20 } {
	set_bjerrum $tmp_bjerrum
	set_lj_eps $tmp_lj_eps
	set_lj_cut $tmp_lj_cut
	inter ljforcecap 0
    }

    frame .star.sim.integ
    label .star.sim.integ.title -text "Integ:"
    label .star.sim.integ.time -text "time = *"
    label .star.sim.integ.re   -text "re   = *"
    pack .star.sim.integ.title .star.sim.integ.time .star.sim.integ.re -in .star.sim.integ
    pack .star.sim.integ -in .star.sim
    
    analyze set chains 0 $n_arms $l_arms
    while { 1 > 0 } {
	integrate 20
	imd positions
	set time [setmd time]
	set re [lindex [analyze re] 0]
	.star.sim.integ.time   conf -text "time = $time"
	.star.sim.integ.re     conf -text "re   = $re*"
	update
    }

}


#############################################################
#      Helpers                                              #
#############################################################

proc onExitHandler {} {
    
}

proc popup_message {text} {
    toplevel .message
    button .message.btn -text "$text" -justify left -command "destroy .message"
    pack .message.btn
    wm geometry .message +60+0
    tkwait window .message
    update
}