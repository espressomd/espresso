#############################################################
#                                                           #
#  Demo of a polyelectrolyte star                           #
#                                                           #
#############################################################
#
# Copyright (C) 2010,2011,2012,2013 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
#  Max-Planck-Institute for Polymer Research, Theory Group
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

#############################################################
#  Initialization                                           #
#############################################################
require_feature "TK"
require_feature "LENNARD_JONES"
require_feature "ELECTROSTATICS"
require_feature "CONSTRAINTS"
require_feature "PARTIAL_PERIODIC"

source ../../scripts/bundle.tcl

global name level
set name ""
set level 0

#############################################################
#  Setup GUI                                                #
#############################################################

bind . <Destroy> onExitHandler

wm geometry . +0+0

# star
#############################################################
frame .star -relief raised -border 2 
pack .star -expand 1 -fill both -in .

# title
frame .star.title
label .star.title.l -text "Polyelectrolyte Star Simulation" -bg lightblue -height 2

canvas .star.title.logo -width 234 -height 72 -bg black
set esp_logo [image create photo -file "logo.gif"]
.star.title.logo create image 117 36 -image $esp_logo

pack .star.title.l .star.title.logo -fill both  -in .star.title
pack .star.title -in .star

# menu
frame .star.menu -relief raised -border 1
label  .star.menu.t1 -text "Main Menu" -bg lightyellow
button .star.menu.b1 -text "Setup System"     -command SetupStar
button .star.menu.b2 -text "Read in System"   -command ReadStar
button .star.menu.b3 -text "Test System"   -command TestStar
button .star.menu.b4 -text "Exit"   -command exit

pack .star.menu.t1 .star.menu.b1 .star.menu.b2 .star.menu.b3 .star.menu.b4 \
    -expand 1 -fill both -in .star.menu
pack .star.menu -expand 1 -fill both -in .star

#############################################################
#     Main Menu                                             #
#############################################################

proc SetupStar  {} {
    global name level
    if { $name != "" } { popup_message "System exists already. Exit and start new!"; return }
    if { $level > 0 } { popup_message "Close other submenu first!"; return }
    set level 1
    frame .star.setup -relief raised -border 1 -width 20c
    label .star.setup.title -text "Setup System" -bg lightyellow 

    frame .star.setup.inp1
    label .star.setup.inp1.t -anchor n -text "Name: "
    entry .star.setup.inp1.i 
    pack  .star.setup.inp1.t -fill both -side left  -in .star.setup.inp1
    pack  .star.setup.inp1.i  -in .star.setup.inp1

    frame .star.setup.inp2
    label .star.setup.inp2.t -anchor n -text "Arm Number: "
    entry .star.setup.inp2.i
    pack  .star.setup.inp2.t -fill both -side left -in .star.setup.inp2
    pack  .star.setup.inp2.i -expand 1 -in .star.setup.inp2

    frame .star.setup.inp3
    label .star.setup.inp3.t -anchor n  -text "Arm Length:      "
    entry .star.setup.inp3.i
    pack  .star.setup.inp3.t -fill both -side left -in .star.setup.inp3
    pack  .star.setup.inp3.i -in .star.setup.inp3

    frame .star.setup.inp4
    label .star.setup.inp4.t -anchor n  -text "Charge Distance: "
    entry .star.setup.inp4.i
    pack  .star.setup.inp4.t -fill both -side left -in .star.setup.inp4
    pack  .star.setup.inp4.i -in .star.setup.inp4

    frame .star.setup.inp5
    label .star.setup.inp5.t -anchor n  -text "Density:         "
    entry .star.setup.inp5.i
    pack  .star.setup.inp5.t -fill both -side left -in .star.setup.inp5
    pack  .star.setup.inp5.i -in .star.setup.inp5

    button .star.setup.b1 -text "OK" -command "get_Setup .star.setup"
    button .star.setup.b2 -text "Close" -command "destroy .star.setup; set level 0"

    pack  .star.setup.title .star.setup.inp1 .star.setup.inp2             \
	.star.setup.inp3 .star.setup.inp4 .star.setup.inp5 .star.setup.b1 .star.setup.b2 \
	-expand 1 -fill both -in .star.setup
    pack .star.setup -expand 1 -in .star
    tkwait window .star.setup
    update
}


proc ReadStar {} {
    global name 
    if { $name != "" } { popup_message "System exists already. Exit and start new"; return }
    puts "Not ready"
}

proc TestStar {} {
    global name n_arms l_arms c_dist density level
    if { $name != "" } { popup_message "System exists already. Exit and start new"; return }
    if { $level > 0 } { popup_message "Close other submenu first!"; return }
    set level 1
    set name    test1
    set n_arms  5
    set l_arms  20
    set c_dist  3
    set density 0.001

    CreateSystem
}


#############################################################
#     Create System                                         #
#############################################################

proc get_Setup {win} {
    global name n_arms l_arms c_dist density
    set name [$win.inp1.i get]
    if {$name == ""} { popup_message "System needs a name!"; return }
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

#############################################################
#     Create System                                         #
#############################################################

proc CreateSystem {} {
    global name n_arms l_arms c_dist density rad

    # Calculate number of particles
    set n_ci   [expr $n_arms*floor($l_arms/$c_dist)]
    set n_part [expr $n_arms*$l_arms + 1 + $n_ci]

    # Volume
    set vol [expr $n_part/$density]
    set rad [expr pow((3.0*$vol)/(4.0*[PI]),1.0/3.0)]
    set box [expr 4.0*$rad+20.0]
    set cen [expr $box/2.0]

    # Simulation Box
    setmd box_l $box $box $box
    constraint sphere center $cen $cen $cen radius $rad type 4
    setmd periodic  0 0 0

    # Integrator
    set_time_step -2
    set_skin 0.5
    thermostat langevin 1.0 1.0

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

    # Window appearance
    frame .star.create -relief raised -border 1
    label .star.create.title -text "Create System" -bg lightyellow -height 2
    label .star.create.t1 -justify left -text "System: $name 
Star with $n_arms arms of length $l_arms
Put $n_part Particles, $n_ci Counterions
in Sphere with radius $rad"	
    button .star.create.int_but -text "Change Integrator settings" -command change_integrator
    button .star.create.inter_but -text "Change Interactions" -command change_inter
    button .star.create.start_but -text "Start Simulation" -command simulation

    pack .star.create.title .star.create.t1                                   \
	.star.create.int_but .star.create.inter_but .star.create.start_but   \
	-expand 1 -fill both -in .star.create
    pack .star.create  -expand 1  -in .star

}


#############################################################
#     Change Parameters                                     #
#############################################################

proc change_integrator {} {
    global time_step skin gamma temp
    toplevel .integrator
    label .integrator.title -text "Change Integrator parameters" \
	-bg lightyellow -height 2 -relief raised -border 1

    frame .integrator.slider1 -relief raised -border 1
    label .integrator.slider1.t -text "log(time_step)"
    scale .integrator.slider1.s -orient h -from -3.5 -to -0.5 -resolution [expr 1/100.] -command set_time_step
    .integrator.slider1.s set $time_step
    pack .integrator.slider1.t .integrator.slider1.s -side left  -in .integrator.slider1

    frame .integrator.slider2 -relief raised -border 1
    label .integrator.slider2.t -text "Skin"
    scale .integrator.slider2.s -orient h -from 0 -to 4 -resolution [expr 1/10.] -command set_skin
    .integrator.slider2.s set $skin
    pack .integrator.slider2.t .integrator.slider2.s -expand 1 -side left  -in .integrator.slider2
 
    button .integrator.but1 -text "OK" -command "destroy .integrator"
    pack .integrator.title .integrator.slider1 .integrator.slider2 .integrator.but1 \
	-expand 1 -fill both -in .integrator

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
    label .inter.title -text "Change Interaction parameters"  \
	-bg lightyellow -height 2 -relief raised -border 1

    frame .inter.slider1 -relief raised -border 1
    label .inter.slider1.t -text "Solvent Quality (lj_eps)"
    scale .inter.slider1.s -orient h -from 0 -to 2 -resolution [expr 1/100.] -command set_lj_eps
    .inter.slider1.s set $lj_eps
    pack .inter.slider1.t .inter.slider1.s -side left -in .inter.slider1

    frame .inter.slider2 -relief raised -border 1
    label .inter.slider2.t -text "Poor Solvent range (lj_cut)"
    scale .inter.slider2.s -orient h -from 1.12246 -to 3.0 -resolution [expr 1/100.] -command set_lj_cut
    .inter.slider2.s set $lj_cut
    pack .inter.slider2.t .inter.slider2.s -side left -in .inter.slider2

    frame .inter.slider3 -relief raised -border 1
    label .inter.slider3.t -text "Electrostatics (Bjerrum length)"
    scale .inter.slider3.s -orient h -from 0 -to 10 -resolution [expr 1/10.] -command set_bjerrum
    .inter.slider3.s set $bjerrum
    pack .inter.slider3.t .inter.slider3.s -side left -in .inter.slider3

    button .inter.but1 -text "OK" -command "destroy .inter"

    pack .inter.title .inter.slider1 .inter.slider2 .inter.slider3 .inter.but1  \
	-expand 1 -fill both -in .inter

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

#############################################################
#     Simulation                                            #
#############################################################


proc simulation {} {
    global name bjerrum lj_cut lj_eps n_arms l_arms vmd_on
    global cobs_name
    set cobs_name "jhdjuyuw"
    set vmd_on 0

    frame .star.sim -relief raised -border 1
    label .star.sim.title -text "Simulation" -bg lightyellow -height 2 

    frame .star.sim.warmup -relief raised -border 1
    label .star.sim.warmup.title -text ""  -height 2 
    label .star.sim.warmup.md   -text ""
    label .star.sim.warmup.time -text ""
    pack .star.sim.warmup.title .star.sim.warmup.md .star.sim.warmup.time -in .star.sim.warmup

    frame .star.sim.integ -relief raised -border 1
    label .star.sim.integ.title -text "Integ:" -height 2 
    label .star.sim.integ.time -text "time = *"
    label .star.sim.integ.re   -text "re   = *"
    pack .star.sim.integ.title .star.sim.integ.time .star.sim.integ.re  -in .star.sim.integ

    
    frame .star.sim.graphic -relief raised -border 1
    label .star.sim.graphic.title -text "Visualisation"  -bg lightyellow -height 2
    button .star.sim.graphic.vmd -text "VMD" -command "start_vmd"
    button .star.sim.graphic.re -text "RE Graph" -command "create_obs_window re t RE 0 200 400 500 300"
    button .star.sim.graphic.et -text "Energy Graph" -command "create_obs_window totenergy t E_t 0 200 200 500 300"
    button .star.sim.graphic.cus -text "Custom Graph" -command "customize_obs_window"

    pack .star.sim.graphic.title .star.sim.graphic.vmd .star.sim.graphic.re .star.sim.graphic.et .star.sim.graphic.cus -expand 1 -fill both -in .star.sim.graphic

    pack .star.sim.title .star.sim.warmup .star.sim.integ .star.sim.graphic \
	-expand 1 -fill both -in .star.sim
    pack .star.sim -expand 1 -fill both -in .star


    # Warmup
    ########################################
    set warmup_dist 0.9
    set cap 20  
    set act_min_dist [analyze mindist]
    if { $act_min_dist < $warmup_dist } {
	.star.sim.warmup.title conf -text "Warmup:" 
    
	set tmp_bjerrum $bjerrum
	set tmp_lj_cut $lj_cut
	set tmp_lj_eps $lj_eps
	set_bjerrum 0
	set_lj_eps 1.0
	set_lj_cut 1.12246
	inter forcecap $cap
    }
    while { $act_min_dist < $warmup_dist } {
	integrate 20
	if { $vmd_on == 1 } { imd positions }
	set act_min_dist [analyze mindist]
	set time [setmd time]
	.star.sim.warmup.md   conf -text "mindist: $act_min_dist"
	.star.sim.warmup.time conf -text "time:  : $time"
	update
	set cap [expr $cap+5]
	inter forcecap $cap  
    }
    if { $cap > 20 } {
	set_bjerrum $tmp_bjerrum
	set_lj_eps $tmp_lj_eps
	set_lj_cut $tmp_lj_cut
	inter forcecap 0
    }

 
    # Integration
    ########################################
    analyze set chains 0 $n_arms $l_arms
    set step 0
    while { 1 > 0 } {
	.star.sim.integ.title conf -text "Integration:" 
	integrate 20

	if { $vmd_on == 1 } { imd positions }

	set time [setmd time]
	set re [lindex [analyze re] 0]
	.star.sim.integ.time   conf -text "time = $time"
	.star.sim.integ.re     conf -text "re   = $re*"
	if { [winfo exists .re] && [expr $step%5]==0 } {
	    add_data_point re $time $re
	}
	if { [winfo exists .totenergy] && [expr $step%10]==0 } {
	    set energy [analyze energy]
	    set tot_e [lindex [lindex $energy 0] 1]
	    add_data_point totenergy $time $tot_e
	}
	if { [winfo exists .$cobs_name] } {
	    get_data_point
	}
	update
	incr step
    }

}

#############################################################
#      VMD Connection                                       #
#############################################################


proc start_vmd {} {
    global vmd_on name
    if { $vmd_on == 1 } { return }
    set vmd_on 1
    prepare_vmd_connection "$name" 0 1
}

#############################################################
#      Observable Graphics                                  #
#############################################################

proc create_obs_window {obs xlabel ylabel xrange yrange max_points width height } {
    if { [winfo exists .$obs] } { return }

    # Global variables associated with the window:
    # plot region
    global xpix_$obs ypix_$obs 
    set xpix_$obs $width
    set ypix_$obs $height
    # Margins
    global lmarg_$obs rmarg_$obs tmarg_$obs bmarg_$obs
    set lmarg_$obs 60
    set rmarg_$obs 25
    set tmarg_$obs 30
    set bmarg_$obs 40
    # data range
    global xrng_$obs yrng_$obs pmax_$obs
    global ymin_$obs ymax_$obs
    set xrng_$obs $xrange
    set yrng_$obs $yrange
    set pmax_$obs $max_points
    # data
    global data_$obs show_$obs

    # Window elements
    toplevel .$obs
    canvas .$obs.c -width [expr \$lmarg_$obs + \$xpix_$obs + \$rmarg_$obs]   \
	-height [expr \$tmarg_$obs + \$ypix_$obs + \$bmarg_$obs] -bg white
    frame .$obs.b
    button .$obs.b1 -text "PS-File" -command "obs_window_to_ps $obs"
    button .$obs.b2 -text "Options" -command "puts \"Sorry not ready yet!\""
#    button .$obs.b2 -text "Options" -command "obs_window_options $obs"
    button .$obs.b3 -text "Close"   -command \
	"destroy .$obs; unset data_$obs show_$obs ymin_$obs ymax_$obs"
    pack .$obs.b1 .$obs.b2 .$obs.b3 -side left -fill both -in .$obs.b
    pack .$obs.c .$obs.b -in .$obs

    # Draw axes
    .$obs.c create line                                                          \
	[expr \$lmarg_$obs]  [expr \$tmarg_$obs - 10]                            \
	[expr \$lmarg_$obs]  [expr \$tmarg_$obs + \$ypix_$obs]                   \
	[expr \$lmarg_$obs + \$xpix_$obs +10] [expr \$tmarg_$obs + \$ypix_$obs]  \
	-arrow both -tag "axes"
    .$obs.c create text [expr \$lmarg_$obs] [expr \$tmarg_$obs-20] \
	-text "$ylabel" 
    .$obs.c create text [expr \$lmarg_$obs + \$xpix_$obs +15]      \
	[expr \$tmarg_$obs + \$ypix_$obs]                          \
	-text "$xlabel"
}

proc obs_window_to_ps {obs} {
    if {[ winfo exists .ps_$obs]} { return }
    toplevel .ps_$obs
    frame .ps_$obs.f1
    label .ps_$obs.f1.t -text "File name: "
    entry .ps_$obs.f1.i
    .ps_$obs.f1.i insert 0 "$obs.ps"
    button .ps_$obs.f1.b -text "OK" -command \
	"update; .$obs.c postscript -file [.ps_$obs.f1.i get] ; destroy .ps_$obs"
    pack .ps_$obs.f1.t .ps_$obs.f1.i .ps_$obs.f1.b -side left -in .ps_$obs.f1
    pack .ps_$obs.f1 -in .ps_$obs
}
proc obs_window_options {obs} {
    if {[ winfo exists .opt_$obs]} { return }
    global xrng_$obs yrng_$obs pmax_$obs xpix_$obs ypix_$obs 

    toplevel .opt_$obs

    frame .opt_$obs.s1
    label .opt_$obs.s1.t -text "Width (pixel):"
    scale .opt_$obs.s1.s -orient h                              \
	-from [expr \$ypix_$obs/2.0] -to [expr \$ypix_$obs*2.0] \
	-resolution 1 -command "change_var $obs xpix_$obs"
    .opt_$obs.s1.s set [expr \$ypix_$obs]
    pack .opt_$obs.s1.t .opt_$obs.s1.s -side left -in .opt_$obs.s1




    button .opt_$obs.b -text "OK" -command "destroy .opt_$obs"
    pack .opt_$obs.s1 ..opt_$obs.b  -in .opt_$obs
}

proc change_var {obs var val} {
    set var $val
} 


proc customize_obs_window {} {
    global ycmd yres
    set cobs_name ""
    set ycmd "set dummy 0 "
    

    toplevel .cobs

    label .cobs.title -text "Customize new observable window" -height 2 -bg lightyellow

    frame .cobs.name
    label .cobs.name.t -text "Name:         "
    entry .cobs.name.i  -width 30
    .cobs.name.i insert 0 "custom" 
    pack .cobs.name.t .cobs.name.i -side left -fill both -in .cobs.name

    frame .cobs.ycmd
    label .cobs.ycmd.t -text "Y-Command:"
    entry .cobs.ycmd.i -width 30
    .cobs.ycmd.i insert 0 "analyze" 
    pack .cobs.ycmd.t .cobs.ycmd.i -side left -fill both -in .cobs.ycmd
    
    frame .cobs.yfct
    label .cobs.yfct.t -text "Y-Function:  "
    entry .cobs.yfct.i -width 30
    .cobs.yfct.i insert 0 "\$y" 
    pack .cobs.yfct.t .cobs.yfct.i  -side left -fill both -in .cobs.yfct

    frame .cobs.size
    label .cobs.size.t -text "Size (w h):  "
    entry .cobs.size.i -width 30
    .cobs.size.i insert 0 "300 200" 
    pack .cobs.size.t .cobs.size.i  -side left -fill both -in .cobs.size

    frame .cobs.num
    label .cobs.num.t -text "Data Points:"
    entry .cobs.num.i -width 30
    .cobs.num.i insert 0 "100" 
    pack .cobs.num.t .cobs.num.i  -side left -fill both -in .cobs.num

    button .cobs.b1 -text "OK" -command "get_ycmd"

    frame .cobs.yres 
    label .cobs.yres.t0  -text "Instructions:"
    label .cobs.yres.t1  -text ""
    label .cobs.yres.t2  -text ""
    label .cobs.yres.t3  -text ""
    pack .cobs.yres.t0 .cobs.yres.t1 .cobs.yres.t2 .cobs.yres.t3  -in .cobs.yres

    pack .cobs.title  .cobs.name .cobs.ycmd .cobs.yfct .cobs.size .cobs.num .cobs.b1 .cobs.yres -expand 1 -fill both -in  .cobs
}

proc get_ycmd {} {
    global ycmd yfct cobs_name 
    update
    set cobs_name [ .cobs.name.i get ]
    set ycmd [ .cobs.ycmd.i get ]
    set yfct [ .cobs.yfct.i get ]

    if { $cobs_name == "" || [winfo exists .$cobs_name] } {
	label .cobs.yres.t1  -text ""
	label .cobs.yres.t2  -text "Need another name!"
	label .cobs.yres.t3  -text ""
	return
    }
    set size [ .cobs.size.i get ]
    if { [llength $size] < 1 } {
	set size_x [lindex $size 0]
	set size_y [lindex $size 1]
	if { ![string is integer -strict $size_x] || ![string is integer -strict $size_y] ||
	     $size_x < 20 || $size_y < 20 } {
	    label .cobs.yres.t1  -text ""
	    label .cobs.yres.t2  -text "Wrong size! (must be larger than 20 20)"
	    label .cobs.yres.t3  -text ""
	    return
	}
    }
    set size_x [lindex $size 0]
    set size_y [lindex $size 1]
    set data_num [ .cobs.num.i get ]
    if { ![string is integer -strict $data_num] || $data_num < 1} {
	    label .cobs.yres.t1  -text ""
	    label .cobs.yres.t2  -text "Need at least 1 data point!"
	    label .cobs.yres.t3  -text ""
	    return
    }
    if { $ycmd == "analyze" } {
	    label .cobs.yres.t1  -text ""
	    label .cobs.yres.t2  -text "You have to modify the Y-command!"
	    label .cobs.yres.t3  -text ""
	    return
    }
    if { $ycmd != "analyze" } {
	set yres [eval $ycmd] 
	puts $yres
	if { [llength $yres] > 1 } { 
	    .cobs.yres.t1 conf -text "Observable structure:\nElement   Structure"
	    set tmp "0\t[lindex $yres 0]"
	    for { set i 1 } { $i < [llength $yres] } {incr i} {
		set tmp "$tmp\n$i\t[lindex $yres $i]"
	    }
	    .cobs.yres.t2 conf -justify left -text "$tmp"
	    .cobs.yres.t3 conf -text "Choose an element!"
	    .cobs.ycmd.i insert 0 "lindex \["
	    .cobs.ycmd.i insert end "\] 0"
	    return
	}
	if { [string is double -strict [eval $ycmd]] } { 
	    set y [eval $ycmd]
	    puts $y
	    if { [string is double -strict [expr $yfct]] } { 
		.cobs.yres.t0 conf -text ""
		.cobs.yres.t1 conf -text ""
		.cobs.yres.t2 conf -text "Start Graphic"
		.cobs.yres.t3 conf -text ""
		.cobs.b1 conf -text "Close Window" -command "destroy .cobs"
		create_obs_window $cobs_name t $cobs_name 0 100 $data_num $size_x $size_y
	    } else {
		.cobs.yres.t1 conf -text "Y-Function does not result a double"
		.cobs.yres.t2 conf -text "Change Y-Function entry!"
		.cobs.yres.t3 conf -text ""
	    } 
	} else {
	    puts [eval $ycmd]
	    .cobs.yres.t1 conf -text "Y-Command does not result a double"
	    .cobs.yres.t2 conf -text "Change Y-Command entry!"
	    .cobs.yres.t3 conf -text ""
	}

    }
}

proc get_data_point {} {
    global ycmd yfct cobs_name
    set time [setmd time]
    set y [eval $ycmd]
    set yval [expr $yfct]
    add_data_point $cobs_name $time $yval
}

proc add_data_point {obs xdat ydat} {
    global data_$obs show_$obs

    lappend data_$obs $xdat $ydat
    lappend show_$obs $xdat $ydat

    if { [llength [expr \$data_$obs]] > 2 } {
	change_xrange $obs
	change_yrange $obs
        if { [llength [expr \$data_$obs]] > 4 } { .$obs.c delete "data_line" }
	.$obs.c create line [expr \$show_$obs] -width 3 -fill red -tag "data_line"
	update
    }
}

# Time coordinate:
# Sucessive coordinates are assumed to be increasing
proc change_xrange {obs} {
    global data_$obs show_$obs
    global xrng_$obs pmax_$obs xpix_$obs
    global lmarg_$obs

    # case 1: xrange given
    if { [expr \$xrng_$obs] > 0 } {
	set xrange [expr \$xrng_$obs + 1]
	while { $xrange > [expr \$xrng_$obs] } {
	    set dat_len [llength [expr \$data_$obs]]
	    set xmin [lindex [expr \$data_$obs] 0]
	    set xmax [lindex [expr \$data_$obs] [expr $dat_len - 2]]
	    set xrange [expr $xmax - $xmin]
	    if { $xrange > [expr \$xrng_$obs] } {
		set data_$obs [lrange [expr \$data_$obs] 2 end]
		set show_$obs [lrange [expr \$show_$obs] 2 end]
	    }
	}
	set xfac [expr \$xpix_$obs/(1.0*$xrange)]
	set xoff [expr \$lmarg_$obs - ($xmin*$xfac)]
	for { set i 0 } { $i < $dat_len } { incr i 2 } {
	    lset show_$obs $i [expr [lindex [expr \$data_$obs] $i]*$xfac+$xoff]
	}
	draw_xtics $obs $xmin $xmax $xrange $xfac $xoff
	return
    }
    # case 2: number of data points to shown given
    while { [llength [expr \$data_$obs]] > [expr 2*\$pmax_$obs] } {
		set data_$obs [lrange [expr \$data_$obs] 2 end]
		set show_$obs [lrange [expr \$show_$obs] 2 end]
    }
    set dat_len [llength [expr \$data_$obs]]
    set xmin [lindex [expr \$data_$obs] 0]
    set xmax [lindex [expr \$data_$obs] [expr $dat_len - 2]]
    set xrange [expr $xmax - $xmin]
    set xfac [expr \$xpix_$obs/(1.0*$xrange)]
    set xoff [expr \$lmarg_$obs - ($xmin*$xfac)]
    for { set i 0 } { $i < $dat_len } { incr i 2 } {
	lset show_$obs $i [expr [lindex [expr \$data_$obs] $i]*$xfac+$xoff]
    }
    draw_xtics $obs $xmin $xmax $xrange $xfac $xoff
}

# This is a scattering value
proc change_yrange {obs} {
    global data_$obs show_$obs
    global yrng_$obs pmax_$obs ypix_$obs
    global ymin_$obs ymax_$obs
    global tmarg_$obs

    set dat_len [llength [expr \$data_$obs]]
    # init
    if { $dat_len == 4 } {
        set ymin_$obs [lindex [expr \$data_$obs] 1]
	set ymax_$obs [lindex [expr \$data_$obs] 1]
    }
    # check yrange change
    incr dat_len -1
    if { [lindex [expr \$data_$obs] $dat_len] < [expr \$ymin_$obs] } {
	set ymin_$obs [lindex [expr \$data_$obs] $dat_len]
    }
    if { [lindex [expr \$data_$obs] $dat_len] > [expr \$ymax_$obs] } {
	set ymax_$obs [lindex [expr \$data_$obs] $dat_len]
    }
    set yrange [expr \$ymax_$obs - \$ymin_$obs]
    if { $yrange == 0 } { set yrange 1 }
    set yfac [expr \$ypix_$obs/(1.0*$yrange)]
    set yoff [expr \$tmarg_$obs + (\$ymax_$obs*$yfac)]
    for { set i 1 } { $i <= $dat_len } { incr i 2 } {
	lset show_$obs $i [expr $yoff - [lindex [expr \$data_$obs] $i]*$yfac]
    }
    draw_ytics $obs [expr \$ymin_$obs] [expr \$ymax_$obs] $yrange $yfac $yoff

}

proc draw_xtics { obs xmin xmax xrange xfac xoff } {
    global ypix_$obs
    global lmarg_$obs tmarg_$obs
    .$obs.c delete "xtics"
    if { $xrange > 0.0 } {
	set xexp [expr floor(log10($xrange))]
	set xpow [expr pow(10,$xexp)]
	set bas_range [expr $xrange/(1.0*$xpow)]
	set dist [expr 2.0*$xpow]
	if {  $bas_range < 6.0 } { set dist [expr 1.0*$xpow] }
	if {  $bas_range < 3.0 } { set dist [expr 0.5*$xpow] }
	if {  $bas_range < 1.7 } { set dist [expr 0.25*$xpow] }
	
	set xtics [list [expr $dist*ceil($xmin/(1.0*$dist))] ]
	.$obs.c create text                                                             \
	    [expr [lindex $xtics 0]*$xfac+$xoff] [expr 15 + \$tmarg_$obs + \$ypix_$obs] \
	    -text "[lindex $xtics 0]" -tag "xtics"
	.$obs.c create line                                                            \
	    [expr [lindex $xtics 0]*$xfac+$xoff] [expr \$tmarg_$obs + \$ypix_$obs]     \
	    [expr [lindex $xtics 0]*$xfac+$xoff] [expr 5 + \$tmarg_$obs + \$ypix_$obs] \
	    -width 2 -tag "xtics"
	set i 0
	while { [expr [lindex $xtics $i]+$dist] <= $xmax } {
	    lappend xtics [expr [lindex $xtics $i]+$dist]
	    incr i
	    .$obs.c create text                                                          \
		[expr [lindex $xtics $i]*$xfac+$xoff] [expr 15+\$tmarg_$obs+\$ypix_$obs] \
		-text "[lindex $xtics $i]" -tag "xtics"
	    .$obs.c create line                                                         \
		[expr [lindex $xtics $i]*$xfac+$xoff] [expr \$tmarg_$obs+\$ypix_$obs]   \
		[expr [lindex $xtics $i]*$xfac+$xoff] [expr 5+\$tmarg_$obs+\$ypix_$obs] \
		-width 2 -tag "xtics"
	}  
    }
}

proc draw_ytics { obs ymin ymax yrange yfac yoff } {
    global lmarg_$obs 
    .$obs.c delete "ytics"
    if { $yrange > 0.0 } {
	set yexp [expr floor(log10($yrange))]
	set ypow [expr pow(10,$yexp)]
	set bas_range [expr $yrange/(1.0*$ypow)]
	set dist [expr 2.0*$ypow]
	if {  $bas_range < 6.0 } { set dist [expr 1.0*$ypow] }
	if {  $bas_range < 3.0 } { set dist [expr 0.5*$ypow] }
	if {  $bas_range < 1.7 } { set dist [expr 0.25*$ypow] }

	set ytics [list [expr $dist*ceil($ymin/(1.0*$dist))] ]
	.$obs.c create text                                               \
	    [expr \$lmarg_$obs - 20] [expr $yoff-[lindex $ytics 0]*$yfac] \
	    -text "[lindex $ytics 0]" -tag "ytics"
	.$obs.c create line                                               \
	    [expr \$lmarg_$obs - 5] [expr $yoff-[lindex $ytics 0]*$yfac]  \
	    [expr \$lmarg_$obs] [expr $yoff-[lindex $ytics 0]*$yfac]      \
	    -width 2 -tag "ytics"
	set i 0
	while { [expr [lindex $ytics $i]+$dist] <= $ymax } {
	    lappend ytics [expr [lindex $ytics $i]+$dist]
	    incr i
	    .$obs.c create text                                                \
		[expr \$lmarg_$obs - 20] [expr $yoff-[lindex $ytics $i]*$yfac] \
		-text "[lindex $ytics $i]" -tag "ytics"
	    .$obs.c create line                                                \
		[expr \$lmarg_$obs - 5] [expr $yoff-[lindex $ytics $i]*$yfac]  \
		[expr \$lmarg_$obs] [expr $yoff-[lindex $ytics $i]*$yfac]      \
		-width 2 -tag "ytics"
	}  
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

