#############################################################
#                                                           #
#  Pong for ESPResSo                                        #
#                                                           #
#############################################################
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

set box_x 30
set box_y 30
set box_z 50

set paddle_pos_x 12
set paddle_pos_y 12
set paddle_pos_z 1
set paddle_spheres_per_row 9
set paddle_step_size 3

set wall_parts_x 15
set wall_parts_y 15
set wall_parts_z 30

set ball_stopped 0

#settings for ball-paddle-interaction
set epsilon 0.8
set sigma 4
set cutoff 2
set shift [expr sqrt(sqrt(6)) * $sigma]
set offset 0

#settings for ball-constraints-interaction
set c_epsilon 1
set c_sigma 0.5
set c_cutoff 2
set c_shift 1
set c_offset 0

set skin 0.01
set timestep 0.0001
set steps 10

setmd periodicity 0 0 0

setmd box_l $box_x $box_y $box_z


#tk-stuff
bind . <Destroy> onExitHandler
wm geometry . -0+0

proc onExitHandler {} {
    # terminate vmd
    eval "exec rm -f [glob -nocomplain .lock*]"    
}

frame .case2 -relief raised -border 2
pack .case2 -expand 1 -fill both -in .
frame .case2.title
label .case2.title.l -justify left -text "Pong for ESPResSo"
pack .case2.title.l -fill both -side left -in .case2.title
pack .case2.title -pady {0 10} -fill x -side top -in .case2

focus .case2
update

#game-setup

constraint wall normal -1 0 0 dist [expr 1 - $box_x] type 3
constraint wall normal 1 0 0 dist 1 type 3
constraint wall normal 0 -1 0 dist [expr 1 - $box_y] type 3
constraint wall normal 0 1 0 dist 1 type 3

#set pid [expr $paddle_spheres_per_row * $paddle_spheres_per_row + 1]
#for {set i 0} {$i < $wall_parts_x} {incr i} {
#    for {set j 0} {$j < $wall_parts_z} {incr j} {
#	part $pid pos [expr $i * ($box_x / $wall_parts_x)] 0 [expr $j * ($box_z / $wall_parts_z)] fix type 4
#	incr pid
#	part $pid pos [expr $i * ($box_x / $wall_parts_x)] $box_y [expr $j * ($box_z / $wall_parts_z)] fix type 4
#	incr pid
#    }
#}
#for {set i 0} {$i < $wall_parts_y} {incr i} {
#    for {set j 0} {$j < $wall_parts_z} {incr j} {
#	part $pid pos 0 [expr $i * ($box_y / $wall_parts_y)] [expr $j * ($box_z / $wall_parts_z)] fix type 4
#	incr pid
#	part $pid pos $box_x [expr $i * ($box_y / $wall_parts_y)] [expr $j * ($box_z / $wall_parts_z)] fix type 4
#	incr pid
#    }
#}

inter 0 fene 30 5

part 0 pos [expr $box_x / 2] [expr $box_y / 2] [expr $box_z / 2] type 2 mass 0.1 ext_force 0 0 -15

proc del_paddle {} {
    global paddle_spheres_per_row
    # delete old paddle
    set pid 1
    for {set i 0} {$i < $paddle_spheres_per_row} {incr i} {
	for {set j 0} {$j < $paddle_spheres_per_row} {incr j} {
	    part $pid delete
	    incr pid
	}
    }    
}

proc set_paddle {} {
    global paddle_pos_x paddle_pos_y paddle_pos_z paddle_spheres_per_row
    #create new paddle
    set pid 1
    for {set i 0} {$i < $paddle_spheres_per_row} {incr i} {
	for {set j 0} {$j < $paddle_spheres_per_row} {incr j} {
	    set pid [expr $i * $paddle_spheres_per_row + $j + 1]
	    part $pid pos [expr $paddle_pos_x + $i] [expr $paddle_pos_y + $j] $paddle_pos_z type 1 mass 3.0
	    if {$i==0 || $i==[expr $paddle_spheres_per_row - 1] || $j==0 || $j==[expr $paddle_spheres_per_row - 1]} {
	        part $pid fix
	    }
	}
    }
    #create bonds
    for {set i 0} {$i < $paddle_spheres_per_row} {incr i} {
	for {set j 2} {$j <= $paddle_spheres_per_row} {incr j} {
	    set pid [expr $i * $paddle_spheres_per_row + $j]
	    part $pid bond 0 [expr $pid - 1]
	}
   }
    for {set i 1} {$i < $paddle_spheres_per_row} {incr i} {
	for {set j 1} {$j <= $paddle_spheres_per_row} {incr j} {
	    set pid [expr $i * $paddle_spheres_per_row + $j]
	    part $pid bond 0 [expr $pid - $paddle_spheres_per_row]
	}
   }
}

#paddle moves
proc paddle_left {} {
    global paddle_pos_x box_x paddle_spheres_per_row paddle_step_size
    if {$paddle_pos_x >= $paddle_step_size} {
	for {set i 1} {$i <= [expr $paddle_spheres_per_row * $paddle_spheres_per_row]} {incr i} {
	    scan [part $i print pos] "%f %f %f" px py pz
	    part $i pos [expr $px - $paddle_step_size] $py $pz
	}
	set paddle_pos_x [expr $paddle_pos_x - $paddle_step_size]
    }
}

proc paddle_right {} {
    global paddle_pos_x box_x paddle_spheres_per_row paddle_step_size
    if {$paddle_pos_x <= [expr $box_x - $paddle_step_size - $paddle_spheres_per_row]} {
	for {set i 1} {$i <= [expr $paddle_spheres_per_row * $paddle_spheres_per_row]} {incr i} {
	    scan [part $i print pos] "%f %f %f" px py pz
	    part $i pos [expr $px + $paddle_step_size] $py $pz
 	}
	set paddle_pos_x [expr $paddle_pos_x + $paddle_step_size]
    }
}

proc paddle_up {} {
    global paddle_pos_y box_y paddle_spheres_per_row paddle_step_size
    if {$paddle_pos_y <= [expr $box_y - $paddle_step_size - $paddle_spheres_per_row]} {
	for {set i 1} {$i <= [expr $paddle_spheres_per_row * $paddle_spheres_per_row]} {incr i} {
	    scan [part $i print pos] "%f %f %f" px py pz
	    part $i pos $px [expr $py + $paddle_step_size] $pz
	}
	set paddle_pos_y [expr $paddle_pos_y + $paddle_step_size]
    }
}

proc paddle_down {} {
    global paddle_pos_y box_y paddle_spheres_per_row paddle_step_size
    if {$paddle_pos_y >= $paddle_step_size} {
	for {set i 1} {$i <= [expr $paddle_spheres_per_row * $paddle_spheres_per_row]} {incr i} {
	    scan [part $i print pos] "%f %f %f" px py pz
	    part $i pos $px [expr $py - $paddle_step_size] $pz
	}
	set paddle_pos_y [expr $paddle_pos_y - $paddle_step_size]
    }
}

bind all <KeyPress-a> {paddle_left}
bind all <KeyPress-d> {paddle_right}
bind all <KeyPress-w> {paddle_up}
bind all <KeyPress-s> {paddle_down}

proc checkball {} {
    global ball_stopped
    scan [part 0 print pos] "%f %f %f" px py pz
    if {$pz > 40 && $ball_stopped == 0} {
	part 0 v 0 0 0
	set ball_stopped 1
    }
    if {$pz < 40} {
	set ball_stopped 0
    }
    if {$pz < 1} {
	exit;
    }
}

#vmd
proc prepare_vmd_connection_special { {filename "vmd"} {wait "0"} {start "1" } } {
    global box_x box_y box_z
    writepsf "$filename.psf"
    writepdb "$filename.pdb"
    for {set port 10000} { $port < 65000 } { incr port } {
        catch {imd connect $port} res
        if {$res == ""} break
    }
    set HOSTNAME [exec hostname]
    set vmdout_file [open "vmd_start.script" "w"]
    puts $vmdout_file "mol load psf $filename.psf pdb $filename.pdb"
    puts $vmdout_file "rotate z by 180"
    puts $vmdout_file "rotate stop"
    puts $vmdout_file "scale 0.7"
    puts $vmdout_file "display resize 400 400"
    puts $vmdout_file "display distance 4"
#    puts $vmdout_file "logfile vmd.log"
#    puts $vmdout_file "rotate x by 0"
#    puts $vmdout_file "rotate stop"
#    puts $vmdout_file "axes location off"
    puts $vmdout_file "logfile off"

    puts $vmdout_file "mol selection {segname T001}"
    puts $vmdout_file "mol representation CPK"
    puts $vmdout_file "mol addrep 0"
    puts $vmdout_file "mol color SegName      "

    puts $vmdout_file "mol selection {segname T002}"
    puts $vmdout_file "mol representation CPK"
    puts $vmdout_file "mol addrep 0"
    puts $vmdout_file "mol color SegName      "

    puts $vmdout_file "mol selection {segname T004}"
    puts $vmdout_file "mol representation Lines"
    puts $vmdout_file "mol addrep 0"
    puts $vmdout_file "mol color SegName      "

#    puts $vmdout_file "mol selection {segname T001}     "
#    puts $vmdout_file "mol material Opaque      "
#    puts $vmdout_file "mol addrep top      "
#    puts $vmdout_file "mol selupdate 0 top 0    "
#    puts $vmdout_file "mol colupdate 0 top 0    "
#    puts $vmdout_file "mol scaleminmax top 0 0.000000 0.000000   "
#    puts $vmdout_file "mol smoothrep top 0 0    "
#    puts $vmdout_file "mol drawframes top 0 {now}    "
#    puts $vmdout_file "mol representation CPK 50.0 0.300000 20.0 1.000000  "
#    puts $vmdout_file "mol color SegName      "
#    puts $vmdout_file "mol selection {segname T000}     "
#    puts $vmdout_file "mol material Opaque      "
#    puts $vmdout_file "mol addrep top      "
#    puts $vmdout_file "mol selupdate 1 top 0    "
#    puts $vmdout_file "mol colupdate 1 top 0    "
#    puts $vmdout_file "mol scaleminmax top 1 0.000000 0.000000   "
#    puts $vmdout_file "mol smoothrep top 1 0    "
#    puts $vmdout_file "mol drawframes top 1 {now}    "
#    puts $vmdout_file "color Segname {T001} red"
#    puts $vmdout_file "color Segname {T000} blue"
#    puts $vmdout_file "color Display {Background} black"
    puts $vmdout_file ""
    puts $vmdout_file "graphics 0 materials on"
    puts $vmdout_file "material add wall"
    puts $vmdout_file "material change ambient wall 0.5"
    puts $vmdout_file "material change shininess wall 0"
    puts $vmdout_file "material change specular wall 1.0"
    puts $vmdout_file "material change diffuse wall 0"

    puts $vmdout_file "graphics 0 material wall"


    puts $vmdout_file "graphics 0 triangle {0 0 0} {$box_x 0 0} {$box_x 0 $box_z}"
    puts $vmdout_file "graphics 0 triangle {0 0 0} {0 0 $box_z} {$box_x 0 $box_z}"

    puts $vmdout_file "graphics 0 triangle {0 $box_y 0} {$box_x $box_y 0} {$box_x $box_y $box_z}"
    puts $vmdout_file "graphics 0 triangle {0 $box_y 0} {0 $box_y $box_z} {$box_x $box_y $box_z}"

    puts $vmdout_file "graphics 0 triangle {0 0 0} {0 $box_y 0} {0 $box_y $box_z}"
    puts $vmdout_file "graphics 0 triangle {0 0 0} {0 0 $box_z} {0 $box_y $box_z}"

    puts $vmdout_file "graphics 0 triangle {$box_x 0 0} {$box_x $box_y 0} {$box_x $box_y $box_z}"
    puts $vmdout_file "graphics 0 triangle {$box_x 0 0} {$box_x 0 $box_z} {$box_x $box_y $box_z}"

    puts $vmdout_file ""
    puts $vmdout_file ""
    puts $vmdout_file "imd connect $HOSTNAME $port"
    puts $vmdout_file "imd transfer 1"
    puts $vmdout_file "imd keep 1"
     close $vmdout_file
    if { $start == 0 } {
        puts "Start VMD in the same directory on the machine you with :"
        puts "vmd -e vmd_start.script &"
        imd listen $wait
    } else {
        exec vmd -e vmd_start.script &
    }
}


#start game

set_paddle

#ball-paddle-interaction
inter 1 2 lennard-jones $epsilon $sigma $cutoff $shift $offset

#ball-constraints-interaction
inter 3 2 lennard-jones $c_epsilon $c_sigma $c_cutoff $c_shift $c_offset

thermostat set langevin 1.0 0.5
setmd time_step $timestep
setmd skin $skin

prepare_vmd_connection_special

after 2000

while {1} {
    integrate $steps
    imd positions
    checkball
    update
    focus .case2
}
