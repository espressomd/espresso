#############################################################
#                                                           #
# auxiliary.tcl                                             #
# =======                                                   #
#                                                           #
# Several additional auxiliary functions for Espresso.      #
#                                                           #
#############################################################
#
# Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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

# Deprecation warning
proc warn_deprecated { fname version } {
    puts stderr "WARNING: The function $fname is deprecated since version $version"
    puts stderr "         and will be removed in some future version."
}


# timeStamp
# ---------
# 
# Squeezes a prefix into a given filename and adds the 
# current date as postfix before the suffix.
#
# Input:
# - complete path 'destination'
# - 'prefix' to be squeezed between path- and file-name
# - 'postfix' to be added before 'suffix'; 
#   if it is '-1', the current date will be used
# 
# 
#############################################################

proc timeStamp { destination prefix postfix suffix } {

    # Customize 'destination' to include informations about the current time and the Deserno-wishlist used
    puts -nonewline "    Building a customized destination file name... "
    flush stdout

    # Squeeze a timestamp just before the 'suffix'
    if { [string compare $postfix -1]==0 } {
	set tmp_var [clock format [clock seconds] -format %y%m%d]
    } else {
	set tmp_var $postfix
    }
    set destination "[string range $destination 0 [string last .$suffix $destination]]$tmp_var.$suffix"

    # Squeeze the 'prefix' in between path and filename
    set tmp_var [expr [string first [lindex [split $destination \"/\"] end] $destination]-1]
    set destination "[string range $destination 0 $tmp_var]$prefix\_[lindex [split $destination \"/\"] end]"
    puts "Done."

    return $destination
}

proc analysisInit { stat stat_out N_P MPC simtime { noted "na" } { notedD "na" } } {
    if {[llength $stat]>0} {
	puts -nonewline "    Analyzing at t=$simtime: "; flush stdout
	puts -nonewline "Init all ([analyze set chains 0 $N_P $MPC]): "; flush stdout
	puts -nonewline $stat_out "t "; set tmp_stat "$simtime "
	foreach i $stat {
	    if { $i == "g123" } { puts -nonewline "init g123 ([analyze g123 -init 0 $N_P $MPC ]), "; flush stdout }
	    set tmp_var [analyze $i]
	    puts -nonewline $stat_out "$i "; set tmp_stat "$tmp_stat $tmp_var"
	    puts -nonewline "$i=$tmp_var; "; flush stdout
	}
	if { $noted != "na" } { puts -nonewline $stat_out "$noted "; set tmp_stat "$tmp_stat $notedD" }
	puts $stat_out "\n$tmp_stat"; flush $stat_out
	puts "Done."
    }
}

proc analysis { stat stat_out N_P MPC simtime { noted "na" } } {
    puts -nonewline $stat_out "$simtime "
    foreach i $stat {
	set tmp_stat [analyze $i]
	puts -nonewline ", $i=$tmp_stat"; flush stdout
	puts -nonewline $stat_out "$tmp_stat "
    }
    if { $noted != "na" } { puts -nonewline $stat_out "$noted " }
    puts $stat_out " "; flush $stat_out
}

#
# calc_aniso_box
# --------------
#
# Calculate box dimensions of an unisotropic box given the volume
# and the box ratios { box_x/box_y box_y/box_z } (HJL).
#
#############################################################
proc calc_aniso_box { volume box_ratio } {
    set box_x 1.0
    set box_y [expr 1.0/[lindex $box_ratio 0]]
    set box_z [expr 1.0/([lindex $box_ratio 0]*[lindex $box_ratio 1])]
    set temp [expr $box_x * $box_y * $box_z]
    set temp [expr pow($volume/$temp, 1.0/3.0)]
    set box_x [expr $box_x*$temp]
    set box_y [expr $box_y*$temp]
    set box_z [expr $box_z*$temp]
    set box "$box_x $box_y $box_z"
    return $box
}

#
# tune_cells
# ----------
# 
# Tunes the skin and cell-grid to balance verlet_reuse vs. 
# organisational overhead.
#
#############################################################

proc tune_cells { { int_steps 1000 } { min_cells 1 } { max_cells 30 } { tol_cells 0.099 } } {
    set msg_string "Tuning current system with [setmd n_part] particles and initial skin=[setmd skin] and cell_grid=[setmd cell_grid] to:"
    array set int_times [list "0 0 0" 0.0]
    proc eval_cells { i num_cells } {
	upvar int_steps int_steps  int_times int_times  msg_string msg_string

	set maxcells [expr int(pow($num_cells,3))]
	# calculate minimal number of cells
	if { [setmd n_part] >= 1000 } {
	    set mincells [expr int(5e-7*pow([setmd n_part]/[setmd n_nodes],2))]
	} {
	    # for small particle numbers, the boundary is wrong, and even one cells does not
	    # mean really excessive computation time
	    set mincells 1
	}
	# minimal number of cells due to node grid
	set phys_cell_min 8
	foreach dim [setmd node_grid] {
	    if {$dim > 1} { set phys_cell_min [expr $phys_cell_min/2] }
	}
	if { $mincells < $phys_cell_min } { set mincells $phys_cell_min }
	if { $maxcells <= $mincells } { set maxcells [expr 8*$mincells] }

	if { [setmd min_num_cells] < $mincells } {
	    # minimum is larger than before, increase maximum before to be safe
	    setmd max_num_cells $maxcells
	    setmd min_num_cells $mincells
	} {
	    # minimum is not larger than before, changing it first is safe
	    setmd min_num_cells $mincells
	    setmd max_num_cells $maxcells
	}
	setmd skin 0.0
	set msg_string "$msg_string\n    run [expr $i+1] (t=[setmd time]; max_cells=[format %.3f $num_cells],"
	set msg_string "$msg_string skin=[format %.3g [setmd skin [expr 0.999*[setmd max_skin]]]],"
	set msg_string "$msg_string grid=[setmd cell_grid], cell_size=[format %.3f [lindex [set tmp_cell [setmd cell_size]] 0]]"
	set msg_string "$msg_string [format %.3f [lindex $tmp_cell 1]] [format %.3f [lindex $tmp_cell 2]], max_cut=[format %.3f [setmd max_cut]],"
	set msg_string "$msg_string max_range=[format %.3f [setmd max_range]])"
	if { [llength [set int_time [array get int_times "[setmd cell_grid]"]]]==2 } {
	    set integ_time [lindex $int_time 1]
	} else {
	    set integ_time [lindex [time {
		integrate $int_steps
	    } ] 0]
	    array set int_times [list "[setmd cell_grid]" $integ_time]
	}
	set msg_string "$msg_string => [expr $integ_time/1.0e6]s (vr=[format %.3f [setmd verlet_reuse]])"
	return $integ_time
    }
    set i 0; set ax $min_cells; set bx [expr 0.5*($min_cells+$max_cells)]; set cx $max_cells
    set GR 0.61803399; set GC [expr 1.0-$GR]; set x0 $ax; set x3 $cx
    if { abs($cx-$bx) > abs($bx-$ax) } { set x1 $bx; set x2 [expr $bx + $GC*($cx-$bx)] } else { set x2 $bx; set x1 [expr $bx - $GC*($bx-$ax)] }
    set f1 [eval_cells $i $x1]; incr i 
    set f2 [eval_cells $i $x2]; incr i
    while { abs($x3-$x0) > $tol_cells*(abs($x1) + abs($x2))} {
	if { $f2 < $f1 } {
	    set x0 $x1; set x1 $x2; set x2 [expr $GR*$x1 + $GC*$x3]
	    set f1 $f2; set f2 [eval_cells $i $x2]; incr i
	} else {
	    set x3 $x2; set x2 $x1; set x1 [expr $GR*$x2 + $GC*$x0]
	    set f2 $f1; set f1 [eval_cells $i $x1]; incr i
	}
    }
    set msg_string "$msg_string\nFinal set: "
    if { $f1 < $f2 } { set fs [eval_cells $i $x1] } else { set fs [eval_cells $i $x2] }
    return $msg_string
}

#
# stop_particles
# -------------
# 
# Sets the velocities and forces of all particles currently
# stored in Espresso to zero.
#
#############################################################

proc stop_particles { } { 
  warn_deprecated stop_particles 3.1.2
  puts stderr "         Use the function kill_particle_motion instead."
  kill_particle_motion
}

#
# stopParticles
# -------------
# 
# Sets the velocities and forces of all particles currently
# stored in Espresso to zero.
#
#############################################################



proc stopParticles { } {
  warn_deprecated stopParticles 3.1.2
  puts stderr "         Use the function kill_particle_motion and kill_particle_forces instead."
  puts -nonewline "        Setting all particles' velocities and forces to zero... "; flush stdout
  kill_particle_motion
  kill_particle_forces
  puts ".,. Done (stopped [setmd n_part] particles and as many forces)."; flush stdout
}

#
# center of mass
# --------------
#
# Calculate the center of mass of the system stored in Espresso
#
#############################################################

proc system_com { } {
  puts stderr "         Use the function system_CMS instead."
  return [system_CMS]
}

#
# center of mass motion
# ---------------------
#
# Calculate the center of mass velocity of the system stored in Espresso
#
#############################################################

proc system_com_vel {} {
  puts stderr "         Use the function system_CMS_velocity instead."
  return [system_CMS_velocity]
}

#
# galileiTransformParticles 
# -------------------------
# 
# Perform a Galilei Transformation on all Particles stored 
# in Espresso such that the center of mass motion of the 
# system is zero. Returns old center of mass motion.
#
#############################################################

proc galileiTransformParticles {} {
  puts stderr "         Use the function galilei_transform instead."
  set com_vel [system_CMS_velocity]
  
  galilei_transform
  
  return $com_vel
}

#
# Prepare Connection to VMD
# -------------------------
# 
#
#############################################################

proc prepare_vmd_connection { args } {
  # default values, as before
  set filename "vmd"
  set wait 0
  set start 1
  set hostname [info hostname]
  set draw_constraints 0

  # parse off filename, both for old and new style
  if {[llength $args] > 0} {
    set filename [lindex $args 0]
    set args [lrange $args 1 end]
  }
	
  if {[llength $args] > 0 && [string is integer [lindex $args 0]]} {
    # old, fixed position format of other args: wait start draw_constraints
    # detect simply if a number instead of keyword follows
    set wait [lindex $args 0]
    if {[llength $args] > 1} { set start [lindex $args 1] }
    if {[llength $args] > 2} { set draw_constraints [lindex $args 2] }
    # old format did not allow to specify arguments for writevsf
    set vsf_args ""
  } {
    # new format with keyword arguments
    while {[llength $args] > 0} {
      switch [lindex $args 0] {
        "wait" {
          set wait [lindex $args 1]
          set args [lrange $args 2 end]
        }
        "start" {
          set start 1
          # if we start VMD, it is necessarily on localhost
          set hostname "localhost"
          set args [lrange $args 1 end]
        }
        "constraints" {
          set draw_constraints 1
          set args [lrange $args 1 end]
        }
        "localhost" {
          set hostname "localhost"
          set args [lrange $args 1 end]
        }
        default {
          # unknown keyword, assume it is for writevsf
          break
        }
      }
    }
    # remaining arguments we assume to be for writevsf
    set vsf_args $args
  }
  
  # structure information only
  set f [open "$filename.vsf" "w"]
  eval writevsf $f $vsf_args
  close $f
  
  for {set port 10000} { $port < 65000 } { incr port } {
      catch {imd connect $port} res
      if {$res == ""} {
          break
      }
  }
  
  set vmdout_file [open "${filename}.vmd_start.script" "w"]
  
  puts $vmdout_file "logfile vmd.log"
  puts $vmdout_file "mol load vsf $filename.vsf"
  puts $vmdout_file "rotate stop"
  puts $vmdout_file "mol modstyle 0 0 CPK 1.800000 0.300000 8.000000 6.000000"
  puts $vmdout_file "mol modcolor 0 0 Name"
  puts $vmdout_file "imd connect $hostname $port"
  puts $vmdout_file "imd transfer 1"
  puts $vmdout_file "imd keep 1"
  
  # draw constraints  
  if {$draw_constraints != "0"} {
    foreach c [ constraint ] {
      set id [ lindex $c 0 ]
      set type [ lindex $c 1 ]
      
      puts $vmdout_file "draw color gray"
      
      if {$type == "sphere"} {
        set centerx [ lindex $c 3 ]
        set centery [ lindex $c 4 ]
        set centerz [ lindex $c 5 ]
        set radius [ lindex $c 7 ]
        
        puts $vmdout_file "draw sphere \{ $centerx $centery $centerz \} radius $radius"
      } elseif {$type == "rhomboid"} {
        set p_x [lindex $c 3]
        set p_y [lindex $c 4]
        set p_z [lindex $c 5]
        
        set a_x [lindex $c 7]
        set a_y [lindex $c 8]
        set a_z [lindex $c 9]
        
        set b_x [lindex $c 11]
        set b_y [lindex $c 12]
        set b_z [lindex $c 13]
        
        set c_x [lindex $c 15]
        set c_y [lindex $c 16]
        set c_z [lindex $c 17]
        
        puts $vmdout_file "draw triangle \{$p_x $p_y $p_z\} \{[expr $p_x+$a_x] [expr $p_y+$a_y] [expr $p_z+$a_z]\} \{[expr $p_x+$b_x] [expr $p_y+$b_y] [expr $p_z+$b_z]\}"
        puts $vmdout_file "draw triangle \{$p_x $p_y $p_z\} \{[expr $p_x+$b_x] [expr $p_y+$b_y] [expr $p_z+$b_z]\} \{[expr $p_x+$c_x] [expr $p_y+$c_y] [expr $p_z+$c_z]\}"
        puts $vmdout_file "draw triangle \{$p_x $p_y $p_z\} \{[expr $p_x+$a_x] [expr $p_y+$a_y] [expr $p_z+$a_z]\} \{[expr $p_x+$c_x] [expr $p_y+$c_y] [expr $p_z+$c_z]\}"
        puts $vmdout_file "draw triangle \{[expr $p_x+$a_x+$b_x+$c_x] [expr $p_y+$a_y+$b_y+$c_y] [expr $p_z+$a_z+$b_z+$c_z]\} \{[expr $p_x+$a_x+$c_x] [expr $p_y+$a_y+$c_y] [expr $p_z+$a_z+$c_z]\} \{[expr $p_x+$b_x+$c_x] [expr $p_y+$b_y+$c_y] [expr $p_z+$b_z+$c_z]\}"
        puts $vmdout_file "draw triangle \{[expr $p_x+$a_x+$b_x+$c_x] [expr $p_y+$a_y+$b_y+$c_y] [expr $p_z+$a_z+$b_z+$c_z]\} \{[expr $p_x+$a_x+$c_x] [expr $p_y+$a_y+$c_y] [expr $p_z+$a_z+$c_z]\} \{[expr $p_x+$a_x+$b_x] [expr $p_y+$a_y+$b_y] [expr $p_z+$a_z+$b_z]\}"
        puts $vmdout_file "draw triangle \{[expr $p_x+$a_x+$b_x+$c_x] [expr $p_y+$a_y+$b_y+$c_y] [expr $p_z+$a_z+$b_z+$c_z]\} \{[expr $p_x+$b_x+$c_x] [expr $p_y+$b_y+$c_y] [expr $p_z+$b_z+$c_z]\} \{[expr $p_x+$a_x+$b_x] [expr $p_y+$a_y+$b_y] [expr $p_z+$a_z+$b_z]\}"
        puts $vmdout_file "draw triangle \{[expr $p_x+$c_x] [expr $p_y+$c_y] [expr $p_z+$c_z]\} \{[expr $p_x+$a_x+$c_x] [expr $p_y+$a_y+$c_y] [expr $p_z+$a_z+$c_z]\} \{[expr $p_x+$a_x] [expr $p_y+$a_y] [expr $p_z+$a_z]\}"
        puts $vmdout_file "draw triangle \{[expr $p_x+$c_x] [expr $p_y+$c_y] [expr $p_z+$c_z]\} \{[expr $p_x+$a_x+$c_x] [expr $p_y+$a_y+$c_y] [expr $p_z+$a_z+$c_z]\} \{[expr $p_x+$b_x+$c_x] [expr $p_y+$b_y+$c_y] [expr $p_z+$b_z+$c_z]\}"
        puts $vmdout_file "draw triangle \{[expr $p_x+$a_x] [expr $p_y+$a_y] [expr $p_z+$a_z]\} \{[expr $p_x+$b_x] [expr $p_y+$b_y] [expr $p_z+$b_z]\} \{[expr $p_x+$a_x+$b_x] [expr $p_y+$a_y+$b_y] [expr $p_z+$a_z+$b_z]\}"
        puts $vmdout_file "draw triangle \{[expr $p_x+$b_x+$c_x] [expr $p_y+$b_y+$c_y] [expr $p_z+$b_z+$c_z]\} \{[expr $p_x+$b_x] [expr $p_y+$b_y] [expr $p_z+$b_z]\} \{[expr $p_x+$a_x+$b_x] [expr $p_y+$a_y+$b_y] [expr $p_z+$a_z+$b_z]\}"
        puts $vmdout_file "draw triangle \{[expr $p_x+$a_x] [expr $p_y+$a_y] [expr $p_z+$a_z]\} \{[expr $p_x+$a_x+$b_x] [expr $p_y+$a_y+$b_y] [expr $p_z+$a_z+$b_z]\} \{[expr $p_x+$a_x+$c_x] [expr $p_y+$a_y+$c_y] [expr $p_z+$a_z+$c_z]\}"
        puts $vmdout_file "draw triangle \{[expr $p_x+$b_x] [expr $p_y+$b_y] [expr $p_z+$b_z]\} \{[expr $p_x+$c_x] [expr $p_y+$c_y] [expr $p_z+$c_z]\} \{[expr $p_x+$b_x+$c_x] [expr $p_y+$b_y+$c_y] [expr $p_z+$b_z+$c_z]\}"
      } elseif {$type == "cylinder"} {
        set c_x [lindex $c 3]
        set c_y [lindex $c 4]
        set c_z [lindex $c 5]
        
        set a_x [lindex $c 7]
        set a_y [lindex $c 8]
        set a_z [lindex $c 9]
        
        set r [lindex $c 11]
        
        set l [lindex $c 13]
        
        puts $vmdout_file "draw cylinder \{[expr $l*($c_x/$l -$a_x)] [expr $l*($c_y/$l -$a_y)] [expr $l*($c_z/$l -$a_z)]\} \{[expr $l*($c_x/$l +$a_x)] [expr $l*($c_y/$l +$a_y)] [expr $l*($c_z/$l +$a_z)]\} radius $r resolution 36"
      }
    }
  }
 
  close $vmdout_file
 
  if { $start == 0 } {
    puts "Start VMD in the same directory with:"
    puts "vmd -e ${filename}.vmd_start.script &"
  } else {
    exec vmd -e "${filename}.vmd_start.script" &
  }
  
  imd listen $wait
}

#
# has_feature
# -----------
# 
# Returns, if all given features have been compiled into the Espresso
# kernel.
# 
#############################################################

proc has_feature { args } {
    foreach feature $args {
	if { ! [regexp $feature [code_info]] } then { return 0 }
    }
    return 1
}

#
# require_feature 
# ---------------
# 
# Check, if the given features have been compiled into the Espresso
# kernel. Exits with an error if it is not.
# 
# Example:
# require_feature "ELECTROSTATICS" "LJCOS"
# 
#############################################################

proc require_feature { args } {
    if { ! [eval has_feature $args]} {
	puts "not compiled in: $args"
	exit -42
    }
}

#
# Create box for psf/pdb files
#
# Do not run integration with that, it will crash!!!
proc vmd_create_box { {type "0"} } {
    set pid [expr [setmd max_part] +1]; set res $pid
    set box [setmd box_l]
    part $pid pos 0 0 0 type $type ; incr pid
    part $pid pos 0 0 [lindex $box 2] type $type bond 0 [expr $pid-1]; incr pid
    part $pid pos 0 [lindex $box 1] [lindex $box 2] type $type bond 0 [expr $pid-1]; incr pid
    part $pid pos 0 [lindex $box 1] 0 type $type bond 0 [expr $pid-1]
    part $pid bond 0 [expr $pid-3]; incr pid
    part $pid pos [lindex $box 0] 0 0 type $type
    part $pid bond 0 [expr $pid-4]; incr pid
    part $pid pos [lindex $box 0] 0 [lindex $box 2] type $type bond 0 [expr $pid-1]
    part $pid bond 0 [expr $pid-4]; incr pid
    part $pid pos [lindex $box 0] [lindex $box 1] [lindex $box 2] type $type bond 0 [expr $pid-1]
    part $pid bond 0 [expr $pid-4]; incr pid
    part $pid pos [lindex $box 0] [lindex $box 1] 0 type $type bond 0 [expr $pid-1]
    part $pid bond 0 [expr $pid-4];    part $pid bond 0 [expr $pid-3]; incr pid
    return $res
}

# Delete the box for vmd 
proc vmd_delete_box { start } {
    for { set i $start } { $i < [expr $start+8] } { incr i } { part $i del } 
}
	
# Return the degrees of freedom of the system
# This is important, if you calculate the temperature
# of your system from the kinetic energy:
#   temp = 2 * E_kin / ( dgf * N )
proc degrees_of_freedom {} {
 if { [regexp "ROTATION" [code_info]] } {
   set dgf 6
 } else {
     set dgf 3
 }
 return $dgf
}

#
# copy_particles 
# ---------------
# 
# Copy a group of particles including their bonds. Positions can be
# shifted by an offset.  The particles can be given as a combination
# of lists or a ranges. The new particles in any case obtain
# consecutive identities after the largest current identity. Bonds and
# exclusions within the particle set are copied, but not bonds with
# particles outside the list. That is, if the particle set corresponds
# to a molecule, intramolecular bonds are preserved, but not
# intermolecular ones.
# 
# Example:
# copy_particles set {1 2 3 4} shift 0.0 0.0 0.0
# copy_particles set {1 2} set {3 4}
# copy_particles range 1 4
# 
#############################################################

proc copy_particles { args } {
    set shift {0 0 0}
    set parts {}

    # parse arguments
    while { $args != "" } {
	switch [lindex $args 0] {
	    "set" {
		set parts [concat $parts [lindex $args 1]]
		set args [lrange $args 2 end]
	    }
	    "range" {
		for {set id [lindex $args 1]} {$id <= [lindex $args 2]} {incr id} {
		    lappend parts $id
		}
		set args [lrange $args 3 end]
	    }
	    "shift" {
		set shift [lrange $args 1 3]
		set args [lrange $args 4 end]
	    }
	    default {
		error "copy_particles: wrong argument [lindex $args 0]"
	    }
	}
    }
    set parts [lsort -unique $parts]

    # copy loop.
    # step 1: all except bonds and exclusions
    set newid [setmd max_part]
    foreach id $parts {
	incr newid
	set newids($id) $newid

	# new shifted position
        set pos [vecadd [part $id print pos] $shift]

	# parameters of the particle for creating a copy
	set partcmd [lrange [part $id] 1 end]
	# look for the position and fix it up
	set p [lsearch $partcmd "pos"]
	set partcmd [eval lreplace {$partcmd} [expr $p + 1] [expr $p + 3] $pos]

	# remove bonds
	set p [lsearch $partcmd "bond"]
	if {$p != -1} { set partcmd [lreplace $partcmd $p [expr $p + 1]] }

	# remove exclusions
	set p [lsearch $partcmd "exclude"]
	if {$p != -1} {
            set pend [expr $p + 1]
            while {1} {
                set element [lindex $partcmd $pend]
                if {$element == "" || ![string is integer $element]} { break }
                incr pend
            }
            set partcmd [lreplace $partcmd $p $pend]
        }

        # and set the new particle
	eval part $newid $partcmd
    }

    # step 2: bonds. This makes sure all bonding partners exist
    foreach id $parts {
	set newid $newids($id)

	# new bond list
	set bonds {}
	foreach bond [lindex [part $id print bonds] 0] {
	    set newbond [lindex $bond 0]
	    # check whether partners are in the molecule and translate
	    foreach partner [lrange $bond 1 end] {
		if {[array names newids -exact $partner] != ""} {
		    lappend newbond $newids($partner)
		} {
		    # this bond has a partner outside, ignore
		    set newbond {}; break
		}
	    }
	    if {$newbond != {}} {
		eval part $newid bond $newbond
	    }
	}

	# new exclusions list
	set exclusions {}
        # check whether exclusions are in the molecule and translate
	foreach exclude [part $id print exclusions] {
            if {[array names newids -exact $exclude] != ""} {
                lappend exclusions $newids($exclude)
            }
        }
        if {$exclusions != {}} {
            eval part $newid exclude $exclusions
        }
    }
}
