#!/bin/sh
# tricking... the line after a these comments are interpreted as standard shell script \
    PLATFORM=`uname -s`; if [ "$1" != "" ]; then NP=$1; else NP=2; fi
# OSF1 \
    if test $PLATFORM = OSF1; then  exec dmpirun -np $NP $TCLMD_SOURCE/$PLATFORM/tcl_md $0 $*
# AIX \
    elif test $PLATFORM = AIX; then exec poe $TCLMD_SOURCE/$PLATFORM/tcl_md $0 $* -procs $NP
# Linux \
    else export EF_ALLOW_MALLOC_0=1; lamboot; exec mpirun -np $NP -nsigs $TCLMD_SOURCE/$PLATFORM/tcl_md $0 $*;
# \
    fi;

#############################################################
#                                                           #
#  Test System 1: NaCl crystal - Madelung constant         #
#                                                           #
#                                                           #
#  Created:       18.03.2003 by HL                          #
#  Last modified: 18.03.2003 by HL                          #
#                                                           #
#############################################################

puts " "
puts "==================================================="
puts "=                 madelung.tcl                    ="
puts "==================================================="
puts " "

#############################################################
#  Parameters                                               #
#############################################################

# The Madelung constant
set madelung_nacl 1.747564594633182190636212035544397403481

# System identification: 
set name  "madelung"
set ident "_t1"

# System parameters
#############################################################

# grid constant for NaCl crystal
set nacl_const 2.0
# replications of the unit cell in each direcrtion
set nacl_repl  4

set charge 1.0

# for further checks shift the grid in the simulation box
set x0 0.0
set y0 0.0
set z0 0.0

# Interaction parameters coulomb p3m)
#############################################################

set bjerrum   1.0
set p3m_alpha 0.408905
set p3m_rcut  7.95 
set p3m_mesh  16
set p3m_cao   7

# Integration parameters
#############################################################
setmd time_step 0.01
setmd skin      0.03
setmd gamma     1.0
setmd temp      0.0

set int_steps   1


# Other parameters
#############################################################
set tcl_precision 12

#############################################################
#  Setup System                                             #
#############################################################

set box_l [expr 2.0*$nacl_const*$nacl_repl]
setmd box_l $box_l $box_l $box_l

# Interaction setup
#############################################################
inter coulomb $bjerrum p3m $p3m_rcut $p3m_mesh $p3m_cao $p3m_alpha 0.0
#inter coulomb mesh_off 0.0 0.0 0.0

# Particle setup
#############################################################

set part_id 0
set nc $nacl_const
for {set i 0} { $i < $nacl_repl } {incr i} {
    set xoff [expr $x0+($i * 2.0*$nacl_const)]
    for {set j 0} { $j < $nacl_repl } {incr j} {    
	set yoff [expr $y0+($j * 2.0*$nacl_const)]
	for {set k 0} { $k < $nacl_repl } {incr k} {
	    set zoff [expr $z0+($k * 2.0*$nacl_const)]
	    
	    part $part_id pos [expr $xoff] [expr $yoff] [expr $zoff] type 0 q $charge
	    incr part_id
	    part $part_id pos [expr $xoff+$nc] [expr $yoff] [expr $zoff] type 1 q -$charge
	    incr part_id
	    part $part_id pos [expr $xoff] [expr $yoff+$nc] [expr $zoff] type 1 q -$charge
	    incr part_id
	    part $part_id pos [expr $xoff+$nc] [expr $yoff+$nc] [expr $zoff] type 0 q $charge
	    incr part_id
	    part $part_id pos [expr $xoff] [expr $yoff] [expr $zoff+$nc] type 1 q -$charge
	    incr part_id
	    part $part_id pos [expr $xoff+$nc] [expr $yoff] [expr $zoff+$nc] type 0 q $charge
	    incr part_id
	    part $part_id pos [expr $xoff] [expr $yoff+$nc] [expr $zoff+$nc] type 0 q $charge
	    incr part_id
	    part $part_id pos [expr $xoff+$nc] [expr $yoff+$nc] [expr $zoff+$nc] type 1 q -$charge
	    incr part_id

	}
    }
}

# print system specification
set n_part [setmd n_part]
puts "NaCl crystal with $n_part particle (grid $nc, replic $nacl_repl)"
puts "in cubic box of size $box_l. Crystal offset: $x0 $y0 $z0"
puts "Interactions: \n{ [inter] }"

writepdb "$name$ident.pdb"

puts "\nIntegrate $int_steps step(s)"
integrate $int_steps

# calculate Madelung constant
set energy "[analyze energy]"
puts "\nEnergy: $energy"
set calc_mad_const [expr -2.0*$nacl_const*[lindex [lindex $energy 0] 1] / $n_part]
puts "\nMadelung constant:      $calc_mad_const"
puts "relative error:         [expr ($calc_mad_const-$madelung_nacl)/$madelung_nacl]"

# check forces
set max_force 0.0
for {set i 0} { $i < $n_part } {incr i} {
    set force [part $i print f]
    set fx [lindex $force 0]
    if { $fx > $max_force } { set max_force $fx }
    set fy [lindex $force 1]
    if { $fy > $max_force } { set max_force $fy }
    set fz [lindex $force 2]
    if { $fz > $max_force } { set max_force $fz }
}
puts "Maximal force component: $max_force"