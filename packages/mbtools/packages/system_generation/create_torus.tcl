#
# Routines for placing lipids uniformly on a torus
#
#

namespace eval ::system_generation {}

proc ::system_generation::find_deltafudge { r1 r2 delta_in nlipids } {
    set pi 3.141592
    set deltafudge 0.01
    set fudgechange 0.001
    set ncount 0
    set changecount 0
    while { $ncount != $nlipids && $changecount < 500 } {
	set delta [expr $delta_in + $deltafudge]    
	set mtheta [expr int(floor(2*$pi*$r2/sqrt($delta)))]
	set dtheta [expr 2*$pi*$r2/(1.0*$mtheta)]
	set dphi [expr $delta/$dtheta]
	set ncount 0
	for { set m 0 } { $m < $mtheta } { incr m } {
	    set theta [expr 2*$pi*$m/($mtheta*1.0)]
	    set mphi [expr int(floor(2*$pi*($r1-$r2*cos($theta))/$dphi))]
	    for { set n 0 } { $n < $mphi } { incr n } {
#		set phi [expr 2*$pi*$n/(1.0*$mphi)]
		incr ncount
	    }
	}

	if { $ncount > $nlipids && $fudgechange < 0 } {
	    set fudgechange [expr -$fudgechange/(2.0) ]
	    incr changecount
	}
	if { $ncount < $nlipids && $fudgechange > 0 } {
	    set fudgechange [expr -$fudgechange/(2.0) ]
	    incr changecount
	}

	set deltafudge [expr $deltafudge + $fudgechange ]

#	puts "$ncount $fudgechange $deltafudge"
    }
    set deltafudge [expr $deltafudge - $fudgechange ]
    return $deltafudge 
 
}

proc ::system_generation::create_torus {n_lipids box_l beads_per_lipid area_lipid ratio {bond_l "1.0"} } { 
    puts "Creating Torus: A/Lipid: $area_lipid ratio: $ratio bond_l: $bond_l"
    set pi 3.141592
    set del [expr $beads_per_lipid/2.0]

    set rad2_tot [expr sqrt($n_lipids*$area_lipid/(8.0*$pi*$pi*$ratio)) ]
    set rad1 [expr $rad2_tot*$ratio]

    puts "r1: $rad1 r2: $rad2_tot"

    set n_lipids_outer [expr 4*$pi*$pi*$rad1*($rad2_tot+$del)/$area_lipid]
    set n_lipids_inner [expr 4*$pi*$pi*$rad1*($rad2_tot-$del)/$area_lipid]

    set n_lipids_outer [expr int(($n_lipids_outer))]
    set n_lipids_inner [expr $n_lipids - $n_lipids_outer]
    puts "N+: $n_lipids_outer N-:$n_lipids_inner"


    # Determine r1 and r2 for the outer layer
    set rad2 [expr $rad2_tot + $del ]
    puts "r1: $rad1 r2: $rad2"
 
    set deltafudge [ find_deltafudge $rad1 $rad2 $area_lipid $n_lipids_outer ]
    set delta [expr $area_lipid + $deltafudge]
    
    puts "Old Area_lipid: $area_lipid New Area_lipid: $delta"

    # Now go through and place all lipids in the outer layer
    set mtheta [expr int((2*$pi*$rad2/sqrt($delta)))]
    set dtheta [expr 2*$pi*$rad2/(1.0*$mtheta)]
    set dphi [expr $delta/$dtheta]
    set ncount 0
    for { set m 0 } { $m < $mtheta } { incr m } {
	set theta [expr 2*$pi*$m/($mtheta*1.0)]
	set mphi [expr int((2*$pi*($rad1-$rad2*cos($theta))/$dphi))]
	for { set n 0 } { $n < $mphi } { incr n } {
	    set phi [expr 2*$pi*$n/(1.0*$mphi)]
	    set rx [expr ($rad1 - $rad2*cos($theta))*cos($phi) + [lindex $box_l 0]*0.5]
	    set ry [expr ($rad1 - $rad2*cos($theta))*sin($phi) + [lindex $box_l 1]*0.5]
	    set rz [expr $rad2*sin($theta) + [lindex $box_l 2]*0.5 ]

	    # Find the normal vector
	    set nx [expr -cos($theta)*cos($phi)]
	    set ny [expr -cos($theta)*sin($phi)]
	    set nz [expr sin($theta)]

	    # Place the lipid
	    set mol 0

	    for { set b 0 } { $b < $beads_per_lipid } {incr b } {
		set partnum [expr $ncount*$beads_per_lipid + $b ]
		lappend mol $partnum
#		set bond_l 1.0
		set posx [expr $rx+(1.0-$b)*$bond_l*$nx]
		set posy [expr $ry+(1.0-$b)*$bond_l*$ny]
		set posz [expr $rz+(1.0-$b)*$bond_l*$nz]
		if { $b == 0 } { set ptype 0 } else { set ptype 1 }
		part $partnum pos $posx $posy $posz type $ptype fix 1 1 1

		#Bonding potentials
		if { $b > 0 } {
		    part [expr $partnum - 1] bond 0 $partnum
		}

		# Pseudo bending potentials
		if { $b > 1 } {
		    part [expr $partnum -2] bond 1 $partnum	
		}

	    }
	    lappend topo $mol
	    incr ncount
	}

    }
    puts "Placed $ncount lipids in outer layer"

    #------------------ Inner ---------#
    # Determine r1 and r2 for the inner layer
    set rad2 [expr $rad2_tot - $del ]
    puts "r1: $rad1 r2: $rad2"
 
    set deltafudge [ find_deltafudge $rad1 $rad2 $area_lipid $n_lipids_inner ]
    set delta [expr $area_lipid + $deltafudge]

    # Now go through and place all lipids in the outer layer
    set mtheta [expr int((2*$pi*$rad2/sqrt($delta)))]
    set dtheta [expr 2*$pi*$rad2/(1.0*$mtheta)]
    set dphi [expr $delta/$dtheta]
    set innercount 0
    for { set m 0 } { $m < $mtheta } { incr m } {
	set theta [expr 2*$pi*$m/($mtheta*1.0)]
	set mphi [expr int((2*$pi*($rad1-$rad2*cos($theta))/$dphi))]
	for { set n 0 } { $n < $mphi } { incr n } {
	    set phi [expr 2*$pi*$n/(1.0*$mphi)]
	    set rx [expr ($rad1 - $rad2*cos($theta))*cos($phi) + [lindex $box_l 0]*0.5]
	    set ry [expr ($rad1 - $rad2*cos($theta))*sin($phi) + [lindex $box_l 1]*0.5]
	    set rz [expr $rad2*sin($theta) + [lindex $box_l 2]*0.5 ]
	    
	    # Find the normal vector pointing inwards
	    set nx [expr cos($theta)*cos($phi)]
	    set ny [expr cos($theta)*sin($phi)]
	    set nz [expr -sin($theta)]
	    
	    # Place the lipid
	    set mol 1
	    
	    for { set b 0 } { $b < $beads_per_lipid } {incr b } {
		set partnum [expr $ncount*$beads_per_lipid + $b ]
		lappend mol $partnum
#		set bond_l 1.0
		set posx [expr $rx+(1.0-$b)*$bond_l*$nx]
		set posy [expr $ry+(1.0-$b)*$bond_l*$ny]
		set posz [expr $rz+(1.0-$b)*$bond_l*$nz]
		if { $b == 0 } { set ptype 0 } else { set ptype 1 }
		part $partnum pos $posx $posy $posz type $ptype fix 1 1 1
		
		#Bonding potentials
		if { $b > 0 } {
		    part [expr $partnum - 1] bond 0 $partnum
		}
		
		# Pseudo bending potentials
		if { $b > 1 } {
		    part [expr $partnum -2] bond 1 $partnum	
		}
	
	    }
	    lappend topo $mol
	    incr ncount
	    incr innercount
	}
    }
    puts "Placed $innercount lipids in inner layer"
    #    exit
    
    puts "Toroidal vesicle created" 
    flush stdout
    
    return $topo
}


