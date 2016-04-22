proc h5md_init { data_path } {
	if { [setmd n_part] == 0 } {
		puts "Please set your particles before h5md initialisation\n"; flush stdout
		return
	}
	if { [file exists ${data_path}] == 1 } {
		h5mdfile H5Fopen "${data_path}"
		return
	} elseif { [file exists ${data_path}] == 0 } {
		# Create hdf5 file
		h5mdfile H5Fcreate "${data_path}"
	}
	# Create data groups
	h5mdfile H5Gcreate2 "particles"
	h5mdfile H5Gcreate2 "particles/atoms"
	h5mdfile H5Gcreate2 "particles/atoms/box"
	h5mdfile H5Gcreate2 "particles/atoms/mass"
	h5mdfile H5Gcreate2 "particles/atoms/position"
	h5mdfile H5Gcreate2 "particles/atoms/velocity"
	h5mdfile H5Gcreate2 "particles/atoms/force"
	h5mdfile H5Gcreate2 "parameters"
	h5mdfile H5Gcreate2 "parameters/vmd_structure"
	h5mdfile H5Gcreate2 "observables" 
	h5mdfile H5Gcreate2 "files"
	h5mdfile H5Gcreate2 "files/scripts"
	# Create datasets
	
	#box
	h5mdfile H5Screate_simple type double dims 3 3
	h5mdfile H5Pset_chunk dims 3 3
	h5mdfile H5Dcreate2 "particles/atoms/box/edges"
	h5mdfile H5Dopen2 "particles/atoms/box/edges"
	#assuming cuboid box
	for {set i 0 } {$i<3} {incr i} {
		h5mdfile H5_write_value value [lindex [setmd box_l] $i] index $i $i
	}
	h5mdfile H5Dwrite
	
	#mass	
	h5mdfile H5Screate_simple type double dims [setmd n_part] 
	h5mdfile H5Pset_chunk dims [setmd n_part] 
	h5mdfile H5Dcreate2 "particles/atoms/mass/value"
	h5mdfile H5Dopen2 "particles/atoms/mass/value"
	if { [string first MASS [code_info]] != -1 } {
		for  { set i 0 } { $i < [setmd n_part] } { incr i } {
			h5mdfile H5_write_value value [part $i print mass] index $i
		}
		h5mdfile H5Dwrite
	}
	# position
	h5mdfile H5Screate_simple type int dims 0 
	h5mdfile H5Pset_chunk dims 1 
	h5mdfile H5Dcreate2 "particles/atoms/position/step"
	h5mdfile H5Screate_simple type double dims 0 
	h5mdfile H5Pset_chunk dims 1 
	h5mdfile H5Dcreate2 "particles/atoms/position/time"
	h5mdfile H5Screate_simple type double dims 0 [setmd n_part] 3
	h5mdfile H5Pset_chunk dims 5 [setmd n_part] 3
	h5mdfile H5Dcreate2 "particles/atoms/position/value"
	# velocity
	h5mdfile H5Screate_simple type int dims 0 
	h5mdfile H5Pset_chunk dims 1 
	h5mdfile H5Dcreate2 "particles/atoms/velocity/step"
	h5mdfile H5Screate_simple type double dims 0 
	h5mdfile H5Pset_chunk dims 1 
	h5mdfile H5Dcreate2 "particles/atoms/velocity/time"
	h5mdfile H5Screate_simple type double dims 0 [setmd n_part] 3
	h5mdfile H5Pset_chunk dims 5 [setmd n_part] 3
	h5mdfile H5Dcreate2 "particles/atoms/velocity/value"
	# force
	h5mdfile H5Screate_simple type int dims 0 
	h5mdfile H5Pset_chunk dims 1 
	h5mdfile H5Dcreate2 "particles/atoms/force/step"
	h5mdfile H5Screate_simple type double dims 0 
	h5mdfile H5Pset_chunk dims 1 
	h5mdfile H5Dcreate2 "particles/atoms/force/time"
	h5mdfile H5Screate_simple type double dims 0 [setmd n_part] 3
	h5mdfile H5Pset_chunk dims 5 [setmd n_part] 3
	h5mdfile H5Dcreate2 "particles/atoms/force/value"
	# species
	h5mdfile H5Screate_simple type int dims [setmd n_part] 
	h5mdfile H5Pset_chunk dims [setmd n_part] 
	h5mdfile H5Dcreate2 "particles/atoms/species"

	h5md_write_species
	
	##################
	# vmd
	##################
    	set types [h5mdutil_get_types]
    	set names "X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg"
    	
    	#index of species
	h5mdfile H5Screate_simple type int dims [llength $types]
	h5mdfile H5Pset_chunk dims [llength $types]
	h5mdfile H5Dcreate2 "parameters/vmd_structure/indexOfSpecies"
	h5mdfile H5Dopen2 "parameters/vmd_structure/indexOfSpecies"
	set i 0
	foreach type $types {
		h5mdfile H5_write_value value $type index $i
		incr i
	}
	h5mdfile H5Dwrite
	
	#name
	h5mdfile H5Screate_simple type str dims [llength $types]
	h5mdfile H5Pset_chunk dims [llength $types]
	h5mdfile H5Dcreate2 "parameters/vmd_structure/name"
	h5mdfile H5Dopen2 "parameters/vmd_structure/name"
	set i 0
	foreach type $types {
		h5mdfile H5_write_value value [lindex $names [expr $i%112]] index $i
		incr i
	}
	h5mdfile H5Dwrite
	
	#type
	h5mdfile H5Screate_simple type str dims [llength $types]
	h5mdfile H5Pset_chunk dims [llength $types]
	h5mdfile H5Dcreate2 "parameters/vmd_structure/type"
	h5mdfile H5Dopen2 "parameters/vmd_structure/type"
	for {set i 0 } {$i < [llength $types] } {incr i} {
		h5mdfile H5_write_value value [lindex $names [expr $i%112]] index $i
	}
	h5mdfile H5Dwrite
	
	#write charge for vmd_structure
	h5mdfile H5Screate_simple type double dims [setmd n_part]
	h5mdfile H5Pset_chunk dims [setmd n_part]
	h5mdfile H5Dcreate2 "parameters/vmd_structure/charge"
	h5mdfile H5Dopen2 "parameters/vmd_structure/charge"
	for {set i 0 } {$i < [setmd n_part] } {incr i} {
		h5mdfile H5_write_value value [part $i pr q] index $i
	}
	h5mdfile H5Dwrite
	
	#bonds
	set froms ""
	set tos ""
	for { set from 0 } { $from < [setmd n_part] } { incr from } {
		if { [part $from] != "na" } then {
		    set bonds [lindex [part $from print bond] 0]
		    if { [llength $bonds] >0 } {
			    for { set i 0 } { $i < [llength $bonds] } { incr i } {
				set to [lindex $bonds $i 1]
				lappend froms [expr $from +1]; #VMD counts from 1 and not from 0 therefore +1
				lappend tos [expr $to +1]
			    }
		    }
		}
	}
	
	if { [llength $froms] >0 } {
		h5mdfile H5Screate_simple type int dims [llength $froms]
		h5mdfile H5Pset_chunk dims [llength $froms]
		h5mdfile H5Dcreate2 "parameters/vmd_structure/bond_from"
		h5mdfile H5Dopen2 "parameters/vmd_structure/bond_from"
		for {set i 0 } {$i < [llength $froms] } {incr i} {
			h5mdfile H5_write_value value [lindex $froms $i] index $i
		}
		h5mdfile H5Dwrite

		h5mdfile H5Screate_simple type int dims [llength $tos]
		h5mdfile H5Pset_chunk dims [llength $tos]
		h5mdfile H5Dcreate2 "parameters/vmd_structure/bond_to"
		h5mdfile H5Dopen2 "parameters/vmd_structure/bond_to"
		for {set i 0 } {$i < [llength $tos] } {incr i} {
			h5mdfile H5_write_value value [lindex $tos $i] index $i
		}
		h5mdfile H5Dwrite
	}
		
	
}

proc h5md_write_positions { args } {
	h5mdfile H5Dopen2 "particles/atoms/position/value"
	set offset [h5mdfile get_dataset_dims]
	# Write positions of all particles
	h5mdfile H5Dextend dims [expr [lindex $offset 0]+1 ] [setmd n_part] 3
	h5mdfile H5Sselect_hyperslab offset [expr [lindex $offset 0]] 0 0
	h5mdfile H5Screate_simple type double dims 1 [setmd n_part] 3
	for { set i 0 } { $i < [setmd n_part] } { incr i } {
		if { [lindex $args 0 ] == "folded" } {
			set pos [part $i print pos folded]
		} else {
			set pos [part $i print pos]
		}
		h5mdfile H5_write_value value [lindex $pos 0] index 0 $i 0		
		h5mdfile H5_write_value value [lindex $pos 1] index 0 $i 1		
		h5mdfile H5_write_value value [lindex $pos 2] index 0 $i 2
	}
	h5mdfile H5Dwrite
	# Write simulation step (assumes that time_step hasnt changed)
	h5mdfile H5Dopen2 "particles/atoms/position/step"
	set offset [h5mdfile get_dataset_dims]
	h5mdfile H5Dextend dims [expr [lindex $offset 0]+1] 
	h5mdfile H5Sselect_hyperslab offset [lindex $offset 0] 
	h5mdfile H5Screate_simple type int dims 1 
	h5mdfile H5_write_value value [expr int([setmd time]/[setmd time_step])] index 0 
	h5mdfile H5Dwrite
	# Write simulation time
	h5mdfile H5Dopen2 "particles/atoms/position/time"
	set offset [h5mdfile get_dataset_dims]
	h5mdfile H5Dextend dims [expr [lindex $offset 0]+1] 
	h5mdfile H5Sselect_hyperslab offset [lindex $offset 0] 
	h5mdfile H5Screate_simple type double dims 1 
	h5mdfile H5_write_value value [setmd time] index 0 
	h5mdfile H5Dwrite
	h5mdfile H5_Fflush
}


proc h5md_write_velocities {} {
	h5mdfile H5Dopen2 "particles/atoms/velocity/value"
	set offset [h5mdfile get_dataset_dims]
	# Write velocities of all particles
	h5mdfile H5Dextend dims [expr [lindex $offset 0]+1] [setmd n_part] 3
	h5mdfile H5Sselect_hyperslab offset [expr [lindex $offset 0]] 0 0
	h5mdfile H5Screate_simple type double dims 1 [setmd n_part] 3
	for { set i 0 } { $i < [setmd n_part] } { incr i } {
		set vel [part $i print v]
		h5mdfile H5_write_value value [lindex $vel 0] index 0 $i 0
		h5mdfile H5_write_value value [lindex $vel 1] index 0 $i 1
		h5mdfile H5_write_value value [lindex $vel 2] index 0 $i 2
	}
	h5mdfile H5Dwrite
	# Write simulation step (assumes that time_step hasnt changed)
	h5mdfile H5Dopen2 "particles/atoms/velocity/step"
	set offset [h5mdfile get_dataset_dims]
	h5mdfile H5Dextend dims [expr [lindex $offset 0]+1] 
	h5mdfile H5Sselect_hyperslab offset [lindex $offset 0] 
	h5mdfile H5Screate_simple type int dims 1 
	h5mdfile H5_write_value value [expr int([setmd time]/[setmd time_step])] index 0 
	h5mdfile H5Dwrite
	# Write simulation time
	h5mdfile H5Dopen2 "particles/atoms/velocity/time"
	set offset [h5mdfile get_dataset_dims]
	h5mdfile H5Dextend dims [expr [lindex $offset 0]+1] 
	h5mdfile H5Sselect_hyperslab offset [lindex $offset 0] 
	h5mdfile H5Screate_simple type double dims 1 
	h5mdfile H5_write_value value [setmd time] index 0 
	h5mdfile H5Dwrite
	h5mdfile H5_Fflush
}


proc h5md_write_forces {} {
	h5mdfile H5Dopen2 "particles/atoms/force/value"
	set offset [h5mdfile get_dataset_dims]
	# Write forces of all particles
	h5mdfile H5Dextend dims [expr [lindex $offset 0]+1] [setmd n_part] 3
	h5mdfile H5Sselect_hyperslab offset [expr [lindex $offset 0]] 0 0
	h5mdfile H5Screate_simple type double dims 1 [setmd n_part] 3
	for { set i 0 } { $i < [setmd n_part] } { incr i } {
		set force [part $i print f]
		h5mdfile H5_write_value value [lindex $force 0] index 0 $i 0
		h5mdfile H5_write_value value [lindex $force 1] index 0 $i 1
		h5mdfile H5_write_value value [lindex $force 2] index 0 $i 2
	}
	h5mdfile H5Dwrite
	# Write simulation step (assumes that time_step hasnt changed)
	h5mdfile H5Dopen2 "particles/atoms/force/step"
	set offset [h5mdfile get_dataset_dims]
	h5mdfile H5Dextend dims [expr [lindex $offset 0]+1] 
	h5mdfile H5Sselect_hyperslab offset [lindex $offset 0] 
	h5mdfile H5Screate_simple type int dims 1 
	h5mdfile H5_write_value value [expr int([setmd time]/[setmd time_step])] index 0 
	h5mdfile H5Dwrite
	# Write simulation time
	h5mdfile H5Dopen2 "particles/atoms/force/time"
	set offset [h5mdfile get_dataset_dims]
	h5mdfile H5Dextend dims [expr [lindex $offset 0]+1] 
	h5mdfile H5Sselect_hyperslab offset [lindex $offset 0] 
	h5mdfile H5Screate_simple type double dims 1 
	h5mdfile H5_write_value value [setmd time] index 0 
	h5mdfile H5Dwrite
	h5mdfile H5_Fflush
}

proc h5md_write_species {} {
	# Write particle type of all particles in the system
	h5mdfile H5Dopen2 "particles/atoms/species"
	h5mdfile H5Screate_simple type double dims [setmd n_part] 
	for { set i 0 } { $i < [setmd n_part] } { incr i } {
		h5mdfile H5_write_value value [part $i print type] index $i 
	}
	h5mdfile H5Dwrite
	h5mdfile H5_Fflush
}

# Initialize a user defined one dimensional observable, first argument is name of observable
proc h5md_observable1D_init { args } {
	h5mdfile H5Gcreate2 "observables/[lindex $args 0]"
	h5mdfile H5Screate_simple type int dims 0 
	h5mdfile H5Pset_chunk dims 1 
	h5mdfile H5Dcreate2 "observables/[lindex $args 0]/step"
	h5mdfile H5Screate_simple type double dims 0 
	h5mdfile H5Pset_chunk dims 1 
	h5mdfile H5Dcreate2 "observables/[lindex $args 0]/time"
	h5mdfile H5Screate_simple type double dims 0 
	h5mdfile H5Pset_chunk dims 1
	h5mdfile H5Dcreate2 "observables/[lindex $args 0]/value"
}

# Writes to a user defined zero dimensional but timedependent observable dataset
proc h5md_observable1D_write { args } {
	# Write simulation step (assumes that time_step hasnt changed)
	h5mdfile H5Dopen2 "observables/[lindex $args 0]/step"
	set offset [h5mdfile get_dataset_dims]
	h5mdfile H5Dextend dims [expr [lindex $offset 0]+1] 
	h5mdfile H5Sselect_hyperslab offset [lindex $offset 0] 
	h5mdfile H5Screate_simple type int dims 1 
	h5mdfile H5_write_value value [expr int([setmd time]/[setmd time_step])] index 0 
	h5mdfile H5Dwrite
	# Write simulation time
	h5mdfile H5Dopen2 "observables/[lindex $args 0]/time"
	set offset [h5mdfile get_dataset_dims]
	h5mdfile H5Dextend dims [expr [lindex $offset 0]+1] 
	h5mdfile H5Sselect_hyperslab offset [lindex $offset 0] 
	h5mdfile H5Screate_simple type double dims 1 
	h5mdfile H5_write_value value [setmd time] index 0 
	h5mdfile H5Dwrite
	h5mdfile H5Dopen2 "observables/[lindex $args 0]/value"
	set offset [h5mdfile get_dataset_dims]
	h5mdfile H5Dextend dims [expr [lindex $offset 0]+1]
	h5mdfile H5Sselect_hyperslab offset [lindex $offset 0]
	h5mdfile H5Screate_simple type double dims 1
	h5mdfile H5_write_value value [lindex $args 1] index 0
	h5mdfile H5Dwrite
	h5mdfile H5_Fflush
}

# Writes to a user defined 1 dimensional but timedependent observable dataset
proc h5md_observable2D_init_new { data_path name p_ids } {
	set num_particles [llength $p_ids]
	if { [setmd n_part] == 0 || $num_particles <1 } {
		puts "Please set your particles before h5md initialisation\n"; flush stdout
		return
	}
	if { [file exists ${data_path}] == 1 } {
		h5mdfile H5Fopen "${data_path}"
	} elseif { [file exists ${data_path}] == 0 } {
		# Create hdf5 file
		h5mdfile H5Fcreate "${data_path}"
	}

	h5mdfile H5Gcreate2 "observables" ; #XXX only if it not already exists
	h5mdfile H5Gcreate2 "observables/$name"
	
        h5mdfile H5Screate_simple type int dims 0 
	h5mdfile H5Pset_chunk dims 1 
	h5mdfile H5Dcreate2 "observables/$name/step"
	h5mdfile H5Screate_simple type double dims 0 
	h5mdfile H5Pset_chunk dims 1 
	h5mdfile H5Dcreate2 "observables/$name/time"
	h5mdfile H5Screate_simple type double dims 0 $num_particles 3
	h5mdfile H5Pset_chunk dims 5 $num_particles 3
	h5mdfile H5Dcreate2 "observables/$name/value"

}

proc h5md_observable2D_write_new { name p_ids } {
	h5mdfile H5Dopen2 "observables/$name/value"
	set offset [h5mdfile get_dataset_dims]
	puts $offset
	# Write observable of all particles
	set num_particles [llength $p_ids]
	h5mdfile H5Dextend dims [expr [lindex $offset 0]+1] $num_particles 3
	h5mdfile H5Sselect_hyperslab offset [expr [lindex $offset 0]] 0 0
	h5mdfile H5Screate_simple type double dims 1 $num_particles 3
	for { set j 0 } { $j < $num_particles } { incr j } {
		set p_id [lindex $p_ids $j]
		set property [part $p_id print $name]
		h5mdfile H5_write_value value [lindex $property 0] index 0 $j 0
		h5mdfile H5_write_value value [lindex $property 1] index 0 $j 1
		h5mdfile H5_write_value value [lindex $property 2] index 0 $j 2
	}
	h5mdfile H5Dwrite
	# Write simulation step (assumes that time_step hasnt changed)
	h5mdfile H5Dopen2 "observables/$name/step"
	set offset [h5mdfile get_dataset_dims]
	h5mdfile H5Dextend dims [expr [lindex $offset 0]+1] 
	h5mdfile H5Sselect_hyperslab offset [lindex $offset 0] 
	h5mdfile H5Screate_simple type int dims 1 
	h5mdfile H5_write_value value [expr int([setmd time]/[setmd time_step])] index 0 
	h5mdfile H5Dwrite
	# Write simulation time
	h5mdfile H5Dopen2 "observables/$name/time"
	set offset [h5mdfile get_dataset_dims]
	h5mdfile H5Dextend dims [expr [lindex $offset 0]+1] 
	h5mdfile H5Sselect_hyperslab offset [lindex $offset 0] 
	h5mdfile H5Screate_simple type double dims 1 
	h5mdfile H5_write_value value [setmd time] index 0 
	h5mdfile H5Dwrite
	h5mdfile H5_Fflush

}



# Close all h5md groups and datasets and free memory at the end
proc h5md_close {} {
	h5mdfile H5Dclose
	h5mdfile H5Fclose
	h5mdfile H5_free_memory
}

proc h5mdutil_get_types {} {
	set type_list ""
	for { set i 0 } { $i < [setmd n_part] } { incr i } {
		if { [expr [lsearch $type_list [part $i print type]] == -1] } {
			lappend type_list [part $i print type]
		}
	}
	return $type_list
}
