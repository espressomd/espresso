proc h5md_init { data_path } {
	if { [file exists ${data_path}] == 1 } {
		h5md_open ${data_path}
	} else {
		# create hdf5 file
		h5mdfile H5Fcreate "${data_path}"
	}
	if { [setmd n_part] == 0 } {
		puts "Please set your particles before h5md initialisation\n"; flush stdout
		return
	}
	# create data groups
	h5mdfile H5Gcreate2 "particles"
	h5mdfile H5Gcreate2 "particles/atoms"
	h5mdfile H5Gcreate2 "particles/atoms/box"
	h5mdfile H5Gcreate2 "particles/atoms/mass"
	h5mdfile H5Gcreate2 "particles/atoms/position"
	h5mdfile H5Gcreate2 "particles/atoms/velocity"
	h5mdfile H5Gcreate2 "parameters"
	h5mdfile H5Gcreate2 "parameters/vmd_structure"
	h5mdfile H5Gcreate2 "observables" 
	h5mdfile H5Gcreate2 "files"
	h5mdfile H5Gcreate2 "files/scripts"
	# create datasets
	h5mdfile H5Screate_simple type double dims 3 3
	h5mdfile H5Pset_chunk dims 3 3
	h5mdfile H5Dcreate2 "particles/atoms/box/edges"
	h5mdfile H5Screate_simple type double dims [setmd n_part] 
	h5mdfile H5Pset_chunk dims [setmd n_part] 
	h5mdfile H5Dcreate2 "particles/atoms/mass/value"
	h5mdfile H5Dopen2 "particles/atoms/mass/value"
	if { [string first MASS [code_info]] != -1 } {
		for  { set i 0 } { $i < [setmd n_part] } { incr i } {
			h5mdfile H5_write_value value [part $i print mass]
		}
		h5mdfile H5Dwrite
	}
	h5mdfile H5Screate_simple type int dims 0 
	h5mdfile H5Pset_chunk dims 1 
	h5mdfile H5Dcreate2 "particles/atoms/position/step"
	h5mdfile H5Screate_simple type double dims 0 
	h5mdfile H5Pset_chunk dims 1 
	h5mdfile H5Dcreate2 "particles/atoms/position/time"
	h5mdfile H5Screate_simple type double dims 0 [setmd n_part] 3
	h5mdfile H5Pset_chunk dims 5 [setmd n_part] 3
	h5mdfile H5Dcreate2 "particles/atoms/position/value"
	h5mdfile H5Screate_simple type int dims [setmd n_part] 
	h5mdfile H5Pset_chunk dims [setmd n_part] 
	h5mdfile H5Dcreate2 "particles/atoms/species"
	h5mdfile H5Screate_simple type int dims 0 
	h5mdfile H5Pset_chunk dims 1 
	h5mdfile H5Dcreate2 "particles/atoms/velocity/step"
	h5mdfile H5Screate_simple type double dims 0 
	h5mdfile H5Pset_chunk dims 1 
	h5mdfile H5Dcreate2 "particles/atoms/velocity/time"
	h5mdfile H5Screate_simple type double dims 0 [setmd n_part] 3
	h5mdfile H5Pset_chunk dims 5 [setmd n_part] 3
	h5mdfile H5Dcreate2 "particles/atoms/velocity/value"
	h5mdfile H5Screate_simple type int dims [setmd n_part_types]
	h5mdfile H5Pset_chunk dims [setmd n_part_types]
	h5mdfile H5Dcreate2 "parameters/vmd_structure/indexOfSpecies"
	h5mdfile H5Dopen2 "parameters/vmd_structure/indexOfSpecies"
	for { set i 0 } { $i < [setmd n_part_types] } { incr i } {
		h5mdfile H5_write_value value [expr [lindex [h5mdutil_get_types] $i]] index $i
	}
	h5mdfile H5Dwrite
	h5md_write_species
}

proc h5md_write_positions { args } {
	h5mdfile H5Dopen2 "particles/atoms/position/value"
	set offset [h5mdfile get_dataset_dims]
	# write positions of all particles
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
	# write simulation step (assumes that time_step hasnt changed)
	h5mdfile H5Dopen2 "particles/atoms/position/step"
	set offset [h5mdfile get_dataset_dims]
	h5mdfile H5Dextend dims [expr [lindex $offset 0]+1] 
	h5mdfile H5Sselect_hyperslab offset [lindex $offset 0] 
	h5mdfile H5Screate_simple type int dims 1 
	h5mdfile H5_write_value value [expr int([setmd time]/[setmd time_step])] index 0 
	h5mdfile H5Dwrite
	# write simulation time
	h5mdfile H5Dopen2 "particles/atoms/position/time"
	set offset [h5mdfile get_dataset_dims]
	h5mdfile H5Dextend dims [expr [lindex $offset 0]+1] 
	h5mdfile H5Sselect_hyperslab offset [lindex $offset 0] 
	h5mdfile H5Screate_simple type double dims 1 
	h5mdfile H5_write_value value [setmd time] index 0 
	h5mdfile H5Dwrite
}

proc h5md_write_velocities {} {
	h5mdfile H5Dopen2 "particles/atoms/velocity/value"
	set offset [h5mdfile get_dataset_dims]
	# write velocities of all particles
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
	# write simulation step (assumes that time_step hasnt changed)
	h5mdfile H5Dopen2 "particles/atoms/velocity/step"
	set offset [h5mdfile get_dataset_dims]
	h5mdfile H5Dextend dims [expr [lindex $offset 0]+1] 
	h5mdfile H5Sselect_hyperslab offset [lindex $offset 0] 
	h5mdfile H5Screate_simple type int dims 1 
	h5mdfile H5_write_value value [expr int([setmd time]/[setmd time_step])] index 0 
	h5mdfile H5Dwrite
	# write simulation time
	h5mdfile H5Dopen2 "particles/atoms/velocity/time"
	set offset [h5mdfile get_dataset_dims]
	h5mdfile H5Dextend dims [expr [lindex $offset 0]+1] 
	h5mdfile H5Sselect_hyperslab offset [lindex $offset 0] 
	h5mdfile H5Screate_simple type double dims 1 
	h5mdfile H5_write_value value [setmd time] index 0 
	h5mdfile H5Dwrite
}

proc h5md_write_species {} {
	# write particle type of all particles in the system
	h5mdfile H5Dopen2 "particles/atoms/species"
	h5mdfile H5Screate_simple type double dims [setmd n_part] 
	for { set i 0 } { $i < [setmd n_part] } { incr i } {
		h5mdfile H5_write_value value [part $i print type] index $i 
	}
	h5mdfile H5Dwrite
}

# initialize a user defined zero dimensional observable, first argument is name of observable
proc h5md_observable0D_init { args } {
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

# writes to a user defined zero dimensional observable dataset
proc h5md_write_observable0D { args } {
	# write simulation step (assumes that time_step hasnt changed)
	h5mdfile H5Dopen2 "observables/[lindex $args 0]/step"
	set offset [h5mdfile get_dataset_dims]
	h5mdfile H5Dextend dims [expr [lindex $offset 0]+1] 
	h5mdfile H5Sselect_hyperslab offset [lindex $offset 0] 
	h5mdfile H5Screate_simple type int dims 1 
	h5mdfile H5_write_value value [expr int([setmd time]/[setmd time_step])] index 0 
	h5mdfile H5Dwrite
	# write simulation time
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
}

proc h5md_open { args } {
	h5mdfile H5Fopen "[lindex $args 0]"
}

# close all h5md groups and datasets and free memory at the end
proc h5md_close {} {
	#h5mdfile H5Pclose
	h5mdfile H5Dclose
	#h5mdfile H5Sclose
	#h5mdfile H5Gclose
	h5mdfile H5Fclose
	h5mdfile H5_free_memory
}

proc h5mdutil_get_types {} {
	set type_list ""
	for { set i 0 } { $i < [setmd n_part] } { incr i } {
		if { [expr [lsearch $type_list [part $i print type]] == -1] } {
			puts "adding [part $i print type]"; flush stdout
			lappend type_list [part $i print type]
		}
	}
	return $type_list
}
