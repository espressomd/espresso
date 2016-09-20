set h5md_num_part 0
set h5md_p_ids ""

proc h5md_init { data_path {_h5md_p_ids ""} } {
    global h5md_num_part h5md_p_ids
    global __h5md_filename
    global __h5md_filename_tmp
    set __h5md_filename $data_path
    set h5md_p_ids ""
    if { $_h5md_p_ids eq ""} {
        for {set i 0} {$i<[setmd n_part]} {incr i} {
            lappend h5md_p_ids $i
        }
    } else {
        set h5md_p_ids $_h5md_p_ids
    }
    set h5md_num_part [llength $h5md_p_ids] 
    if { [setmd n_part] == 0 } {
        puts "Please set your particles before h5md initialisation\n"; flush stdout
        return
    }
    set datadir [file dirname $data_path]
    set filenamebasepath [file rootname $data_path]
    set filenamebase [file tail $filenamebasepath]
    set datadirfiles [glob -type f -nocomplain -dir $datadir *]
    if { [expr [file exists "${filenamebasepath}.h5"] && [file exists "${filenamebasepath}.h5_tmp"]] } {
        puts "It seems like either your simulation crashed before or you didnt call the close function."
        exit 111
    } elseif { [file exists "${filenamebasepath}.h5"] } {
        puts "Found existing h5 file matching given pattern."; flush stdout
        set __h5md_filename_tmp "${filenamebasepath}.h5_tmp"
        file copy "${__h5md_filename}" "${__h5md_filename_tmp}"
        h5mdfile H5Fopen "${__h5md_filename_tmp}"
        return
    } else {
        set __h5md_filename_tmp "${filenamebasepath}.h5_tmp"
        # Create hdf5 file
        h5mdfile H5Fcreate "${__h5md_filename_tmp}"
    }
    # Create data groups
    h5mdfile H5Gcreate2 "particles"
    h5mdfile H5Gclose
    h5mdfile H5Gcreate2 "particles/atoms"
    h5mdfile H5Gclose
    h5mdfile H5Gcreate2 "particles/atoms/box"
    h5mdfile H5Gclose
    h5mdfile H5Gcreate2 "particles/atoms/mass"
    h5mdfile H5Gclose
    h5mdfile H5Gcreate2 "particles/atoms/position"
    h5mdfile H5Gclose
    h5mdfile H5Gcreate2 "particles/atoms/velocity"
    h5mdfile H5Gclose
    h5mdfile H5Gcreate2 "particles/atoms/force"
    h5mdfile H5Gclose
    h5mdfile H5Gcreate2 "parameters"
    h5mdfile H5Gclose
    h5mdfile H5Gcreate2 "parameters/vmd_structure"
    h5mdfile H5Gclose
    h5mdfile H5Gcreate2 "parameters/files"
    h5mdfile H5Gclose
    h5mdfile H5Gcreate2 "observables" 
    h5mdfile H5Gclose
    # Create datasets
    #box
    h5mdfile H5Screate_simple type double dims 3 3
    h5mdfile H5Pset_chunk dims 3 3
    h5mdfile H5Dcreate2 "particles/atoms/box/edges"
    h5mdfile H5Dopen2 "particles/atoms/box/edges"
    ##assuming cuboid box
    for {set i 0 } {$i<3} {incr i} {
        h5mdfile H5_write_value value [lindex [setmd box_l] $i] index $i $i
    }
    h5mdfile H5Dwrite
    h5mdfile H5Dclose
    h5mdfile H5Sclose_simple
    #mass   
    h5mdfile H5Screate_simple type double dims $h5md_num_part
    h5mdfile H5Pset_chunk dims $h5md_num_part
    h5mdfile H5Dcreate2 "particles/atoms/mass/value"
    h5mdfile H5Dopen2 "particles/atoms/mass/value"
    if { [string first MASS [code_info]] != -1 } {
        set i 0
        foreach p_id $h5md_p_ids {
            h5mdfile H5_write_value value [part $p_id print mass] index $i
            incr i
        }
    }
    h5mdfile H5Dwrite
    h5mdfile H5Sclose_simple
    h5mdfile H5Dclose
    ## position
    h5mdfile H5Screate_simple type int dims 0 
    h5mdfile H5Pset_chunk dims 5 
    h5mdfile H5Dcreate2 "particles/atoms/position/step"
    h5mdfile H5Sclose_simple
    h5mdfile H5_Tclose
    h5mdfile H5Screate_simple type double dims 0 
    h5mdfile H5Pset_chunk dims 5 
    h5mdfile H5Dcreate2 "particles/atoms/position/time"
    h5mdfile H5Sclose_simple
    h5mdfile H5_Tclose
    h5mdfile H5Screate_simple type double dims 0 $h5md_num_part 3
    h5mdfile H5Pset_chunk dims 5 $h5md_num_part 3
    h5mdfile H5Dcreate2 "particles/atoms/position/value"
    h5mdfile H5Sclose_simple
    h5mdfile H5_Tclose
    ## velocity
    h5mdfile H5Screate_simple type int dims 0 
    h5mdfile H5Pset_chunk dims 5 
    h5mdfile H5Dcreate2 "particles/atoms/velocity/step"
    h5mdfile H5Sclose_simple
    h5mdfile H5_Tclose
    h5mdfile H5Screate_simple type double dims 0 
    h5mdfile H5Pset_chunk dims 5 
    h5mdfile H5Dcreate2 "particles/atoms/velocity/time"
    h5mdfile H5Sclose_simple
    h5mdfile H5_Tclose
    h5mdfile H5Screate_simple type double dims 0 $h5md_num_part 3
    h5mdfile H5Pset_chunk dims 5 $h5md_num_part 3
    h5mdfile H5Dcreate2 "particles/atoms/velocity/value"
    h5mdfile H5Sclose_simple
    h5mdfile H5_Tclose
    ## force
    h5mdfile H5Screate_simple type int dims 0 
    h5mdfile H5Pset_chunk dims 1 
    h5mdfile H5Dcreate2 "particles/atoms/force/step"
    h5mdfile H5Sclose_simple
    h5mdfile H5_Tclose
    h5mdfile H5Screate_simple type double dims 0 
    h5mdfile H5Pset_chunk dims 1 
    h5mdfile H5Dcreate2 "particles/atoms/force/time"
    h5mdfile H5Sclose_simple
    h5mdfile H5_Tclose
    h5mdfile H5Screate_simple type double dims 0 $h5md_num_part 3
    h5mdfile H5Pset_chunk dims 5 $h5md_num_part 3
    h5mdfile H5Dcreate2 "particles/atoms/force/value"
    h5mdfile H5Sclose_simple
    h5mdfile H5_Tclose
    ## species
    h5mdfile H5Screate_simple type int dims $h5md_num_part 
    h5mdfile H5Pset_chunk dims $h5md_num_part
    h5mdfile H5Dcreate2 "particles/atoms/species"
    h5mdfile H5Dopen2 "particles/atoms/species"
    set i 0
    foreach p_id $h5md_p_ids {
        h5mdfile H5_write_value value [part $p_id print type] index $i 
        incr i
    }
    h5mdfile H5Dwrite
    h5mdfile H5Sclose_simple
    h5mdfile H5Dclose
    ## Tcl script
    dump_script_to_h5md 
    ## VMD visualization
    set types [h5mdutil_get_types]
    set names "X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg"; # for coloring in VMD only
    ##index of species
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
    h5mdfile H5Sclose_simple
    h5mdfile H5Dclose
    ##name
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
    h5mdfile H5Sclose_simple
    h5mdfile H5Dclose
    ##type
    h5mdfile H5Screate_simple type str dims [llength $types]
    h5mdfile H5Pset_chunk dims [llength $types]
    h5mdfile H5Dcreate2 "parameters/vmd_structure/type"
    h5mdfile H5Dopen2 "parameters/vmd_structure/type"
    for {set i 0 } {$i < [llength $types] } {incr i} {
        h5mdfile H5_write_value value [lindex $names [expr $i%112]] index $i
    }
    h5mdfile H5Dwrite
    h5mdfile H5Sclose_simple
    h5mdfile H5Dclose
    ##write charge for vmd_structure
    h5mdfile H5Screate_simple type double dims $h5md_num_part
    h5mdfile H5Pset_chunk dims $h5md_num_part
    h5mdfile H5Dcreate2 "particles/atoms/charge"
    h5mdfile H5Dopen2 "particles/atoms/charge"
    set i 0
    foreach p_id $h5md_p_ids {
        h5mdfile H5_write_value value [part $p_id pr q] index $i
        incr i
    }
    h5mdfile H5Dwrite
    h5mdfile H5Sclose_simple
    h5mdfile H5Dclose
    ##bonds
    set froms ""
    set tos ""
    foreach from $h5md_p_ids {
        if { [part $from] != "na" } then {
            set bonds [lindex [part $from print bond] 0]
            if { [llength $bonds] >0 } {
                for { set i 0 } { $i < [llength $bonds] } { incr i } {
                set to [lindex $bonds $i 1]
                if { [expr [lsearch $h5md_p_ids $to] != -1] } {
                    lappend froms [expr $from +1]; #VMD counts from 1 and not from 0 therefore +1
                    lappend tos [expr $to +1]
                }
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
        h5mdfile H5Sclose_simple
        h5mdfile H5Dclose

        h5mdfile H5Screate_simple type int dims [llength $tos]
        h5mdfile H5Pset_chunk dims [llength $tos]
        h5mdfile H5Dcreate2 "parameters/vmd_structure/bond_to"
        h5mdfile H5Dopen2 "parameters/vmd_structure/bond_to"
        for {set i 0 } {$i < [llength $tos] } {incr i} {
            h5mdfile H5_write_value value [lindex $tos $i] index $i
        }
        h5mdfile H5Dwrite
        h5mdfile H5Sclose_simple
        h5mdfile H5Dclose
    }
    h5mdfile H5_Fflush
}

proc h5md_write_positions { args } {
    global h5md_num_part h5md_p_ids
    h5mdfile H5Dopen2 "particles/atoms/position/value"
    set offset [h5mdfile get_dataset_dims]
    # Write positions of all particles
    h5mdfile H5Dextend dims [expr [lindex $offset 0]+1 ] $h5md_num_part 3
    h5mdfile H5Sselect_hyperslab offset [expr [lindex $offset 0]] 0 0
    h5mdfile H5Screate_simple type double dims 1 $h5md_num_part 3
    set i 0
    foreach p_id $h5md_p_ids {
        if { [lindex $args 0 ] == "folded" } {
            set pos [part $p_id print pos folded]
        } else {
            set pos [part $p_id print pos]
        }
        h5mdfile H5_write_value value [lindex $pos 0] index 0 $i 0      
        h5mdfile H5_write_value value [lindex $pos 1] index 0 $i 1      
        h5mdfile H5_write_value value [lindex $pos 2] index 0 $i 2
        incr i
    }
    h5mdfile H5Dwrite
    h5mdfile H5Dclose
    h5mdfile H5Sclose_simple
    ### Write simulation step (assumes that time_step hasnt changed)
    h5mdfile H5Dopen2 "particles/atoms/position/step"
    set offset [h5mdfile get_dataset_dims]
    h5mdfile H5Dextend dims [expr [lindex $offset 0]+1] 
    h5mdfile H5Sselect_hyperslab offset [lindex $offset 0] 
    h5mdfile H5Screate_simple type int dims 1 
    h5mdfile H5_write_value value [expr int([setmd time]/[setmd time_step])] index 0 
    h5mdfile H5Dwrite
    h5mdfile H5Dclose
    h5mdfile H5Sclose_simple
    # Write simulation time
    h5mdfile H5Dopen2 "particles/atoms/position/time"
    set offset [h5mdfile get_dataset_dims]
    h5mdfile H5Dextend dims [expr [lindex $offset 0]+1] 
    h5mdfile H5Sselect_hyperslab offset [lindex $offset 0] 
    h5mdfile H5Screate_simple type double dims 1 
    h5mdfile H5_write_value value [setmd time] index 0 
    h5mdfile H5Dwrite
    h5mdfile H5Dclose
    h5mdfile H5Sclose_simple
}


proc h5md_write_velocities {} {
    global h5md_num_part h5md_p_ids
    h5mdfile H5Dopen2 "particles/atoms/velocity/value"
    set offset [h5mdfile get_dataset_dims]
    # Write velocities of all particles
    h5mdfile H5Dextend dims [expr [lindex $offset 0]+1] $h5md_num_part 3
    h5mdfile H5Sselect_hyperslab offset [expr [lindex $offset 0]] 0 0
    h5mdfile H5Screate_simple type double dims 1 $h5md_num_part 3
    set i 0
    foreach p_id $h5md_p_ids {
        set vel [part $p_id print v]
        h5mdfile H5_write_value value [lindex $vel 0] index 0 $i 0
        h5mdfile H5_write_value value [lindex $vel 1] index 0 $i 1
        h5mdfile H5_write_value value [lindex $vel 2] index 0 $i 2
        incr i
    }
    h5mdfile H5Dwrite
    h5mdfile H5Dclose
    h5mdfile H5Sclose_simple
    # Write simulation step (assumes that time_step hasnt changed)
    h5mdfile H5Dopen2 "particles/atoms/velocity/step"
    set offset [h5mdfile get_dataset_dims]
    h5mdfile H5Dextend dims [expr [lindex $offset 0]+1] 
    h5mdfile H5Sselect_hyperslab offset [lindex $offset 0] 
    h5mdfile H5Screate_simple type int dims 1 
    h5mdfile H5_write_value value [expr int([setmd time]/[setmd time_step])] index 0 
    h5mdfile H5Dwrite
    h5mdfile H5Dclose
    h5mdfile H5Sclose_simple
    # Write simulation time
    h5mdfile H5Dopen2 "particles/atoms/velocity/time"
    set offset [h5mdfile get_dataset_dims]
    h5mdfile H5Dextend dims [expr [lindex $offset 0]+1] 
    h5mdfile H5Sselect_hyperslab offset [lindex $offset 0] 
    h5mdfile H5Screate_simple type double dims 1 
    h5mdfile H5_write_value value [setmd time] index 0 
    h5mdfile H5Dwrite
    h5mdfile H5Dclose
    h5mdfile H5Sclose_simple
}


proc h5md_write_forces {} {
    global h5md_num_part h5md_p_ids
    h5mdfile H5Dopen2 "particles/atoms/force/value"
    set offset [h5mdfile get_dataset_dims]
    # Write forces of all particles
    h5mdfile H5Dextend dims [expr [lindex $offset 0]+1] $h5md_num_part 3
    h5mdfile H5Sselect_hyperslab offset [expr [lindex $offset 0]] 0 0
    h5mdfile H5Screate_simple type double dims 1 $h5md_num_part 3
    set i 0
    foreach p_id $h5md_p_ids {
        set force [part $p_id print f]
        h5mdfile H5_write_value value [lindex $force 0] index 0 $i 0
        h5mdfile H5_write_value value [lindex $force 1] index 0 $i 1
        h5mdfile H5_write_value value [lindex $force 2] index 0 $i 2
        incr i
    }
    h5mdfile H5Dwrite
    h5mdfile H5_free_memory
    h5mdfile H5Dclose
    h5mdfile H5Sclose_simple
    # Write simulation step (assumes that time_step hasnt changed)
    h5mdfile H5Dopen2 "particles/atoms/force/step"
    set offset [h5mdfile get_dataset_dims]
    h5mdfile H5Dextend dims [expr [lindex $offset 0]+1] 
    h5mdfile H5Sselect_hyperslab offset [lindex $offset 0] 
    h5mdfile H5Screate_simple type int dims 1 
    h5mdfile H5_write_value value [expr int([setmd time]/[setmd time_step])] index 0 
    h5mdfile H5Dwrite
    h5mdfile H5_free_memory
    h5mdfile H5Dclose
    h5mdfile H5Sclose_simple
    # Write simulation time
    h5mdfile H5Dopen2 "particles/atoms/force/time"
    set offset [h5mdfile get_dataset_dims]
    h5mdfile H5Dextend dims [expr [lindex $offset 0]+1] 
    h5mdfile H5Sselect_hyperslab offset [lindex $offset 0] 
    h5mdfile H5Screate_simple type double dims 1 
    h5mdfile H5_write_value value [setmd time] index 0 
    h5mdfile H5Dwrite
    h5mdfile H5_free_memory
    h5mdfile H5Dclose
    h5mdfile H5Sclose_simple
}


# Initialize a user defined one dimensional observable, first argument is name of observable
proc h5md_observable1D_init { args } {
    h5mdfile H5Gcreate2 "observables/[lindex $args 0]"
    h5mdfile H5Screate_simple type int dims 0 
    h5mdfile H5Pset_chunk dims 1 
    h5mdfile H5Dcreate2 "observables/[lindex $args 0]/step"
    h5mdfile H5Dclose
    h5mdfile H5Sclose_simple
    h5mdfile H5Screate_simple type double dims 0 
    h5mdfile H5Pset_chunk dims 1 
    h5mdfile H5Dcreate2 "observables/[lindex $args 0]/time"
    h5mdfile H5Dclose
    h5mdfile H5Sclose_simple
    h5mdfile H5Screate_simple type double dims 0 
    h5mdfile H5Pset_chunk dims 1
    h5mdfile H5Dcreate2 "observables/[lindex $args 0]/value"
    h5mdfile H5Dclose
    h5mdfile H5Sclose_simple
}

# Writes to a user defined zero dimensional but timedependent observable dataset
proc h5md_observable1D_write { args } {
    # Write simulation step (assumes that time_step has not changed)
    h5mdfile H5Dopen2 "observables/[lindex $args 0]/step"
    set offset [h5mdfile get_dataset_dims]
    h5mdfile H5Dextend dims [expr [lindex $offset 0]+1] 
    h5mdfile H5Sselect_hyperslab offset [lindex $offset 0] 
    h5mdfile H5Screate_simple type int dims 1 
    h5mdfile H5_write_value value [expr int([setmd time]/[setmd time_step])] index 0 
    h5mdfile H5Dwrite
    h5mdfile H5_free_memory
    h5mdfile H5Dclose
    h5mdfile H5Sclose_simple
    # Write simulation time
    h5mdfile H5Dopen2 "observables/[lindex $args 0]/time"
    set offset [h5mdfile get_dataset_dims]
    h5mdfile H5Dextend dims [expr [lindex $offset 0]+1] 
    h5mdfile H5Sselect_hyperslab offset [lindex $offset 0]
    h5mdfile H5Screate_simple type double dims 1 
    h5mdfile H5_write_value value [setmd time] index 0 
    h5mdfile H5Dwrite
    h5mdfile H5_free_memory
    h5mdfile H5Dclose
    h5mdfile H5Sclose_simple
    # Write actual data
    h5mdfile H5Dopen2 "observables/[lindex $args 0]/value"
    set offset [h5mdfile get_dataset_dims]
    h5mdfile H5Dextend dims [expr [lindex $offset 0]+1]
    h5mdfile H5Sselect_hyperslab offset [lindex $offset 0]
    h5mdfile H5Screate_simple type double dims 1
    h5mdfile H5_write_value value [lindex $args 1] index 0
    h5mdfile H5Dwrite
    h5mdfile H5_free_memory
    h5mdfile H5Dclose
    h5mdfile H5Sclose_simple
}


# Writes to a user defined 1 dimensional but timedependent observable dataset
proc h5md_observable2D_init { name } {
    global h5md_num_part h5md_p_ids
    h5mdfile H5Gcreate2 "observables/$name"
    h5mdfile H5Screate_simple type int dims 0 
    h5mdfile H5Pset_chunk dims 1 
    h5mdfile H5Dcreate2 "observables/$name/step"
    h5mdfile H5Dclose
    h5mdfile H5Sclose_simple
    h5mdfile H5Screate_simple type double dims 0 
    h5mdfile H5Pset_chunk dims 1 
    h5mdfile H5Dcreate2 "observables/$name/time"
    h5mdfile H5Dclose
    h5mdfile H5Sclose_simple
    h5mdfile H5Screate_simple type double dims 0 $h5md_num_part
    h5mdfile H5Pset_chunk dims 5 $h5md_num_part
    h5mdfile H5Dcreate2 "observables/$name/value"
    h5mdfile H5Dclose
    h5mdfile H5Sclose_simple
}


proc h5md_observable2D_write { name } {
    global h5md_num_part h5md_p_ids
    h5mdfile H5Dopen2 "observables/$name/value"
    set offset [h5mdfile get_dataset_dims]
    # Write observable of all particles
    h5mdfile H5Dextend dims [expr [lindex $offset 0]+1] $h5md_num_part
    h5mdfile H5Sselect_hyperslab offset [expr [lindex $offset 0]] 0
    h5mdfile H5Screate_simple type double dims 1 $h5md_num_part
    set j 0
    foreach p_id $h5md_p_ids {
        set property [part $p_id print $name]
        h5mdfile H5_write_value value $property index 0 $j
        incr j
    }
    h5mdfile H5Dwrite
    h5mdfile H5_free_memory
    h5mdfile H5Dclose
    h5mdfile H5Sclose_simple
    # Write simulation step (assumes that time_step hasnt changed)
    h5mdfile H5Dopen2 "observables/$name/step"
    set offset [h5mdfile get_dataset_dims]
    h5mdfile H5Dextend dims [expr [lindex $offset 0]+1] 
    h5mdfile H5Sselect_hyperslab offset [lindex $offset 0] 
    h5mdfile H5Screate_simple type int dims 1 
    h5mdfile H5_write_value value [expr int([setmd time]/[setmd time_step])] index 0 
    h5mdfile H5Dwrite
    h5mdfile H5_free_memory
    h5mdfile H5Dclose
    h5mdfile H5Sclose_simple
    # Write simulation time
    h5mdfile H5Dopen2 "observables/$name/time"
    set offset [h5mdfile get_dataset_dims]
    h5mdfile H5Dextend dims [expr [lindex $offset 0]+1] 
    h5mdfile H5Sselect_hyperslab offset [lindex $offset 0] 
    h5mdfile H5Screate_simple type double dims 1 
    h5mdfile H5_write_value value [setmd time] index 0 
    h5mdfile H5Dwrite
    h5mdfile H5_free_memory
    h5mdfile H5Dclose
    h5mdfile H5Sclose_simple
}



# Close all h5md groups and datasets and free memory at the end
proc h5md_close {} {
    global __h5md_filename __h5md_filename_tmp
    h5mdfile H5_free_memory
    h5mdfile H5Fclose
    file delete ${__h5md_filename}
    file copy ${__h5md_filename_tmp} ${__h5md_filename}
    file delete ${__h5md_filename_tmp}
}

proc h5mdutil_get_types {} {
    global h5md_p_ids
    set type_list ""
    foreach p_id $h5md_p_ids {
        if { [expr [lsearch $type_list [part $p_id print type]] == -1] } {
            lappend type_list [part $p_id print type]
        }
    }
    return $type_list
}


proc dump_script_to_h5md {} {
    set script ""
    set fp [open [file normalize [info script]] "r"]
    set file_data [read $fp]
    close $fp
    set data [split $file_data "\n"]
    foreach line $data {
        lappend script $line
    }
    ### write each line to h5md
    h5mdfile H5Screate_simple type str dims [llength $script]
    h5mdfile H5Pset_chunk dims 1
    h5mdfile H5Dcreate2 "parameters/files/script"
    h5mdfile H5Sclose_simple
    h5mdfile H5Dopen2 "parameters/files/script"
    for { set index 0 } { $index < [llength $script] } { incr index } {
        h5mdfile H5_write_value value [lindex $script $index] index $index
    }
    h5mdfile H5Dwrite
    h5mdfile H5_free_memory
    h5mdfile H5Dclose
    ### write command line arguments to seperate dataset
    if { $::argc > 0 } {
        h5mdfile H5Screate_simple type str dims $::argc
        h5mdfile H5Pset_chunk dims 1
        h5mdfile H5Dcreate2 "parameters/files/arguments"
        h5mdfile H5Sclose_simple
        h5mdfile H5Dopen2 "parameters/files/arguments"
        for { set index 0 } { $index < $::argc } { incr index } {
            h5mdfile H5_write_value value [lindex $::argv $index] index $index
        }
        h5mdfile H5Dwrite
        h5mdfile H5_free_memory
        h5mdfile H5Dclose
    }
    set which [info globals]
    foreach ixx {quiet tcl_version argv argv0 tcl_interactive auto_oldpath errorCode  
        auto_path errorInfo tcl_pkgPath tcl_patchLevel argc tcl_libPath
        tcl_library tk_strictMotif tk_library tk_version tk_patchLevel tcl_platform env} {
        set ind [lsearch -exact $which $ixx]
        set which [lreplace $which $ind $ind]
    }
    if { [llength $which] > 0 } {
        h5mdfile H5Screate_simple type str dims [llength $which] 2
        h5mdfile H5Pset_chunk dims [llength $which] 1
        h5mdfile H5Dcreate2 "parameters/files/tclvariables"
        h5mdfile H5Sclose_simple
        h5mdfile H5Dopen2 "parameters/files/tclvariables"
        for { set i 0 } { $i < [llength $which] } { incr i } {
            if {[array exists ::[lindex $which $i] ]==0 } {
                    h5mdfile H5_write_value value [lindex $which $i] index $i 0 
                    h5mdfile H5_write_value value [set ::[lindex $which $i]] index $i 1
            } else {
                puts "Did not write tcl array variable [lindex $which $i] to h5md dataset"    
            }
        }
        h5mdfile H5Dwrite
        h5mdfile H5_free_memory
        h5mdfile H5Dclose
    }
}
