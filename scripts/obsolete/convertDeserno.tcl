# Copyright (C) 2010,2012,2013 The ESPResSo project
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
#############################################################
#                                                           #
# convertDeserno.tcl                                        #
# ==================                                        #
#                                                           #
# convertDeserno2MD                                         #
# -----------------                                         #
# Reads in an (initial) particle configuration given in     #
# Deserno-format (e.g. created by 'polygelSetup') and       #
# converts it for 'Espresso', setting up global variables   #
# corresponding to the parameters and broadcasting the      #
# appropriate values to the 'Espresso' program, including   #
# particle configuration and interactions.                  #
#                                                           #
# Input:                                                    #
# - A file $origin containing the configuration             #
#   in Deserno-format (e.g. as created by 'polygelSetup')   #
#   supplied as parameter containing the full path.         #
#                                                           #
# Output:                                                   #
# - A gzipped file $destination ready for use with          # 
#   'Espresso' supplied in AxA's-blockfile-format;          #
#   the output file is suppressed if{$destination=="-1"}.   #
# - Global variables for use with any Espresso-script, i.e.:#
#   + lists 'conf' & 'corr', which contain the parameter's  #
#     names in Deserno's and A&L's world, respectively      #
#   + bunch of parameters having the same names as in       #
#     'polygelSetup.c'                                      #
#                                                           #
#                                                           #
# convertMD2Deserno                                         #
# -----------------                                         #
# The same as 'convertDeserno2MD', just vice versa,         #
# i.e. it reads in a configuration file in AxA's new        #
# 'Espresso'-blockfileformat and creates a Deserno-compatible # 
# style-file $destination ( if{$origin=="-1"} no input-file #
# is used, informations are taken from 'Espresso' instead). #
# Meanwhile, the blockfile-format does also provide info's  #
# about interactions etc., therefore if $origin is provided #
# as an input file, additional informations may not be      #
# required to be taken from 'Espresso'.                     #
# Since Deserno needs to know the chain-structure of the    #
# (sometimes crosslinked) polymers, this script does also   #
# have to analyze the topology of the given network/melt.   #
# Note that a polymer solution is always identified,        #
# whereas a cross-linked network can only be reconstructed  #
# if the chains are aligned consecutively such as           #
# 0-1-...-(MPC-1) MPC-...-(2*MPC-1) ...                     #
#                                                           #
#                                                           #
# ! =============                                         ! #
# ! = IMPORTANT =                                         ! #
# ! =============                                         ! #
# !                                                       ! #
# ! - All additional parameters are defined in            ! #
# !   'initConversion', so check that procedure regularly!! #
# ! - Make sure all 'global'ly defined variables are      ! #
# !   loaded into every other procedure as well!          ! #
# ! - It is possible to bypass initialization by directly ! #
# !   accessing the 'XXXmain'-scripts, which comes handy  ! #
# !   in particular when using 'MD2Deserno' on course of  ! #
# !   a simulation run with 'Espresso', because everything! #
# !   needed for a Deserno-compatible output is available ! #
# !   in some program's variable by then - all missing    ! #
# !   informations not extractable from $origin or        ! #
# !   'Espresso' itself (such as 'saveConfig' etc.) may be! #
# !   manually added (e.g. by setting 'saveConfig' in the ! #
# !   main script before calling the conversion-script)   ! #
# !   and will be included in the Deserno-file.           ! #
# !   However, use with extreme caution!                  ! #
# !   You have to know what you're doing, otherwise       ! # 
# !   unexpected behaviour may occur!                     ! #
# !   A recommendable good compromise would be:           ! #
# !   + manually invoke 'initConversion'                  ! #
# !   + change all the parameters you know but 'Espresso' ! #
# !     doesn't (such as 'saveConfig', 'seed', etc.)      ! #
# !   + continue with 'conversionMD2DesernoMain'          ! #
#                                                           #
#############################################################




#############################################################
#  initConversion                                           #
#############################################################

proc initConversion {} {
    puts "    Initializing conversion script..."

    # Include needed functions
    puts -nonewline "        Checking if all auxiliary functions have been loaded... "
    flush stdout
    if { [lsearch [info procs] "countBonds"] == -1 } {
	puts -nonewline "\n            Function 'countBonds' is missing - trying to reload... "
	flush stdout
	if { [catch [source countBonds.tcl]]!=0 } {
	    if { [catch [source ./scripts/countBonds.tcl]]!=0 && [catch [source $ESPRESSO_SCRIPTS/countBonds.tcl]]!=0 } {
		puts "Failed.\n"
		puts "----------------"
		puts "--> Warning! <--"
		puts "----------------"
		puts "--> Could not find 'countBonds.tcl' which is required for execution!"
		puts "Aborting..."
		exit
	    }
	}
	puts "Done."
    } elseif { [lsearch [info procs] "polyBlockWrite"] == -1 } {
	puts -nonewline "\n            Function 'polyBlockWrite' is missing - trying to reload... "
	flush stdout
	if { [catch [source auxiliary.tcl]]!=0 } {
	    if { [catch [source ./scripts/auxiliary.tcl]]!=0 && [catch [source $ESPRESSO_SCRIPTS/auxiliary.tcl]]!=0 } {
		puts "Failed.\n"
		puts "----------------"
		puts "--> Warning! <--"
		puts "----------------"
		puts "--> Could not find 'auxiliary.tcl' which is required for execution!"
		puts "Aborting..."
		exit
	    }
	}
	puts "Done."
    } else {
	puts "Done (all required functions accounted for)."
    }

    # The internal names of the Deserno-variables from 'polygelSetup' will be stored in 'conf',
    # the corresponding counterparts for 'Espresso' in 'corr' at the same relative position within the lists;
    # 'spec' specifies any special properties one has to take into account when using the parameters
    puts -nonewline "        Correlate parameter names... "
    flush stdout
    global conf corr spec iptq
    set conf [concat { prefix postfix seed startTime endTime deltaT 
	integrationSteps saveResults saveConfig N_P N_CPP chargeDist 
	val_CI N_Salt val_pS val_nS N_T boxLen subbox_1D 
	skin rcut_P3M k_FENE r_FENE eps_LJ rcut_LJ alpha mesh ip 
	Bjerrum Temp Gamma N_CR } ]
    # ...and correlate it with the ones from 'Espresso'
    # ('NA' = there's no implementation in 'Espresso')
    set corr [concat { NA NA NA NA NA time_step 
        NA NA NA NA NA NA 
        NA NA NA NA max_part box_l NA
        skin p3m_r_cut NA NA NA NA p3m_alpha p3m_mesh NA
        bjerrum temp gamma NA } ]
    # Note the special properties some 'Espresso' variables have:
    # 'ro' = read only; '3d' = 3D input required; 'OK' = nothing special
    set spec [concat { OK OK OK OK OK OK
	OK OK OK OK OK OK
	OK OK OK OK ro 3d OK
	OK OK OK OK OK OK OK 3d OK
	OK OK OK OK } ]
    set iptq "id pos type q"
    # sets which informations (out of pos|type|q|v|f) on the particles should be saved to disk ('Deserno2MD' only)
    puts "Done."
    
    puts -nonewline "        Preparing global variables... "
    flush stdout
    # Determine which type-number should be assigned to monomers, counter-ions, salt-ions, ...
    global type_P type_CI type_S 
    set type_P 0
    set type_CI 1
    set type_S 2
    # Determine which type-number should be assigned to the interactions
    global type_FENE type_bend
    set type_FENE 0
    set type_bend 1
    # Specify additional parameters for the potentials
    global pot_bend pot_LJ_sigma pot_LJ_shift pot_LJ_offset pot_ramp_cut pot_ramp_fmax
    set pot_bend 1
    set pot_LJ_sigma 1.0
    set pot_LJ_shift 0
    set pot_LJ_offset 0
    set pot_ramp_cut 1.3
    set pot_ramp_fmax 100
    # Define additional parameters which are not supplied by Deserno but have to be calculated later on
    global MPC N_CI N_pS N_nS N_S
       # MPC = [expr ($N_CPP-1)*$chargeDist+1] is the number of monomers per polymer chain
       # N_CI = [expr $N_P*$N_CPP/$val_CI] is the number of counterions
       # N_pS = [expr $N_Salt*$val_pS] is the number of positively charged salt ions
       # N_nS = [expr $N_Salt*$val_nS] is the number of negatively charged salt ions
       # N_S = [expr $N_pS+$N_nS] is the total number of salt ions 
       #       (not to be confused with N_Salt, which is the number of salt _molecules_)
    # Define auxiliary parameters
    global step foldCoord
       # step = amount of integration steps taken so far
       #        (only used for MD2Deserno; has to be supplied by user)
    set foldCoord 0
       # foldCoord = decides whether to periodically fold the particle's coordinates before working with them; 
       #             this is important, since 'polygelSetup' by default creates an unfolded system of polymers,
       #             centered in a box with boxlenght L, but extending outside, too.
    puts "Done."

    # Mark this initialization procedure as being completed
    global convertInitDone
    set convertInitDone 1
    puts "    Initialization completed."
}



#############################################################
#  convertDeserno2MD                                        #
#############################################################

proc convertDeserno2MD {origin destination} {
    # Initialize necessary parameters
    initConversion

    # Continue
    convertDeserno2MDmain "$origin" "$destination"
}


#############################################################
#  convertDeserno2MDmain                                    #
#############################################################

proc convertDeserno2MDmain {origin destination} {
    # Check if initialization was done
    global convertInitDone
    if { $convertInitDone!=1 } {
	puts "    !!! Warning !!! You're bypassing the initialization sequence! Continuing at your own risk..."
    }

    # Note that tcl requires to have the name of used global variables repeated here
    # => cross-check with 'initConversion' to make sure everything there matches an entry here
    global conf corr spec iptq
    global type_P type_CI type_S
    global type_FENE type_bend
    global pot_bend pot_LJ_sigma pot_LJ_shift pot_LJ_offset pot_ramp_cut pot_ramp_fmax
    global MPC N_CI N_pS N_nS N_S
    global foldCoord


    # Check the supplied file-names for existence and correctness
    if { ![ file exists $origin ] } {
	puts "Failed.\n "
	puts "----------------"
	puts "--> Warning! <--"
	puts "----------------"
	puts "--> The file containing the Deserno-configuration could not be found!"
	puts "--> (You supplied: $origin)"
	puts "Aborting...\n"
	exit
    }
    if { $destination == "-1" } {
	puts "    No output-file was given, everything will be submitted without saving to 'Espresso' only..."
    }


    # Open the input-file which must have the Deserno-format, e.g. created by 'polygelSetup'
    # looking like this:
    # Line defining a parameter             --> '# Description............: <value>'
    # Line defining coordinates             --> '<x> <y> <z> <vx> <vy> <vz> <charge>'
    # Line defining crosslinking parameters --> '# Description......: var = <value>'
    # Line defining the crosslinks          --> '<crosslink origin> <crosslink destination>'
    puts -nonewline "    Preparing Deserno-configuration to be read from '[lindex [split $origin \"/\"] end]'... "
    flush stdout
    set in [open $origin "r"]
    puts "Done."


#  Import entries in Deserno-file
#############################################################

    # Now get the file's entries and assign them to the corresponding variable
    puts "    Conversion in progress..."
    flush stdout
    set i 0
    set j 0
    while { ![eof $in] } {
	# Consecutively load in each line of the input-file
	# and assign it to the names in 'conf'
	# (call it 'dummy' if you've run out of names)
	if {$i<[llength $conf]} {
	    set tmp_var [lindex $conf $i]
	} else {
	    set tmp_var "dummy"
	}
	eval global $tmp_var T$tmp_var
	eval gets $in $tmp_var
	# This trick allows to access the current variable's name by '$tmp_var'
	# while using the variable itself by 'eval ... $$tmp_var'

#  Submit Parameters
#############################################################

	# If the current line starts with '#' it contains a parameter
	if { [eval string compare -length 1 \"$$tmp_var\" \"\#\"]==0 } then {
	    if {$j==0} { set sep ":" } else { set sep "=" }

	    # Check for consistency
	    if { ![expr $i<[llength $conf]] } {
		puts "Failed.\n"
		puts "----------------"
		puts "--> Warning! <--"
		puts "----------------"
		puts "--> There are more parameters found in the file than specified in 'conf'!"
		puts "--> Supply more tags in 'conf', 'corr' and 'spec'!"
		puts "Aborting...\n"
		exit
	    }

	    # Split the line into its components:
	    # Left part is the parameter's description -> e.g. 'Tprefix' (with 'T'='Text')
	    eval set T$tmp_var \"[lindex [eval split \"$$tmp_var\" $sep] 0]: \"
	    # Right part is the value of the parameter -> e.g. 'prefix' (as stored in 'conf')
	    eval set $tmp_var \"[string trim [lindex [eval split \"$$tmp_var\" $sep] 1] ]\"
	    
	    # If possible, sent the value to 'Espresso':
	    if { [string compare "[lindex $corr $i]" "NA"]!=0 } {
		# Intercept 'read-only' parameters until they are available in 'Espresso'
		if { [string compare "[lindex $spec $i]" "ro"]!=0 } {
		    puts -nonewline "        ---> Espresso <--- Submitting "
		    eval puts -nonewline \"'$[concat "T$tmp_var"] [lindex $conf $i]' \"
		    puts -nonewline "\t to   "
                    eval puts \"'[lindex $corr $i]' (= $$tmp_var).\"
		    
		    # some parameters require 3D, although Deserno specified only one
		    if { [string compare "[lindex $spec $i]" "3d"]==0 } {
			setmd [lindex $corr $i] [eval concat \"$$tmp_var\"] \
			    [eval concat \"$$tmp_var\"] [eval concat \"$$tmp_var\"]
		    } else {
			# regular parameters don't make trouble  ;-/
			setmd [lindex $corr $i] [eval concat \"$$tmp_var\"]
		    }
		} else {
		    # If read-only, the parameter is not submitted, but stored in a global variable
		    # for further use in this or any other script
		    puts -nonewline "        Setting read-only variable  "
		    eval puts \"'$[concat "T$tmp_var"] [lindex $conf $i]' (= $$tmp_var).\"
		}
	    } else {
		# Even if they're not sent to 'Espresso' the parameters are set up globally
		# for further use in this or any other script
		puts -nonewline "        Assigning global variable   "
		eval puts \"'$[concat "T$tmp_var"] [lindex $conf $i]' (= $$tmp_var).\"
	    }

	    # Goto next parameter in 'conf'
	    incr i

	# Otherwise (i.e. if the current line read from $origin does not start with '#') 
	# it's probably either the particles' coordinates or the crosslinks' configurations
	} else {
	    # Switch from indirect variables (needed for the parameters) to direct access
	    eval set tmp_var $$tmp_var
	    
	    # Distinguish between particles (second block) and crosslinks (fourth block)
	    if {$j==0} {
		# Derive some additional parameters
		set MPC  [expr ($N_CPP-1)*$chargeDist+1]
		set N_CI [expr $N_P*$N_CPP/$val_CI]
		set N_pS [expr $N_Salt*$val_pS]
		set N_nS [expr $N_Salt*$val_nS]
		set N_S  [expr $N_pS+$N_nS] 

#  Submit Interactions
#############################################################

		# Create the bonded interactions
		puts "        Preparing bonded interactions: "
		puts -nonewline "            Lennard-Jones (eps=[format %4.2f $eps_LJ], "
		puts -nonewline "sigma=[format %4.2f $pot_LJ_sigma], rcut=[format %4.2f $rcut_LJ], "
		puts "shift=[format %4.2f $pot_LJ_shift], offset=[format %4.2f $pot_LJ_offset]) between: "
		puts -nonewline "                Polymer-Polymer... "
		flush stdout
		inter $type_P $type_P lennard-jones $eps_LJ $pot_LJ_sigma $rcut_LJ $pot_LJ_shift $pot_LJ_offset
		puts -nonewline "Polymer-Counterion... "
		flush stdout
		inter $type_P $type_CI lennard-jones $eps_LJ $pot_LJ_sigma $rcut_LJ $pot_LJ_shift $pot_LJ_offset
		puts -nonewline "Polymer-Salt... "
		flush stdout
		inter $type_P $type_S lennard-jones $eps_LJ $pot_LJ_sigma $rcut_LJ $pot_LJ_shift $pot_LJ_offset
		puts " "
		puts -nonewline "                Counterion-Counterion... "
		flush stdout
		inter $type_CI $type_CI lennard-jones $eps_LJ $pot_LJ_sigma $rcut_LJ $pot_LJ_shift $pot_LJ_offset
		puts -nonewline "Counterion-Salt... "
		flush stdout
		inter $type_CI $type_S lennard-jones $eps_LJ $pot_LJ_sigma $rcut_LJ $pot_LJ_shift $pot_LJ_offset
		puts -nonewline "Salt-Salt... "
		flush stdout
		inter $type_S $type_S lennard-jones $eps_LJ $pot_LJ_sigma $rcut_LJ $pot_LJ_shift $pot_LJ_offset
		puts " "

		# Create the non-bonded interactions
		puts "        Preparing non-bonded interactions: "
		puts "            Type $type_FENE: FENE (k=[format %4.2f $k_FENE], r=[format %4.2f $r_FENE])... "
		inter $type_FENE fene $k_FENE $r_FENE
		puts "            Type $type_bend: Cosine (pot_bend=$pot_bend)... "
		inter $type_bend angle $pot_bend

		puts "        Interactions completed."

#  Submit Polymers
#############################################################

		# Setting up the polymers
		puts -nonewline "        Submitting $N_P polymers with $MPC monomers (of type $type_P) each... "
		flush stdout
		set old [expr 10*($N_P*$MPC)]
		set k 0
		for {set j 0} {$j<[expr $N_P*$MPC]} {incr j} {
		    # Draw a status bar for every 10% of progress
		    if { [expr $j*100]>=$old} {
			set old [expr $old + 10*($N_P*$MPC)]
			puts -nonewline ".,."
			flush stdout
		    } 
		    # upon entering this part of the program we already have an unprocessed line stored in 'tmp_var'
		    # (which we used to detect that it does not contain a parameter), hence we only load additional
		    # lines once we've processed this one, i.e. for $j>0
		    if { $j>0 } { gets $in tmp_var }

		    # Submit the data to 'Espresso'
		    if { $foldCoord == 1 } {
			part $j pos [expr [lindex $tmp_var 0] - floor([lindex $tmp_var 0] / $boxLen)*$boxLen] \
			    [expr [lindex $tmp_var 1] - floor([lindex $tmp_var 1] / $boxLen)*$boxLen] \
			    [expr [lindex $tmp_var 2] - floor([lindex $tmp_var 2] / $boxLen)*$boxLen]
		    } else {
			part $j pos [lindex $tmp_var 0] [lindex $tmp_var 1] [lindex $tmp_var 2]
		    }
		    part $j v   [lindex $tmp_var 3] [lindex $tmp_var 4] [lindex $tmp_var 5]
		    part $j q   [lindex $tmp_var 6]
		    part $j type $type_P

		    # Setting a bond between all $MPC monomers of the same polymer
		    # => omit monomer 0, 1*$MPC, 2*$MPC,..., ($N_P-1)*$MPC
		    if { [expr $j%$MPC]!=0 } {
			part $j bond $type_FENE [expr $j-1]
			# Note that Espresso stores bonds only with the particle of lower particle-number;
			# hence, a bond between 'a' and 'b' for 'a<b' is stored as '{0 b}' with particle 'a',
			# whereas particle 'b' doesn't seem to have any bond...
			# => use the function 'countBonds' to get a tcl-list 'bonds' 
			# which has '{0 a}' at particle 'b' as well!
			incr k
		    }
		}
		puts ".,. $j particles total, $k of which bounded (by type $type_FENE)."

#  Submit remaining Particles
#############################################################

		# Setting up some counterions
		puts -nonewline "        Submitting $N_CI counterions (of type $type_CI)... "
		flush stdout
		set old [expr 10*$N_CI+100*($N_P*$MPC)]
		for {set j $j} {$j<[expr $N_CI+($N_P*$MPC)] } {incr j} {
		    # Draw a status bar for every 10% of progress
		    if { [expr $j*100]>=$old } {
			set old [expr $old + 10*$N_CI]
			puts -nonewline ".,."
			flush stdout
		    } 
		    gets $in tmp_var
		    # Submit the data to 'Espresso'
		    if { $foldCoord == 1 } {
			part $j pos [expr [lindex $tmp_var 0] - floor([lindex $tmp_var 0] / $boxLen)*$boxLen] \
			    [expr [lindex $tmp_var 1] - floor([lindex $tmp_var 1] / $boxLen)*$boxLen] \
			    [expr [lindex $tmp_var 2] - floor([lindex $tmp_var 2] / $boxLen)*$boxLen]
		    } else {
			part $j pos [lindex $tmp_var 0] [lindex $tmp_var 1] [lindex $tmp_var 2]
		    }
		    part $j v   [lindex $tmp_var 3] [lindex $tmp_var 4] [lindex $tmp_var 5]
		    part $j q   [lindex $tmp_var 6]
		    part $j type $type_CI
		}
		puts ".,. [expr $j-($N_P*$MPC)] particles with valency $val_CI total."

		# Setting up the salt-molecules
		puts -nonewline "        Submitting $N_nS "
		for {set k 0} {$k<$val_nS} {incr k} {puts -nonewline "-"}
		puts -nonewline " and $N_pS "
		for {set k 0} {$k<$val_pS} {incr k} {puts -nonewline "+"}
		puts -nonewline " salt ions (of type $type_S)... "
		flush stdout
		set old [expr 10*$N_S+100*($N_CI+$N_P*$MPC)]
		for {set j $j} {$j<[expr $N_S+($N_CI+$N_P*$MPC)] } {incr j} {
		    # Draw a status bar for every 10% of progress
		    if { [expr $j*100]>=$old } {
			set old [expr $old + 10*$N_S]
			puts -nonewline ".,."
			flush stdout
		    } 
		    gets $in tmp_var
		    # Submit the data to 'Espresso'
		    if { $foldCoord == 1 } {
			part $j pos [expr [lindex $tmp_var 0] - floor([lindex $tmp_var 0] / $boxLen)*$boxLen] \
			    [expr [lindex $tmp_var 1] - floor([lindex $tmp_var 1] / $boxLen)*$boxLen] \
			    [expr [lindex $tmp_var 2] - floor([lindex $tmp_var 2] / $boxLen)*$boxLen]
		    } else {
			part $j pos [lindex $tmp_var 0] [lindex $tmp_var 1] [lindex $tmp_var 2]
		    }
		    part $j v   [lindex $tmp_var 3] [lindex $tmp_var 4] [lindex $tmp_var 5]
		    part $j q   [lindex $tmp_var 6]
		    part $j type $type_S
		}
		puts ".,. [expr $j-($N_CI+$N_P*$MPC)] particles in $N_Salt molecules total."

#  Submit Crosslinks
#############################################################

	    # Deal with crosslinks between the polymer chains
	    } elseif {$j==[expr ($N_CI+$N_P*$MPC)+$N_S] && [string length $tmp_var]>0} {
		# Connecting cross-linked monomers
		puts -nonewline "        Crosslinking [expr 2*$N_CR] monomers... "
		flush stdout
		set old [expr 10*$N_CR+100*($N_S+$N_CI+$N_P*$MPC)]
		for {set j $j} {$j<[expr $N_CR+($N_S+$N_CI+$N_P*$MPC)]} {incr j} {
		    # Draw a status bar for every 10% of progress
		    if { [expr $j*100]>=$old } {
			set old [expr $old + 10*$N_CR]
			puts -nonewline ".,."
			flush stdout
		    }
		    # upon entering this part of the program we already have an unprocessed line stored in 'tmp_var'
		    # (which we used to detect that it does not contain a parameter), hence we only load additional
		    # lines once we've processed this one, i.e. for $j>($N_S+$N_CI+$N_P*$MPC)
		    if { $j>[expr $N_S+$N_CI+$N_P*$MPC] } { gets $in tmp_var }

		    # Crosslink particles
		    part [lindex $tmp_var 0] bond $type_FENE [lindex $tmp_var 1]
		    # Note that Espresso stores bonds only with the particle of lower particle-number;
		    # hence, a bond between 'a' and 'b' for 'a<b' is stored as '{0 b}' with particle 'a',
		    # whereas particle 'b' doesn't seem to have any bond...
		    # => use the function 'countBonds' to get a tcl-list 'bonds' 
		    # which has '{0 a}' at particle 'b' as well!
		}
		puts ".,. [expr $j-($N_S+$N_CI+$N_P*$MPC)] additional bonds (of type $type_FENE) created."
	    } elseif {[string length $N_CR]==0 } {
		set N_CR 0
	    }
	}
    }
    # Finish up the conversion part by closing the input-file
    close $in
    puts "    Conversion completed."


#  Write output file
#############################################################

    # Save all converted configurations, parameters, and interactions to $destination
    # using AxA's blockfile-format for 'Espresso'; skip if user specified "-1" as file-name
    if { $destination != "-1" } {
	# The converted configurations should be compressed on-the-fly using 
	# 'set f [open "|gzip -c - >$destination" w]' for output, allowing later usage with
	# 'set f [open "|gzip -cd $destination" r]' to read it in
	# This is done in 'polyBlockWrite' which writes all we've done so far to a 'Espresso'-compatible file
	puts -nonewline "    Saving all configurations in 'Espresso'-format to '[lindex [split $destination \"/\"] end]'... "
	flush stdout
	foreach j $corr {
	    if { $j != "NA" } { lappend tmp_corr $j }
	}
	polyBlockWrite $destination $tmp_corr $iptq
	puts "Done "

    } else {
	puts "    Skipping to save converted configurations: Nothing will be written to disk, only 'Espresso' will have them!"
    }


    # End this script
    puts "    Function successfully completed. Returning control to calling script..."
}




#############################################################
#  convertMD2Deserno                                        #
#############################################################

proc convertMD2Deserno {origin destination} {
    # Initialize necessary parameters
    initConversion

    # Set parameters which cannot be determined from 'Espresso' or from $origin
    global conf
    foreach i $conf { global $i }
    global step
    set prefix AA0000
    set postfix 0
    set seed -1
    set startTime -1
    set endTime -1
    set integrationSteps -1
    set saveResults -1
    set saveConfig -1
    set subbox_1D -1
    set ip -1
    set step -1

    # Continue
    convertMD2DesernoMain "$origin" "$destination"
}


#############################################################
#  convertMD2DesernoMain                                    #
#############################################################

proc convertMD2DesernoMain {origin destination} {
    # Check if initialization was done
    global convertInitDone
    if { $convertInitDone!=1 } {
	puts "    !!! Warning !!! You're bypassing the initialization sequence! Continuing at your own risk..."
    }

    # Note that tcl requires to have the name of used global variables repeated here
    # => cross-check with 'initConversion' to make sure everything there matches an entry here
    global conf corr spec iptq
    foreach i $conf { global $i }
    global type_P type_CI type_S
    global type_FENE type_bend
    global pot_bend pot_LJ_sigma pot_LJ_shift pot_LJ_offset pot_ramp_cut pot_ramp_fmax
    global MPC N_CI N_pS N_nS N_S
    global step


#  Import entries from file
#############################################################

    # If the user specified '-1' for $origin, then there is no input-file and everything is directly taken from Espresso
    if { [string compare "$origin" "-1"]==0 } {
	puts "    No input-file was specified, let's hope that all required informations are accessible from 'Espresso'!"
    } else {
	# Open the (compressed) input-file which is expected to contain a configuration in 'Espresso'-blockfile-format
	puts -nonewline "    Preparing 'Espresso'-configuration to be read from '[lindex [split $origin \"/\"] end]'... "
	flush stdout
	if { [string compare [lindex [split $origin "."] end] "gz"]==0 } {
	    set in [open "|gzip -cd $origin" "r"]
	} else {
	    set in [open "$origin" "r"]
	}
	puts "Done."

	# Now get the file's entries
	puts "    Reading content of '[lindex [split $origin \"/\"] end]' (this may take a while)... "
	while {[blockfile $in read auto] != "eof" } {}
	puts "        Got coordinates for [expr [setmd max_part] + 1] particles total..."
	puts "        Got all the bonds between particles..."
	puts "        Got all the interactions..."	
	close $in
	puts "    ... and closed the file."
    }


#  Get the rest from 'tcl-md'
#############################################################

    # Now everything should either have been loaded from the blockfile $origin into 'Espresso'
    # or (in case $origin==-1) it was already there
    puts "    Gathering informations from 'Espresso'... "
    flush stdout

    # Load everything useful out of 'Espresso'
    puts -nonewline "        Getting"
    foreach i $corr j $conf {
	if { [string compare "$i" "NA"]!=0 } {
	    puts -nonewline " $i/"
	    flush stdout
	    set $j [lindex [setmd $i] 0]
	    puts -nonewline "$j ."
	}
    }

    # Load interactions parameters:
    puts "..\n        Load interactions parameters..."
    flush stdout
    set tmp_FENE 0
    set tmp_LJ 0
    set tmp_int [inter]
    for {set i 0} {$i<[llength $tmp_int]} {incr i} {
	set tmp_var [lindex $tmp_int $i]
	#  FENE,...
	if { [string compare [lindex $tmp_var 1] "FENE"]==0 } {
	    if { $tmp_FENE==0 } {
		puts -nonewline "            Setting FENE-parameters "
		set k_FENE [lindex $tmp_var 2]
		puts -nonewline " (k_FENE = $k_FENE, "
		flush stdout
		set r_FENE [lindex $tmp_var 3]
		puts -nonewline "r_FENE = $r_FENE)... "
		set tmp_FENE 1
		puts "Done."
	    } elseif { $k_FENE!=[lindex $tmp_var 2] || $r_FENE!=[lindex $tmp_var 3] } {
		puts "Failed.\n"
		puts "----------------"
		puts "--> Warning! <--"
		puts "----------------"
		puts "--> FENE-interaction parameters must be constant for all bonds!"
		puts "--> (Found $k_FENE=k_FENE=[lindex $tmp_var 2] and $r_FENE=r_FENE=[lindex $tmp_var 3])"
		puts "Aborting...\n"
		exit
	    }
	# ...and Lennard-Jones
	} elseif {[string compare [lindex $tmp_var 2] "lennard-jones"]==0} {
	    if { $tmp_LJ==0 } {
		puts -nonewline "            Setting Lennard-Jones-parameters "
		set eps_LJ [lindex $tmp_var 3]
		puts -nonewline " (eps_LJ = $eps_LJ, "
		set rcut_LJ [lindex $tmp_var 5]
		puts -nonewline "rcut_LJ = $rcut_LJ)... "
		set tmp_LJ 1
		puts "Done."
	    } elseif { $eps_LJ!=[lindex $tmp_var 3] || $rcut_LJ!=[lindex $tmp_var 5] } {
		puts "Failed.\n"
		puts "----------------"
		puts "--> Warning! <--"
		puts "----------------"
		puts "--> Lennard-Jones parameters must be constant for all particles!"
		puts "--> (Found $eps_LJ=eps_LJ=[lindex $tmp_var 3] and $rcut_LJ=rcut_LJ=[lindex $tmp_var 5])"
		puts "Aborting...\n"
		exit
	    }
	}
    }
    puts "        Done with interactions."

    # Load all particle data
    puts -nonewline "        Getting informations on any particle... "
    flush stdout
    set part_all [part]
    set N_T [llength $part_all]
    puts "Done (got $N_T particles)."

    # Everything should be back in this function now
    puts "    Gathering completed."


#  Re-build network structure
#############################################################

    # Now create a tcl-list with all bonds in it
    set bonds [countBonds $part_all]


    # Get an idea of the network structure
    puts "    Analyzing network structure... "

    # In a polymer melt the monomers on a chain have usually two bonds while their ends have just one;
    # in a end-to-end-linked network the ends are cross-linked, too, so that their new partners have three bonds
    puts -nonewline "        Checking number of bonds on $N_T particles... "
    flush stdout
    set part_end {}
    set part_lin {}
    set part_crs {}
    set old [expr 10*$N_T]
    for {set i 0} {$i<$N_T} {incr i} {
	# Draw a status bar for every 10% of progress
	if { [expr $i*100]>=$old } {
	    set old [expr $old + 10*$N_T]
	    puts -nonewline ".,."
	    flush stdout
	} 
	# Check if the current particle has any cross-links (e.g. counter-ions have not).
	# Note, that $bonds contains an entry for every particle (i.e. the particle_number),
	# hence the number of bonds is -1 smaller than the length of that entry.
	if { [llength [lindex $bonds $i]]>1 } {
	    # Check the number of cross-links:
	    switch [llength [lindex $bonds $i]] {
		2 { # 1 bond:  particle is the end of a polymer chain
		    lappend part_end $i }
		3 { # 2 bonds: particle is monomer in a chain
		    lappend part_lin $i }
		4 { # 3 bonds: particle also has a cross-link to another chain
		    lappend part_crs $i }
		default {
		    # 4+ bonds??? => hmm... too much?!
		    puts "Failed.\n"
		    puts "----------------"
		    puts "--> Warning! <--"
		    puts "----------------"
		    puts "--> Current particle [lindex [lindex $bonds $i] 0] (found at index $i) has too many bonds!"
		    puts "--> (Got: [expr [llength [lindex $bonds $i]]-1]; expected: 1-3)"
		    puts "Aborting...\n"
		    exit
		}
	    }
	} else {
	    # no bonds, no cross-links => must be something else (counter-ion, etc.), but must _not_ be a monomer
	    if { [lindex [lindex $part_all $i] [findPropPos [lindex $part_all $i] type]]!=$type_P } {
		lappend part_else $i
	    } else {
		puts "Failed.\n"
		puts "----------------"
		puts "--> Warning! <--"
		puts "----------------"
		puts -nonewline "--> Current particle [lindex [lindex $bonds $i] 0] (found at index $i) "
		puts "is a loose monomer (type $type_P) without any bonds!"
		puts "--> (Stats: [lindex $part_all $i])"
		puts "Aborting...\n"
		exit
	    }
	}
    }
    puts ".,. found [llength $part_end] ends, [llength $part_lin] linear bonds, [llength $part_crs] cross-links."


#  Identify Polymer Melt
#############################################################

    # If 'part_crs' is empty, there's no network but only a melt; hence, the polymer chains are easy to identify
    if { [llength $part_crs]==0 } {
	puts "        => polymer melt/solution detected."
	set N_CR 0
	puts -nonewline "        Identifying polymer chains... "
	flush stdout
	set tmp_CPP 0
	set tmp_cD 0
	set tmp_NP 0
	set tmp_MPC 0
	set old_i 0
	set old_e [llength $part_end]
	set old [expr 10*$old_e]

	# loop all particles found having only one bond, because these are the chains' ends
	foreach i $part_end {
	    # draw a status bar for every 10% of progress
	    if { [expr $old_i*100]>=$old } {
		set old [expr $old + 10*$old_e]
		puts -nonewline ".,."
		flush stdout
	    }
	    incr old_i

	    # make sure that current end is not the end of a previously reconstructed chain
	    set tmp_var 0
	    if { $tmp_NP>0} {
		for {set j 0} {$j<[llength $polymer_chains]} {incr j} {
		    if { [expr [lindex [lindex $polymer_chains $j] 0]==$i] || \
			     [expr [lindex [lindex $polymer_chains $j] end]==$i] } {
			set tmp_var 1
		    }
		}
	    }
            if { $tmp_var==1 } continue

	    # build the chain following the bond information stored with each particle
	    set chains [lindex [lindex $bonds $i] 0]
	    set tmp_CPP_k 0
	    set tmp_cD_k -1
	    set tmp_MPC_k 0
            while { $tmp_MPC_k<$N_T } {
		set j [lindex $chains end]
		set j_bond [lrange [lindex $bonds $j] 1 end]
		set j_bond_nr 0
		# is current particle charged? => keep track of the amount of charges and the distances between them
		if { [lindex [lindex $part_all $j] [findPropPos [lindex $part_all $j] q]]!=0 } {
		    if { [expr $tmp_cD_k+1]==$tmp_cD || $tmp_cD_k==-1 } {
			incr tmp_CPP_k
			set tmp_cD_k 0
		    } elseif { $tmp_cD==0 } {
			incr tmp_CPP_k
			set tmp_cD [expr $tmp_cD_k+1]
			set tmp_cD_k 0
		    } else {
			puts "Failed.\n"
			puts "----------------"
			puts "--> Warning! <--"
			puts "----------------"
			puts "--> The distance between charges on a chain must be constant!"
			puts "--> (Found a chain with distances $tmp_cD and $tmp_cD_k)"
			puts "Aborting...\n"
			exit
		    }
		} elseif { $tmp_cD_k>-1 } { incr tmp_cD_k }
		# agglomerate the monomers
		if { [llength $j_bond]==1 } {
		    if { $tmp_MPC_k>0 } {
			# last particle -> wrap up evaluation
			incr tmp_MPC_k
			break
		    } 
		} else {
		    # check which of the bonds lead "forward" on the chain
		    if { [lindex [lindex $j_bond 0] 1]==[lindex $chains end-1] } {
			set j_bond_nr 1
		    }
		}
		lappend chains [lindex [lindex $j_bond $j_bond_nr] 1]
		incr tmp_MPC_k
	    }

            # check for consistency
	    if { $tmp_CPP_k!=$tmp_CPP } {
		if { $tmp_CPP==0 } {
		    set tmp_CPP $tmp_CPP_k
		} else {
		    puts "Failed.\n"
		    puts "----------------"
		    puts "--> Warning! <--"
		    puts "----------------"
		    puts "--> The amount of charges on a chain must be constant!"
		    puts "--> (Found chains with $tmp_CPP and $tmp_CPP_k)"
		    puts "Aborting...\n"
		    exit
		}
	    }		
	    if { $tmp_MPC_k==$tmp_MPC } {
		lappend polymer_chains $chains
		incr tmp_NP
	    } elseif { [lindex $part_end 0]==$i} {
		set tmp_MPC $tmp_MPC_k
		lappend polymer_chains $chains
		incr tmp_NP
	    } else {
		puts "Failed.\n"
		puts "----------------"
		puts "--> Warning! <--"
		puts "----------------"
		puts "--> The number of monomers on all chains must be constant!"
		puts "--> (Found chains with $tmp_MPC and $tmp_MPC_k)"
		puts "Aborting...\n"
		exit
	    }
        }
	set N_CPP $tmp_CPP
	set chargeDist $tmp_cD
	set N_P $tmp_NP
	set MPC $tmp_MPC
        puts -nonewline ".,. reconstructed $N_P chains with $MPC monomers each "
	puts -nonewline "($N_CPP charged with distance $chargeDist"
	if { $chargeDist>0 } {
	    puts -nonewline ", like  "
	    for {set j 0} {$j<3} {incr j} {
		puts -nonewline "c"
		for {set i 0} {$i<[expr $chargeDist-1]} {incr i} { puts -nonewline "." }
	    }   
            puts -nonewline "c  etc."
	}
	puts ")."


#  Identify cross-linked Network
#############################################################

    # If 'part_end' is empty, this is an end-to-end cross-linked network, 
    # where any of the three partners in 'part_end' may be one of the chain's end;
    # therefore, one has to have a closer look
    } elseif { [llength $part_end]==0 } {
	puts "        => end-to-end crosslinked network detected."
	puts "        (Note that the reconstruction will only be successful if the chains are aligned consecutively!)"
	puts -nonewline "        Attempting to identify polymer chains... "
	flush stdout
	set tmp_CPP 0
	set tmp_cD 0
	set tmp_NP 0
	set tmp_MPC 0
	set tmp_NCR 0
	set old_i 0
	set old_e [expr [llength $part_crs]+[llength $part_lin]]
	set old [expr 10*$old_e]

	# loop all particles with at least two bonds
	for {set i 0} {$i<$old_e} {incr i} {
	    # draw a status bar for every 10% of progress
	    if { [expr $old_i*100]>=$old } {
		set old [expr $old + 10*$old_e]
		puts -nonewline ".,."
		flush stdout
	    }
	    incr old_i

	    # get the bonding informations on the current particle & check if they comply with the rules
	    set tmp_bond [lindex $bonds $i]
	    set tmp_current [lindex $tmp_bond 0]
	    if {[llength $tmp_bond]<2} {
		puts "Failed.\n"
		puts "----------------"
		puts "--> Warning! <--"
		puts "----------------"
		puts "--> An end-to-end crosslinked network may only be reconstructed"
		puts "--> if the chains are aligned consecutively!"
		puts -nonewline "--> (Expected $old_e particles with bonds, but found particle $tmp_current "
		puts "(at index $i) having none)"
		puts "Aborting...\n"
		exit
	    }
	    set tmp_bond [lrange $tmp_bond 1 end]
	    
	    # build the chain following the bond information stored with each particle
	    # and making use of the requirement that all chains have to be aligned
	    set tmp_end -1
	    for {set j 0} {$j<[llength $tmp_bond]} {incr j} {
		set tmp_partner [lindex [lindex $tmp_bond $j] 1]
		if { $tmp_partner==[expr $tmp_current-1] } {
		    # bonding partner is left neighbour => append $i to the chain
		    lappend chains $tmp_current
		    incr tmp_MPC_k
		    # is current particle charged? => keep track of the amount of charges and the distances between them
		    if { [lindex [lindex $part_all $i] [findPropPos [lindex $part_all $i] q]]!=0 } {
			if { [expr $tmp_cD_k+1]==$tmp_cD || $tmp_cD_k==-1 } {
			    incr tmp_CPP_k
			    set tmp_cD_k 0
			} elseif { $tmp_cD==0 } {
			    incr tmp_CPP_k
			    set tmp_cD [expr $tmp_cD_k+1]
			    set tmp_cD_k 0
			} else {
			    puts "Failed.\n"
			    puts "----------------"
			    puts "--> Warning! <--"
			    puts "----------------"
			    puts "--> The distance between charges on a chain must be constant!"
			    puts "--> (Found a chain with distances $tmp_cD and $tmp_cD_k)"
			    puts "Aborting...\n"
			    exit
			}
		    } elseif { $tmp_cD_k>-1 } { incr tmp_cD_k }
		} elseif { $tmp_partner!=[expr $tmp_current+1] } {
		    # bonding partner is not on this chain => append to crosslinks
		    # if the current particle is at the end of a chain
		    if { [llength $tmp_bond]==2 } {
			set tmp_end $j
			if {$tmp_partner < $tmp_current} {
			    lappend cross "{$tmp_partner $tmp_current}"
			} else {
			    lappend cross "{$tmp_current $tmp_partner}"
			}
			incr tmp_NCR
		    }
		}
	    }
	    # if current particle is at the end of a chain,
	    # check if it's the left one => in that case open a new chain
	    if {$tmp_end > -1} {
		set tmp_partner [lindex [lindex $tmp_bond [expr $tmp_end*(-1)+1]] 1]
		if { $tmp_partner==[expr $tmp_current+1] } {
		    # left end => new chain
		    set chains $tmp_current
		    set tmp_MPC_k 1
		    # is current particle charged? => keep track of the amount of charges and the distances between them
		    if { [lindex [lindex $part_all $i] [findPropPos [lindex $part_all $i] q]]!=0 } {
			set tmp_CPP_k 1
			set tmp_cD_k 0
		    } else {
			set tmp_CPP_k 0
			set tmp_cD_k -1
		    }
		} else {
		    # right end => wrap up evaluation of this chain by checking for consistency
		    if { $tmp_CPP_k!=$tmp_CPP } {
			if { $tmp_CPP==0 } {
			    set tmp_CPP $tmp_CPP_k
			} else {
			    puts "Failed.\n"
			    puts "----------------"
			    puts "--> Warning! <--"
			    puts "----------------"
			    puts "--> The amount of charges on a chain must be constant!"
			    puts "--> (Found chains with $tmp_CPP and $tmp_CPP_k)"
			    puts "Aborting...\n"
			    exit
			}
		    }
		    if { $tmp_MPC_k==$tmp_MPC } {
			lappend polymer_chains $chains
			incr tmp_NP
		    } elseif { $tmp_MPC==0 && $tmp_NP==0 } {
			set tmp_MPC $tmp_MPC_k
			lappend polymer_chains $chains
			incr tmp_NP
		    } else {
			puts "Failed.\n"
			puts "----------------"
			puts "--> Warning! <--"
			puts "----------------"
			puts "--> The number of monomers on all chains must be constant!"
			puts "--> (Found chains with $tmp_MPC and $tmp_MPC_k)"
			puts "Aborting...\n"
			exit
		    }
		}
	    }
	}
	set N_CPP $tmp_CPP
	set chargeDist $tmp_cD
	set N_P $tmp_NP
	set MPC $tmp_MPC
	set N_CR $tmp_NCR
        puts -nonewline ".,. reconstructed $N_P chains with $MPC monomers each ($N_CPP charged with distance $chargeDist"
	if { $chargeDist>0 } {
	    puts -nonewline ", like  "
	    for {set j 0} {$j<3} {incr j} {
		puts -nonewline "c"
		for {set i 0} {$i<[expr $chargeDist-1]} {incr i} { puts -nonewline "." }
	    }   
            puts -nonewline "c  etc."
	}
	puts "), found $N_CR crosslinks."
    }
    puts "    Analysis completed."


#  Derive remaining particle informations
#############################################################

    # Determine N_CPP, chargeDist, val_CI, N_Salt, val_pS, val_nS
    puts "    Analyzing particle data..."

    # Determine val_CI, N_Salt, val_pS, val_nS
    puts -nonewline "        Determining valency of counter-ions and salt-molecules... "
    flush stdout
    set tmp_CI 0
    set tmp_valCI 0
    set tmp_Salt 0
    set tmp_pS 0
    set tmp_nS 0
    set old_i 0
    set old_e [llength $part_else]
    set old [expr 10*$old_e]

    # loop all remaining particles with no bonds
    foreach i $part_else {
	# draw a status bar for every 10% of progress
	if { [expr $old_i*100]>=$old } {
	    set old [expr $old + 10*$old_e]
	    puts -nonewline ".,."
	    flush stdout
	}
	incr old_i

	set tmp_if [lindex [lindex $part_all $i] [findPropPos [lindex $part_all $i] type]]
	if { $tmp_if==$type_CI } {
	    # it's a counter-ion
	    set tmp_var [lindex [lindex $part_all $i] [findPropPos [lindex $part_all $i] q]]
	    if { $tmp_CI==0} {
		set tmp_valCI $tmp_var
	    } elseif { $tmp_valCI!=$tmp_var } {
		puts "Failed.\n"
		puts "----------------"
		puts "--> Warning! <--"
		puts "----------------"
		puts "--> The valency of the counter-ions must be constant!"
		puts "--> (Found some with $tmp_valCI and others with $tmp_var)"
		puts "Aborting...\n"
		exit
	    }
	    incr tmp_CI 
	} elseif { $tmp_if==$type_S } {
	    # it's a salt-molecule
	    set tmp_var [lindex [lindex $part_all $i] [findPropPos [lindex $part_all $i] q]]
	    if { $tmp_var>0 } { set tmp_S tmp_pS } else { set tmp_S tmp_nS }
	    if { [eval expr $$tmp_S==0] } {
		eval set $tmp_S $tmp_var
	    } elseif { [eval expr $$tmp_S!=$tmp_var] } {
		puts "Failed.\n"
		puts "----------------"
		puts "--> Warning! <--"
		puts "----------------"
		puts "--> The valency of the salt ions must be constant!"
		eval puts "--> (Found some with $$tmp_S and others with $tmp_var)"
		puts "Aborting...\n"
		exit
	    }
	    incr tmp_Salt 
	} else {
	    # neither salt nor counter-ion but having no bonds at all??? => hmm... something's wrong?!
	    puts "Failed.\n"
	    puts "----------------"
	    puts "--> Warning! <--"
	    puts "----------------"
	    puts -nonewline "--> Found unbound particle $i which is neither counter-ion ($type_CI) "
	    puts "nor salt ($type_S), but of type $tmp_if!"
	    puts "--> Confused: What shall I do with it?"
	    puts "--> (Stats: [lindex $part_all $i])"
	    puts "Aborting...\n"
	    exit
	}
    }
    set N_CI $tmp_CI
    set val_CI $tmp_valCI
    set N_Salt $tmp_Salt
    set val_pS $tmp_pS
    set val_nS [expr (-1)*$tmp_nS]
    puts -nonewline ".,. found $N_CI counter-ions with valency $val_CI and $N_Salt salt-molecules"
    if { $val_nS>0 || $val_pS>0} {
	puts -nonewline " with valency "
	for {set k 0} {$k<$val_nS} {incr k} {puts -nonewline "-"}
	puts -nonewline " or "
	for {set k 0} {$k<$val_pS} {incr k} {puts -nonewline "+"}
    }
    puts "."


    # So much for the particle data
    puts "    Analysis completed."


#  Write output file
#############################################################

    # Open output-file
    puts -nonewline "    Preparing converted configurations to be written to '[lindex [split $destination \"/\"] end]'... "
    flush stdout
    set out [open $destination "w"]
    puts "Done."

    # Now write all we've done so far to a Deserno-compatible file
    puts "    Saving all configurations in Deserno-format to disk... "
    puts -nonewline "        Writing file header... "
    flush stdout
    puts $out "# Praefix-String......................: $prefix"
    puts $out "# Postfix-Nummer......................: $postfix"
    puts $out "# Zufallszahlen-Seed..................: $seed"
    puts $out "# Physikalische Zeit..................: [expr $startTime+$step*$deltaT]"
    puts $out "# Physikalische Endzeit...............: $endTime"
    puts $out "# Diskretisierungs-Schritt............: $deltaT"
    puts $out "# Anzahl der Integrationsschritte.....: $integrationSteps"
    puts $out "# Messergebnisse herausschreiben......: $saveResults"
    puts $out "# Konfigurationen herausschreiben.....: $saveConfig"
    puts $out "# Anzahl der Polymere.................: $N_P"
    puts $out "# Ladungen pro Polymer................: $N_CPP"
    puts $out "# Ladungsabstand......................: $chargeDist"
    puts $out "# Valenz der Polymer-Gegenionen.......: $val_CI"
    puts $out "# Anzahl der Salzmolekuele............: $N_Salt"
    puts $out "# Valenz der positiven Salzionen......: $val_pS"
    puts $out "# Valenz der negativen Salzionen......: $val_nS"
    puts $out "# Teilchenzahl........................: $N_T"
    puts $out "# Boxlaenge...........................: $boxLen"
    puts $out "# Subboxen pro Boxlaenge..............: $subbox_1D"
    puts $out "# Skin................................: $skin"
    puts $out "# Ortsraum-Cutoff.....................: $rcut_P3M"
    puts $out "# k-FENE..............................: $k_FENE"
    puts $out "# r-FENE..............................: $r_FENE"
    puts $out "# epsilon-LJ..........................: $eps_LJ"
    puts $out "# LJ-Cutoff...........................: $rcut_LJ"
    puts $out "# alpha...............................: $alpha"
    puts $out "# FFT-Mesh............................: $mesh"
    puts $out "# Gitter-ip...........................: $ip"
    puts $out "# Bjerrum-Laenge......................: $Bjerrum"
    puts $out "# Temperatur..........................: $Temp"
    puts $out "# Gamma...............................: $Gamma"
    puts "Done."

    # Note the format requirements:
    # 1. Monomers               MMMM MMMM MMMM MMMM | ++++ ++++ ++++ | ++++ ++++ ++++ | ---- ---- ---- 
    # 2. Counter-ions           First, N_P polymer  |     Second,    |     Third,     |      Last,     
    # 3. Positive Salt            chains with MPC   |      N_CI      |  val_pS*N_Salt |  val_nS*N_Salt 
    # 4. Negative Salt             monomers each    |  counter-ions  |  positive salt |  negative salt 

    # Write out the polymer chains
    set old_i 0
    set old_e [expr $N_P*$MPC]
    set old [expr 10*$old_e]
    puts -nonewline "        Writing $N_P polymer chains with $MPC monomers and $N_CPP charges each... "
    flush stdout
    foreach i $polymer_chains {
	# The Deserno-file-format requires that first the monomer's data has to be written in the order of the chains
	foreach j $i {
	    # draw a status bar for every 10% of progress
	    if { [expr $old_i*100]>=$old } {
		set old [expr $old + 10*$old_e]
		puts -nonewline ".,."
		flush stdout
	    }
	    incr old_i

	    # write monomer
	    set tmp_part [lindex $part_all $j]
	    set tmp_var [findPropPos $tmp_part pos]
	    puts -nonewline $out "[lindex $tmp_part $tmp_var]\t"
	    puts -nonewline $out "[lindex $tmp_part [expr $tmp_var+1]]\t"
	    puts -nonewline $out "[lindex $tmp_part [expr $tmp_var+2]]\t"
	    set tmp_var [findPropPos $tmp_part v]
	    puts -nonewline $out "[lindex $tmp_part $tmp_var]\t"
	    puts -nonewline $out "[lindex $tmp_part [expr $tmp_var+1]]\t"
	    puts -nonewline $out "[lindex $tmp_part [expr $tmp_var+2]]\t"
	    set tmp_var [findPropPos $tmp_part q]
	    puts -nonewline $out "[lindex $tmp_part $tmp_var]\n"
	}
    }
    puts ".,. completed."

    # Write out any remaining counter-ions or salt molecules
    puts -nonewline "        Writing [expr $N_CPP*$N_P/$val_CI] counter-ions and $N_Salt salt molecules... "
    flush stdout
    set old_i 0
    set old_e [expr 3*[llength $part_else]]
    set old [expr 10*$old_e]
    for {set j 0} {$j<3} {incr j} {
	# The Deserno-file-format requires that after the monomer's datas
	# first the counter-ions, then the positive salt ions, and last the negative salt ions
	# have to be written; hence we'll pass this loop three times... *sigh*
	foreach i $part_else {
	    # draw a status bar for every 10% of progress
	    if { [expr $old_i*100]>=$old } {
		set old [expr $old + 10*$old_e]
		puts -nonewline ".,."
		flush stdout
	    }
	    incr old_i
	    
	    # write particles
	    set tmp_part [lindex $part_all $i]
	    set tmp_if [lindex $tmp_part [findPropPos $tmp_part type]]
	    set tmp_fi [lindex $tmp_part [findPropPos $tmp_part q]]
	    set tmp_d 0
	    if { $j==0 && $tmp_if==$type_CI } { set tmp_d 1 }
	    if { $j==1 && $tmp_if==$type_S && $tmp_fi>0 } { set tmp_d 1 }
	    if { $j==2 && $tmp_if==$type_S && $tmp_fi<0 } { set tmp_d 1 }
	    if { $tmp_d==1 } {
		set tmp_var [findPropPos $tmp_part pos]
		puts -nonewline $out "[lindex $tmp_part $tmp_var]\t"
		puts -nonewline $out "[lindex $tmp_part [expr $tmp_var+1]]\t"
		puts -nonewline $out "[lindex $tmp_part [expr $tmp_var+2]]\t"
		set tmp_var [findPropPos $tmp_part v]
		puts -nonewline $out "[lindex $tmp_part $tmp_var]\t"
		puts -nonewline $out "[lindex $tmp_part [expr $tmp_var+1]]\t"
		puts -nonewline $out "[lindex $tmp_part [expr $tmp_var+2]]\t"
		set tmp_var [findPropPos $tmp_part q]
		puts -nonewline $out "[lindex $tmp_part $tmp_var]\n"
	    }
	}
    }
    puts ".,. completed."

    # Write out cross-links (if any)
    if { $N_CR>0 } {
	puts -nonewline "        Writing $N_CR cross-links... "
	flush stdout
	puts $out "# Number of Crosslinks.........: N_CR = $N_CR"
	set old_i 0
	set old_e $N_CR
	set old [expr 10*$old_e]
	foreach i $cross {
	    # draw a status bar for every 10% of progress
	    if { [expr $old_i*100]>=$old } {
		set old [expr $old + 10*$old_e]
		puts -nonewline ".,."
		flush stdout
	    }
	    incr old_i
	    
	    # write cross-links
	    puts $out "[lindex $i 0]\t[lindex $i 1]"
	}
	puts ".,. completed."
    }
	
    # Done
    puts "    Saving completed."


    # Now finish up by closing the file
    puts -nonewline "    Closing output file... "
    flush stdout
    close $out
    puts "Done."
    puts "    Function successfully completed. Returning control to calling script..."
}
