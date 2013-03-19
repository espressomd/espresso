#############################################################
#                                                           #
# pe_network.tcl                                            #
# ==============                                            #
#                                                           #
# Sets up a polymer network either from scratch or from a   #
# previously created polymer configuration, relaxing it to  #
# crosslinked equilibrium, integrating afterwards.          #
#                                                           #
#############################################################
#
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

puts " "
puts " "
puts "============================================================="
puts "=                         Welcome!                          ="
puts "============================================================="
puts " "
puts "Starting 'pe_network.tcl'..."
puts " "
set section 1
puts "($section) Initializing:"
set prg_step 1



#############################################################
#  Parameters                                               #
#############################################################

puts -nonewline "    Applying parameters... "


# System parameters
#############################################################

set N_P 20
   # number of polymers
set MPC 106
   # number of monomers per polymer chain (= its length)
set bond_l 1.0
   # length of the bonds between monomers
set cM_dist 2
   # distance between charged monomers
set val_cM 1.0
   # valency of charged monomers

set val_CI -1
   # valency of counterions
set N_Salt 0
   # number of salt molecules
set val_pS 2
set val_nS -3
   # valency of the positive (pS) and negative (nS) salt-ions, respectively

set density 0.5
   # particle density of the system (from which the box_l will be derived)


# System setup options
#############################################################

set poly_mode SAW
set shield 0.1
   # $poly_mode=SAW uses a Self-Avoiding Random Walk to place the monomers, observing an excluded volume $shield around each particle 
   # $poly_mode=RW  uses a simple random walk without any additional fancy stuff
set max_try 30000
   # how many trials to attempt to place a particle, crosslink the network, etc.

set v_max 0
   # determines the maximum absolute velocity assigned to the particles initially

set r_catch 1.9
   # the radius around each end monomer in which possible binding partners have to be before crosslinking could occur
# if set to '0', no crosslinking is attempted (note that the equilibration integration will nevertheless take place)
set distLink 2
   # to prevent the polymer chains to become too stiff, the difference between two crosslinked monomers' indices has to be >= $distLink
set distChain $MPC
   # same as $distLink, but for the minimum distance between the chain's end and a crosslinked monomer
   # if >=$MPC, no end monomer may bind to a monomer on the same chain


# Interaction parameters
#############################################################

set LJ_eps       1.0
set LJ_sigma     1.0
set LJ_cut       1.122462048
set LJ_shift     "auto"
set LJ_offset    0

set fene_k       7.0
set fene_r       2.0

set bjerrum      1.0
set p3m_accuracy 1e-3


# Integration parameters
#############################################################

# parameters for initial warm-up process with capped LJ-interactions
set ljcp_dt 0.010
   # time step to be used during initial warm-up
set ljcp_cap1 10.0
set ljcp_incr 2.0
set ljcp_step 50
set ljcp_loop 100
   # force magnitude at which the LJ-potential is to be capped (1) initially, increasing it by $ljcp_incr for $ljcp_loop cycles of $ljcp_step integrations
set ljcp_dist 0.95
   # desired minimum distance [analyze mindist] during the warm-up process (in units of LJ_cut); if reached, the warm-up is considered done
set ljcp_temp 1.0
set ljcp_gamma 5.0
   # temperature friction coefficient to be used during the warm-up (may be even up to 1e4 if desired)

# parameters for the equilibration process using full LJ-interactions (and no electrostatics)
set equl_dt 0.010
   # time-step to be used during equilibration
set equl_step 50
set equl_loop 100
   # specifies how many more integration-steps will be derived with full LJ-interactions to equilibrate the 'ljcp-ed' configuration above
set equl_zero 1
   # if '1' all particles' velocities and forces will be set to zero before the equilibration process commences
set equl_stop 1
   # if '1' the equilibration integration will be aborted once all 2*$N_P end monomers have been successfully crosslinked
set equl_temp 1.0
set equl_gamma 1.0
   # temperature and friction coefficient to be used during the equilibration process

# parameters for the main integration process using all interactions (including electrostatics)
set main_dt 0.010
   # time-step to be used during the main integration
set main_step 50
set main_loop 100
   # specifies how many more integration-steps will be derived with all interactions on the equilibrated configuration from above
set main_zero 1
   # if '1' all particles' velocities and forces will be set to zero before the main integration process commences
set main_temp 1.0
set main_gamma 1.0
   # temperature and friction coefficient to be used during the equilibration process


# Analysis options & I/0-options
#############################################################

set stat {mindist re rg rh g123}
   # sets which mode of analysis should be performed (analysis skipped if 'stat' is empty):
   # 'mindist' returns the minimum distance between all particles
   # 're'      returns the end-to-end-distance averaged over all polymers
   # 'rg'      returns the radius of gyration averaged over all chains
   # 'rh'      returns the hydrodynamic radius 
   # 'g123'    returns the mean-square displacement g1(t) of a monomer, the mean-square displacement g2(t) of in the center of gravity of the chain itself,
   #                   the motion of the center of mass g3(t)as a tcl-list {g1(t) g2(t) g3(t)}

set vmd_output yes
   # 'yes' writes out .pdf- and .psf-files to which vmd can be connected

set write movie
   # 'finish' only writes out the final state as a tcl-blockfile
   # 'yes' also writes each configuration during integration after '$write_loop_xxxx' cycles
   # --> Requires directory '$output_path' to store everything!
   # 'movie' takes those configurations and converts them to a 'poly2pdb'-movie for 'vmd'
set write_loop_ljcp 1
set write_loop_equl 1
set write_loop_main 1
   # after this many integration-cycles (out of $ljcp_loop/$equl_loop) the current configuration may be saved if{$write=="yes" || $write=="movie"}
set write_param all
   # sets which tcl-parameters (out of { box_l cell_grid cell_size gamma local_box_l max_cut max_num_cells max_part max_range max_skin n_nodes n_part n_part_types node_grid periodicity skin temperature time time_step transfer_rate verlet_flag verlet_reuse } or "all") should be saved to disk
set write_part "id pos type q v f"
   # sets which informations (out of pos|type|q|v|f) on the particles should be saved to disk
   # note that the newer versions of the blockfile-format do NOT allow to save bonds here, 
   # but rather require to do that separately by calling 'blockfile_write_bonds'


# Filenames & paths
#############################################################

set input_path  ./configs/
set input_file    simulation_input1-WARM5000

set output_path ./configs/
set output_prfx   simulation_output

set vmd_path    ./movie/
set vmd_prfx      simulation_vmd

set movie_path  ./movie/
set movie_prfx    simulation_movies

set stat_path   ./configs/
set stat_prfx     simulation_stats


# Espresso parameters
#############################################################

setmd skin 1.333333
   # specifies the skin around each particle
# setmd node_grid x y z
   # number of nodes in all 3D
   # 3D node grid for real space domain decomposition (optional, if unset an optimal set is chosen automatically).
# setmd max_cut r
   # max. real space interaction cutoff
# setmd periodic 1 1 1
   # sets period boundaries in all 3D
   # READ-ONLY???
setmd max_num_cells 512
   # maximal number of cells for the link cell algorithm
   # reasonable values are between 125 and 1000, or for some problems (n_total_particles/n_nodes)

set type_FENE 0
   # Determine which type-number should be assigned to the interactions
set type_nM 0
set type_cM 1
set type_CI 2
set type_pS 3
set type_nS 4
set types [list $type_nM $type_cM $type_CI $type_pS $type_nS]
   # Determine which type-number should be assigned to neutral/charged monomers, counter-ions, salt-ions, ...
set max_type [llength $types]
   # specifies how many particle types are present (e.g. 5 when using neutral/charged monomers, counterions, and positive/negative salt-molecules)
   # should correspond to '[setmd n_part_types]'


# Other parameters
#############################################################

set tcl_precision 5
   # tells tcl which precision to use for formating the output
set random_seeds {  }
   # seeds to be used to initialize the random number generator on every node (give at least 'n_nodes' seeds, exceeding ones are ignored)
   # (if empty, the default seed 'seed = (10*this_node+1)*1103515245 + 12345; seed = (seed/65536) % 32768;' is used)


puts "Done."



#############################################################
#  Preparations                                             #
#############################################################

puts "    Preparing environement... "

# Check in/out directories
#############################################################

proc makedirexist {path} {
    if { ![file exists $path] } {
	file mkdir $path
    } {
	if { ![file isdir $path]} {
	    puts "$path has to be a directory. Remove it or change the path setup"
	    exit
	}
    }
}

makedirexist $output_path
makedirexist $movie_path
makedirexist $vmd_path


# Clean up fragments of previous runs
#############################################################

# clean up directories, if something is to be written later on
if { $write  == "yes" || $write == "finish" || $write == "movie" } {
    puts -nonewline "        Cleaning saved configurations '$output_prfx\*' out of $output_path... "; flush stdout
    foreach i [glob -nocomplain "$output_path$output_prfx*"] { eval exec rm -f $i }
    puts "Done."
    puts -nonewline "        Cleaning saved vmd-files '$vmd_prfx\*' out of $vmd_path... "; flush stdout
    foreach i [glob -nocomplain "$vmd_path$vmd_prfx*"] { eval exec rm -f $i }
    puts "Done."
    puts -nonewline "        Cleaning saved movies '$movie_prfx\*' out of $movie_path... "; flush stdout
    foreach i [glob -nocomplain "$movie_path$movie_prfx*"] { eval exec rm -f $i }
    puts "Done."
}

# clean up old statistics, prepare new ones
if { [llength $stat]!=0 } {
    puts -nonewline "        Cleaning saved statistics '$stat_prfx\*' out of $stat_path... "; flush stdout
    foreach i [glob -nocomplain "$stat_path$stat_prfx*"] { eval exec rm -f $i }
    puts "Done."
    puts -nonewline "        Open file $stat_prfx.dat for statistical output... "; flush stdout
    set stat_out [open "$stat_path$stat_prfx.dat" "w"]
    puts "  Done."
}


# Internal variables
#############################################################

puts -nonewline "        Setting up internal variables... "; flush stdout
set simtime 0
   # how much time has been simulated so far
set polyConfAux   [list $output_path $output_prfx \{$write_param\} \{$write_part\} $movie_path $movie_prfx]
   # auxiliary variables abbreviating often needed parameters for file-I/O ('$N_P $MPC $N_CI $N_pS $N_nS' are added further down)
puts "Done."

puts "    Preparation completed."



#############################################################
#  Setup the system                                         #
#############################################################

puts " "
incr section
puts "($section) Setting up the system:"


# Complete parameter set
#############################################################

puts -nonewline "    Processing parameter set... "
flush stdout
set N_M    [expr ($N_P * $MPC)]
   # total number of monomers
set N_CPP  [expr int(($MPC + $cM_dist -1) / $cM_dist)]
   # charges per polymer chain
set N_CI   [expr -int($N_CPP * $N_P / $val_CI)]
   # number of counter-ions
set N_pS   [expr $N_Salt*$val_pS] 
set N_nS   [expr $N_Salt*(-1)*$val_nS] 
   # number of positively and negatively charged salt ions, respectively
set N_tS   [expr $N_pS+$N_nS] 
   # total number of salt ions 
   # (not to be confused with N_Salt, which is the number of salt _molecules_)
set N_T    [expr $N_M + $N_CI + $N_tS]
   # total number of particles
set polyConfAuxL   $polyConfAux;   lappend polyConfAuxL   $N_P $MPC $N_CI $N_pS $N_nS

set box_l  [expr pow(($N_T / $density),(1.0/3.0))]
setmd box_l $box_l $box_l $box_l
   # size of the simulation box

puts "Done (expecting $N_T particles of density $density in a box of length $box_l)."


# Consistency checks
#############################################################

puts -nonewline "    Checking your input data for consistency... "
flush stdout
set tmp_chk 0
if { [expr $N_CPP % $val_CI]!=0 } {
    puts "Failed.\n    The number of charges on a chain ($N_CPP) has to be a multiple of the counterions' valency ($val_CI)!\nAborting...\n"
    exit }
puts "Passed."

# Random number generator setup
#############################################################

# Initialize random generators
if { [llength $random_seeds] > 0 } {
    puts -nonewline "    Setting seed of random number generator on every node to $random_seeds... "; flush stdout
    puts "Done ([eval t_random seed $random_seeds])."
} else { puts "        Using default seeds for the random number generator on every node ([t_random seed])." }


# Interaction setup
#############################################################

puts "    Creating interactions..."

# fene
inter $type_FENE fene $fene_k $fene_r
puts "        [inter $type_FENE]"

# pairwise lennard_jones for all particles
foreach i $types {
    foreach j $types {
	inter $i $j lennard-jones $LJ_eps $LJ_sigma $LJ_cut $LJ_shift $LJ_offset
	puts "        [inter $i $j]"
    }
}
puts "    Interaction setup complete."


# Particle setup
#############################################################

# Open start configuration or create new one
if { [file exists "$input_path$input_file" ] } {
    set inp [open "$input_path$input_file" r]; set warmup "$input_file"
} elseif { [file exists "$input_path$input_file.gz" ] } {
    set inp [open "| gzip -cd $input_path$input_file.gz" r]; set warmup "$input_file.gz"
} else { set warmup 1 }
if { $warmup != 1 } {
    puts -nonewline "    Found file '$warmup' (currently reading it... "; flush stdout
    while { [blockfile $inp read auto] != "eof" } {}
    close $inp
    puts -nonewline "done) with [setmd n_part] particles; "; flush stdout
    if { $N_T != [setmd n_part] } {
	puts "however, configuration differs from this setup (expected $N_T particles).\n    Creating new start configuration now..."
	puts -nonewline "        First, erasing all previously stored particles (i.e. those just read) from Espresso... "; flush stdout
	puts "Done ([part deleteall])."
	set warmup 1
    } else {
	puts "will continue using these, hence skipping warm-up."
	set warmup 0
    }
} else {
    puts "    No previously saved start configuration found; creating a new one now..."
    set warmup 1
}
if { $warmup==1 } {
    puts -nonewline "        Creating $N_P polymers with $MPC monomers (of which $N_CPP have charge $val_cM) each... "; flush stdout
    set tmp_out [polymer $N_P $MPC $bond_l mode $poly_mode $shield $max_try charge $val_cM distance $cM_dist types $type_nM $type_cM FENE $type_FENE]
    puts "Done (retried placing a monomer no more than $tmp_out times)."
    puts -nonewline "        Creating $N_CI counterions with charge $val_CI each... "; flush stdout
    puts "Done (retried placing a counterion no more than [counterions $N_CI mode $poly_mode $shield $max_try charge $val_CI type $type_CI] times)."
    puts -nonewline "        Creating $N_pS positive ($val_pS) and $N_nS negative ($val_nS) salt molecules... "; flush stdout
    set tmp_out [salt $N_pS $N_nS mode $poly_mode $shield $max_try charges $val_pS $val_nS types $type_pS $type_nS]
    puts "Done (retried placing a salt particle no more than $tmp_out times)."
    puts -nonewline "        Setting all $N_T particles' initial velocity to a random vector of length \[0, $v_max\]... "; flush stdout
    puts "Done (average length [expr [velocities $v_max]/$N_T])."
    puts "    New start configuration successfully created."
}

# Do some initial analysis
analysisInit $stat $stat_out $N_P $MPC $simtime


#  Prepare IMD real-time visualization connection
# (refer to http://www.mpip-mainz.mpg.de/www/theory/computers/visualization/IMD.html for further informations)
#############################################################

if { $vmd_output=="yes" } {
    puts -nonewline "        Preparing connection of real-time visualization tool 'imd'... "; flush stdout
    writepsf $vmd_path$vmd_prfx.psf $N_P $MPC $N_CI $N_pS $N_nS
    writepdb $vmd_path$vmd_prfx.pdb
    puts -nonewline "Output created, establishing link... "; flush stdout
    for {set port 10000} { $port < 65000 } { incr port } {
	catch {imd connect $port} res
	if {$res == ""} break
    }
    if { $port==65000 } { puts "Failed." } else { puts "Done (now listening at port $port)." 
	puts "        What you have to do now for a VMD connection:"
	puts "        (1) Start vmd in directory $vmd_path (best before running the script)."
	puts "        (2) Enter on vmd command line: 'mol load psf $vmd_prfx.psf pdb $vmd_prfx.pdb'"
	set HOSTNAME [exec hostname]
	puts "        (3) Enter on vmd command line: 'imd connect $HOSTNAME $port'"
	puts "        (4) To have the chains coloured individually, set 'Coloring-Method' to 'ResName' in the 'Graphics'-menu"
	imd listen 0
    }
}

puts "   System setup completed."



#############################################################
#  Warm-up Integration                                      #
#############################################################

puts " "
incr section
puts "($section) Warm-up integration:"

if {$max_type!=[setmd n_part_types]} {
    puts "Failed.\n\n----------------\n--> Warning! <--\n----------------"
    puts "--> The number of particle types you specified differs from the one known to Espresso!"
    puts "--> (You said: $max_type, but Espresso has [setmd n_part_types])\nAborting..."; exit
}
if { $warmup == 1 } {
    puts -nonewline "    Prepare environement for the warm-up process ("; flush stdout
    thermostat langevin $ljcp_temp $ljcp_gamma
    puts -nonewline "temperature T=[setmd temp], friction gamma=[setmd gamma], "; flush stdout
    setmd time_step $ljcp_dt; puts -nonewline "dt=[setmd time_step]"; flush stdout
    puts ")... "
    puts -nonewline "        Capping all LJ-interactions... "; flush stdout
    set tmp_cap $ljcp_cap1; inter forcecap $tmp_cap; puts "Done (cap initially at $tmp_cap)."
    puts -nonewline "        Making sure that electrostatics are disabled... "; flush stdout
    inter coulomb 0; puts "Done ([inter coulomb])."
    puts "    Preparation completed."

    polyConfMovWriteInit $write "1-WARM" $polyConfAuxL

    puts "    Initiating [expr $ljcp_loop*$ljcp_step] warm-up integrations in $ljcp_loop cycles with capped LJ-interactions and no electrostatics..."
    set tmp_digit [expr int(log10($ljcp_loop*$ljcp_step))+1]
    for {set i 0} { $i < $ljcp_loop } { incr i} {
	integrate $ljcp_step
	set simtime [expr $simtime + [setmd time_step]]; set tmp_dist [analyze mindist]
	puts -nonewline "        Step [expr ($i+1)*$ljcp_step]/[expr $ljcp_loop*$ljcp_step] (t=$simtime): "; flush stdout
	puts -nonewline "LJ's cap = $tmp_cap"; flush stdout

	analysis $stat $stat_out $N_P $MPC $simtime
	if { [expr ($i+1) % $write_loop_ljcp == 0] } {
	    polyConfMovWrite $write "1-WARM" $tmp_digit [expr ($i+1)*$ljcp_step] $polyConfAux
	    puts -nonewline "...\r"; flush stdout } else { puts -nonewline "...                \r"; flush stdout }
	if { $vmd_output=="yes" } { imd positions }

	if { $tmp_dist >= $ljcp_dist } { break }
	inter forcecap $tmp_cap; set tmp_cap [expr $tmp_cap + $ljcp_incr]
    }
    puts -nonewline "\n    Initial relaxation completed"
    if { [expr $i % $write_loop_ljcp != 0] || ($tmp_dist >= $ljcp_dist) } {
        polyConfMovWrite $write "1-WARM" $tmp_digit [expr $i*$ljcp_step] $polyConfAux
    } else { puts -nonewline " (warmed up)" }
    puts "."
} else { puts "    Skipped." }



#############################################################
#  Equilibrating & Crosslinking                             #
#############################################################

puts " "
incr section
puts "($section) Equilibrating & Crosslinking:"

if { $equl_loop > 0 } {
    puts -nonewline "    Prepare environement for the equilibration process ("; flush stdout
    thermostat langevin $equl_temp $equl_gamma
    puts -nonewline "temperature T=[setmd temp], friction gamma=[setmd gamma], "; flush stdout
    setmd time_step $equl_dt; puts -nonewline "dt=[setmd time_step]"; flush stdout
    puts ")... "
    puts -nonewline "        Remove capping of LJ-interactions... "; flush stdout
    inter forcecap 0; puts "Done ([inter forcecap])."
    puts -nonewline "        Making sure that electrostatics are disabled... "; flush stdout
    inter coulomb 0; puts "Done ([inter coulomb])."

    if { $equl_zero==1 } {
      kill_particle_motion
      kill_particle_forces
    } else { 
      puts "        As requested, all particles' velocities and forces remain at their current values!"
    }
    
    # Attempt to crosslink the polymers (if requested)
    set tmp_cl1 0; set tmp_cl2 [expr 2*$N_P]
    if { $r_catch > 0 } {
	puts -nonewline "        Connecting chains ends to random monomers <$r_catch away (if >$distChain apart and >$distLink from next link)... "
	set tmp_cl1 [crosslink $N_P $MPC catch $r_catch distLink $distLink distChain $distChain FENE $type_FENE trials $max_try]
	puts "Done (now $tmp_cl1 of $tmp_cl2 ends are crosslinked)."
    } else { puts "        As requested, crosslinking has been disabled." }
    puts "    Preparation completed."
    
    polyConfMovWriteInit $write "2-EQUL" $polyConfAuxL
    
    puts "    Equilibrating for [expr $equl_loop*$equl_step] integrations in $equl_loop cycles with full LJ-interactions but no electrostatics..."
    set tmp_digit [expr int(log10($equl_loop*$equl_step))+1]; 
    for {set i 0} { $i < $equl_loop } { incr i} {
	integrate $equl_step
	set simtime [expr $simtime + [setmd time_step]]
	puts -nonewline "        Step [expr ($i+1)*$equl_step]/[expr $equl_loop*$equl_step] (t=$simtime): "; flush stdout
	puts -nonewline "Currently $tmp_cl1/$tmp_cl2 ends crosslinked"; flush stdout
	
	analysis $stat $stat_out $N_P $MPC $simtime
	if { [expr ($i+1) % $write_loop_equl == 0] } {
	    polyConfMovWrite $write "2-EQUL" $tmp_digit [expr ($i+1)*$equl_step] $polyConfAux
	    puts -nonewline "...\r"; flush stdout } else { puts -nonewline "...                \r"; flush stdout }
	if { $vmd_output=="yes" } { imd positions }
	
	if { ($tmp_cl1==$tmp_cl2) && ($equl_stop==1) } { break }
	if { ($r_catch > 0) && ($tmp_cl1 < $tmp_cl2) } { 
	    set tmp_cl1 [crosslink $N_P $MPC catch $r_catch distLink $distLink distChain $distChain FENE $type_FENE trials $max_try] }
    }
    puts -nonewline "\n    Equilibration completed"
    if { [expr $i % $write_loop_equl != 0] } { polyConfMovWrite $write "2-EQUL" $tmp_digit [expr $i*$equl_step] $polyConfAux
    } else { puts -nonewline " (done)" }; puts "."
    if { ($r_catch > 0) && ($tmp_cl1 < $tmp_cl2) } { 
	puts "    WARNING! Even after [expr $equl_loop*$equl_step] integration steps the system is not yet fully crosslinked ($tmp_cl1 of $tmp_cl2 ends)!"
	puts "             Recommend more a longer equilibrium process, a higher density than $density, and/or a bigger catching radius than $r_catch!"
    }
} else { puts "    Skipped." }



#############################################################
#  Main Integration                                         #
#############################################################

puts " "
incr section
puts "($section) Main Integration:"

if { $main_loop > 0 } {
    puts -nonewline "    Prepare environement for the main integration process ("; flush stdout
    thermostat $main_temp $main_gamma
    puts -nonewline "temperature T=[setmd temp], friction gamma=[setmd gamma], "; flush stdout
    setmd time_step $main_dt; puts -nonewline "dt=[setmd time_step]"; flush stdout
    puts ")... Done."
    puts -nonewline "        Making sure that capping of the LJ-interactions is removed... "; flush stdout; 
    inter forcecap 0; puts "Done ([inter forcecap])."
    if { $bjerrum > 0 } {
	puts -nonewline "        Activating p3m-subsystem & enabling full electrostatics... "; flush stdout
	inter coulomb $bjerrum p3m accuracy $p3m_accuracy; inter coulomb epsilon $p3m_epsilon
	puts "Done ([inter coulomb])."
    } else { inter coulomb $bjerrum; puts "        As requested, electrostatics have been disabled ([inter coulomb])!" }

    if { $main_zero==1 } {
      kill_particle_motion
      kill_particle_forces
    } else { 
      puts "        As requested, all particles' velocities and forces remain at their current values!"
    }

    puts "    Preparation completed."
    
    polyConfMovWriteInit $write "3-MAIN" $polyConfAuxL
    
    puts "    Main integrating for [expr $main_loop*$main_step] steps in $main_loop cycles with full interactions..."
    set tmp_digit [expr int(log10($main_loop*$main_step))+1]; 
    for {set i 0} { $i < $main_loop } { incr i} {
	integrate $main_step
	set simtime [expr $simtime + [setmd time_step]]
	puts -nonewline "        Step [expr ($i+1)*$main_step]/[expr $main_loop*$main_step] (t=$simtime)"; flush stdout
	
	analysis $stat $stat_out $N_P $MPC $simtime
	if { [expr ($i+1) % $write_loop_main == 0] } {
	    polyConfMovWrite $write "3-MAIN" $tmp_digit [expr ($i+1)*$main_step] $polyConfAux
	    puts -nonewline "...\r"; flush stdout } else { puts -nonewline "...                \r"; flush stdout }
	if { $vmd_output=="yes" } { imd positions }
    }
    puts -nonewline "\n    Main integration completed"
    if { [expr $i % $write_loop_main != 0] } { polyConfMovWrite $write "3-MAIN" $tmp_digit [expr $i*$main_step] $polyConfAux
    } else { puts -nonewline " (done)" }; puts "."
} else { puts "    Skipped." }



#############################################################
#  Some final thoughts                                      #
#############################################################

puts " "
incr section
puts "($section) Some final thoughts:"

puts "    Wrapping up current proceedings..."
if { $write == "finish" } { set $write yes }
puts -nonewline "    "; polyConfMovWriteInit $write "4-FINE" $polyConfAuxL
if { [llength $stat]!=0 } {
    puts -nonewline "        Closing analysis file $stat_prfx.dat... "; flush stdout; flush $stat_out; close $stat_out; puts "Done." 
}
if { $vmd_output=="yes" } {
    puts -nonewline  "        Closing connection to 'imd'... "; flush stdout
    if { "[imd disconnect]"=="no connection" } { puts "Done (no 'imd' was listening)." 
    } else { puts "Done (successfully disconnected 'imd')." }
}
puts "    Wrapped."
puts "\nThis simulation run is now complete. Thank you for using 'Espresso'!"
puts " "
puts "============================================================="
puts "=                         Finished.                         ="
puts "============================================================="
puts " "
