#############################################################
#                                                           #
# aux.tcl                                                   #
# =======                                                   #
#                                                           #
# Several additional auxiliary functions for tcl_md.        #
#                                                           #
# Created:       01.10.2002 by BAM                          #
# Last modified: 08.11.2002 by BAM                          #
#                                                           #
#############################################################




#
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
# Created:       01.10.2002 by BAM
# Last modified: 01.10.2002 by BAM
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




#
# polyBlockWrite
# --------------
# 
# Writes current polymer configuration to disk,
# including all bonds and interactions there are,
# using Axa's blockfile-format.
#
# Input:
# - complete path 'destination';
#   if the filename ends with '.gz', the file will be compressed
# - a list of 'tcl_md' parameters to be saved (out of node_grid|box_l|niatypes|time_step|skin|gamma|bjerrum|...
#   ...p3m_alpha|p3m_r_cut|p3m_mesh|p3m_cao|p3m_epsilon|p3m_mesh_offset|max_num_cells|periodicity);
#   if an empty list '{}' is supplied, no parameters are written
# - a string containing which informations (out of pos|type|q|v|f) on the particles should be saved to disk;
#   if an empty string ' "" 'is provided, no particles, no bonds, and no interactions are written
# 
# Created:       08.11.2002 by BAM
# Last modified: 08.11.2002 by BAM
# 
#############################################################

proc polyBlockWrite { destination write_param write_part } {
    
    # Open output-file - compressed, if desired
    if { [string compare [lindex [split $destination "."] end] "gz"]==0 } {
	set f [open "|gzip -c - >$destination" w]
    } else {
	set f [open "$destination" "w"]
    }
    
    # Write parameters, if desired
    if { "$write_param" != "{}" } {
	foreach j $write_param {
	    blockfile $f write variable $j
	}
    }
    
    # Write particles, bonds, and interactions, if desired
    if { "$write_part" != "" } {
	blockfile $f write particles $write_part all
	blockfile $f write bonds all
	blockfile $f write interactions
    }
    
    # Close file
    close $f
}
