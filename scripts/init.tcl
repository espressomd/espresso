#############################################################
#                                                           #
# init.tcl                                                  #
# ========                                                  #
#                                                           #
# Used for initializing tcl-scripts powered by the 'tcl_md' #
# MPI parallel molecular dynamics program package.          #
#                                                           #
# It is a very good idea to include all custom tcl-scripts  #
# here which may enhance other programs running 'tcl_md'    #
# by providing additional tcl-commands/-functions.          #
#                                                           #
# Created:       20.09.2002 by AxA                          #
# Last modified: 25.09.2002 by BAM                          #
#                                                           #
#############################################################


# here do everything you want to do upon initialization.
# e. g.
puts ""
puts "*******************************************************"
puts "*                                                     *"
puts "*                    - tcl_md -                       *"
puts "*                      ======                         *"
puts "*      A MPI Parallel Molecular Dynamics Program      *"
puts "*                                                     *"
puts "*                                                     *"
puts "* (c) 2002                                            *"
puts "* Max-Planck-Institute for Polymer Research           *"
puts "* Mainz, Germany                                      *"
puts "*                                                     *"
puts "*******************************************************"
puts ""


# Include useful tcl-scripts providing new functions etc.

source convertDeserno.tcl
# adds 'convertDeserno2MD' & 'convertMD2Deserno' for converting particle configurations
#       stored in Deserno-file-format into tcl_md format and vice versa
# adds 'convertDeserno2MDmain' & 'convertMD2DesernoMain' for directly accessing the conversion scripts
#       therefore bypassing and/or overriding the initialization procedure

source countBonds.tcl
# adds 'countBonds' returning a tcl-script with complete bonding-informations on any particle
# adds 'findPropPos' to determine the index-position of a certain property in the output of 'part'
# adds 'findBondPos' to do the same for property 'bonds'

source blockfile_support.tcl