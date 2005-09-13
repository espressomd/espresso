# mbtools::utils -- 
#
# This package provides a set of routines for basic bookkeeping tasks
# that typically arise in a membrane simulation #
# Author: Ira
# 

package require ::mmsg 1.0.0
package provide ::mbtools::utils 1.0.0

namespace eval ::mbtools::utils {
    variable verbosity 0
    variable firstconfignum 0

    # The identity number for written out warmup configurations
    variable warmcfg 0

    variable maxnbinteractiontype 0

 


    

}

source [file join [file dirname [info script]] warmup.tcl]
source [file join [file dirname [info script]] topo.tcl]
source [file join [file dirname [info script]] math.tcl]
source [file join [file dirname [info script]] misc.tcl]
source [file join [file dirname [info script]] setup.tcl]



