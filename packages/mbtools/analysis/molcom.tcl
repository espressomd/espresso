# ::mbtools::analysis::analyze_molforce --
#
# Extract the centre of mass of a molecule from Espresso
# 

namespace eval ::mbtools::analysis {}

namespace eval ::mbtools::analysis::molcom {
    variable av_molcom
    variable av_molcom_i
    variable f_tvsmf
    variable verbose
    variable mollist

    namespace export setup_molcom
    namespace export analyze_molcom
    namespace export printav_molcom
    namespace export resetav_molcom
}

proc ::mbtools::analysis::molcom::resetav_molcom { } {
    variable av_molcom
    variable av_molcom_i
    variable mollist
    
    set av_molcom ""
    foreach mol $mollist {
       lappend av_molcom { 0.0 0.0 0.0 }
    }
    set av_molcom_i 0
}

proc ::mbtools::analysis::molcom::printav_molcom { } {
    variable av_molcom_i
    variable av_molcom
    variable f_tvsmf
    global ::mbtools::analysis::time
    variable mollist
    
    if { $av_molcom_i > 0 } {
	set i 0;
	puts $f_tvsmf "$time " nonewline
	foreach mol $mollist {
	    set this_av_molcom [lindex $av_molcom $i]
	    set avmfx [expr [lindex $this_av_molcom 0]/($av_molcom_i*1.0)]
	    set avmfy [expr [lindex $this_av_molcom 1]/($av_molcom_i*1.0)]
	    set avmfz [expr [lindex $this_av_molcom 2]/($av_molcom_i*1.0)]
	    puts $f_tvsmf "$avmfx $avmfy $avmfz " nonewline
	    incr i
	}
	puts $f_tvsmf ""
    } else {
	::mmsg::warn [namespace current] "can't print average molcom"
	flush stdout
    }
    flush $f_tvsmf
}

proc ::mbtools::analysis::molcom::setup_molcom { args } {
    global  ::mbtools::analysis::outputdir
    variable f_tvsmf
    global ::mbtools::analysis::iotype
    global ::mbtools::analysis::suffix
    variable verbose
    variable av_molcom
    variable av_molcom_i
    variable mollist
    
    ::mmsg::debug [namespace current] "opening $outputdir/time_vs_molcom$suffix "
    
    set options {
	{verbose "print out lots of stuff" }
	{mollist.arg "" "find com for which molecules"}
    }
    set usage "Usage: setup_molcom mollist:verbose"
    array set params [::cmdline::getoptions args $options $usage]

    set verbose $params(verbose)
    set mollist $params(mollist)

    if { [file exists "$outputdir/time_vs_molcom$suffix"] } {
	set newfile 0
    } else { 
	set newfile 1
    }
    set f_tvsmf [open "$outputdir/time_vs_molcom$suffix" $iotype]
    if { $newfile || $iotype == "w"} {
	puts $f_tvsmf "\#      " nonewline
	foreach mol $mollist {
	    puts $f_tvsmf "$mol           " nonewline
	}
	puts $f_tvsmf ""
	puts $f_tvsmf "\# Time " nonewline
	foreach mol $mollist {
	    puts $f_tvsmf "x   y   z   " nonewline
	}
	puts $f_tvsmf ""
    }
    foreach mol $mollist {
	lappend av_molcom { 0.0 0.0 0.0 }
    }
    set av_molcom_i 0
}

proc ::mbtools::analysis::molcom::analyze_molcom {  } {
    ::mmsg::send [namespace current] "analyzing molcom"
    variable av_molcom_i
    variable av_molcom
    variable verbose
    variable mollist

    set i 0;
    set av_molcom_new ""
    incr av_molcom_i
    foreach mol $mollist {
	set inst_molcom [analyze mol com $mol]
	set this_av_molcom [lindex $av_molcom $i]
	set this_av_molcom "[expr [lindex $this_av_molcom 0] + [lindex $inst_molcom 0] ]  \
                            [expr [lindex $this_av_molcom 1] + [lindex $inst_molcom 1] ]  \
                            [expr [lindex $this_av_molcom 2] + [lindex $inst_molcom 2] ]"
	lappend av_molcom_new $this_av_molcom
	if { $verbose } {
	    set avmfx [expr [lindex $this_av_molcom 0]/($av_molcom_i*1.0)]
	    set avmfy [expr [lindex $this_av_molcom 1]/($av_molcom_i*1.0)]
	    set avmfz [expr [lindex $this_av_molcom 2]/($av_molcom_i*1.0)]
	    ::mmsg::send [namespace current]  "For mol $mol - L: [lindex $inst_molcom 0] [lindex $inst_molcom 1] [lindex $inst_molcom 2] :: <L> $avmfx $avmfy $avmfz"
	    flush stdout
	}

	incr i
    }   
    set av_molcom $av_molcom_new

    ::mmsg::debug [namespace current] "done"
}
