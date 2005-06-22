# Messages
#
#  Utility for printing messages of various types
#

package provide mmsg 0.1.0
package require cmdline
namespace eval mmsg {
    namespace export send
    namespace export err
    namespace export warn
    namespace export debug

    # The following variables determine what mmsg will print.  To
    # disable printing from a particular module delete it from the
    # list of allowable namespaces or vice versa if you want to
    # enable.  To enable printing of error warnings etc then set the
    # value of the appropriate variable to 1 otherwise set it to zero.
    variable allowablenamespaces { :: ::system_generation ::setup_utilities ::std_analysis }
    variable enablesend 1
    variable enableerr 1
    variable enablewarn 1
    variable enabledebug 1

    # For keeping track of whether we are printing again on the
    # sameline or not.
    variable sameline 0

}

# ::mmsg::setnamespaces -- 
# 
# Proceedure for specifying which namespaces to allow output from.
# This might not work on all systems (eg macintosh) if multiple copies
# of mmsg get loaded for some reason.
#
proc ::mmsg::setnamespaces { namespacelist } {

    variable allowablenamespaces 
    set allowablenamespaces $namespacelist

}

# ::mmsg::send --
#
# Wrapper for the ::mmsg::print command which prints messages with no prefix
#
proc ::mmsg::send { namespace string { newline "yes" } } {
    variable enablesend
    if { $enablesend } {
	::mmsg::print $namespace $string "" $newline
    }
}

# ::mmsg::err --
#
# Wrapper for the ::mmsg::print command which prints messages with the
# prefix Error and which calls exit
#
proc ::mmsg::err { namespace string { newline "yes" } } {
    variable enableerr
    if { $enableerr } {
	::mmsg::print $namespace $string "Error: " $newline -force
    }
    exit
}

# ::mmsg::warn --
#
# Wrapper for the ::mmsg::print command which prints messages with the
# prefix Warning 
#
proc ::mmsg::warn { namespace string { newline "yes" } } {
    variable enablewarn
    if { $enablewarn } {
	::mmsg::print $namespace $string "Warning: " $newline
    }
}

# ::mmsg::debug --
#
# Wrapper for the ::mmsg::print command which prints debug messages
#
proc ::mmsg::debug { namespace string { newline "yes" } } {
    variable enabledebug
    if { $enabledebug } {
	::mmsg::print $namespace $string "Debug: " $newline
    }
}

# ::mmsg::checknamespace --
#
# Check the namespace provided against the list of allowable
# namespaces return 1 if there is a match or 0 if not
#
proc ::mmsg::checknamespace { namespace allowable } {
    foreach name $allowable {
	if { $name == $namespace } {
	    return 1
	}
    }
    return 0
}

# ::mmsg::print --
#
# Print a message provided it is from an allowable namespace
#
#
proc ::mmsg::print {namespace string prefix newline args } {
    variable allowablenamespaces
    variable sameline
    set options {
	{force "override namespace restrictions"}
    }
    set usage "Usage: print \[force] "
    array set params [::cmdline::getoptions args $options $usage]

    if { $params(force) } {
	set namespace_OK 1
    } else {
	set namespace_OK [checknamespace $namespace $allowablenamespaces]
    }

    if { $namespace_OK } {
	if { $newline != "yes" } {
	    if { $sameline } {
		puts -nonewline  "$string"
	    } else {
		puts -nonewline  "$namespace > $prefix $string"
		set sameline 1
	    }
	} else {
	    if { $sameline } {
		puts  "$string"
		set sameline 0
	    } else {
		puts  "$namespace > $prefix $string"
	    }
	}
    }
}


