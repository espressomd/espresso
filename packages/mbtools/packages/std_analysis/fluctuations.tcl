# ::std_analysis::analyze_fluctuations --
#
# Routine for calculating the power spectrum of fluctuations for a
# flat bilayer sheet.  Uses the "modes_2d" routine in espresso to
# calculate the height function and perform the fft.  
#
#
#
# A summary of what happens is as follows
#
# (1) do_analysis calls analyze_fluctuations which in turn calls
#  <analyze modes2d> in espresso.

# (2) The resulting raw mode output (raw output from fftw) is sorted
# into increasing mode (q) order and some correction factors are
# applied. This is done via the function <arrange_modes

# (3) The result is stored in av_pow 

# (4) When print_averages is called the values stored in av_pow are
# printed by the (perhaps poorly named) power_analysis
# function. Generally for a mode analysis we would only call
# print_averages right at the end so that the entire simulation is
# used as an average.  For this reason this analysis is best done on
# stored configuration files rather than during simulation.
#
# Calculating kappa
#
# In order to calculate the bending stiffness from the results of this
# analysis one needs to carefully choose the fitting range so that
# only the q^4 regime is used.  An example fitting gnuplot script is
# as follows:
#
# boxl = Set the average tensionless boxlength of your system
# gridsize = Set the value of mgrid you used for the analysis
# 
#const = (boxl/gridsize)**2
#
#
#f(x) = a+b*x**2
#
#fit [0.02:0.2] f(x) 'powav.dat' u 1:($2*$1**2) via a,b
#
#set log
#
#plot 'powav.dat' u 1:2,\
#     f(x)/x**2,\
#     a/x**2
#
#print "kappa = ", 1.0/a/const
# 
#
# See the our membrane paper for more details behind this.
#
#

namespace eval ::std_analysis {}

proc ::std_analysis::analyze_fluctuations { printflag } {
    variable this
    variable rawmodes
    variable av_pow_i
    variable av_pow
    variable mgrid 
    
    # Obtain unprocessed results from espresso
    set errnum [ catch { set rawmodes "[analyze modes2d]" } ]
    if { !$errnum } {
	set boxlen [setmd box_l]
	# Assume a square membrane in xy plane
	set Lx [lindex $boxlen 0]
	
	#Calculate the power spectrum for this configuration
	set pow_tmp [arrange_modes $rawmodes $Lx $mgrid ]
	# Add this set of results to the accumulated data
	for { set k 0 } { $k < [llength $av_pow] } {incr k } {
	    set tmpval "[expr [lindex $av_pow $k] + [lindex $pow_tmp $k]]"
	    lset av_pow $k $tmpval
	}	
	incr av_pow_i
    } else {
	mmsg::warn $this "fluctuation analysis failed"
    }    
}

# ::std_analysis::arrange_modes --
#
# Arranges rawmodes into a power spectrum
#
proc ::std_analysis::arrange_modes { rawmodes box_l mgrid } {
     # The arrangement is done in two stages 
    set round1 [modes_analysis_half $rawmodes $box_l $mgrid]
    set final [expand_powspec $round1 $box_l $mgrid]
    return $final
}


# ::std_analysis::modes_analysis_half --
#
# Performs magic corrections and sorting.  Bloody painful
#
proc ::std_analysis::modes_analysis_half { rawmodes box_l mgrid } {
	
    set norm [expr 1.0/[expr $mgrid*1.0]]
    for { set i  0 } { $i < $mgrid } { incr i } {
	set qi $i
	if { $qi > [expr $mgrid/2 - 1] } { set qi [expr $i - $mgrid] }
	set sincx [sinc [expr $qi/($mgrid*1.0)] ]
	set normx [expr 1.0/$sincx]
	set normx [expr $normx*$normx ]
	
	for { set j  0 } { $j < [expr $mgrid/2 + 1] } { incr j } {  
	    set qj $j
	    if { $qj > [expr $mgrid/2 - 1] } { set qj [expr $j - $mgrid] }
	    set sincy [sinc [expr $qj/($mgrid*1.0)]]
	    set normy [expr 1.0/$sincy]
	    set normy [expr $normy*$normy]
	    
	    set re [expr [lindex [lindex $rawmodes 1 $i $j] 0]*$norm]
	    set im [expr [lindex [lindex $rawmodes 1 $i $j] 1]*$norm]
	    
	    # Add the Deserno correction factor
	    set pow [expr ($re*$re + $im*$im)*$normy*$normx] 
	    
	    lappend powspec $pow
	}
    }
    return $powspec
}

# ::std_analysis::expand_powspec --
#
# More of the nastiness
#
proc ::std_analysis::expand_powspec { rawmodes box_l mgrid } {	
    for { set i  0 } { $i < $mgrid } { incr i } {
	for { set j  0 } { $j < [expr $mgrid] } { incr j } {  
	    set xindex $i
	    set yindex $j
	    if { $j > [expr $mgrid/2] } { 
		if { $i !=  0 } {
		    set xindex [expr $mgrid - $i ]
		    set yindex [expr $mgrid - $j ]
		} else {
		    set yindex [expr $mgrid - $j ]
		}
	    }
	    set index [expr $xindex*($mgrid/2+1) + $yindex]
	    lappend powspec [lindex $rawmodes $index]
	}
    }
    return $powspec
}


# ::std_analysis::sinc --
#
# Calculate the sinc(x) function required in order to apply the magic
# Deserno correction
#
proc ::std_analysis::sinc { x } {	
    set retval 0
    if { $x != 0.0 } {
	set retval [expr sin($x*3.141592654)/($x*3.141592654) ]
    } else {
	set retval 1
    }
    return $retval
}
    
# Basically just creates a file of the squared modulus of the
# wavevector vs the squared power.  This is essentially just a sorting
# and printing function
proc ::std_analysis::power_analysis { } {
    variable mgrid 
    variable outputdir
    variable av_pow
    variable av_pow_i
    set file "$outputdir/powav.dat"
    set Lx [lindex [setmd box_l] 0]

    set f_powav [open $file w]
    puts $f_powav "\# q^2 <|h(q)^2|>"
    set wavevectors [modes_q2 $Lx $mgrid]
    set k 0
    foreach val $av_pow {
	puts $f_powav "[lindex $wavevectors $k] [ expr $val/($av_pow_i*1.0) ]"
	incr k
    }
    close $f_powav
    
}

## -------------------------------------------------- ####
## Stuff we don't use anymore and which might not work ###
## -------------------------------------------------- ####



proc modes_q2 { box_l mgrid } {

    for { set i  0 } { $i < $mgrid } { incr i } {
	for { set j  0 } { $j < [expr $mgrid] } { incr j } { 
	    set qi $i
	    set qj $j
	    #		puts -nonewline "raw $qi $qj "	   
	    if { $qi > [expr $mgrid/2 - 1] } { set qi [expr $i - $mgrid] }
	    if { $qj > [expr $mgrid/2 - 1] } { set qj [expr $j - $mgrid] }
	    
	    set qx [expr $qi*2*3.141592/$box_l]
	    set qy [expr $qj*2*3.141592/$box_l]
	    #		puts -nonewline "$qi $qj "	    
	    set q2 [expr $qx*$qx + $qy*$qy]
	    lappend wavecs $q2
	    #		puts "$q2 "
	}
    }
    
    return $wavecs
} 


