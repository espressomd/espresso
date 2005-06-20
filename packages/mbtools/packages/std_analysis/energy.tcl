# ::std_analysis::analyze_energy --
#
# Calculate the total energy of the system and break it into components
#   

namespace eval ::std_analysis {}

proc ::std_analysis::analyze_energy { printflag } {
    variable this
    mmsg::send $this "analyzing energy"
    variable av_components_en
    variable av_total_en
    variable av_kin_en
    variable av_fene_en
    variable av_harm_en
    variable av_nb_en
    variable av_en_i
    set energy_all [analyze energy]
    #	puts $energy_all
    set nb_en 0
    set nbcount 0
    for { set i 0 } { $i < [llength $energy_all ] } { incr i } {
	set tmp [lindex $energy_all $i]
	set ntmp [llength $tmp]	
	if { [ regexp "energy" $tmp ] } {
	    set total_en [lindex $tmp [expr $ntmp -1]]
	}
	if { [ regexp "kinetic" $tmp ] } {
	    set kin_en [lindex $tmp [expr $ntmp -1]]
	}
	if { [ regexp "FENE" $tmp ] } {
	    set fene_en [lindex $tmp [expr $ntmp -1]]
	}
	if { [ regexp "HARMONIC" $tmp ] } {
	    set harm_en [lindex $tmp [expr $ntmp -1]]
	}
	if { [ regexp "nonbonded" $tmp ] } {
	    set nb_en [expr $nb_en + [lindex $tmp [expr $ntmp -1]]]
	    lset av_components_en $nbcount [expr [lindex $av_components_en $nbcount] + [lindex $tmp [expr $ntmp -1]] ]

	    incr nbcount
	}
	
    }
    incr av_en_i
    set av_total_en  [expr $av_total_en + $total_en/(1.0)]
    set av_kin_en  [expr $av_kin_en + $kin_en/(1.0)]
    set av_fene_en  [expr $av_fene_en + $fene_en/(1.0)]
    set av_harm_en  [expr $av_harm_en + $harm_en/(1.0)]
    set av_nb_en  [expr $av_nb_en + $nb_en/(1.0)]
    if { $printflag } {
	mmsg::send $this "energy: [expr $av_total_en/($av_en_i*1.0)] kinetic: [expr $av_kin_en/($av_en_i*1.0)] FENE: [expr $av_fene_en/($av_en_i*1.0)] HARMONIC: [expr $av_harm_en/($av_en_i*1.0)] nonbonded: [expr $av_nb_en/($av_en_i*1.0)]"
    }
    
    mmsg::debug $this "done"
}


