#To call first:


proc add_drude_to_core { bondId_drude id_core id_drude type_drude polarization initialSpread } {

	#PARAMETERS:
	#--------------------------------------------------------------------------------
	#bondId_drude:	setted by 'inter $bondId_drude drude [...]' 
	#id_core:	    particle ID of existing core particle 
	#id_drude:	    free particle ID for drude particle
	#type_drude:	particle type for drude particle
    #polarization:  in [length]^-3 gets unitless with particle volume (*sigma_core^-3)

	#SCRIPT-DESCRIPTION:
	#--------------------------------------------------------------------------------
    #       -User creates drude bond: inter $bondId_drude drude $temp_core $gamma_core $temp_drude $gamma_drude $k_drude $mass_drude $r_cut
	#		-Adds drude particle to existing particle 'id_core'
	#		-Disables thermostat for core and drude via LPP and temp=gamma=0
	#		-Adds drude-bond between core and drude: Langevin T* on distance vector and T on center of mass + Subtract electrostatic core <-> drude + Harmonic bond
    #		-Subtracts drude mass and charge from core (always uses negative charge on drude particle)

	#REQUIREMENTS:
	#--------------------------------------------------------------------------------
	#		Features: LANGEVIN_PER_PARTICLE, ELECTROSTATICS, MASS
	#		Existing charged core particle

	set warnings ""

    if {$polarization <= 0} {
        return "ERROR: Polarization must be a positive number."
    } elseif { [part $id_core print] == "na" } {
        return "ERROR: Can't find core particle with id $id_core."
    } elseif { [inter $bondId_drude] == "unknown bonded interaction number $bondId_drude" } {
        return "ERROR: Can't find drude bond with bond id $bondId_drude. Create a drude-bond interaction type first."
    } elseif { [part $id_drude] != "na" } {
        set warnings "WARNING: Particle with id $id_drude already exists."
    }

	set q_tot [part $id_core print q]
	set mass_tot [part $id_core print mass]
	
    set k_drude [lindex [inter $bondId_drude] 7]
	set mass_drude [lindex [inter $bondId_drude] 8]

    set sign_q -1
#    if {$q_tot > 0} {
#        set sign_q 1
#    }

    set q_drude [expr $sign_q * pow($k_drude*$polarization, 0.5)]
    set q_core [expr $q_tot - $q_drude]
    set mass_core [expr $mass_tot-$mass_drude]
    
    set dx [expr 2.0*([t_random]-0.5)*$initialSpread]
    set dy [expr 2.0*([t_random]-0.5)*$initialSpread]
    set dz [expr 2.0*([t_random]-0.5)*$initialSpread]

	set drude_px [expr [lindex [part $id_core print pos] 0] + $dx]
	set drude_py [expr [lindex [part $id_core print pos] 1] + $dy]
	set drude_pz [expr [lindex [part $id_core print pos] 2] + $dz]

    part $id_drude pos $drude_px $drude_py $drude_pz v 0 0 0 q $q_drude type $type_drude mass $mass_drude temp 0 gamma 0 exclude $id_core
	part $id_core mass $mass_core q $q_core temp 0 gamma 0 bond $bondId_drude $id_drude exclude $id_drude

	return "M Created Drude particle:\n\
		id: $id_drude | charge: [format %.3f $q_drude] | mass: [format %.3f $mass_drude]\n\
		New core charge: [format %.3f $q_core] | New core mass: [format %.3f $mass_core]\n\
		$warnings"
}
