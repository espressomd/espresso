#!/bin/sh

TESTCASES="nve_pe.tcl npt.tcl madelung.tcl \
    lj.tcl lj-cos.tcl harm.tcl fene.tcl dh.tcl tabulated.tcl mmm1d.tcl gb.tcl rotation.tcl constraints.tcl \
    kinetic.tcl thermostat.tcl \
    intpbc.tcl intppbc.tcl layered.tcl nsquare.tcl \
    comforce.tcl comfixed.tcl analysis.tcl"

# List of testcases to be done (and people responsible for them):
#################################################################
#
# - madelung: (HL / Status: Done.)
#   Calculating the Madelung constant in a NaCl-crystal.
# - lj/lj-cos: (AxA / Status: Done.)
#   Testing forces, energies, pressures of the LJ/LJ-cos-interaction.
# - harm/fene: (BAM / Status: Done.)
#   Testing forces, energies, pressures of the harmonic-/FENE-interaction.
# - p3m/dh: (HL / Status: Done.)
#   Testing electrostatic interactions.
# - tabulated (Ira / Status: Done. )
#   Testing forces, energies, pressures of the tabulated interaction.
# - mmm1d: (AxA / Status: Done.)
#   Testing forces of the mmm1D-interaction.
# - gb: (DmA / Status: Done.)
#   Testing forces, energies, pressures of the Gay-Berne-Potential.
# - rotation: (DmA / Status: Done.)
#   Testing system with rotational degrees of freedom.
# - constraints: (HL / Status: Done.)
#   Check if constraints are working.
# - kinetic: (BAM / Status: Done.)
#   Small check of velocities, forces and kinetic energy for a two-particle-electrostatic system
# - thermostat: (FRM / Status: Done.)
#   10000 timesteps only thermostat - does the given temperature (i.e. Ekin) remain constant?
# - Int-PBC/Int-PPBC: (AxA / Status: Done.)
#   (Partial) Periodic Boundary Integration.
# - layered/nsquare: (AxA / Status: Done.)
#   Testcases for the layered- and N^2-particle structures available since v1.5.
# - Analysis: (BAM / Status: Done.)
#   Checking checkpoints and analysis routines.
# - comfixed: (MS / Status: Done.)
#   Testing forces, energies, pressures of the comfixed.
# - comforce: (MS / Status: Done.)
#   Testing forces, energies, pressures of the comfixed.
# - nve_pe.tcl: (MS / Status: Done.)
#   Testing energy conservation with PE chain 
#
#################################################################

errf=_test.sh_error.$$

# list of tests that are not supported by this version
blacklist=
# and what is missing for what test
missing=
for np in 1 2 3 4 6 8; do
    for f in $TESTCASES; do
	ignore=0
	for ft in $blacklist; do
	    if test "$ft" = "$f"; then
		ignore=1;
		break;
	    fi
	done
	# only not blacklisted tests
	if test $ignore -eq 0; then
	    # this is removed if the script runs through
	    echo "execution of script failed at unexpected point" > $errf
	    ./$f $np $errf
	    if test -f $errf; then
		if grep -q -e "^not compiled in:" $errf; then
		    missing="$missing$f: `cat $errf`\n"
		    blacklist="$blacklist $f"
		    rm -f $errf
		else
		    echo "Testcase $f failed for $np nodes with error `cat $errf`."
		    rm -f $errf
		    exit -666
		fi
	    fi
	else
	    echo "Test $f is blacklisted, ignoring..."
	fi
    done
done

if test "$missing" != ""; then
    echo -e "\n\n===============================================\n\n"
    echo -e "                Tests not done:\n                "
    echo -e $missing
fi

echo -e "\n\n===============================================\n\n"
echo -e "   Congratulations! ESPResSo seems to be ok."
echo -e "\n\n===============================================\n\n"
