#!/bin/sh
TESTCASES="lj.tcl tabulated.tcl madelung.tcl kinetic.tcl lj.tcl lj-cos.tcl harm.tcl fene.tcl dh.tcl mmm1d.tcl gb.tcl rotation.tcl constraints.tcl thermostat.tcl intpbc.tcl intppbc.tcl analysis.tcl"
# 



# List of testcases to be done (and people responsible for them):
#################################################################
#
# - lj: (AxA / Status: Done.)
#   Testing forces, energies, pressures of the LJ-interaction.
# - tabulated (Ira / Status: Done. )
#   Testing forces, energies, pressures of the tabulated interaction.
# - lj-cos: (AxA / Status: Done.)
#   Testing forces, energies, pressures of the LJ-cos-interaction.
# - fene/harm: (BAM / Status: Done.)
#   Testing forces, energies, pressures of the FENE-/harmonic-interaction.
# - p3m/dh: (HL / Status: Done.)
#   Testing electrostatic interactions.
# - mmm1d: (AxA / Status: Done.)
#   Testing forces of the mmm1D-interaction.
# - gb: (DmA / Status: Done.)
#   Testing forces, energies, pressures of the Gay-Berne-Potential.
# - rotation: (DmA / Status: Done.)
#   Testing system with rotational degrees of freedom.
# - const: (HL / Status: Done.)
#   Check if constraints are working.
# - Ekin: (BAM / Status: Done.)
#   Small check of velocities, forces and kinetic energy for a two-particle-electrostatic system
# - thermo: (FRM / Status: Done.)
#   10000 timesteps only thermostat - does the given temperature (i.e. Ekin) remain constant?
# - Int-PBC: (AxA / Status: Done.)
#   Periodic Boundary Integration.
# - Int-PPBC: (AxA / Status: Done.)
#   Partial Periodic Boundary Integration.
# - Analysis: (BAM / Status: Done.)
#   Checking checkpoints and analysis routines.
# - Vol: (BAM / Status: Planned.)
#   Testing systems with resizing boxes.
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
