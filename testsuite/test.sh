#!/bin/sh
TESTCASES="lj.tcl lj-cos.tcl harm.tcl fene.tcl mmm1d.tcl dh.tcl madelung.tcl"
# lj-cos.tcl FENE.tcl harmonic.tcl... constraints thermosim energy pressure rotation gay-berne
# 
# List of testcases to be done (and people responsible for them):
#################################################################
#
# - lj: (AxA / Status: Done.)
#   Testing forces, energies, pressures of the LJ-interaction.
# - lj-cos: (AxA / Status: Done.)
#   Testing forces, energies, pressures of the LJ-cos-interaction.
# - fene/harm: (BAM / Status: Finishing.)
#   Testing forces, energies, pressures of the FENE-/harmonic-interaction.
# - p3m/dh: (HL)
#   Testing electrostatic interactions.
# - mmm1d: (AxA / Status: Done.)
#   Testing forces of the mmm1D-interaction.
# - GB: (DmA / Status: Scheduled.)
#   Testing forces, energies, pressures of the Gay-Berne-Potential.
# - Rot: (DmA / Status: Scheduled.)
#   Testing system with rotational degrees of freedom.
# - const: (HL / Status: Scheduled.)
#   Check if constraints are working.
# - Ekin: (BAM / Status: Scheduled.)
#   Small check of velocities, forces and kinetic energy for a two-particle-electrostatic system
# - thermo: (FRM / Status: Scheduled.)
#   10000 timesteps only thermostat - does the given temperature (i.e. Ekin) remain constant?
# - Int-PBC: (AxA / Status: Scheduled.)
#   Periodic Boundary Integration.
# - Int-PPBC: (AxA / Status: Scheduled.)
#   Partial Periodic Boundary Integration.
# - Anal: (BAM / Status: Scheduled.)
#   Checking checkpoints and anlysis routines.
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
	    if test "$ft" == "$f"; then
		ignore=1;
		break;
	    fi
	done
	# only not blacklisted tests
	if test $ignore -eq 0; then
	    if ./$f $np $errf; then
		echo "Test $f done"
	    else
		echo "execution of script failed" > $errf
	    fi
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
    echo -e "\n\n*******************************************\n\n"
    echo -e "           Tests not done:\n"
    echo -e $missing
fi

echo -e "\n\n********************************************\n\n"
echo -e "   Gratulations! Espresso seems to be ok."
