#!/bin/sh
TESTCASES="lj.tcl mmm1d.tcl dh.tcl madelung.tcl"
# lj-cos.tcl FENE.tcl harmonic.tcl... constraints thermosim energy pressure rotation gay-berne
# 
# List of testcases to be done (and people responsible for them):
#################################################################
#
# - lj: (AxA / Status: Done/Scheduled.)
#   Testing forces, energies, pressures of the LJ-interaction.
#   Enhancements scheduled: test of lj-cos, shift, offset, different sizes.
# - fene/harm: (BAM / Status: Finishing.)
#   Testing forces, energies, pressures of the FENE-/harmonic-interaction.
# - p3m/dh: (HL)
#   Testing electrostatic interactions.
# - mmm1D: (AxA / Status: Finishing.)
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

for np in 1 2 3 4 6 8; do
    for f in $TESTCASES; do
	if ! ./$f $np; then
	    set error 1
	    break
	fi
    done
    if test "$error" = "1"; then
	break;
    fi	
done

if test "$error" == "1"; then
    echo "Testcase $f failed for $np nodes."
fi
