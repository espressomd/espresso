# Copyright (C) 2012 The ESPResSo project
# Copyright (C) 2007,2008,2009,2010,2011 Olaf Lenz, Axel Arnold
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
SUBDIRS = configs .

# alphabetically sorted list of test scripts
tests = \
	analysis.tcl \
	collision-detection-centers.tcl \
	collision-detection-poc.tcl \
	comforce.tcl \
	comfixed.tcl \
	command_syntax.tcl \
	constraints.tcl \
	constraints_reflecting.tcl \
	correlation.tcl \
	constraints_rhomboid.tcl \
	dh.tcl \
	dielectric_test.tcl \
	el2d.tcl \
	el2d_die.tcl \
	el2d_nonneutral.tcl \
	exclusions.tcl \
	fene.tcl \
	gb.tcl \
	harm.tcl \
	iccp3m.tcl \
	intpbc.tcl \
	intppbc.tcl \
	kinetic.tcl \
	layered.tcl \
	lb.tcl \
	lb_stokes_sphere.tcl \
	lb_gpu.tcl \
	lj.tcl \
	lj-cos.tcl \
	lj-generic.tcl \
	madelung.tcl \
	maggs.tcl \
	magnetic-field.tcl \
	mass.tcl \
	mass-and-rinertia.tcl \
	mdlc.tcl \
	mmm1d.tcl \
	npt.tcl \
	nsquare.tcl \
	nve_pe.tcl \
	observable.tcl \
	p3m.tcl \
	p3m_magnetostatics.tcl \
	p3m_magnetostatics2.tcl \
	p3m_simple_noncubic.tcl \
	p3m_wall.tcl \
	rotation.tcl \
	tabulated.tcl \
	thermostat.tcl \
        tunable_slip.tcl \
	virtual-sites.tcl
# please keep the alphabetic ordering of the above list!

deactivated_tests =

# add data files for the tests here
EXTRA_DIST = $(tests) $(deactivated_tests) \
	configs \
	tests_common.tcl \
	analysis_system.data.00.gz analysis_system.data.01.gz \
	analysis_system.data.02.gz analysis_system.data.03.gz \
	analysis_system.data.04.gz analysis_system.data.05.gz \
	analysis_system.data.06.gz analysis_system.data.07.gz \
	analysis_system.data.08.gz analysis_system.data.09.gz \
	analysis_system.data.10.gz analysis_system.data.chk \
	comforce_system.data comfixed_system.data \
	constraints_system.data \
        p3m_magnetostatics.data \
        p3m_magnetostatics2_system.data p3m_magnetostatics2_expected.data \
	p3m_wall_system.data \
	p3m_system.data el2d_system.data el2d_system_die.data \
	mdlc_system.data mdlc_expected_energy.data mdlc_expected_force_torque.data \
	mmm1d_system.data dh_system.data \
	lj_system.data lj-cos_system.data \
	gb_system.data \
	tabulated_system.data lj1.tab lj2.tab lj3.tab \
	fene_system.data.gz harm_system.data.gz \
	intpbc_system.data.gz intppbc.tcl intppbc_system.data.gz \
	npt_lj_system.data \
	thermostat.data thermostat_rot.data \
	mass_system.data \
	lb_system.data \
	\
	pe_micelle.tcl pe_micelle.data \
	uwerr_test.tcl uwerr_test.data \
	\
	gen_fene.tcl gen_harm.tcl \
	\
	maggs_correct_rdf.data

# run the testsuite
check-local: runtest.sh
	if test -n "$(processors)"; then \
	  @builddir@/runtest.sh -p "${processors}" $(tests); \
	else \
	  @builddir@/runtest.sh $(tests); \
	fi

CLEANFILES = runtest.log

DISTCLEANFILES = \
	fftw3_1d_wisdom_*.file \
	dfftw3_1d_wisdom_*.file
