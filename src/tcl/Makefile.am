# Copyright (C) 2013,2014 Olaf Lenz
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

EXTRA_DIST =
CLEANFILES =
AM_CPPFLAGS = -I$(top_srcdir)/src/core -I$(top_srcdir)/src -I$(builddir)/../core
if MPI_FAKE
AM_CPPFLAGS += -I$(top_srcdir)/src/core/mpifake
endif
AM_DEFAULT_SOURCE_EXT = .cpp

noinst_LTLIBRARIES = libEspressoTcl.la
libEspressoTcl_la_SOURCES = \
	TclOutputHelper.hpp \
	bin_tcl.cpp \
	binary_file_tcl.cpp binary_file_tcl.hpp \
	blockfile_tcl.cpp \
	cells_tcl.cpp \
	cells_tcl.hpp \
	channels_tcl.cpp \
	collision_tcl.cpp \
	comfixed_tcl.cpp comfixed_tcl.hpp \
	comforce_tcl.cpp comforce_tcl.hpp \
	config_tcl.cpp \
	constraint_tcl.cpp constraint_tcl.hpp \
	domain_decomposition_tcl.cpp domain_decomposition_tcl.hpp \
	electrokinetics_tcl.cpp electrokinetics_tcl.hpp \
	external_potential_tcl.cpp external_potential_tcl.hpp \
	energy_tcl.cpp \
	galilei_tcl.cpp galilei_tcl.hpp \
	global_tcl.cpp global_tcl.hpp \
	grid_tcl.cpp grid_tcl.hpp \
	iccp3m_tcl.cpp iccp3m_tcl.hpp \
	imd_tcl.cpp \
	initialize_interpreter.cpp \
	integrate_tcl.cpp integrate_tcl.hpp \
	interaction_data_tcl.cpp interaction_data_tcl.hpp \
	lb-boundaries_tcl.cpp lb-boundaries_tcl.hpp \
	lb_tcl.cpp lb_tcl.hpp \
	metadynamics_tcl.cpp metadynamics_tcl.hpp \
	nemd_tcl.cpp \
	mol_cut_tcl.cpp mol_cut_tcl.hpp \
	parser.cpp parser.hpp \
	particle_data_tcl.cpp \
	polymer_tcl.cpp polymer_tcl.hpp \
	pressure_tcl.cpp pressure_tcl.hpp \
	random_tcl.cpp random_tcl.hpp \
	reaction_tcl.cpp reaction_tcl.hpp \
	rattle_tcl.cpp rattle_tcl.hpp \
	statistics_chain_tcl.cpp statistics_chain_tcl.hpp \
	statistics_cluster_tcl.cpp statistics_cluster_tcl.hpp \
	statistics_correlation_tcl.cpp statistics_correlation_tcl.hpp \
	statistics_fluid_tcl.cpp statistics_fluid_tcl.hpp \
	statistics_observable_tcl.cpp statistics_observable_tcl.hpp \
	statistics_wallstuff_tcl.cpp statistics_wallstuff_tcl.hpp \
	statistics_tcl.cpp statistics_tcl.hpp \
	thermostat_tcl.cpp thermostat_tcl.hpp \
	topology_tcl.cpp \
	tuning_tcl.cpp \
	uwerr_tcl.cpp \
	virtual_sites_com_tcl.cpp virtual_sites_com_tcl.hpp \
	ghmc_tcl.cpp ghmc_tcl.hpp

# nonbonded potentials and forces
libEspressoTcl_la_SOURCES += \
	bmhtf-nacl_tcl.cpp bmhtf-nacl_tcl.hpp \
	buckingham_tcl.cpp buckingham_tcl.hpp \
	dpd_tcl.cpp dpd_tcl.hpp \
	forcecap_tcl.cpp forcecap_tcl.hpp \
	gaussian_tcl.cpp gaussian_tcl.hpp \
	gb_tcl.cpp gb_tcl.hpp \
	hat_tcl.cpp hat_tcl.hpp \
	hertzian_tcl.cpp hertzian_tcl.hpp \
	lj_tcl.cpp lj_tcl.hpp \
	ljangle_tcl.cpp ljangle_tcl.hpp \
	ljcos2_tcl.cpp ljcos2_tcl.hpp \
	ljcos_tcl.cpp ljcos_tcl.hpp \
	ljgen_tcl.cpp ljgen_tcl.hpp \
	morse_tcl.cpp morse_tcl.hpp \
	soft_sphere_tcl.cpp soft_sphere_tcl.hpp \
	steppot_tcl.cpp steppot_tcl.hpp \
	tab_tcl.cpp tab_tcl.hpp \
	tunable_slip_tcl.cpp tunable_slip_tcl.hpp

# bonded potentials and forces
libEspressoTcl_la_SOURCES += \
	angle_tcl.cpp angle_tcl.hpp \
	angle_harmonic_tcl.cpp angle_harmonic_tcl.hpp \
	angle_cosine_tcl.cpp angle_cosine_tcl.hpp \
	angle_cossquare_tcl.cpp angle_cossquare_tcl.hpp \
	angledist_tcl.cpp angledist_tcl.hpp \
	dihedral_tcl.cpp dihedral_tcl.hpp \
	endangledist_tcl.cpp endangledist_tcl.hpp \
	fene_tcl.cpp fene_tcl.hpp \
	harmonic_tcl.cpp harmonic_tcl.hpp \
	overlap_tcl.cpp overlap_tcl.hpp \
	subt_lj_tcl.cpp subt_lj_tcl.hpp \
	object-in-fluid/area_force_local_tcl.cpp \
	object-in-fluid/area_force_local_tcl.hpp \
	object-in-fluid/area_force_global_tcl.cpp \
	object-in-fluid/area_force_global_tcl.hpp \
	object-in-fluid/bending_force_tcl.cpp \
	object-in-fluid/bending_force_tcl.hpp \
	object-in-fluid/stretching_force_tcl.cpp \
	object-in-fluid/stretching_force_tcl.hpp \
	object-in-fluid/stretchlin_force_tcl.cpp \
	object-in-fluid/stretchlin_force_tcl.hpp \
	object-in-fluid/volume_force_tcl.cpp \
	object-in-fluid/volume_force_tcl.hpp

# Coulomb methods
libEspressoTcl_la_SOURCES += \
	debye_hueckel_tcl.cpp debye_hueckel_tcl.hpp \
	elc_tcl.cpp elc_tcl.hpp \
	magnetic_non_p3m_methods_tcl.cpp magnetic_non_p3m_methods_tcl.hpp \
	maggs_tcl.cpp maggs_tcl.hpp \
	mmm1d_tcl.cpp mmm1d_tcl.hpp \
	mmm2d_tcl.cpp mmm2d_tcl.hpp \
	p3m-dipolar_tcl.cpp p3m-dipolar_tcl.hpp \
	p3m_tcl.cpp p3m_tcl.hpp \
	reaction_field_tcl.cpp reaction_field_tcl.hpp \
	mdlc_correction_tcl.cpp mdlc_correction_tcl.hpp	

# Generic actors
libEspressoTcl_la_SOURCES += \
	actor/HarmonicWell_tcl.cpp actor/HarmonicWell_tcl.hpp

if CUDA
libEspressoTcl_la_SOURCES += cuda_init_tcl.cpp
endif

#################################################################
# Compile the (TCL) main program
#################################################################

# Two binaries are generated: Espresso for the build dir,
# Espresso.install for the installation dir. 
# The ".install" suffix is removed upon installation.
noinst_PROGRAMS = Espresso
Espresso_CPPFLAGS = -D ESPRESSO_SCRIPTS_DEFAULT=\"$(buildscriptsdir)\" $(AM_CPPFLAGS)
Espresso_SOURCES = scriptsdir.cpp main.cpp
Espresso_LDADD = libEspressoTcl.la ../core/libEspresso.la

bin_PROGRAMS = Espresso.install
Espresso_install_CPPFLAGS = -D ESPRESSO_SCRIPTS_DEFAULT=\"$(scriptsdir)\" $(AM_CPPFLAGS)
Espresso_install_SOURCES = scriptsdir.cpp main.cpp
Espresso_install_LDADD = libEspressoTcl.la ../core/libEspresso.la 

ESPRESSO = `echo Espresso | sed '$(transform)'`$(EXEEXT)
ESPRESSO_INSTALL = `echo Espresso.install | sed '$(transform)'`$(EXEEXT)
# rename Espresso after installation
install-exec-hook:
	$(am__mv) \
		$(DESTDIR)$(bindir)/$(ESPRESSO_INSTALL) \
		$(DESTDIR)$(bindir)/$(ESPRESSO)

uninstall-local:
	-rm -f $(DESTDIR)$(bindir)/$(ESPRESSO)

#################################################################
# Build an empty CUDA object
#################################################################
# this is necessary since Xcode 4.3 otherwise builds binaries
# with a broken dyld relocation table
if CUDA

SUFFIXES=.cu

cuda_verbose = $(cuda_verbose_@AM_V@)
cuda_verbose_ = $(cuda_verbose_@AM_DEFAULT_V@)
cuda_verbose_0 = @echo "  NVCC     $@";

# nvcc does not allow for option MF for dependency generation, 
# therefore make a separate run to generate dependencies
# putting them into DEPDIR ensures configure creates empties dummies
cuda_workaround.o: cuda_workaround.cu
	$(cuda_verbose)$(NVCC) \
		$(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \
		$(AM_NVCCFLAGS) $(NVCCFLAGS) -c -o $@ $<

Espresso_LDADD += cuda_workaround.o
Espresso_install_LDADD += cuda_workaround.o
EXTRA_DIST += cuda_workaround.cu

endif
