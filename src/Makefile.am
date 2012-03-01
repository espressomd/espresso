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

# List headers that are not used by the program here
noinst_HEADERS = \
	myconfig-default.h

noinst_LIBRARIES = libEspresso.a

# List the sources of the Espresso binary here
libEspresso_a_SOURCES = \
	adresso.c adresso.h \
	cells.c cells.h \
	collision.c collision.h \
	communication.c communication.h \
	comfixed.c comfixed.h \
	comforce.h comforce.c \
	config.h \
	constraint.c constraint.h \
	debug.c debug.h \
	domain_decomposition.c domain_decomposition.h \
	energy.c energy.h \
	errorhandling.c	errorhandling.h \
	debye_hueckel.h \
	fft.c fft.h \
	fft-common.c fft-common.h \
	fft-dipolar.c fft-dipolar.h \
	forces.c forces.h \
	ghosts.c ghosts.h \
	global.c global.h \
	grid.c grid.h \
	halo.c halo.h \
	iccp3m.c iccp3m.h \
	imd.c imd.h \
	initialize.c initialize.h \
	integrate.c integrate.h \
	interaction_data.c interaction_data.h \
	lattice.c lattice.h \
	layered.c layered.h \
	lb.c lb.h \
	lb-boundaries.c lb-boundaries.h \
	lb-d3q18.h lb-d3q19.h \
	lbgpu_cfile.c \
	metadynamics.c metadynamics.h \
	modes.c	modes.h \
	molforces.c molforces.h \
	mol_cut.c mol_cut.h \
	nemd.c nemd.h \
	npt.h \
	nsquare.c nsquare.h \
	particle_data.c	particle_data.h \
	polymer.c polymer.h \
	polynom.h \
	pressure.c pressure.h \
	random.c random.h \
	rattle.c rattle.h \
	rotation.c rotation.h \
	specfunc.c specfunc.h \
	statistics.c statistics.h \
	statistics_chain.c statistics_chain.h \
	statistics_cluster.c statistics_cluster.h \
	statistics_correlation.c statistics_correlation.h \
	statistics_fluid.c statistics_fluid.h \
	statistics_molecule.c statistics_molecule.h \
	statistics_observable.c statistics_observable.h \
	thermostat.c thermostat.h \
	topology.c topology.h \
	tuning.c tuning.h \
	utils.h \
	uwerr.c uwerr.h \
	verlet.c verlet.h \
	virtual_sites.c virtual_sites.h \
	virtual_sites_com.c virtual_sites_com.h \
	virtual_sites_relative.c virtual_sites_relative.h \
	vmdsock.c vmdsock.h

# nonbonded potentials and forces
libEspresso_a_SOURCES += \
	bmhtf-nacl.c bmhtf-nacl.h \
	buckingham.c buckingham.h \
	dpd.c dpd.h \
	gb.c gb.h \
	hertzian.c hertzian.h \
	lj.c lj.h \
	ljangle.c ljangle.h \
	ljcos.c ljcos.h \
	ljcos2.c ljcos2.h \
	ljgen.c ljgen.h \
	morse.c morse.h \
	soft_sphere.c soft_sphere.h \
	steppot.c steppot.h \
	tab.c tab.h \
	tunable_slip.c tunable_slip.h

# bonded potentials and forces
libEspresso_a_SOURCES += \
	angle.c angle.h \
	angledist.c angledist.h \
	dihedral.c dihedral.h \
	endangledist.c endangledist.h \
	fene.c fene.h \
	harmonic.c harmonic.h \
	overlap.c overlap.h \
	subt_lj.c subt_lj.h

# Coulomb methods
libEspresso_a_SOURCES += \
	debye_hueckel.c debye_hueckel.h \
	elc.c elc.h \
	ewald.c ewald.h \
	magnetic_non_p3m_methods.c magnetic_non_p3m_methods.h \
	mdlc_correction.c  mdlc_correction.h \
	maggs.c maggs.h \
	mmm1d.c mmm1d.h \
	mmm2d.c mmm2d.h \
	mmm-common.c mmm-common.h \
	p3m.c p3m.h \
	p3m-common.c p3m-common.h \
	p3m-dipolar.c p3m-dipolar.h \
	reaction_field.c reaction_field.h

if TCL

libEspresso_a_SOURCES += \
	tcl/adresso_tcl.c tcl/adresso_tcl.h \
	tcl/bin_tcl.c \
	tcl/binary_file_tcl.c tcl/binary_file_tcl.h \
	tcl/blockfile_tcl.c \
	tcl/cells_tcl.c \
	tcl/channels_tcl.c \
	tcl/collision_tcl.c \
	tcl/comfixed_tcl.c tcl/comfixed_tcl.h \
	tcl/comforce_tcl.c tcl/comforce_tcl.h \
	tcl/config_tcl.c \
	tcl/constraint_tcl.c tcl/constraint_tcl.h \
	tcl/domain_decomposition_tcl.c tcl/domain_decomposition_tcl.h \
	tcl/energy_tcl.c \
	tcl/global_tcl.c tcl/global_tcl.h \
	tcl/grid_tcl.c tcl/grid_tcl.h \
	tcl/iccp3m_tcl.c tcl/iccp3m_tcl.h \
	tcl/imd_tcl.c \
	tcl/initialize_interpreter.c \
	tcl/integrate_tcl.c tcl/integrate_tcl.h \
	tcl/interaction_data_tcl.c tcl/interaction_data_tcl.h \
	tcl/lb-boundaries_tcl.c \
	tcl/lb_tcl.c tcl/lb_tcl.h \
	tcl/metadynamics_tcl.c tcl/metadynamics_tcl.h \
	tcl/nemd_tcl.c \
	tcl/mol_cut_tcl.c tcl/mol_cut_tcl.h \
	tcl/parser.c tcl/parser.h \
	tcl/particle_data_tcl.c \
	tcl/polymer_tcl.c tcl/polymer_tcl.h \
	tcl/pressure_tcl.c tcl/pressure_tcl.h \
	tcl/random_tcl.c tcl/random_tcl.h \
	tcl/rattle_tcl.c tcl/rattle_tcl.h \
	tcl/statistics_chain_tcl.c tcl/statistics_chain_tcl.h \
	tcl/statistics_cluster_tcl.c tcl/statistics_cluster_tcl.h \
	tcl/statistics_correlation_tcl.c tcl/statistics_correlation_tcl.h \
	tcl/statistics_fluid_tcl.c tcl/statistics_fluid_tcl.h \
	tcl/statistics_observable_tcl.c tcl/statistics_observable_tcl.h \
	tcl/statistics_tcl.c tcl/statistics_tcl.h \
	tcl/thermostat_tcl.c tcl/thermostat_tcl.h \
	tcl/topology_tcl.c \
	tcl/tuning_tcl.c \
	tcl/uwerr_tcl.c \
	tcl/virtual_sites_com_tcl.c tcl/virtual_sites_com_tcl.h

# nonbonded potentials and forces
libEspresso_a_SOURCES += \
	tcl/bmhtf-nacl_tcl.c tcl/bmhtf-nacl_tcl.h \
	tcl/buckingham_tcl.c tcl/buckingham_tcl.h \
	tcl/dpd_tcl.c tcl/dpd_tcl.h \
	tcl/gb_tcl.c tcl/gb_tcl.h \
	tcl/hertzian_tcl.c tcl/hertzian_tcl.h \
	tcl/lj_tcl.c tcl/lj_tcl.h \
	tcl/ljangle_tcl.c tcl/ljangle_tcl.h \
	tcl/ljcos2_tcl.c tcl/ljcos2_tcl.h \
	tcl/ljcos_tcl.c tcl/ljcos_tcl.h \
	tcl/ljgen_tcl.c tcl/ljgen_tcl.h \
	tcl/morse_tcl.c tcl/morse_tcl.h \
	tcl/soft_sphere_tcl.c tcl/soft_sphere_tcl.h \
	tcl/steppot_tcl.c tcl/steppot_tcl.h \
	tcl/tab_tcl.c tcl/tab_tcl.h \
	tcl/tunable_slip_tcl.c tcl/tunable_slip_tcl.h

# bonded potentials and forces
libEspresso_a_SOURCES += \
	tcl/angle_tcl.c tcl/angle_tcl.h \
	tcl/angledist_tcl.c tcl/angledist_tcl.h \
	tcl/dihedral_tcl.c tcl/dihedral_tcl.h \
	tcl/endangledist_tcl.c tcl/endangledist_tcl.h \
	tcl/fene_tcl.c tcl/fene_tcl.h \
	tcl/harmonic_tcl.c tcl/harmonic_tcl.h \
	tcl/overlap_tcl.c tcl/overlap_tcl.h \
	tcl/subt_lj_tcl.c tcl/subt_lj_tcl.h

# Coulomb methods
libEspresso_a_SOURCES += \
	tcl/debye_hueckel_tcl.c tcl/debye_hueckel_tcl.h \
	tcl/elc_tcl.c tcl/elc_tcl.h \
	tcl/ewald_tcl.c tcl/ewald_tcl.h \
	tcl/magnetic_non_p3m_methods_tcl.c tcl/magnetic_non_p3m_methods_tcl.h \
	tcl/maggs_tcl.c tcl/maggs_tcl.h \
	tcl/mmm1d_tcl.c tcl/mmm1d_tcl.h \
	tcl/mmm2d_tcl.c tcl/mmm2d_tcl.h \
	tcl/p3m-dipolar_tcl.c tcl/p3m-dipolar_tcl.h \
	tcl/p3m_tcl.c tcl/p3m_tcl.h \
	tcl/reaction_field_tcl.c tcl/reaction_field_tcl.h \
	tcl/mdlc_correction_tcl.c tcl/mdlc_correction_tcl.h	

endif

if TCL

#################################################################
# Compile the (TCL) main program
#################################################################

# Two binaries are generated: Espresso for the build dir,
# Espresso.install for the installation dir. 
# The ".install" suffix is removed upon installation.
noinst_PROGRAMS = Espresso
Espresso_CPPFLAGS = -D ESPRESSO_SCRIPTS_DEFAULT=\"$(buildscriptsdir)\"
Espresso_SOURCES = tcl/scriptsdir.c tcl/main.c
Espresso_LDADD = libEspresso.a

bin_PROGRAMS = Espresso.install
Espresso_install_CPPFLAGS = -D ESPRESSO_SCRIPTS_DEFAULT=\"$(scriptsdir)\"
Espresso_install_SOURCES = tcl/scriptsdir.c tcl/main.c
Espresso_install_LDADD = libEspresso.a

# rename Espresso after installation
install-exec-hook:
	$(am__mv) \
		$(DESTDIR)$(bindir)/Espresso.install \
		$(DESTDIR)$(bindir)/Espresso

uninstall-local:
	-rm -f $(DESTDIR)$(bindir)/Espresso

endif

#################################################################
# Handling of myconfig.h
#################################################################
libEspresso_a_SOURCES += myconfig-final.h
BUILT_SOURCES = myconfig-final.h
CLEANFILES = myconfig-final.h

myconfig-final.h: FORCE
	@echo -n "Determining myconfig file..."
	@config_files="\
	  $(top_builddir)/$(myconfig) \
	  $(top_srcdir)/$(myconfig)"; \
	if test $(myconfig) = myconfig.h; then \
	  config_files="$$config_files $(top_srcdir)/src/myconfig-default.h"; \
	fi; \
	for file in $$config_files; do \
	  if test -e $$file; then \
	    echo " $$file"; \
	    myconfig_found=1; \
	    if test -e $@; then \
	      cmp -s $$file $@ || cp $$file $@; \
	    else \
	       cp $$file $@; \
	    fi; \
	    break; \
	  fi; \
	done; \
	if test ! $$myconfig_found; then \
	  echo "ERROR: Couldn't find $(myconfig) anywhere"; \
	  exit 1; \
	fi

#################################################################
# Handling of the features
#################################################################
libEspresso_a_SOURCES += config-features.c
EXTRA_DIST = \
	features.def \
	config-features.h config-features.c \
	config-doxygen.h

if HAVE_PYTHON

BUILT_SOURCES +=  config-features.h config-features.c
CLEANFILES += config-features.h config-features.c config-doxygen.h
config-features.h config-features.c: features.def $(top_srcdir)/config/gen_featureconfig.py
	$(PYTHON) $(top_srcdir)/config/gen_featureconfig.py $(srcdir)/features.def config-features.h config-features.c

config-doxygen.h: features.def $(top_srcdir)/config/gen_doxyconfig.py
	$(PYTHON) $(top_srcdir)/config/gen_doxyconfig.py $(srcdir)/features.def config-doxygen.h

else

config-doxygen.h config-features.h config-features.c: features.def
	@echo "Python was not found in your PATH."
	@echo "If you change the set of features, you need Python!"
	@echo "Install it and rerun configure."

endif

doc doxygen: config-doxygen.h

#################################################################
# Handling of the version
#################################################################
libEspresso_a_SOURCES += config-version.c 
BUILT_SOURCES += config-version.c
CLEANFILES += config-version.c
config-version.c: FORCE
	VERSION=`cd $(top_srcdir); sh config/genversion.sh -c`; \
	if ! test -e "config-version.c" || \
		test "$$VERSION" != "`cat config-version.c`"; then \
		echo "$$VERSION" > config-version.c; \
	fi

#################################################################
# Handling of the MPI fake implementation
#################################################################
if MPI_FAKE
libEspresso_a_SOURCES += mpifake/mpi.h mpifake/mpi.c
# mpifake should come before any system includes
AM_CPPFLAGS = -I$(srcdir)/mpifake
if TCL
Espresso_CPPFLAGS += $(AM_CPPFLAGS)
Espresso_install_CPPFLAGS += $(AM_CPPFLAGS)
endif
endif

##################################################
# CUDA rules
##################################################
if CUDA
SUFFIXES=.cu
.cu.o:
	$(NVCC) -c $(NVCCFLAGS) $(DEFAULT_INCLUDES) $(CPPFLAGS) $(INCLUDES) $(DEFS) -o $@ $<

CUDA_SOURCES = \
	cuda_init.cu cuda_init.h \
	lbgpu.cu lbgpu.h

if TCL

CUDA_SOURCES += \
	tcl/cuda_init_tcl.c

endif

libEspresso_a_SOURCES += $(CUDA_SOURCES)

EXTRA_DIST += $(CUDA_SOURCES)
endif

.PHONY: FORCE
FORCE:
