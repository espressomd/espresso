# Copyright (C) 2012,2013 The ESPResSo project
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
	myconfig-default.hpp

noinst_LIBRARIES = libEspresso.a

# List the sources of the Espresso binary here
# config-features.cpp must be at the beginning so that it is compiled first
nodist_libEspresso_a_SOURCES = 
libEspresso_a_SOURCES = \
	config-features.cpp \
	adresso.cpp adresso.hpp \
	cells.cpp cells.hpp \
	collision.cpp collision.hpp \
	communication.cpp communication.hpp \
	comfixed.cpp comfixed.hpp \
	comforce.hpp comforce.cpp \
	config.hpp \
	constraint.cpp constraint.hpp \
	cuda_interface.cpp cuda_interface.hpp\
	cuda_init.hpp \
        debug.cpp debug.hpp \
	domain_decomposition.cpp domain_decomposition.hpp \
	energy.cpp energy.hpp \
	errorhandling.cpp	errorhandling.hpp \
	fft.cpp fft.hpp \
	fft-common.cpp fft-common.hpp \
	fft-dipolar.cpp fft-dipolar.hpp \
	forcecap.cpp forcecap.hpp \
	forces.cpp forces.hpp \
	galilei.cpp galilei.hpp \
	ghosts.cpp ghosts.hpp \
	global.cpp global.hpp \
	grid.cpp grid.hpp \
	halo.cpp halo.hpp \
	iccp3m.cpp iccp3m.hpp \
	imd.cpp imd.hpp \
	initialize.cpp initialize.hpp \
	integrate.cpp integrate.hpp \
	interaction_data.cpp interaction_data.hpp \
	lattice.cpp lattice.hpp \
	layered.cpp layered.hpp \
	lb.cpp lb.hpp \
	lb-boundaries.cpp lb-boundaries.hpp \
	lb-d3q18.hpp lb-d3q19.hpp \
	lbgpu.cpp lbgpu.hpp\
	metadynamics.cpp metadynamics.hpp \
	modes.cpp	modes.hpp \
	molforces.cpp molforces.hpp \
	mol_cut.cpp mol_cut.hpp \
	nemd.cpp nemd.hpp \
	npt.hpp \
	nsquare.cpp nsquare.hpp \
	particle_data.cpp	particle_data.hpp \
	polymer.cpp polymer.hpp \
	polynom.hpp \
	pressure.cpp pressure.hpp \
	random.cpp random.hpp \
	rattle.cpp rattle.hpp \
	reaction.cpp reaction.hpp \
	rotation.cpp rotation.hpp \
	specfunc.cpp specfunc.hpp \
	statistics.cpp statistics.hpp \
	statistics_chain.cpp statistics_chain.hpp \
	statistics_cluster.cpp statistics_cluster.hpp \
	statistics_correlation.cpp statistics_correlation.hpp \
	statistics_fluid.cpp statistics_fluid.hpp \
	statistics_molecule.cpp statistics_molecule.hpp \
	statistics_observable.cpp statistics_observable.hpp \
	statistics_wallstuff.cpp statistics_wallstuff.hpp \
	thermostat.cpp thermostat.hpp \
	topology.cpp topology.hpp \
	tuning.cpp tuning.hpp \
	utils.hpp \
	uwerr.cpp uwerr.hpp \
	verlet.cpp verlet.hpp \
	virtual_sites.cpp virtual_sites.hpp \
	virtual_sites_com.cpp virtual_sites_com.hpp \
	virtual_sites_relative.cpp virtual_sites_relative.hpp \
	vmdsock.cpp vmdsock.hpp \
	ghmc.cpp ghmc.hpp

# nonbonded potentials and forces
libEspresso_a_SOURCES += \
	bmhtf-nacl.cpp bmhtf-nacl.hpp \
	buckingham.cpp buckingham.hpp \
	dpd.cpp dpd.hpp \
	gaussian.cpp gaussian.hpp \
	gb.cpp gb.hpp \
	hat.cpp hat.hpp \
	hertzian.cpp hertzian.hpp \
	lj.cpp lj.hpp \
	ljangle.cpp ljangle.hpp \
	ljcos.cpp ljcos.hpp \
	ljcos2.cpp ljcos2.hpp \
	ljgen.cpp ljgen.hpp \
	morse.cpp morse.hpp \
	soft_sphere.cpp soft_sphere.hpp \
	steppot.cpp steppot.hpp \
	tab.cpp tab.hpp \
	tunable_slip.cpp tunable_slip.hpp 

# bonded potentials and forces
libEspresso_a_SOURCES += \
	angle.cpp angle.hpp \
	angle_harmonic.cpp angle_harmonic.hpp \
	angle_cosine.cpp angle_cosine.hpp \
	angle_cossquare.cpp angle_cossquare.hpp \
	angledist.cpp angledist.hpp \
	dihedral.cpp dihedral.hpp \
	endangledist.cpp endangledist.hpp \
	fene.cpp fene.hpp \
	harmonic.cpp harmonic.hpp \
	overlap.cpp overlap.hpp \
	subt_lj.cpp subt_lj.hpp \
	object-in-fluid/area_force_local.cpp object-in-fluid/area_force_local.hpp \
	object-in-fluid/area_force_global.cpp object-in-fluid/area_force_global.hpp \
	object-in-fluid/bending_force.cpp object-in-fluid/bending_force.hpp \
	object-in-fluid/stretching_force.cpp object-in-fluid/stretching_force.hpp \
	object-in-fluid/stretchlin_force.cpp object-in-fluid/stretchlin_force.hpp \
	object-in-fluid/volume_force.cpp object-in-fluid/volume_force.hpp

# Coulomb methods
libEspresso_a_SOURCES += \
	debye_hueckel.cpp debye_hueckel.hpp \
	elc.cpp elc.hpp \
	magnetic_non_p3m_methods.cpp magnetic_non_p3m_methods.hpp \
	mdlc_correction.cpp  mdlc_correction.hpp \
	maggs.cpp maggs.hpp \
	mmm1d.cpp mmm1d.hpp \
	mmm2d.cpp mmm2d.hpp \
	mmm-common.cpp mmm-common.hpp \
	p3m.cpp p3m.hpp \
	p3m_gpu.hpp \
        p3m-common.cpp p3m-common.hpp \
	p3m-dipolar.cpp p3m-dipolar.hpp \
	reaction_field.cpp reaction_field.hpp

if TCL

libEspresso_a_SOURCES += \
	tcl/adresso_tcl.cpp tcl/adresso_tcl.hpp \
	tcl/bin_tcl.cpp \
	tcl/binary_file_tcl.cpp tcl/binary_file_tcl.hpp \
	tcl/blockfile_tcl.cpp \
	tcl/cells_tcl.cpp \
	tcl/channels_tcl.cpp \
	tcl/collision_tcl.cpp \
	tcl/comfixed_tcl.cpp tcl/comfixed_tcl.hpp \
	tcl/comforce_tcl.cpp tcl/comforce_tcl.hpp \
	tcl/config_tcl.cpp \
	tcl/constraint_tcl.cpp tcl/constraint_tcl.hpp \
	tcl/domain_decomposition_tcl.cpp tcl/domain_decomposition_tcl.hpp \
	tcl/energy_tcl.cpp \
	tcl/galilei_tcl.cpp tcl/galilei_tcl.hpp \
	tcl/global_tcl.cpp tcl/global_tcl.hpp \
	tcl/grid_tcl.cpp tcl/grid_tcl.hpp \
	tcl/iccp3m_tcl.cpp tcl/iccp3m_tcl.hpp \
	tcl/imd_tcl.cpp \
	tcl/initialize_interpreter.cpp \
	tcl/integrate_tcl.cpp tcl/integrate_tcl.hpp \
	tcl/interaction_data_tcl.cpp tcl/interaction_data_tcl.hpp \
	tcl/lb-boundaries_tcl.cpp \
	tcl/lb_tcl.cpp tcl/lb_tcl.hpp \
	tcl/metadynamics_tcl.cpp tcl/metadynamics_tcl.hpp \
	tcl/nemd_tcl.cpp \
	tcl/mol_cut_tcl.cpp tcl/mol_cut_tcl.hpp \
	tcl/parser.cpp tcl/parser.hpp \
	tcl/particle_data_tcl.cpp \
	tcl/polymer_tcl.cpp tcl/polymer_tcl.hpp \
	tcl/pressure_tcl.cpp tcl/pressure_tcl.hpp \
	tcl/random_tcl.cpp tcl/random_tcl.hpp \
	tcl/reaction_tcl.cpp tcl/reaction_tcl.hpp \
	tcl/rattle_tcl.cpp tcl/rattle_tcl.hpp \
	tcl/statistics_chain_tcl.cpp tcl/statistics_chain_tcl.hpp \
	tcl/statistics_cluster_tcl.cpp tcl/statistics_cluster_tcl.hpp \
	tcl/statistics_correlation_tcl.cpp tcl/statistics_correlation_tcl.hpp \
	tcl/statistics_fluid_tcl.cpp tcl/statistics_fluid_tcl.hpp \
	tcl/statistics_observable_tcl.cpp tcl/statistics_observable_tcl.hpp \
	tcl/statistics_wallstuff_tcl.cpp tcl/statistics_wallstuff_tcl.hpp \
	tcl/statistics_tcl.cpp tcl/statistics_tcl.hpp \
	tcl/thermostat_tcl.cpp tcl/thermostat_tcl.hpp \
	tcl/topology_tcl.cpp \
	tcl/tuning_tcl.cpp \
	tcl/uwerr_tcl.cpp \
	tcl/virtual_sites_com_tcl.cpp tcl/virtual_sites_com_tcl.hpp \
	tcl/ghmc_tcl.cpp tcl/ghmc_tcl.hpp

# nonbonded potentials and forces
libEspresso_a_SOURCES += \
	tcl/bmhtf-nacl_tcl.cpp tcl/bmhtf-nacl_tcl.hpp \
	tcl/buckingham_tcl.cpp tcl/buckingham_tcl.hpp \
	tcl/dpd_tcl.cpp tcl/dpd_tcl.hpp \
	tcl/forcecap_tcl.cpp tcl/forcecap_tcl.hpp \
	tcl/gaussian_tcl.cpp tcl/gaussian_tcl.hpp \
	tcl/gb_tcl.cpp tcl/gb_tcl.hpp \
	tcl/hat_tcl.cpp tcl/hat_tcl.hpp \
	tcl/hertzian_tcl.cpp tcl/hertzian_tcl.hpp \
	tcl/lj_tcl.cpp tcl/lj_tcl.hpp \
	tcl/ljangle_tcl.cpp tcl/ljangle_tcl.hpp \
	tcl/ljcos2_tcl.cpp tcl/ljcos2_tcl.hpp \
	tcl/ljcos_tcl.cpp tcl/ljcos_tcl.hpp \
	tcl/ljgen_tcl.cpp tcl/ljgen_tcl.hpp \
	tcl/morse_tcl.cpp tcl/morse_tcl.hpp \
	tcl/soft_sphere_tcl.cpp tcl/soft_sphere_tcl.hpp \
	tcl/steppot_tcl.cpp tcl/steppot_tcl.hpp \
	tcl/tab_tcl.cpp tcl/tab_tcl.hpp \
	tcl/tunable_slip_tcl.cpp tcl/tunable_slip_tcl.hpp

# bonded potentials and forces
libEspresso_a_SOURCES += \
	tcl/angle_tcl.cpp tcl/angle_tcl.hpp \
	tcl/angle_harmonic_tcl.cpp tcl/angle_harmonic_tcl.hpp \
	tcl/angle_cosine_tcl.cpp tcl/angle_cosine_tcl.hpp \
	tcl/angle_cossquare_tcl.cpp tcl/angle_cossquare_tcl.hpp \
	tcl/angledist_tcl.cpp tcl/angledist_tcl.hpp \
	tcl/dihedral_tcl.cpp tcl/dihedral_tcl.hpp \
	tcl/endangledist_tcl.cpp tcl/endangledist_tcl.hpp \
	tcl/fene_tcl.cpp tcl/fene_tcl.hpp \
	tcl/harmonic_tcl.cpp tcl/harmonic_tcl.hpp \
	tcl/overlap_tcl.cpp tcl/overlap_tcl.hpp \
	tcl/subt_lj_tcl.cpp tcl/subt_lj_tcl.hpp \
	tcl/object-in-fluid/area_force_local_tcl.cpp \
	tcl/object-in-fluid/area_force_local_tcl.hpp \
	tcl/object-in-fluid/area_force_global_tcl.cpp \
	tcl/object-in-fluid/area_force_global_tcl.hpp \
	tcl/object-in-fluid/bending_force_tcl.cpp \
	tcl/object-in-fluid/bending_force_tcl.hpp \
	tcl/object-in-fluid/stretching_force_tcl.cpp \
	tcl/object-in-fluid/stretching_force_tcl.hpp \
	tcl/object-in-fluid/stretchlin_force_tcl.cpp \
	tcl/object-in-fluid/stretchlin_force_tcl.hpp \
	tcl/object-in-fluid/volume_force_tcl.cpp \
	tcl/object-in-fluid/volume_force_tcl.hpp

# Coulomb methods
libEspresso_a_SOURCES += \
	tcl/debye_hueckel_tcl.cpp tcl/debye_hueckel_tcl.hpp \
	tcl/elc_tcl.cpp tcl/elc_tcl.hpp \
	tcl/magnetic_non_p3m_methods_tcl.cpp tcl/magnetic_non_p3m_methods_tcl.hpp \
	tcl/maggs_tcl.cpp tcl/maggs_tcl.hpp \
	tcl/mmm1d_tcl.cpp tcl/mmm1d_tcl.hpp \
	tcl/mmm2d_tcl.cpp tcl/mmm2d_tcl.hpp \
	tcl/p3m-dipolar_tcl.cpp tcl/p3m-dipolar_tcl.hpp \
	tcl/p3m_tcl.cpp tcl/p3m_tcl.hpp \
	tcl/reaction_field_tcl.cpp tcl/reaction_field_tcl.hpp \
	tcl/mdlc_correction_tcl.cpp tcl/mdlc_correction_tcl.hpp	

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
Espresso_SOURCES = tcl/scriptsdir.cpp tcl/main.cpp
Espresso_LDADD = libEspresso.a

bin_PROGRAMS = Espresso.install
Espresso_install_CPPFLAGS = -D ESPRESSO_SCRIPTS_DEFAULT=\"$(scriptsdir)\"
Espresso_install_SOURCES = tcl/scriptsdir.cpp tcl/main.cpp
Espresso_install_LDADD = libEspresso.a

ESPRESSO = `echo Espresso | sed '$(transform)'`$(EXEEXT)
ESPRESSO_INSTALL = `echo Espresso.install | sed '$(transform)'`$(EXEEXT)
# rename Espresso after installation
install-exec-hook:
	$(am__mv) \
		$(DESTDIR)$(bindir)/$(ESPRESSO_INSTALL) \
		$(DESTDIR)$(bindir)/$(ESPRESSO)

uninstall-local:
	-rm -f $(DESTDIR)$(bindir)/$(ESPRESSO)

endif

#################################################################
# Handling of myconfig.hpp
#################################################################
nodist_libEspresso_a_SOURCES += myconfig-final.hpp
BUILT_SOURCES = myconfig-final.hpp
CLEANFILES = myconfig-final.hpp

myconfig-final.hpp: FORCE
	@echo -n "Determining myconfig file..."
	@config_files="\
	  $(top_builddir)/$(myconfig) \
	  $(top_srcdir)/$(myconfig)"; \
	if test $(myconfig) = myconfig.hpp; then \
	  config_files="$$config_files $(top_srcdir)/src/myconfig-default.hpp"; \
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
# config-features.cpp must occur at the beginning of the SOURCES so that it is compiled first
#libEspresso_a_SOURCES += config-features.cpp
EXTRA_DIST = \
	features.def \
	config-features.hpp config-features.cpp \
	config-doxygen.hpp

if HAVE_PYTHON

BUILT_SOURCES +=  config-features.hpp config-features.cpp
CLEANFILES += config-features.hpp config-features.cpp config-doxygen.hpp
config-features.hpp config-features.cpp: features.def $(top_srcdir)/config/gen_featureconfig.py
	$(PYTHON) $(top_srcdir)/config/gen_featureconfig.py $(srcdir)/features.def config-features.hpp config-features.cpp

config-doxygen.hpp: features.def $(top_srcdir)/config/gen_doxyconfig.py
	$(PYTHON) $(top_srcdir)/config/gen_doxyconfig.py $(srcdir)/features.def config-doxygen.hpp

else

config-doxygen.hpp config-features.hpp config-features.cpp: features.def
	@echo "Python was not found in your PATH."
	@echo "If you change the set of features, you need Python!"
	@echo "Install it and rerun configure."

endif

doc doxygen: config-doxygen.hpp

#################################################################
# Handling of the version
#################################################################
libEspresso_a_SOURCES += config-version.cpp 
BUILT_SOURCES += config-version.cpp
CLEANFILES += config-version.cpp
config-version.cpp: FORCE
	VERSION=`cd $(top_srcdir); sh config/genversion.sh -c`; \
	if ! test -e "config-version.cpp" || \
		test "$$VERSION" != "`cat config-version.cpp`"; then \
		echo "$$VERSION" > config-version.cpp; \
	fi

#################################################################
# Handling of the MPI fake implementation
#################################################################
if MPI_FAKE
libEspresso_a_SOURCES += mpifake/mpi.h mpifake/mpi.cpp
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
	cuda_init_cuda.cu \
	cuda_common_cuda.cu \
	lbgpu_cuda.cu \
	p3m_gpu_cuda.cu

cuda_interface.o: cuda_interface.hpp myconfig-final.hpp
cuda_init.o: cuda_init.hpp myconfig-final.hpp
lbgpu.o: lbgpu.hpp myconfig-final.hpp
p3m_gpu.o: p3m_gpu.hpp myconfig-final.hpp

if TCL

Espresso_SOURCES += cuda_is_cpp.cpp
Espresso_install_SOURCES += cuda_is_cpp.cpp
CUDA_SOURCES += \
	tcl/cuda_init_tcl.cpp

endif

libEspresso_a_SOURCES += $(CUDA_SOURCES)

EXTRA_DIST += $(CUDA_SOURCES)
endif

.PHONY: FORCE
FORCE:
