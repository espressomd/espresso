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
	myconfig-default.h \
	p3m-charges.c p3m-dipoles.h \
	p3m-dipoles.c p3m-charges.h \
	p3m-assignment.c

# List the sources of the Espresso binary here
libEspresso_a_SOURCES = \
	main.c \
	config.c config.h \
	initialize.c initialize.h \
	global.c global.h \
	communication.c communication.h \
	binary_file.c binary_file.h \
	interaction_data.c interaction_data.h\
	verlet.c verlet.h \
	grid.c grid.h \
	integrate.c integrate.h \
	cells.c cells.h \
	ghosts.c ghosts.h \
	forces.c forces.h \
	rotation.c rotation.h \
	debug.c debug.h \
	particle_data.c	particle_data.h \
	thermostat.c thermostat.h \
	dpd.c dpd.h \
	statistics.c statistics.h \
	statistics_chain.c statistics_chain.h \
	energy.c energy.h \
	pressure.c pressure.h \
	vmdsock.c vmdsock.h \
	imd.c imd.h \
	iccp3m.c iccp3m.h \
	p3m.c p3m.h \
	magnetic_non_p3m_methods.c magnetic_non_p3m_methods.h \
	ewald.c ewald.h \
	fft.c fft.h \
	random.c random.h \
	blockfile.c blockfile.h \
	blockfile_tcl.c blockfile_tcl.h \
	polymer.c polymer.h \
	specfunc.c specfunc.h \
	tuning.c tuning.h \
	uwerr.c	uwerr.h \
	parser.c parser.h \
	domain_decomposition.c domain_decomposition.h \
	nsquare.c nsquare.h \
	layered.c layered.h \
	mmm-common.c mmm-common.h \
	mmm1d.c mmm1d.h \
	mmm2d.c	mmm2d.h \
	modes.c	modes.h \
	topology.c topology.h \
	nemd.c nemd.h \
	statistics_cluster.c statistics_cluster.h \
	elc.c elc.h \
	mdlc_correction.c  mdlc_correction.h \
	statistics_molecule.c statistics_molecule.h \
	errorhandling.c	errorhandling.h \
	constraint.c constraint.h \
	maggs.c maggs.h \
	mol_cut.h \
	rattle.c rattle.h \
	molforces.c molforces.h \
	virtual_sites.c virtual_sites.h \
	virtual_sites_com.c virtual_sites_com.h \
	virtual_sites_relative.c virtual_sites_relative.h \
	metadynamics.c metadynamics.h \
        npt.h \
	lb.c lb.h \
	lb-d3q18.h lb-d3q19.h \
	bin.c bin.h \
	lattice.c lattice.h \
	halo.c halo.h \
	statistics_fluid.c statistics_fluid.h \
	lb-boundaries.c lb-boundaries.h \
	lb_boundaries_gpu.c lb_boundaries_gpu.h \
	lbgpu_cfile.c \
	utils.c utils.h \
	angle.h \
	pwdist.h \
	angledist.h \
	endangledist.h \
	buckingham.h \
	comfixed.h \
	comforce.h \
	debye_hueckel.h \
	reaction_field.h \
	dihedral.h \
	fene.h \
	gb.h \
	harmonic.h \
	hertzian.h \
	imd.h \
	ljcos2.h \
	ljcos.h \
	lj.h \
	ljgen.c ljgen.h \
	steppot.h \
	bmhtf-nacl.h \
	morse.h \
	polynom.h \
	soft_sphere.h \
	subt_lj.h \
	tab.h \
	overlap.h \
	ljangle.h \
	adresso.c adresso.h \
        tunable_slip.h

# All objects are compiled into the library
noinst_LIBRARIES = libEspresso.a

#################################################################
# Handling of myconfig.h
#################################################################
FORCE:
myconfig-final.h: FORCE
	config_files="\
	  $(top_builddir)/$(myconfig) \
	  $(top_srcdir)/$(myconfig)"; \
	if test $(myconfig) = myconfig.h; then \
	  config_files="$$config_files $(top_srcdir)/src/myconfig-default.h"; \
	fi; \
	for file in $$config_files; do \
	  if test -e $$file; then \
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

$(libEspresso_a_OBJECTS): myconfig-final.h

CLEANFILES = myconfig-final.h

#################################################################
# Handling of the MPI fake implementation
#################################################################
if MPI_FAKE
libEspresso_a_SOURCES += mpifake/mpi.h mpifake/mpi.c
# mpifake should come before any system includes
AM_CPPFLAGS = -I$(srcdir)/mpifake
endif

##################################################
# CUDA rules
##################################################
if CUDA
SUFFIXES=.cu
.cu.o:
	$(NVCC) -c $(NVCCFLAGS) $(CPPFLAGS) $(INCLUDES) $(DEFS) -o $@ $<

CUDA_SOURCES= \
	cuda_init.cu \
	cuda_init.h \
	lbgpu.cu lbgpu.h

libEspresso_a_SOURCES += $(CUDA_SOURCES)
endif
