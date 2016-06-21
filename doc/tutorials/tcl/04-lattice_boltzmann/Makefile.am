# Copyright (C) 2012,2013,2014,2015,2016 The ESPResSo project
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

BASENAME=04-lattice_boltzmann
EXTRA_DIST=$(BASENAME).pdf \
    solutions/lb_stokes_force.tcl \
	solutions/polymer_diffusion.tcl \
	solutions/poisseuille.tcl \
	solutions/single_particle_diffusion.tcl \
	src/lb_stokes_force.tcl \
	src/polymer_diffusion.tcl \
	src/poisseuille.tcl \
	src/single_particle_diffusion.tcl

if DEVEL_SRC
tutorial_TEXFILES=$(BASENAME).tex \
				  latex/diffusion.tex \
				  latex/dragforce.tex \
				  latex/interface.tex \
				  latex/intro.tex \
				  latex/lb.tex \
				  latex/poiseuille.tex

tutorial_FIGURES = \
	figures/latticeboltzmann-grid.pdf  \
	figures/latticeboltzmann-momentumexchange.pdf \
	figures/schlitzpore_3d.pdf \
	figures/diffusion/msd.pdf \
	figures/poiseuille.pdf

tutorial_FILES = \
	$(tutorial_TEXFILES) \
	$(tutorial_FIGURES)

.PHONY: doc $(BASENAME).pdf

doc: tutorials
tutorials: $(BASENAME).pdf
$(BASENAME).pdf: $(tutorial_FILES) ../common/common.tex 
if HAVE_LATEX
	sh ../../latexit.sh $(srcdir):$(top_srcdir)/doc/tutorials/common $(BASENAME)
else
	@echo "No complete LaTeX-installation was not found in your PATH."
	@echo "Can't build the tutorial without pdflatex, makeindex and bibtex."
	@echo "Install these and rerun configure."
endif

CLEANFILES = $(BASENAME).pdf
MOSTLYCLEANFILES = \
	*.aux *.aux.bak\
	$(BASENAME).idx $(BASENAME).idx.bak $(BASENAME).ilg $(BASENAME).ind \
	$(BASENAME).log $(BASENAME).out $(BASENAME).toc \
	$(BASENAME).blg $(BASENAME).bbl
endif
