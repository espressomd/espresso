# Copyright (C) 2014 Olaf Lenz
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
BASENAME=03-object_in_fluid
EXTRA_DIST=$(BASENAME).pdf \
	scripts/boundaries.tcl \
	scripts/simulation.tcl \
	scripts/input/cell_nodes.dat \
	scripts/input/cell_triangles.dat \
	scripts/input/sphere_nodes.dat \
	scripts/input/sphere_triangles.dat

if DEVEL_SRC
tutorial_TEXFILES=$(BASENAME).tex

tutorial_FIGURES = \
	figures/1.png \
	figures/2.png \
	figures/3.png \
	figures/channel.png \
	figures/cylinder.png \
	figures/rhomboid.png

tutorial_FILES = \
	$(tutorial_TEXFILES) \
	$(tutorial_FIGURES)

.PHONY: doc $(BASENAME).pdf

doc: tutorials
tutorials: $(BASENAME).pdf
$(BASENAME).pdf: $(tutorial_FILES) ../common/common.tex 
if HAVE_LATEX
	sh ../../latexit.sh $(srcdir):$(srcdir)/../common $(BASENAME)
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
