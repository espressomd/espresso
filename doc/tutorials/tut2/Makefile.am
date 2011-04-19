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
tut2_TEXFILES=\
	tut2.tex

tut2_FIGURES = \
	figures/neutral-rho.pdf \
	figures/nonneutral-rho.pdf \
	figures/nacl-rdf.pdf \
	figures/salt.png

tut2_FILES = \
	$(tut2_TEXFILES) \
	$(tut2_FIGURES)

EXTRA_DIST = \
	$(tut2_FILES) \
	tut2.pdf

.PHONY: doc tut2

doc: tut2
tut2: tut2.pdf
tut2.pdf: $(tut2_FILES) ../common/common.tex 
if HAVE_LATEX
	sh ../../latexit.sh $(srcdir):$(srcdir)/../common tut2
else
	@echo "No complete LaTeX-installation was not found in your PATH."
	@echo "Can't build the tutorial without pdflatex, makeindex and bibtex."
	@echo "Install these and rerun configure."
endif

CLEANFILES = tut2.pdf
MOSTLYCLEANFILES = \
	*.aux *.aux.bak\
	tut2.idx tut2.idx.bak tut2.ilg tut2.ind\
	tut2.log tut2.out tut2.toc
