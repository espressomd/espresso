# Copyright (C) 2012,2013,2014,2015,2016 The ESPResSo project
# Copyright (C) 2011 Olaf Lenz
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
if DEVEL_SRC

dg_TEXFILES = \
	dg.tex

dg_FILES = \
	$(dg_TEXFILES) \
	$(dg_FIGURES) \
	todonotes.sty

# bibliography.bib

.PHONY: doc dg

doc: dg
dg: dg.pdf
dg.pdf: $(dg_FILES)

if HAVE_LATEX
	sh ../latexit.sh $(srcdir) dg;
else
	@echo "No complete LaTeX-installation was not found in your PATH."
	@echo "Can't build the Developer's Guide without pdflatex, makeindex and bibtex."
	@echo "Install these and rerun configure."
endif

CLEANFILES = dg.pdf

MOSTLYCLEANFILES = \
	quickref.inp\
	*.aux *.aux.bak\
	dg.bbl dg.blg \
	dg.idx dg.idx.bak dg.ilg dg.ind\
	dg.log dg.out dg.toc dg.tdo

endif
