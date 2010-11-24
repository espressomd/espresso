# Copyright (C) 2007,2008,2009,2010 Olaf Lenz, Axel Arnold
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
ug_TEXFILES = \
	ug.tex \
	introduction.tex \
	firststeps.tex \
	installation.tex \
	part.tex \
	inter.tex \
	setup.tex \
	run.tex \
	analysis.tex \
	io.tex \
	aux.tex \
	internal.tex \
	contributing.tex \
	quickref.tex \
	features.tex \
	examples.tex \
	deserno.tex \
	mmm.tex \
        mbtools.tex

ug_FIGURES = \
	figures/logo.png \
	figures/nacl-rdf.pdf \
	figures/dihedral-angle.pdf \
	figures/salt.png \
	figures/hbond.pdf

ug_FILES = \
	$(ug_TEXFILES) \
	$(ug_FIGURES) \
	bibliography.bib

EXTRA_DIST = \
	$(ug_FILES) \
	latexit.sh \
	ug.pdf \
	check_consistency.sh assemble_quickref.awk \
	latexmk latexmk.1 latexmkrc

.PHONY: doc ug

doc: ug
ug: ug.pdf
ug.pdf: $(ug_FILES)
if HAVE_LATEX
	cd $(srcdir); $(AWK) -f assemble_quickref.awk $(ug_TEXFILES) > $(abs_builddir)/quickref.inp
	sh latexit.sh
else
	@echo "No complete LaTeX-installation was not found in your PATH."
	@echo "Can't build the User's Guide without pdflatex, makeindex and bibtex."
	@echo "Install these and rerun configure."
endif

CLEANFILES=ug.pdf 
MOSTLYCLEANFILES=\
	quickref.inp\
	*.aux *.aux.bak\
	ug.bbl ug.blg \
	ug.idx ug.idx.bak ug.ilg ug.ind\
	ug.log ug.out ug.toc
