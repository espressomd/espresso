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
ug_TEXFILES = \
	ug-dist.tex \
	ug.tex \
	introduction.tex \
	firststeps.tex \
	installation.tex \
	part.tex \
	inter.tex \
	setup.tex \
	run.tex \
	analysis.tex \
	analysis-core.tex \
	io.tex \
	lb.tex \
	memd.tex \
	aux.tex \
	internal.tex \
	contributing.tex \
	quickref.tex \
	features.tex \
	examples.tex \
	mmm.tex \
        mbtools.tex

ug_FIGURES = \
	figures/nacl-rdf.pdf \
	figures/dihedral-angle.fig \
	figures/dihedral-angle.pdf \
	figures/salt.png \
	figures/hbond.pdf \
	figures/diamond.png \
	figures/elc-errordist.pdf \
	figures/fullerene.png \
	figures/maggs-charge-assignment.pdf \
	figures/maggs-initial-scheme.pdf \
	figures/maggs-rotation.pdf \
	figures/correlator_scheme.fig \
	figures/correlator_scheme.pdf

ug_FILES = \
	$(ug_TEXFILES) \
	$(ug_FIGURES) \
	bibliography.bib \
	todonotes.sty

EXTRA_DIST = \
	$(ug_FILES) \
	check_consistency.sh assemble_quickref.awk \
	latexmk latexmk.1 latexmkrc

.PHONY: doc ug FORCE

doc: ug
ug: ug.pdf

dist-hook: ug-dist.pdf
	cp ug-dist.pdf $(distdir)/ug.pdf
	rm ug-dist.pdf
ug-dist: ug-dist.pdf

FORCE:
ug.pdf: FORCE
if HAVE_LATEX
	@echo "Assembling quickref.inp.."
	@cd $(srcdir); \
		$(AWK) -f assemble_quickref.awk $(ug_TEXFILES) \
			> $(abs_builddir)/quickref.inp
	@echo "Running LaTeX.."
	sh ../latexit.sh $(srcdir) ug
else
	@echo "No complete LaTeX-installation was not found in your PATH."
	@echo "Can't build the User's Guide without pdflatex, makeindex and bibtex."
	@echo "Install these and rerun configure."
endif

ug-dist.pdf: FORCE
if HAVE_LATEX
	@echo "Assembling quickref.inp.."
	@cd $(srcdir); \
		$(AWK) -f assemble_quickref.awk $(ug_TEXFILES) \
			> $(abs_builddir)/quickref.inp
	@echo "Running LaTeX.."
	sh ../latexit.sh $(srcdir) ug-dist 1
endif

CLEANFILES = ug.pdf ug-dist.pdf

MOSTLYCLEANFILES = \
	quickref.inp\
	*.aux *.aux.bak\
	*.brf *.bbl *.blg \
	*.idx *.idx.bak *.ilg *.ind\
	*.log *.out *.toc *.tdo
