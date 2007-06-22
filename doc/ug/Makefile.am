ug_FILES=\
	underscore.sty \
	bibliography.bib \
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
	figures

ug_figures_FILES=\
	logo.jpg \
	nacl-rdf.pdf \
	salt.png

data_DATA=ug.pdf

EXTRA_DIST= $(ug_FILES) \
	ug.pdf \
	latexmk latexmk.1 latexmkrc

.PHONY: doc ug

doc: ug
ug: ug.pdf
ug.pdf: $(ug_FILES)
if HAVE_LATEX
	for file in $(ug_FILES); do \
	  if test ! -e $$file; then $(LN_S) $(srcdir)/$$file .; fi;\
	done
	$(PDFLATEX) ug.tex
	$(MAKEINDEX) ug
	$(BIBTEX) ug
	$(PDFLATEX) ug.tex
	$(PDFLATEX) ug.tex
else
	@echo "No complete LaTeX-installation was not found in your PATH."
	@echo "Can't build the User's Guide without pdflatex, makeindex and bibtex."
	@echo "Install these and rerun configure."
endif

CLEANFILES=ug.pdf
MOSTLYCLEANFILES=\
	*.aux *.aux.bak\
	ug.idx ug.idx.bak ug.ilg ug.ind\
	ug.log ug.out ug.toc

#MOSTLYCLEANFILES=
