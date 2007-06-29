ug_TEXFILES=\
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
	deserno.tex


ug_FILES=$(ug_TEXFILES) \
	underscore.sty \
	bibliography.bib

#data_DATA=ug.pdf

EXTRA_DIST= $(ug_FILES) \
	ug.pdf \
	check_consistency.sh assemble_quickref.awk \
	latexmk latexmk.1 latexmkrc

.PHONY: doc ug

doc: ug
ug: ug.pdf
	cd figures; $(MAKE) --print-directory $@
ug.pdf: $(ug_FILES)
if HAVE_LATEX
	for file in $(ug_FILES); do \
	  test -e $$file || \
	    $(LN_S) $(srcdir)/$$file .; \
	done
	test -d figures || mkdir figures
	fi
	for file in $(srcdir)/figures/*; do \
	  if test ! -e $$file; then \
	    $(LN_S) $$file figures/; \
	  fi \
	done
	$(AWK) -f $(srcdir)/assemble_quickref.awk $(ug_TEXFILES) > quickref.inp
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
	quickref.inp \
	*.aux *.aux.bak\
	ug.idx ug.idx.bak ug.ilg ug.ind\
	ug.log ug.out ug.toc

#MOSTLYCLEANFILES=
