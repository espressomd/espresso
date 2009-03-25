SUBDIRS = figures

tut2_TEXFILES=\
	tut2.tex

tut2_FILES=$(tut2_TEXFILES)

EXTRA_DIST= $(tut2_FILES) \
	tut2.pdf

.PHONY: doc tut2

doc: tut2
tut2: tut2.pdf
tut2.pdf: $(tut2_FILES) ../common/common.tex
if HAVE_LATEX
	cd figures; $(MAKE) --print-directory $@
	for file in $(tut2_FILES); do \
	  test -e $$file || \
	    $(LN_S) $(srcdir)/$$file .; \
	done
	test -e common.tex || $(LN_S) $(srcdir)/../common/common.tex .
	$(PDFLATEX) tut2.tex
	$(PDFLATEX) tut2.tex
else
	@echo "No complete LaTeX-installation was not found in your PATH."
	@echo "Can't build the tutorial without pdflatex, makeindex and bibtex."
	@echo "Install these and rerun configure."
endif

CLEANFILES=tut2.pdf
MOSTLYCLEANFILES=\
	*.aux *.aux.bak\
	tut2.idx tut2.idx.bak tut2.ilg tut2.ind\
	tut2.log tut2.out tut2.toc
