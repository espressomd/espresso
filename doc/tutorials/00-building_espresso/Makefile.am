BASENAME=00-building_espresso
EXTRA_DIST=$(BASENAME).pdf

if DEVEL_SRC
tutorial_TEXFILES=$(BASENAME).tex

tutorial_FIGURES =

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
