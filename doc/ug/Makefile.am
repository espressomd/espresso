ug_SOURCES=ug.tex installation.tex logo.jpg

EXTRA_DIST=$(ug_SOURCES) ug.pdf
.PHONY: doc ug

doc: ug
ug: ug.pdf
ug.pdf: $(ug_SOURCES)
if HAVE_PDFLATEX
	for file in $(ug_SOURCES); do \
	  if test ! -e $$file; then $(LN_S) $(srcdir)/$$file .; fi;\
	done
	$(PDFLATEX) ug.tex
	$(PDFLATEX) ug.tex
else
	@echo "pdflatex was not found in your PATH."
	@echo "Can't build user's guide without pdflatex."
	@echo "Install pdflatex and rerun configure."
endif

CLEANFILES=ug.pdf
MOSTLYCLEANFILES=\
	*.aux *.aux.bak\
	ug.idx ug.idx.bak ug.ilg ug.ind\
	ug.log ug.out ug.toc

#MOSTLYCLEANFILES=
