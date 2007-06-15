ug_SOURCES=\
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

ug_FIGURES=\
	logo.jpg \
	nacl-rdf.pdf \
	salt.png

data_DATA=ug.pdf

EXTRA_DIST=$(ug_SOURCES) figures ug.pdf
.PHONY: doc ug

doc: ug
ug: ug.pdf
ug.pdf: $(ug_SOURCES)
if HAVE_PDFLATEX
	for file in $(ug_SOURCES); do \
	  if test ! -e $$file; then $(LN_S) $(srcdir)/$$file .; fi;\
	done
	mkdir -p figures
	for file in $(ug_FIGURES); do \
	  if test ! -e figures/$$file; then $(LN_S) ../$(srcdir)/figures/$$file figures; fi;\
	done
	$(PDFLATEX) ug.tex
	makeindex ug
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
