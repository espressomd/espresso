ug_FIGURES=logo.png nacl-rdf.pdf salt.png

EXTRA_DIST=$(ug_FIGURES)

ug: FORCE
if HAVE_LATEX
	for file in $(ug_FIGURES); do \
	  test -e $$file || \
	    $(LN_S) $(srcdir)/$$file .; \
	done
endif

FORCE: