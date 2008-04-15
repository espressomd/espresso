ug_FIGURES=logo.png nacl-rdf.pdf dihedral-angle.pdf salt.png hbond.pdf

EXTRA_DIST=$(ug_FIGURES)

.PHONY: ug.pdf

ug.pdf: FORCE
if HAVE_LATEX
	for file in $(ug_FIGURES); do \
	  test -e $$file || \
	    $(LN_S) $(srcdir)/$$file .; \
	done
endif

FORCE:
