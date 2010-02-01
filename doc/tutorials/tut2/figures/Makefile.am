tut2_FIGURES=neutral-rho.pdf nonneutral-rho.pdf nacl-rdf.pdf salt.png

EXTRA_DIST=$(tut2_FIGURES)

.PHONY: tut2.pdf

tut2.pdf: FORCE
if HAVE_LATEX
	for file in $(tut2_FIGURES); do \
	  test -e $$file || \
	    $(LN_S) $(srcdir)/$$file .; \
	done
	test -e logo.png || \
		$(LN_S) $(srcdir)/../../../ug/figures/logo.png .;
endif

FORCE:
