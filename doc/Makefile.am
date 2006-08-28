EXTRA_DIST=\
	background_errors.awk doxygen.cfg \
	text figs

.PHONY: doc doxygen doxygen-html
doc: doxygen
doxygen: doxygen-html
doxygen-html:
	$(AWK) -f $(srcdir)/background_errors.awk $(top_srcdir)/*.c $(top_srcdir)/*.h > background_errors.doc
	SRCDIR=$(top_srcdir) $(DOXYGEN) $(srcdir)/doxygen.cfg | grep -ve "^\(Generating\|Parsing\|Preprocessing\)"
	@echo "***************************************************************************"
	@echo "* The developer's guide is now accessible at"
	@echo "*  file://`pwd`/html/index.html"
	@echo "***************************************************************************"

CLEANFILES=background_errors.doc

clean-local:
	-rm -rf html
