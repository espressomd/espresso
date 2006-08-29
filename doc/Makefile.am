EXTRA_DIST=\
	background_errors.awk doxygen.cfg \
	text figs

.PHONY: doc doxygen doxygen-html
doc: doxygen
doxygen: doxygen-html
doxygen-html:
	SRCDIR=$(srcdir) AWK=$(AWK)\
	  $(SHELL) $(srcdir)/background_errors.sh $(top_srcdir)/*.c $(top_srcdir)/*.h 
	SRCDIR=$(top_srcdir) \
	  $(DOXYGEN) $(srcdir)/doxygen.cfg \
	  | grep -ve "^\(Generating\|Parsing\|Preprocessing\)"
	@echo "***************************************************************************"
	@echo "* The developer's guide is now accessible at"
	@echo "*  file://`pwd`/html/index.html"
	@echo "***************************************************************************"

CLEANFILES=background_errors.doc background_errors.sorted background_errors.unsorted

clean-local:
	-rm -rf html
