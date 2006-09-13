EXTRA_DIST=background_errors.awk doxygen.cfg pages figs

.PHONY: doc dg
doc: dg
dg:
if HAVE_DOXYGEN
	SRCDIR=$(srcdir) AWK=$(AWK)\
	  $(SHELL) $(srcdir)/background_errors.sh $(top_srcdir)/*.c $(top_srcdir)/*.h 
	SRCDIR=$(top_srcdir) PACKAGE="$(PACKAGE_STRING)"\
	  $(DOXYGEN) $(srcdir)/doxygen.cfg \
	  | grep -ve "^\(Generating\|Parsing\|Preprocessing\)"
	@echo "***************************************************************************"
	@echo "* The developer's guide is now accessible at"
	@echo "*  file://`pwd`/html/index.html"
	@echo "***************************************************************************"
else
	@echo "doxygen was not found in your PATH."
	@echo "Can't build developer's guide without doxygen."
	@echo "Install doxygen and rerun configure."

endif

MOSTLYCLEANFILES=background_errors.doc background_errors.sorted background_errors.unsorted
clean-local:
	-rm -rf html
