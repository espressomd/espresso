# generate the file lists via
#  find pages -name CVS -prune -o -printf '\t%p \\\n'
#  find figs -name CVS -prune -o -printf '\t%p \\\n'
EXTRA_DIST=\
	background_errors.sh background_errors.awk \
	doxygen.cfg header.html \
        pages/build.doc \
        pages/cvs.doc \
        pages/doxygen.doc \
        pages/header.doc \
        pages/new_files.doc \
        pages/progguide.doc \
        pages/testsuite.doc \
        pages/ug.doc \
        figs/bond_angle.fig \
        figs/bond_angle.gif \
        figs/datastorage.gif \
        figs/dihedral_angle.fig \
        figs/dihedral_angle.gif \
        figs/dihedral_angle.pdf \
        figs/directions.fig \
        figs/directions.gif \
        figs/elc_errordist.gif \
        figs/ghost_cells.fig \
        figs/ghost_cells.gif \
        figs/ghost_communication.fig \
        figs/ghost_communication.gif \
        figs/linked_cells.fig \
        figs/linked_cells.gif \
        figs/logo.gif \
        figs/move_to_p_buf.fig \
        figs/move_to_p_buf.gif \
        figs/particles.fig \
        figs/particles.gif

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
