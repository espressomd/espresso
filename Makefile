#  This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
#  It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
#  and by which you are legally bound while utilizing this file in any form or way.
#  There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  You should have received a copy of that license along with this program;
#  if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
#  write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
#  Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
# 

########### load platform dependent part

PLATFORM=$(shell config/config.guess)
OUTDIR=obj-$(PLATFORM)
include Makefile.$(PLATFORM)

########### list of fix source files
CSOURCES+= main config initialize global communication binary_file interaction_data \
	  verlet grid integrate cells ghosts forces rotation debug particle_data \
	  thermostat statistics statistics_chain energy pressure vmdsock imd \
	  p3m fft random blockfile blockfile_tcl polymer specfunc mmm1d tuning \
	  uwerr parser domain_decomposition nsquare layered mmm-common mmm2d \
	  modes topology nemd statistics_cluster elc statistics_molecule \
	  errorhandling constraint maggs rattle molforces lb bin

# objects to include in the Espresso library
LIBOBJECTS= c_blockfile.$(OBJEXT)

# flag only used for Espresso, but not the Espresso library
BIN_CFLAGS=-DTCL_FILE_IO

# documentation directories created by doxygen
DOC_RES= doc/html doc/rtf doc/latex doc/man

########### RULES
#################
OBJECTS=$(CSOURCES:%=%.$(OBJEXT))
CFILES=$(CSOURCES:%=%.c)
HEADERS=$(shell ls *.h)
DOCFILES=$(shell ls doc/text/*.doc) doc/text/header.html 

default: all
all: $(OUTDIR)/Espresso_bin$(EXEEXT) $(OUTDIR)/libEspresso.a

install: all
	config/install-sh -c config/config.guess $(pkglibdir)/config.guess
	echo "#!/bin/sh" > __common_Espresso
	echo "export ESPRESSO_SOURCE=$(pkglibdir); export ESPRESSO_SCRIPTS=$(pkglibdir)/scripts" >> __common_Espresso
	echo "if test ! -d \$$ESPRESSO_SCRIPTS; then echo 'ESPRESSO scripts directory \$$ESPRESSO_SCRIPTS is missing, reinstall Espresso' 1>&2; exit -1; fi"  >> __common_Espresso
	echo "if test ! -x $(pkglibdir)/config.guess; then echo 'ESPRESSO configuration guess script $(pkglibdir)/config.guess is missing, reinstall Espresso' 1>&2; exit -1; fi"  >> __common_Espresso
	echo "if test ! -x $(pkglibdir)/obj-\`$(pkglibdir)/config.guess\`/Espresso; then echo 'There is no binary for the current hardware, '\`$(pkglibdir)/config.guess\`', recompile and install Espresso for this platform' 1>&2; exit -1; fi"  >> __common_Espresso
	echo "$(pkglibdir)/obj-\`$(pkglibdir)/config.guess\`/Espresso \$$*" >> __common_Espresso
	config/install-sh -m 755 __common_Espresso $(bindir)/Espresso
	config/install-sh -c $(OUTDIR)/Espresso $(pkglibdir)/$(OUTDIR)/Espresso
	config/install-sh -c $(OUTDIR)/Espresso_bin$(EXEEXT) $(pkglibdir)/$(OUTDIR)/Espresso_bin
	for f in scripts/*.tcl; do config/install-sh -c "$$f" "$(pkglibdir)/$$f"; done

install-doc:
	(cd doc/html; for f in *; do ../../config/install-sh -c "$$f" "$(pkglibdir)/html_doc/$$f"; done)

########### dependencies
ifeq ($(DEPEND),makedepend)
$(OUTDIR)/.depend:
	rm -f $@
	touch $@
	makedepend -f $@ -- $(CFLAGS) -- $(CFILES)
else
ifeq ($(DEPEND),mkdep)
$(OUTDIR)/.depend:
	rm -f $@
	touch $@
	export CC=$(CC); mkdep -f $@ $(CFLAGS) $(CFILES)
else
ifeq ($(DEPEND),gcc)
$(OUTDIR)/.depend:
	rm -f $@
	touch $@
	$(CC) -MM -MF $@ $(CFLAGS) $(CFILES)
else
$(OUTDIR)/.depend:
	rm -f $@
	touch $@
	@echo "***************************************************************************************"
	@echo "dependencies cannot be built with this compiler, install makedepend or mkdep or use gcc"
	@echo "***************************************************************************************"
endif
endif
endif

########### documentation
docu: doc/html/index.html
doc/html/index.html: doxygen_config $(DOCFILES) $(CFILES) $(HEADERS)
################### BACKGROUND_ERROR-CODES
	awk -f ./scripts/background_errors.awk *.c *.h
	sort ./doc/text/background_errors_tmp.doc -o ./doc/text/background_errors_tmp.doc
	echo "/** \\page background_errors background_errors resolved" > ./doc/text/background_errors.doc
	echo "<ul>" >> ./doc/text/background_errors.doc
	cat ./doc/text/background_errors_tmp.doc >> ./doc/text/background_errors.doc
	echo "</ul>" >> ./doc/text/background_errors.doc
	echo "*/" >> ./doc/text/background_errors.doc
	rm ./doc/text/background_errors_tmp.doc
################### END OF BACKGROUND_ERROR-CODES
	doxygen doxygen_config | grep -ve "^\(Generating\|Parsing\|Preprocessing\)"

########### output directory
$(OUTDIR):
	mkdir -p $(OUTDIR)

########### targets
$(OUTDIR)/Espresso_bin$(EXEEXT): $(OBJECTS)
	$(MAKE) $(OUTDIR)
	(cd $(OUTDIR); $(LD) $(LDFLAGS) -o Espresso_bin$(EXEEXT) $(OBJECTS) $(LIBS))

$(OUTDIR)/libEspresso.a: $(LIBOBJECTS)
	$(MAKE) $(OUTDIR)
	(cd $(OUTDIR); ar -crs libEspresso.a $(LIBOBJECTS))

########### dependencies
dep:
	rm -f $(OUTDIR)/.depend
	$(MAKE) $(OUTDIR)/.depend

include $(OUTDIR)/.depend

########## implicit rules
vpath %.$(OBJEXT)  $(OUTDIR)

%.$(OBJEXT): %.c
	$(CC) $(CFLAGS) $(BIN_CFLAGS) -c -o $(OUTDIR)/$@ $<

c_blockfile.$(OBJEXT): blockfile.c
	$(CC) $(CFLAGS) -c -o $(OUTDIR)/$@ $<

########### clean
clean:
	rm -f *~
	(cd $(OUTDIR); rm -f $(OBJECTS))
docclean:
	rm -rf $(DOC_RES:=/*)
mostclean: clean docclean
	for platform in obj-*; do \
		if test -d $$platform; then \
			rm -rf $$platform; \
		fi \
	done

########### dist
TARFILE=Espresso-$(shell date -I).tgz
__EXCLUDES= $(DOC_RES:%=--exclude=%) --exclude=autom4te.cache --exclude=internal \
	--exclude=*.avi --exclude=*Espresso-*.tgz --exclude=*~ --exclude=Makefile.*-*-*\
	--exclude=core --exclude=core.* --exclude=*obj-* --exclude=.\#* --exclude=CVS --exclude=TclTutor \
	$(EXCLUDES)

dist:
	(cd ..; tar -vchzf Espresso/$(TARFILE) $(__EXCLUDES) Espresso)

########## tests
test: $(OUTDIR)/Espresso_bin$(EXEEXT)
	cd testsuite; ./test.sh
