########### load platform dependent part
PLATFORM=$(shell uname -s)
include Makefile.$(PLATFORM)

########### list of source files
CSOURCES= main initialize global communication binary_file interaction_data \
	  verlet grid integrate cells ghosts forces debug particle_data \
	  thermostat statistics vmdsock imd p3m fft random blockfile blockfile_tcl \
	  polymer
CXXSOURCES=

LIBOBJECTS=c_blockfile.o

########### RULES
#################
TCLMD_CFLAGS=-DTCL_FILE_IO
OBJECTS=$(CSOURCES:%=%.o) $(CXXSOURCES:%=%.o)
CFILES=$(CSOURCES:=.c)
CXXFILES=$(CXXSOURCES:=.cc)

default: $(PLATFORM) $(PLATFORM)/tcl_md $(PLATFORM)/libtcl_md.a
all: $(PLATFORM) $(PLATFORM)/tcl_md $(PLATFORM)/libtcl_md.a

########### documentation
doc: doxygen_header $(CFILES) $(CXXFILES)
	doxygen doxygen_config
	# (cd doc/latex; make)

########### output directory
$(PLATFORM):
	mkdir -p $(PLATFORM)

########### final target
$(PLATFORM)/tcl_md: $(OBJECTS)
	(cd $(PLATFORM); $(LINK) $(LDFLAGS) -o tcl_md $(OBJECTS) $(LDLIBS) )

$(PLATFORM)/libtcl_md.a: $(LIBOBJECTS)
	(cd $(PLATFORM); ar -crs libtcl_md.a $(LIBOBJECTS) )

########### clean
clean:
	rm -f *~
	(cd $(PLATFORM); rm -f $(OBJECTS) )
mostclean: clean
	rm -rf $(PLATFORM)

########### dependencies
dep: 
	$(MAKE) $(PLATFORM)
	rm -f $(PLATFORM)/.depend
	$(MAKE) $(PLATFORM)/.depend

$(PLATFORM)/.depend:
	mkdir -p $(PLATFORM)
	rm -f $@
	touch $@
	$(DEPEND) -f $@ -- $(CFLAGS) -- $(CFILES) $(CXXFILES) 2>/dev/null

include $(PLATFORM)/.depend

########## implicit rules
vpath %.o  $(PLATFORM)

%.o: %.c
	$(CC) $(CFLAGS) $(TCLMD_CFLAGS) -c -o $(PLATFORM)/$@ $<

c_blockfile.o: blockfile.c
	$(CC) $(CFLAGS) -c -o $(PLATFORM)/$@ $<

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(TCLMD_CFLAGS) -c -o $(PLATFORM)/$@ $<
