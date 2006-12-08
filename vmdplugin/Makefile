# Location of the VMD include files vmdplugin.h and molfile_plugin.h.
# Typically, these files can be found in the subdir 
# plugins/include of the VMD installation directory.

VMDDIR=/home/users/simbio/software/lib/vmd
VMDINCLUDES=$(VMDDIR)/plugins/include
# comment this line, if zlib is not available
HAVE_ZLIB=1

CPPFLAGS=-I $(VMDINCLUDES)

CC=gcc
CFLAGS=-g -O0
SHLD=$(CC)
SHLDFLAGS=-shared $(LDFLAGS)

# if you use MAX OS X, use the following
#CC=cc
#CFLAGS=-Os -fPIC -dynamic
#SHLD=$(CC)
#SHLDFLAGS=-bundle

ifdef HAVE_ZLIB
# if you want to enable compressed files, use these
CFLAGS += -D HAVE_ZLIB
LDFLAGS += -lz
endif

all: vtftest vtfplugin.so

vtftest: vtfplugin.o
vtfplugin.so: vtfplugin.o
	$(SHLD) $(SHLDFLAGS) $< -o vtfplugin.so

clean:
	-rm vtfplugin.o vtfplugin.so vtftest