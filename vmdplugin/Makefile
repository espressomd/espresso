# Location of the VMD include files vmdplugin.h and molfile_plugin.h.
# Typically, these files can be found in the subdir 
# lib/vmd/plugins/include of the VMD installation directory.
VMDDIR=/usr/local
VMDINCLUDES=$(VMDDIR)/lib/vmd/plugins/include

CPPFLAGS=-I $(VMDDINCLUDES)

CC=gcc
CFLAGS=-g -O0

all: vtftest vtfplugin.so

vtftest: vtfplugin.o
vtfplugin.so: vtfplugin.o
	$(CC) $(LDFLAGS) -shared $< -o vtfplugin.so
