# This makfile tests whether "make" is executed in the object
# directory. 
# If it is, it simply includes the Automakefile Makefile-am.
# Otherwise, it assumes that it is executed in the source
# directory and tries to run make in the objdir.

SUBMAKEFILE := Makefile-am

# Test if we are in the objdir
ifeq ($(shell test -f ./$(SUBMAKEFILE) && echo "yes"),yes)

# if we are, simply use the automakefile
include $(SUBMAKEFILE)

else

# else test if the objdir is already in place and configured
objdir := obj-$(shell config/config.guess)
ifeq ($(shell test -d $(objdir) && test -f $(objdir)/$(SUBMAKEFILE) && echo "yes"),yes)

# if it is, run make in the objdir
all %:
	cd $(objdir); $(MAKE) --print-directory $(MAKEFLAGS) $@

# otherwise try to run configure
else

all:
	@echo "No configured directory found, running configure..."
	@./configure
	@$(MAKE) $(MAKEFLAGS) $@

endif
endif