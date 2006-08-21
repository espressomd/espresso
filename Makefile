objdir := obj-$(shell config/config.guess)

all:
	cd ${objdir}; $(MAKE) --print-directory $(MAKEFLAGS) $@

# dist:
# 	cd ${objdir}; $(MAKE) --print-directory $(MAKEFLAGS) $@

%:
	cd ${objdir}; $(MAKE) --print-directory $(MAKEFLAGS) $@
