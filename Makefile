all:
	cd obj-`config/config.guess`; $(MAKE) --print-directory $(MAKEFLAGS) $@

%:
	cd obj-`config/config.guess`; $(MAKE) --print-directory $(MAKEFLAGS) $@
