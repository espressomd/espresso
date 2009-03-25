SUBDIRS = dg ug tutorials maggs_doc

.PHONY: FORCE doc dg ug tutorials
doc: ug dg tutorials
dg: FORCE
	cd dg; $(MAKE) --print-directory $@
ug: FORCE
	cd ug; $(MAKE) --print-directory $@
tutorials: FORCE
	cd tutorials; $(MAKE) --print-directory $@
FORCE:

