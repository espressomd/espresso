EXTRA_DIST=common/common.tex
SUBDIRS=tut2

.PHONY: FORCE doc tutorials tut2
doc: tutorials
tutorials: tut2
tut2: FORCE
	cd tut2; $(MAKE) --print-directory $@
FORCE:

