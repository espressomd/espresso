# Copyright (C) 2012 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
EXTRA_DIST = \
	gen_featureconfig.py \
	gen_doxyconfig.py \
	gen_sampleconfig.py \
	featuredefs.py \
	myconfig-sample-body.h \
	genversion.sh

CLEANFILES = myconfig-sample-body.h

all-local: myconfig-sample-body.h

if HAVE_PYTHON
myconfig-sample-body.h: FORCE
	$(PYTHON) $(srcdir)/gen_sampleconfig.py $(top_srcdir)/src/features.def myconfig-sample-body.h
else
myconfig-sample-body.h: FORCE
	if test -e $(srcdir)/myconfig-sample-body.h; then \
		if ! test -e myconfig-sample-body.h; then \
			$(LN_S) $(srcdir)/myconfig-sample-body.h .; \
	else \
		echo "ERROR: Python was not found in your PATH."; \
		echo "ERROR: When compiling development code or changing the set of features, you need Python!"; \
		echo "ERROR: Install it and rerun configure."; \
		exit 1; \
		fi; \
	fi
endif

.PHONY: FORCE
FORCE:
