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
	myconfig-sample-body.h

all: myconfig-sample-body.h

if HAVE_PYTHON

myconfig-sample-body.h: $(top_srcdir)/src/features.def
	$(PYTHON) -B $(srcdir)/gen_sampleconfig.py $(top_srcdir)/src/features.def myconfig-sample-body.h

else

myconfig-sample-body.h: $(top_srcdir)/src/features.def
	@echo "Python was not found in your PATH."
	@echo "If you change the set of features, you need Python!"
	@echo "Install it and rerun configure."
	test -e $(top_srcdir)/config/myconfig-sample-body.h && $(LN_S) $(top_srcdir)/config/myconfig-sample-body.h

endif
