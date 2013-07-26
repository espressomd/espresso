# Copyright (C) 2011 Olaf Lenz
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
MYCONFIGS = \
	myconfig-empty.hpp \
	myconfig-nofftw.hpp \
	myconfig-compileonly.hpp \
	\
	myconfig-maxset.hpp \
	myconfig-restcompile1.hpp \
	myconfig-restcompile2.hpp \
	myconfig-LB.hpp \
	myconfig-LBGPU.hpp \
	myconfig-LBGPU-SHANCHEN.hpp \
	myconfig-tunableslip.hpp \
        myconfig-adress.hpp

EXTRA_DIST = $(MYCONFIGS) \
	check_myconfig_complete.py

all-local:
	for file in $(MYCONFIGS); do \
		if ! test -e $$file; then \
			$(LN_S) $(srcdir)/$$file .; \
		fi; \
	done

if HAVE_PYTHON

check-local: $(MYCONFIGS)
	cd $(srcdir); $(PYTHON) check_myconfig_complete.py \
	  $(abs_top_srcdir)/src/features.def $(MYCONFIGS)

else

check:
	@echo "Python was not found in your PATH."
	@echo "Can't check whether test configs are complete."

endif

