# Copyright (C) 2010-2019 The ESPResSo project
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
from . cimport c_profiler


def begin_section(name):
    """
    Start named section in profiler.

    Parameters
    ----------

    name : :obj:`str`
        Name of the section
    """

    c_profiler.begin_section(name)


def end_section(name):
    """
    End named section in profiler.

    Parameters
    ----------

    name : :obj:`str`
        Name of the section
    """

    c_profiler.end_section(name)
