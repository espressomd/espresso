# Copyright (C) 2010-2018 The ESPResSo project
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
#!/usr/bin/env python

"""
Module for methods common to several sample files.
"""

import tempfile


def open(*args, **kwargs):
    """
    Return a temporary file object.

    Parameters
    ----------

    *args: anytype
           All arguments are discarded

    **args: anytype
            All keyword arguments are discarded

    Returns
    -------
    file-like object
    """
    return tempfile.TemporaryFile(mode='w')
