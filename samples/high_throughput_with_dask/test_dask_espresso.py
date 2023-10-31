#
# Copyright (C) 2023 The ESPResSo project
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

"""Unit tests for the ``dask_espresso`` python module."""

import dask_espresso as de
import numpy as np
import sys


def test_encode_decode():
    data = {"a": 1, "b": np.random.random((3, 3))}
    encoded = de.encode_transport_data(data)
    decoded = de.decode_transport_data(encoded)
    for k, v in decoded.items():
        assert np.all(data[k]) == np.all(v)
    assert list(data.keys()) == list(decoded.keys())


def test_espresso_runner():
    data = {"hello": "world"}
    result = de.dask_espresso_task(sys.executable, "echo.py", **data).compute()
    assert result == {"hello": "world", "processed": True}
