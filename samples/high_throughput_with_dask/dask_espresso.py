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

"""Helper functions to use ESPResSo with Dask."""

import pickle
import base64
import sys
import subprocess
import logging
import dask


def encode_transport_data(data):
    """
    Use ``pickle`` and ``base64`` to convert the provided data to a string
    which can be passed safely between the Dask scheduler, worker and ESPResSo.
    """
    return base64.b64encode(pickle.dumps(data)).decode("utf-8")


def decode_transport_data(encoded_data):
    """
    Convert the transport data back to a Python object via ``base64``
    and ``pickle``.
    """
    pickle_data = base64.b64decode(encoded_data)
    return pickle.loads(pickle_data)


def get_data_from_stdin():
    return decode_transport_data(sys.stdin.read())


@dask.delayed
def dask_espresso_task(pypresso, script, **kwargs):
    """
    Run ESPResSo asynchronously as a Dask task.

    pypresso: :obj:`str`
        Path to pypresso
    script: :obj:`str`
        Simulation script to run with pypresso
    kwargs:
        The keyword arguments are passed encoded and sent to
        the standard input of the simulation script.
        Use ``data = get_data_from_stdin()`` to obtain it.
    """

    logger = logging.getLogger(__name__)
    encoded_data = encode_transport_data(kwargs)
    espresso = subprocess.Popen([pypresso, script],
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                text=True)
    espresso.stdin.write(encoded_data)
    out, err = espresso.communicate()
    if err != "":
        logger.warning("STDERR output from ESPResSo\n", err)
    return decode_transport_data(out)    
