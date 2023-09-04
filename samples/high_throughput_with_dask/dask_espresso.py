import pickle
import base64
import sys
import subprocess
import logging
from dask import delayed


def encode_transport_data(data):
    """
    Uses pickle and base64 to convert the provided data to a string which can be passed safely between the dask scheduler, worker and Esprsso
    """
    return base64.b64encode(pickle.dumps(data)).decode("utf-8")


def decode_transport_data(encoded_data):
    """
    Converts the transport data back to a Python object
    via base64 and pickle
    """
    pickle_data = base64.b64decode(encoded_data)
    return pickle.loads(pickle_data)


def get_data_from_stdin():
    return decode_transport_data(sys.stdin.read())


@delayed
def dask_espresso_task(pypresso, script, **kwargs):
    """
    Runs Esprsso asynchronesouly as a Dask task.

    pypresso: string
      Path to pypresso
    script: string
      Simulation script to run with pypresso
    kwargs:
      The keyword arguments are passed encoded and send to the
      standard input of the simulation script. Use 
      `data = get_data_from_Stdin()` to obtain it.
    """

    logger = logging.getLogger(__name__)
    encoded_data = encode_transport_data(kwargs)
    espresso = subprocess.Popen([pypresso,
                                 script],
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                text=True)
    espresso.stdin.write(encoded_data)
    out, err = espresso.communicate()
    if err != "":
        logger.warning("STDERR output from ESPResSo\n", err)
    return decode_transport_data(out)    
