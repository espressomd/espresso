# Introduction

This sample illustrates how to run a large amount of short ESPResSo simulations
with Dask. Dask is a parallel computing library in Python that enables efficient
handling of large datasets and computation tasks.
Note that this sample is not meant to produce meaningful physics results.
The sample consists of the following parts:

- `espresso_dask.py`: contains helper functions that handle running ESPResSo
  within Dask and communicating data between Dask and ESPResSo
- `lj_pressure.py`: simulation script which obtains the average pressure
  for a Lennard-Jones liquid at a given volume fraction
- `run_pv.py`: Uses Dask to run the simulation script at various volume
  fractions and obtain a pressure vs volume fraction curve.
- `test_dask_espresso.py`: corresponding unit tests, to be run with `pytest`
- `echo.py`: Used to mock an ESPResSo simulation for the unit tests

## How to Use

Note: It is not possible to use ESPResSo with `dask.distributed.LocalCluster`.
Instead, follow the procedure described below:

1. Move to the sample directory
   ```bash
   cd samples/high_throughput_with_dask
   ```
1. Open `run_pv.py` in an editor and adapt the `PYPRESSO` variable
   to the correct path to `pypresso`
1. Set the `PYTHONPATH` environment variable such that it includes
   the directory in which `dask_espresso.py` resides:
   ```bash
   export PYTHONPATH="${PYTHONPATH:+$PYTHONPATH:}$(realpath .)"
   ```
1. Start the Dask scheduler
   ```bash
   dask scheduler &
   ```
1. Note the address of the scheduler (e.g., `tcp://127.0.0.1:8786`)
1. Launch a few workers using the correct scheduler address:
   ```bash
   for i in {1..5}; do dask worker SCHEDULER_ADDRESS & done
   ```
1. Run `python3 run_pv.py SCHEDULER_ADDRESS`, again inserting the scheduler address from above
1. Use `fg` and Ctrl-C to shut down the Dask workers and scheduler,
   or use `pkill "dask"` if you don't have any other Dask scheduler
   running in the background.

Note that Dask can also be used on compute clusters with HTCondor and Slurm.

## Technical Notes

- Since currently only one ESPResSo instance can be used in a Python script,
  ESPResSo is run as a separate process. This is accomplished by the
  `dask_espresso_task` function in `dask_espresso.py`.
- Also, the data transfer between Dask and ESPResSo has to be handled such that
  it is safe for inter-process communication. This is achieved via the `pickle`
  and `base64` Python modules. Encoding and decoding functions can be found in
  `dask_espresso.py`
- The communication happens via the standard input and output of the simulation
  script. Therefore, it is essential not to use simple `print()` calls in the
  simulation script. Instead, use the `logging` module for status messages.
  These will go to the standard error stream.
- To use this sample for your own simulations:
   - Use `dask_espresso.py` as is.
   - Adapt `run_pv.py` to run simulations with the parameters you need.
     The keyword arguments passed to `dask_espresso_task()` will be passed
     as a dictionary to the simulation.
   - Use `data = dask_espresso.get_data_from_stdin()` to get the parameters
     at the beginning of the simulation script.
   - Use `print(dask_espresso.encode_transport_data(result))` at the end
     of your simulation to pass the result to Dask.
   - The simulation parameters and results can be any Python object that
     can be safely pickled and do not require additional context. Basic data
     types (int, float, string, list, dict) as well as numpy arrays work,
     whereas objects that require additional context to be valid do not
     (e.g. file objects and ESPResSo particles).
   - To test your simulation script, including the transfer of parameters
     and results outside Dask, you can also use
     the `dask_espresso.dask_espresso_task.py` function.
