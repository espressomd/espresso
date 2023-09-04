# Introduction

This sample illustrates how to run a large amount of short Espresso simulations with Dask. Dask is a parallel computing library in Python that enables efficient handling of large datasets and computation tasks. 
Note that this sample is not meant to produce meaninful physics results.
The sample consists of the following parts:

- `espersso_dask.py`: contains helper functions that handle running Espesso within Dask and communicating data between Dask and ESPResSo
- `lj_pressure.py`: simulation script which obtains the average pressure for a Lennard-Jones liquid at a given volume fraction
- `run_pv.py`: Uses Dask to run the simulation script at various volume fractions and obtain a pressure vs volume fraction curve.
- `test_espresso_dask.py`: corresponding unit tests, to be run with pytest-3
- `echo.py`: Used to mock an ESPResSo simulation for the unit tests

## How to Use

Note: It is not possible to use ESPResSo with `dask.distributed.LocalCluster`. Instead, follow the procedure described below:

1. Set the `PYTHONPATH` environment variable such that it includes the directory in which `dask_espersso.py` resides: 
   ```bash
   export PYTHONPATH=/directory/containing/dask_espresso.py:$PYTHONPATH
   ```
2. Run the dask scheduler 
   ```bash
   dask scheduler &
   ```
3. Note the address of the scheduler (e.g., `tcp://127.0.0.1:8786`)
4. Launch a few workers inserting the correct scheduler address:
   ```bash
   for i in {1..5}; do dask worker SCHEDULER_ADDRESS; done
   ```
5. Open `run_pv.py` in an editor and adapt the `PYPRESSO` variable to the correct path to `pypresso`
6. Run `python3 run_pv.py SCHEDULER_ADDRESS` again inserting the scheduler address from above
7. Use `fg` an CTRL-C co shut down the Dask workers and scheduler.

Note that Dask can also be used on compute clusters with HTCondor and SLURM.


## Technical Notes

• Since currently only one ESPResSo instance can be used in a Python script, ESPResSo is run as a separate process. This is accomplished by the `dask_espresso_task` function in `dask_espersso.py`.
• Also, the data transfer between Dask and ESPResSo has to be handled such that it is safe for inter-process communication. This is achieved via the `pickle` and `base64` Python modules. Encoding and decoding functions can be foud in `dask_espresso.py`
• The communication happens via the standard input and output of the simulation script. Therefore, it is essential not to use simple `print()` calls in the simulation script. Instead, use the `logging` module for status messages. These will go to the standard error stream.
• To use this sample for your own simulations:
i
    • Use `dask_espresso.py` as is.
    • Adapt `run_pv.py` to run simulations with the parameters you need. The keyword arguments passed to `dask_espresso_task()` will be passed as a dictionary to the simulation.
    • Use `data = dask_espresso.get_data_from_stdin()` to get the parameters at the beginning of the simulation script.
    • Use `print(dask_espresso.encode_transport_data(result))` at the end of your simulation to pass the result to Dask.
    • The simulation parameter and result can be any Python object that can be safely pickled and do not require additional context. Basic data types (int, float, string, list, dict) as well as numpy arrays work, whereas objects that require additional context to be valid do not (e.g. file objects and ESPResSo particles).
    * To test your simulation script including the transfer of parameters and results outside Dask, you can also use the `dask_espresso.dask_espresso_task.py` function.
