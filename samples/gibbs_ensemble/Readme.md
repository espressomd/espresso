# Gibbs ensemble simulations using ESPResSo

This sample code shows how to implement the Gibbs ensemble into ESPResSo using
the `multiprocessing` module to create two separate systems and communicate to
them through the main script.

The program consists of 3 files: In the `run_sim.py` script, the client boxes
are started and given instructions. The clients use the `Gibbs_Client` subclass
in `client_system.py` which inherits the necessary methods for communication and
for the trial moves from the `Client` class in `gibbs_ensemble.py`. This subclass
allows for energy correction terms as is seen in this sample. Because ESPResSo
only allows one instance of a system in one process, the system is created using
the `set_up_system()` function and initialized once the `run()` method in the
`Client` class is called.

As a side note: Do not trust the critical point too much as it is quite sensitive to
the init values chosen for the fit.

All equations are taken from Frenkel, Smit: *Understanding Molecular Simulation* 2002,
doi:[10.1016/B978-0-12-267351-1.X5000-7](https://doi.org/10.1016/B978-0-12-267351-1.X5000-7).
