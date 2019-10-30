Tutorials for ESPResSo
======================

Overview
--------


This folder contains tutorials that introduce the use of ESPResSo for different
physical systems. Currently, the following tutorials are available:

* :file:`01-lennard_jones`: Modelling of a single-component and a two-component Lennard-Jones liquid
* :file:`02-charged_system`: Modelling of charged systems such as ionic crystals
* :file:`04-lattice_boltzmann`: Simulations including hydrodynamic interactions using the lattice-Boltzmann method
* :file:`05-raspberry_electrophoresis`: Extended objects in a lattice-Boltzmann fluid, raspberry particles
* :file:`06-active_matter`: Modelling of self-propelling particles
* :file:`07-electrokinetics`: Modelling electrokinetics together with hydrodynamic interactions
* :file:`08-visualization`: Using the online visualizers of ESPResSo
* :file:`10-reaction_ensemble`: Modelling chemical reactions by means of the reaction ensemble
* :file:`11-ferrofluid`: Modelling a monolayer ferrofluid system
* :file:`12-constant_pH`: Modelling the titration of a weak acid using the constant pH method

Using the tutorials
-------------------
For using the tutorials, you need ESPResSo running. For installation
instructions, please see: http://espressomd.org/html/doc/installation.html

Tutorials 1, 2, 4, 5, 8 and 12 are available as IPython notebooks, i.e.
they consist of a `.ipynb` file which contains both the source code
and the corresponding explanations.
They can be viewed, changed and run interactively.


The remaining tutorials consist of a `.pdf`-file containing the explanations and separate `.py`-scripts containing the simulation scripts.

All tutorials can be viewed and the corresponding simulation scripts downloaded
from http://espressomd.org/wordpress/documentation

Using the Jupyter tutorials interactively
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To view the tutorials, either IPython or Jupyter needs to be installed.
To check whether one of them is installed, run:

.. code-block:: bash

    ipython --version
    jupyter --version

If none is found, on Ubuntu and related platforms, IPython can be installed with:

.. code-block:: bash

    sudo apt install ipython-notebook

while Jupyter (the successor of IPython) can be installed with:

.. code-block:: bash

    pip3 install --upgrade jupyter ipython nbconvert

To view the tutorials, first change to the tutorials directory and then run the `ipypresso` script from the directory into which espresso was compiled:

.. code-block:: bash

    cd doc/tutorials
    /path_to_espresso_build/ipypresso notebook

This will launch a web browser in which the notebooks for the tutorials can be viewed and run.
For more details, please see the user guide section on `running ESPResSo
<http://espressomd.org/html/doc/installation.html#running-es>`_, which walks
you through the Jupyter interface.

