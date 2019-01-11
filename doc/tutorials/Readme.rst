Tutorials for ESPResSo
======================

Overview
--------


This folder contains tutorials that introduce the use of ESPResSo for different
physical systems. Currently, the following tutorials are available:

* :file:`01-lennard_jones`: Modelling of a single-component and a two-component Lennard-Jones liquid
* :file:`02-charged_system`: Modelling of charged systems such as ionic crystals
* :file:`04-lattice_boltzmann`: Simulations including hydrodynamic interactions using the Lattice-Boltzmann method
* :file:`05-raspberry_electrophoresis`: Extended objects in a Lattice-Boltzmann fluid, raspberry particles
* :file:`06-active_matter`: Modelling of self-propelling particles
* :file:`07-electrokinetics`: Modelling electrokinetics together with hydrodynamic interactions
* :file:`08-visualization`: Using the online visualizers of ESPResSo
* :file:`09-swimmer_reactions`: Further modelling of self-propelling particles
* :file:`10-reaction_ensemble`: Modelling chemical reactions by means of the reaction ensemble

Using the tutorials
-------------------
For using the tutorials, you need ESPResSo running. For installation
instructions, please see: http://espressomd.org/html/doc/installation.html

Tutorials 1, 2, 4, 5, and 8 are available as IPython notebooks. I.e., they consist of a `.ipynb` file which contains both, the source code and the corresponding explanations.
They can be viewed, changed and run interactively.


The remaining tutorials consist of a `.pdf`-file containing the explanations and separate `.py`-scripts containing the simulation scripts.

All tutorials can be viewed and the corresponding simulation scripts downloaded
from http://espressomd.org/wordpress/documentation

Using the IPython tutorials interactively
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To view the tutorials, IPython needs to be installed.
To check, whether it is installed, run:

.. code-block:: bash

    ipython --version

If it is not found, on Ubuntu and related platforms, it can be installed via:

.. code-block:: bash

    sudo apt install ipython-notebook

To view the tutorials, first change to the tutorials directory and then run the `ipypresso` script from the directory into which espresso was compiled:

.. code-block:: bash

    cd doc/tutorials
    /path_to_espresso_build/ipypresso notebook

This will launch a web browser in which the notebooks for the tutorials can be viewed and run.
For more details, please see: http://jupyter.readthedocs.io/en/latest/running.html
Note that `Jupyter` is the successor of IPython.



