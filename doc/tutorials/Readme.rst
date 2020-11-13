Tutorials for ESPResSo
======================

Overview
--------

This folder contains tutorials that introduce the use of ESPResSo for different
physical systems. Currently, the following tutorials are available:

* ``lennard_jones``: Modelling of a single-component and a two-component Lennard-Jones liquid
* ``charged_system``: Modelling of ion condensation around a charged rod
* ``lattice_boltzmann``: Simulations including hydrodynamic interactions using the lattice-Boltzmann method
* ``raspberry_electrophoresis``: Extended objects in a lattice-Boltzmann fluid, raspberry particles
* ``active_matter``: Modelling of self-propelling particles
* ``electrokinetics``: Modelling electrokinetics together with hydrodynamic interactions
* ``visualization``: Using the online visualizers of ESPResSo
* ``ferrofluid``: Modelling a monolayer ferrofluid system
* ``constant_pH``: Modelling the titration of a weak acid using the constant pH method

Using the tutorials
-------------------

For using the tutorials, you need ESPResSo running. For installation
instructions, please see: http://espressomd.org/html/doc/installation.html

Tutorials are available as Jupyter notebooks, i.e. they consist of a ``.ipynb``
file which contains both the source code and the corresponding explanations.
They can be viewed, changed and run interactively.

All tutorials can be viewed with their solutions at
http://espressomd.org/wordpress/documentation

Running the tutorials interactively
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

When using Jupyter, you will need extra plugins:

.. code-block:: bash

    pip3 install --user 'jupyter_contrib_nbextensions==0.5.1'
    jupyter contrib nbextension install --user
    jupyter nbextension enable rubberband/main
    jupyter nbextension enable exercise2/main

To view the tutorials, first change to the tutorials directory and then run
the `ipypresso` script from the directory into which espresso was compiled:

.. code-block:: bash

    cd doc/tutorials
    ../../ipypresso notebook

This will launch a web browser in which the notebooks for the tutorials can be
viewed and run. For more details, please see the user guide section on `running
ESPResSo <http://espressomd.org/html/doc/installation.html#running-es>`_, which
walks you through the Jupyter interface.
