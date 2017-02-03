Sample scripts
==============

Serveral scripts that can serve as usage examples can be found in the directory ``/samples``

.. todo:: currently resides in ``/samples/python`` should be moved?

* ``billard.py`` needs the Python ``pypopengl`` module

* ``bonds-tst.py``
   Test script that manually creates and deletes different bonds between particles (see :ref:`Bonded interactions`). This shows how to:
  
   * print defined bonded interactions 
   * print bonds on a particle
   * delete a specific (or all) bonds

.. todo:: print to screen the state of bonds, with what the code is doing. 
 

* ``cellsystem_test.py``

.. todo:: remove this file (needed?), or merge with a more complete simulation example. 

* ``coulomb_debye_hueckel.py``,  ``debye_hueckel.py``
    Charged beads with a WCA interaction are simulated using the screened Debye-Huckle potential. See :ref:`Debye-Hückel potential`



* ``ekboundaries.py``

* ``electrophoresis.py``

* ``h5md.py``

* ``lbf.py``

* ``lj-demo.py``

* ``lj_liquid_distribution.py``

* ``lj_liquid.py``
    Simple Lennard-Jones particle liquid. Shows the basic features of how to:

    * set up system parameters, particles and interactions.
    * warm up and integrate. 
    * write parameters, configurations and observables to files
    * handle the connection to VMD.

* ``lj_liquid_structurefactor.py``

* ``load_bonds.py``,  ``store_bonds.py``

* ``load_checkpoint.py``,  ``save_checkpoint.py``
   Basing usage of the checkpointing feature. Shows how to write/load the state of:   
   * custom user variables
   * non bonded interactions
   * particles
   * P3M paremeters
   * thermostat

* ``load_properties.py``,  ``store_properties.py``

* ``minimal-charged-particles.py``
   Simple Lennard-Jones particle liquid where the particles are assigned charges. The P3M method is used to calculate electrostatic interactions. 

* ``minimal-diamond.py``

* ``minimal-polymer.py``
   Sets up a single dilute bead-spring polymer. Shows the basic usage of ``create_polymer``.

* ``minimal_random_number_generator.py``

* ``observables_correlators.py``

* ``p3m.py``
   Simple Lennard-Jones particle liquid where the particles are assigned charges. The P3M method is used to calculate electrostatic interactions. 

* ``slice_input.py``

* ``visualization_bonded.py``

* ``visualization_openGL.py``

* ``visualization.py``


