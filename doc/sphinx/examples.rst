Sample scripts
==============

Serveral scripts that can serve as usage examples can be found in the directory ``/samples``

.. todo:: currently resides in ``/samples/python`` should be moved?

* ``billard.py`` needs the Python ``pypopengl`` module

* ``bonds-tst.py``
   Test script that manually creates and deletes different bonds between particles (see :ref:`Bonded interactions`). This script performs:
  
   * print defined bonded interactions 
   * print bonds on a particle
   * delete bonds by index or name
   * save/load a bond to/from a variable
 

* ``cellsystem_test.py``
    Test script that changes the skin depth parameter.  This should not be seen as a benchmark, but rather as a rough estimate of the effect of the cellsystem.     
    .. todo:: implement the older [tune_cells] call
    .. todo:: add save/load optimal cell parameters from tune_skin()
    

* ``coulomb_debye_hueckel.py``,  ``debye_hueckel.py``
    Charged beads with a WCA interaction are simulated using the screened Debye-Hückel potential. See :ref:`Debye-Hückel potential`


* ``ekboundaries.py``

* ``electrophoresis.py``

* ``h5md.py``

* ``lbf.py``

* ``lj-demo.py``
    Lennard-Jones liquid used for demonstration purposes to showcase |es|.
    Sliders from a MIDI controller can change system variables such as
    temperature and volume. Some thermodynamic observables are analyzed and
    plotted live.

* ``lj_liquid_distribution.py``
    Uses ``analysis.distribution`` (See :ref:`Particle distribution`) to analyze a simple Lennard-Jones liquid.

* ``lj_liquid.py``
    Simple Lennard-Jones particle liquid. Shows the basic features of how to:

    * set up system parameters, particles and interactions.
    * warm up and integrate. 
    * write parameters, configurations and observables to files

* ``lj_liquid_structurefactor.py``
    Uses ``analysis.structure_factor`` (See :ref:`Structure factor`) to analyze a simple Lennard-Jones liquid.


* ``load_bonds.py``,  ``store_bonds.py``
    Uses the Python ``pickle`` module to store and load bond information.

* ``load_checkpoint.py``,  ``save_checkpoint.py``
   Basing usage of the checkpointing feature. Shows how to write/load the state of:   
   * custom user variables
   * non bonded interactions
   * particles
   * P3M paremeters
   * thermostat

* ``load_properties.py``,  ``store_properties.py``
    Uses the Python ``pickle`` module to store and load system information.

* ``MDAnalysisIntegration.py``.
    Shows how to expose configuration to ``MDAnalysis`` at run time. The functions of ``MDAnalysis`` can be used to perform some analysis or 
    convert the frame to other formats (CHARMM, GROMACS, ...)

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
    Uses python array slicing to set and extract various particle properties.

* ``visualization_bonded.py``

* ``visualization_openGL.py``

* ``visualization.py``


