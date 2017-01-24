Sample scripts
==============

In the directory /samples you find several scripts that can serve as
samples how to use .

lj\_liquid.tcl
    Simple Lennard-Jones particle liquid. Shows the basic features of :
    How to set up system parameters, particles and interactions. How to
    warm up and integrate. How to write parameters, configurations and
    observables to files. How to handle the connection to VMD.

pe\_solution.tcl
    Polyelectrolyte solution under poor solvent condition. Test case for
    comparison with data produced by polysim9 from M.Deserno. Note that
    the equilibration of this system takes roughly :math:`15000 \tau`.

pe\_analyze.tcl
    Example for doing the analysis after the actual simulation run
    (offline analysis). Calculates the integrated ion distribution
    :math:`P(r)` for several different time slaps, compares them and
    presents the final result using gnuplot to generate some ps-files.

harmonic\_oscillator.tcl
    A chain of harmonic oscillators. This is a :math:`T=0` simulation to
    test the energy conservation.

espresso\_logo.tcl
    The -logo, the exploding espresso cup, has been created with this
    script. It is a regular simulation of a polyelectrolyte solution. It
    makes use of some nice features of the part command (see section ,
    namely the capability to fix a particle in space and to apply an
    external force.
