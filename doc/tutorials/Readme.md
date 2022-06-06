# Tutorials for ESPResSo

## Overview

This folder contains tutorials that introduce the use of ESPResSo for different
physical systems.

[comment]: # (Begin of tutorials landing page)

### Introductory tutorials

* **Simulate a simple Lennard-Jones liquid**  
  Modelling of a single-component and a two-component Lennard-Jones liquid.  
  [Guide](lennard_jones/lennard_jones.ipynb)
* **Error analysis**  
  Statistical analysis of simulation results  
  Guide
  [Part 1](error_analysis/error_analysis_part1.ipynb) |
  [Part 2](error_analysis/error_analysis_part2.ipynb)
* **Visualization**  
  Using the online visualizers of ESPResSo.  
  [Guide](visualization/visualization.ipynb)

### Intermediate tutorials

* **Charged systems**  
  Modelling of ion condensation around a charged rod.  
  [Guide](charged_system/charged_system.ipynb)
* **Langevin dynamics**  
  Modelling of Brownian motion and measurement of diffusion coefficients.  
  [Guide](langevin_dynamics/langevin_dynamics.ipynb)
* **Ferrofluid**  
  Modelling of a monolayer ferrofluid system.  
  Guide
  [Part 1](ferrofluid/ferrofluid_part1.ipynb) |
  [Part 2](ferrofluid/ferrofluid_part2.ipynb) |
  [Part 3](ferrofluid/ferrofluid_part3.ipynb)
* **Lattice-Boltzmann**  
  Simulations including hydrodynamic interactions using the Lattice-Boltzmann method.  
  Guide
  [Part 1](lattice_boltzmann/lattice_boltzmann_theory.ipynb) |
  [Part 2](lattice_boltzmann/lattice_boltzmann_poiseuille_flow.ipynb)
* **Polymers**  
  Modelling polymers with hydrodynamic interactions.  
  [Guide](polymers/polymers.ipynb)
* **Raspberry electrophoresis**  
  Extended objects in a Lattice-Boltzmann fluid, raspberry particles.  
  [Guide](raspberry_electrophoresis/raspberry_electrophoresis.ipynb)

### Advanced tutorials

* **Active matter**  
  Modelling of self-propelling particles.  
  [Guide](active_matter/active_matter.ipynb) |
  [Exercises](active_matter/exercises) |
  [Solutions](active_matter/solutions)
* **Electrokinetics**  
  Modelling electrokinetics together with hydrodynamic interactions.  
  [Guide](electrokinetics/electrokinetics.ipynb)
* **Constant pH method**  
  Modelling an acid dissociation curve using the constant pH method.  
  [Guide](constant_pH/constant_pH.ipynb)

[comment]: # (End of tutorials landing page)

## Using the tutorials

To run the tutorials, you need ESPResSo and a Jupyter environment.
For installation instructions, please see the user guide sections
[Quick installation](https://espressomd.github.io/doc/installation.html#quick-installation)
and [Setting up a Jupyter environment](https://espressomd.github.io/doc/installation.html#setting-up-a-jupyter-environment).

Tutorials are available as Jupyter notebooks, i.e. they consist of a ``.ipynb``
file which contains both the source code and the corresponding explanations.
They can be viewed, changed and run interactively. To generate the tutorials
in the build folder, do:

```sh
make tutorials
```

The tutorials contain solutions hidden with the ``exercise2`` extension.
Since this extension is only available for Jupyter Notebook, JupyterLab
users need to convert the tutorials:

```sh
for f in doc/tutorials/*/*.ipynb; do
  ./pypresso doc/tutorials/convert.py exercise2 --to-jupyterlab ${f};
done
```

All tutorials can be viewed with their solutions
[online](https://espressomd.github.io/doc/tutorials.html).

### Running the tutorials interactively

To view the tutorials, first change to the tutorials directory and then run
the `ipypresso` script from the directory into which espresso was compiled.

For Jupyter Notebook and IPython users:

```sh
cd doc/tutorials
../../ipypresso notebook
```

For JupyterLab users:

```sh
cd doc/tutorials
../../ipypresso lab
```

This will launch a web browser in which the notebooks for the tutorials can
be viewed and run. For more details, please see the user guide section on
[running notebooks](https://espressomd.github.io/doc/running.html#interactive-notebooks),
which walks you through the Jupyter interface.

## Video lectures

[comment]: # (Begin of videos landing page)

* [Introduction to ESPResSo](https://www.youtube.com/watch?v=aP4jvpD-D1w)
* [How to run ESPResSo inside Visual Studio Code](https://www.youtube.com/watch?v=dlvF1Zk3AAs)
* [Error Estimation in Time-Correlated Data](https://www.youtube.com/watch?v=I-HCxj9dUIU)
* [Electrostatic Algorithms](https://www.youtube.com/watch?v=YPryFf7MQTg)
* [Introduction to Charged Soft Matter](https://www.youtube.com/watch?v=wrnDg-3j2ik)
* [Managing Simulation Data](https://www.youtube.com/watch?v=64rNmTpoS1c)
* [Introduction to Ferrofluids](https://www.youtube.com/watch?v=wbL3EdVCbkI)
* [Molecular Modelling of Polymers](https://www.youtube.com/watch?v=vSF5-eciwms)
* [Simulating Chemical Reactions in ESPResSo](https://www.youtube.com/watch?v=MUG-PSaMFVM)
* [Introduction to Lattice Boltzmann Method](https://www.youtube.com/watch?v=jfk4feD7rFQ)

[comment]: # (End of videos landing page)
