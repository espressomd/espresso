# Hints for Tutors: Tutorial 11 - Ferrofluids

## General Remarks:
Try to put the focus on the first tutorial.
In 1.5h one should theoretically get through the first two tutorials, but key concepts are mainly introduced in the first one, the second one only introducing a external field.
The third tutorial deals mainly with physical concepts, the only new espresso-specific thing is the setup of a 3D system, i.e. no more DLC.

## Learning Objectives (Physics):
 1. Tutorial:
   - Understand role and importance of order parameter lambda ($\lambda \approx 1 \implies$ no large scale strcture, $\lambda > 2 \implies$ formation of chain- and ring-like structures)
   - How to set up the system correctly (especially randomly assigning magnetic moment vectors):
     - Randomly picking vectors would introduce bias
     - Randomly picking angles introduces bias
     - Correct way shown in solution
   - What cluster distribution one expects
   - If more time available: What changes with different lambdas

 2. Tutorial:
   - $\alpha$-parameter comparable to $\lambda$, w.r.t. external field
   - External force/field introduces ordering of moments 
   - Corrections in the analytical models: Due to torque from dipole-dipole interactions


## Learning Objectives (ESPResSo):
 1. Tutorial:
   - How to setup P3M for dipole interactions
   - (Roughly) how to perform cluster analysis
   - Using DLC (analogous to ELC) for quasi 2D system with 3D moment directions

 2. Tutorial:
   - How to setup a external field

## Points to Mention Throughout the Tutorial:
 1. Tutorial:
   - Importance of correct setup for non-biased results
   - Idea of a quasi-2D-System with 3D Moments
 2. Tutorial:
   - Appeal of simulations as "numerical" experiment and benchmark against theoretical models
