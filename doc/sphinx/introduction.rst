.. include:: defs.rst

Introduction ============

|es| is a simulation package designed to perform Molecular Dynamics (MD) and
Monte Carlo (MC) simulations.  It is meant to be a universal tool for
simulations of a variety of soft matter systems.  |es| features a broad range
of interaction potentials which opens up possibilities for performing
simulations using models with different levels of coarse-graining. It also
includes modern and efficient algorithms for treatment of electrostatics (P3M,
MMM-type algorithms, Maggs algorithm, \ldots), hydrodynamic interactions (DPD,
Lattice-Boltzmann), and magnetic interactions. It is designed to exploit the
capabilities of parallel computational environments. The program is being
continuously extended to keep the pace with current developments both in the
algorithms and software.

The kernel of |es| is written in C with computational efficiency in mind.
Interaction between the user and the simulation engine is provided via a Tcl
scripting interface. This enables setup of arbitrarily complex systems which
users might want to simulate in future, as well as modifying simulation
parameters during runtime.


Guiding principles ==================

|es| is a tool for performing computer simulation and this user guide describes
how to use this tool. However, it should be borne in mind that being able to
operate a tool is not sufficient to obtain physically meaningful results. It is
always the responsibility of the user to understand the principles behind the
model, simulation and analysis methods he is using. |es| will *not* do that for
you!

It is expected that the users of |es| and readers of this user guide have a
thorough understanding of simulation methods and algorithms they are planning
to use. They should have passed a basic course on molecular simulations or read
one of the renown textbooks. It is not necessary to understand everything that
is contained in |es|, but it is inevitable to understand all methods that you
want to use.  Using the program as a black box without proper understanding of
the background will most probably result in wasted user and computer time with
no useful output.

To enable future extensions, the functionality of the program is kept as
general as possible.  It is modularized, so that extensions to some parts of
the program (eg implementing a new potential) can be done by modifying or
adding only few files, leaving most of the code untouched.

To facilitate the understanding and the extensibility, much emphasis is put on
readability of the code.  Hard-coded assembler loops are generally avoided in
hope that the overhead in computer time will be more than compensated for by
saving much of the user time while trying to understand what the code is
supposed to do.

Hand-in-hand with the extensibility and readability of the code comes the
flexibility of the whole program.  On the one hand, it is provided by the
generalized functionality of its parts, avoiding highly specialized functions.
An example can be the implementation of the Generic Lennard-Jones potential
described in section~\ref{sec:GenLennardJones} where the user can change all
available parameters. Where possible, default values are avoided, providing the
user with the possibility of choice.  |es| cannot be aware whether your
particles are representing atoms or billiard balls, so it cannot check if the
chosen parameters make sense and it is the user's responsibility to make sure
they do.
