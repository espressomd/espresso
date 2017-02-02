Auxilliary commands
===================

Finding particles and bonds
---------------------------

``countBonds``
~~~~~~~~~~~~~~

countBonds

Returns a Tcl-list of the complete topology described by , which must
have the same format as the output of the command (see section ).

The output list contains only the particle id and the corresponding
bonding information, thus it looks like

106 0 107 107 0 106 0 108 108 0 107 0 109 ... 210 0 209 0 211 211 0 210
212 213 ...

for a single chain of 106 monomers between particle 106 and 211, with
additional loose particles 212, 213, ... ( counter-ions). Note, that the
command stores any bonds only with the particle of lower particle
number, which is why would only return , therefore not revealing the
bond between particle 109 and (the preceding) particle 108, while would
return all bonds particle 109 participates in.

``findPropPos``
~~~~~~~~~~~~~~~

findPropPos

Returns the index of within , which is expected to have the same format
as . If is not found, is returned.

This function is useful to access certain properties of particles
without hard-wiring their index-position, which might change in future
releases of |es|.


This returns the particle type id of particle without fixing where
exactly that information has to be in the output of .

``findBondPos``
~~~~~~~~~~~~~~~

findBondPos

Returns the index of the bonds within , which is expected to have the
same format as ; hence its output is the same as . If the particle does
not have any bonds, is returned.

``timeStamp``
~~~~~~~~~~~~~

timeStamp

Modifies the filename contained within to be preceded by a and having
before the ;

timeStamp ./scripts/config.gz DH863 001 gz

returns . If is :math:`-1`, the current date is used in the format .
This would results in on October 22nd, 2002.

Additional Tcl math-functions
-----------------------------

The following procedures are found in scripts/ABHmath.tcl.

-  CONSTANTS

   -  PI

      returns :math:`\pi` with 16 digits precision.

   -  KBOLTZ

      Returns Boltzmann constant in Joule/Kelvin

   -  ECHARGE

      Returns elementary charge in Coulomb

   -  NAVOGADRO

      Returns Avogadro number

   -  SPEEDOFLIGHT

      Returns speed of light in meter/second

   -  EPSILON0

      Returns dielectric constant of vaccum in Coulomb2/(Joule meter)

   -  ATOMICMASS

      Returns the atomic mass unit u in kilogramms

-  MATHEMATICAL FUNCTIONS

   -  sqr <arg>

      returns the square of .

   -  min <arg1> <arg2>

      returns the minimum of and .

   -  max <arg1> <arg2>

      returns the maximum of and .

   -  sign <arg>

      returns the signum-function of , namely +1 for :math:`>0`, -1 for
      :math:`<0`, and =0 otherwise.

-  RANDOM FUNCTIONS

   -  gauss\_random

      returns random numbers which have a Gaussian distribution

   -  dist\_random <dist> [max]

      returns random numbers in the interval :math:`[0,1]` which have a
      distribution according to the distribution function p(x) which has
      to be given as a tcl list containing equally spaced values of
      p(x). If p(x) contains values larger than 1 (default value of max)
      the maximum or any number larger than that has to be given . This
      routine basically takes the function p(x) and places it into a
      rectangular area ([0,1],[0,max]). Then it uses to random numbers
      to specify a point in this area and checks wether it resides in
      the area under p(x). Attention: Since this is written in tcl it is
      probably not the fastest way to do this!

   -  vec\_random [len]

      returns a random vector of length (uniform distribution on a
      sphere) This is done by chosing 3 uniformly distributed random
      numbers :math:`[-1,1]` If the length of the resulting vector is
      :math:`<= 1.0` the vector is taken and normalized to the desired
      length, otherwise the procedure is repeated until succes. On
      average the procedure needs 5.739 random numbers per vector. (This
      is probably not the most efficient way, but it works!) Ask your
      favorit mathematician for a proof!

   -  phivec\_random <v> <phi> [len]

      return a random vector at angle with and length

-  PARTICLE OPERATIONS

   Operations involving particle positions. The parameters can either
   denote the particle identity (then the particle position is extracted
   with the The part command command) or the particle position directly
   When the optional parameter for minimum image conventions is omited
   the functions use the the command.

   -  bond\_vec <p1> <p2>

      Calculate bond vector pointing from particles to return = (.pos -
      .pos)

   -  bond\_vec\_min <p1> <p2> [box]

      Calculate bond vector pointing from particles to return =
      MinimumImage(.pos - .pos)

   -  bond\_length <p1> <p2>

      Calculate bond length between particles and

   -  bond\_length\_min <p1> <p2> [box]

      Calculate minimum image bond length between particles and

   -  bond\_angle <p1> <p2> <p3> [type]

      Calculate bond angle between particles , and . If is “r” the
      return value is in radiant. If it is “d” the return value is in
      degree. The default for is “r”.

   -  bond\_dihedral <p1> <p2> <p3> <p4> [type]

      Calculate bond dihedral between particles , , and If is “r” the
      return value is in radiant. If it is “d” the return value is in
      degree The default for is “r”.

   -  part\_at\_dist <p> <dist>

      return position of a new particle at distance from with random
      orientation

   -  part\_at\_angle <p1> <p2> <phi> [len]

      return position of a new particle at distance (default=1.0) from
      which builds a bond angle for (, , p-new)

   -  part\_at\_dihedral <p1> <p2> <p3> <theta> [phi] [len]

      return position of a new particle at distance (default=1.0) from
      which builds a bond angle (default=random) for (, , p-new) and a
      dihedral angle for (, , , p-new)

-  INTERACTION RELATED

   Help functions related to interactions implemented in .

   -  calc\_lj\_shift <lj\_sigma> <lj\_cutoff>

      returns the value needed to shift the Lennard Jones potential to
      zero at the cutoff.

-  VECTOR OPERATIONS

   A vector is a tcl list of numbers with an arbitrary length Some
   functions are provided only for three dimensional vectors.
   corresponding functions contain 3d at the end of the name.

   -  veclen <v>

      return the length of a vector

   -  veclensqr <v>

      return the length of a vector squared

   -  vecadd <a> <b>

      add vector to vector : return = (+)

   -  vecsub <a> <b>

      subtract vector from vector : return = (-)

   -  vecscale <s> <v>

      scale vector with factor : return = (\*)

   -  vecdot\_product <a> <b>

      calculate dot product of vectors and : return = (.)

   -  veccross\_product3d <a> <b>

      calculate the cross product of vectors and : return = ( x )

   -  vecnorm <v> [len]

      normalize a vector to length (default 1.0)

   -  unitvec <p1> <p2>

      return unit vector pointing from position to position

   -  orthovec3d <v> [len]

      return orthogonal vector to with length (default 1.0) This vector
      does not have a random orientation in the plane perpendicular to

   -  create\_dihedral\_vec <v1> <v2> <theta> [phi] [len]

      create last vector of a dihedral (, , res) with dihedral angle and
      bond angle (, res) and length (default 1.0). If is ommited or set
      to rnd then is assigned a random value between 0 and 2 Pi.

-  TCL LIST OPERATIONS

   -  average <list>

      Returns the avarage of the provided

   -  list\_add\_value <list> <val>

      Add to each element of

   -  flatten <list>

      flattens a nested

   -  list\_contains <list> <val>

      Checks wether contains . returns the number of occurences of in .

-  REGRESSION

   -  LinRegression <l>

      where is a listof pairs of points ``{ {$x1 $y1} {$x2 $y2} ...} ``.
      ``LinRegression`` returns the least-square linear fit :math:`ax+b`
      and the standard errors :math:`\sigma_a` and :math:`\sigma_b`.

   -  LinRegressionWithSigma <l>

      where is a list of lists of points in the form
      ``{{$x1 $y1 $s1} {$x2 $y2 $s2} ...}`` where ``s`` is the standard
      deviation of ``y``. ``LinRegressionWithSigma`` returns the
      least-square linear fit :math:`ax+b`, the standard errors
      :math:`\sigma_a` and :math:`\sigma_b`, covariance
      :math:`\mathrm{cov}(a,b)` and :math:`\chi`.

``t_random``
~~~~~~~~~~~~

-  Without further arguments (tcl only),

   t\_random

   returns a random double between 0 and 1 using the standard C++
   Mersenne twister random number generator. For drawing random numbers
   in python please use the numpy random number generator.

-  t\_random int <n> (tcl only)

   returns a random integer between 0 and n-1. For python please use the
   numpy random number generator.

-  Note that the best practice in initilizing the random number
   generator is to randomly set its internal state. This is *attempted*
   by seeding the random number generator however the state space of the
   random number generator is typically much bigger than than a single
   random number. Therefore C++ cannot do miracles during the seeding
   with only one random number and cannot come up with more randomness.
   In the python interface Espresso provides the function

   system = espressomd.System() system.set\_random\_state\_PRNG()

   which sets the state of the random number generator with real random
   numbers.

-  Set the state of the random number generator by providing a string
   with a sufficient amount of space separated integers

   t\_random stat “<state(1)> ... <state(n\_nodes\*625)>”

   system = espressomd.System()
   system.random\_number\_generator\_state=[state(1),...
   ,state(n\_nodes\*625)]

   If you want to get the random number generator state, use

   t\_random stat

   system = espressomd.System()
   print(system.random\_number\_generator\_state)

-  Alternatively to setting the full state of the random number
   generator one can also only set a seed and hope for sane
   initilization of the random number generator state. The following
   command sets the seeds to the new values respectively,
   re-initializing the random number generators on each node

   t\_random seed <seed(0)> ... <seed(n\_nodes-1)>

   system = espressomd.System() system.seed=[seed(0), ...
   seed(n\_nodes-1)]

   Note that Espresso automatically comes up with a seed of the random
   number generator, however due to that your simulation will always
   start with the same random sequence on any node *unless you seed your
   random number generator* at the beginning of the simulation.

-  To obtain the seeds of the random number generator which is used in
   the C++ core of Espresso use the command

   t\_random seed

   system = espressomd.System() system.seed

   which returns a list with the seeds of the random number generators
   on each of the ’n\_nodes’ nodes if they were set by the user.

Checking for features of 
-------------------------

In an -Tcl-script, you can get information whether or not one or some of
the features are compiled into the current program with help of the
following Tcl-commands:

-  code\_info

   provides information on the version, compilation status and the debug
   status of the used code. It is highly recommended to store this
   information with your simulation data in order to maintain the
   reproducibility of your results. Exemplaric output:

   ESPRESSO: v1.5.Beta (Neelix), Last Change: 23.01.2004 Compilation
   status PARTIAL\_PERIODIC ELECTROSTATICS EXTERNAL\_FORCES CONSTRAINTS
   TABULATED LENNARD\_JONES BOND\_ANGLE\_COSINE Debug status MPI\_CORE
   FORCE\_CORE

-  has\_feature <feature> ...

   tests, if is compiled into the kernel. A list of possible features
   and their names can be found here.

-  require\_feature <feature> ...

   tests, if is feature is compiled into the kernel, will exit the
   script if it isn’t and return the error code 42. A list of possible
   features and their names can be found here.
