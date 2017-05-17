External package: mbtools
=========================

mbtools [1]_ is a set of tcl packages for setting up, analyzing and
running simulations of lipid membrane systems.

mbtools comes with a basic set of tabulated forces and potentials for
lipid interactions and some example scripts to help explain the syntax
of the commands. If you make use of mbtools or of these potentials
please acknowledge us with a citation to:

\* Cooke, I. R., Kremer, K. and Deserno, M. (2005): Tunable, generic
model for fluid bilayer membranes, Phys. Rev. E. 72 - 011506

Introduction
------------

mbtools is located in the folder .

One of the main features of mbtools is the ability to easily create
initial lipid configurations with interesting geometries. These include
flat membranes, cylinders, spheres, toroids, and randomly distributed
gases. Each of these shapes is referred to as a geometry and any number
of geometries can be combined in a single simulation. Once the geometry
has been chosen the user specifies the molecules which should be placed
in this geometry. For example one could choose sphere as a geometry and
then define two different lipids and/or a protein to be placed on the
sphere. Within reason (e.g. size restrictions) it should be possible to
use any mixture of known molecule types on any geometry. The molecule
types available at present include proteins, lipids of any length, and
spherical colloids.

mbtools includes several miscellaneous utility procedures for performing
tasks such as warmup, setting tabulated interactions, designating
molecules to be trapped and a variety of topology related sorting or
data analysis functions.

The analysis part of the mbtools package is designed to wrap together
all the analysis for a simulation into a single simple interface. At the
beginning of the simulation the user specifies which analyses should be
performed by appending its name and arguments to a variable, . After the
analysis is setup one can then simply call to perform all the specified
proceedures. Analysis will store a data value each time is called. Then
when a call to is made the average of all stored values is printed to a
file and the store of values is reset to nil.

Installing and getting started
------------------------------

Since mbtools is provided as part of the espresso molecular dynamics
simulation package you will need to download and install Espresso before
you can use it. Espresso can be downloaded free from .

Once you have installed espresso you can find mbtools by looking inside
the subdirectory. Inside the directory you will see a directory for each
of the mbtools subpackages as well as an directory. All of the examples
scripts should work out of the box except those involving colloids which
require you to install icover.sh (see documentation for hollowsphere
molecule type). To run the simplebilayer example cd to the examples
directory and then type:

:math:`ESPRESSO_SOURCE/`\ PLATFORM/Espresso scripts/main.tcl
simplebilayer.tcl

The first part of this command is simply the full path to the
appropriate espresso executable on your machine when running on multiple
processors). Obviously you will need to have the and environment
variables set for it to work. After this executable the relative paths
to two tcl scripts are given. The first of these is given as an argument
to espresso and is therefore interpreted first by the espresso tcl
interpreter. The second tcl script is in turn passed as an argument to .

Why separate the tcl commands into two files ?

This is really a matter of preference but if we keep all of the key
commands and complex coding in a single file and delegate simple
parameter setting to a separate file it tends to be much easier to
manage large numbers of jobs with regularly changing requirements.
Regardless of your personal preferences, the important point to note is
that all of the important commands are contained in and you should
probably start there to get an understanding for how mbtools works.

Running the simplebilayer example should produce a directory called
which contains the output from your simulation. To view the results cd
to the simplebilayer directory and look at the contents. The directory
contains many files including:

-  The configurations generated during warmup :

-  pdb files corresponding to warmup configurations :

-  The configurations generated during the main run :

-  pdb files corresponding to main run configs :

-  The most recently generated checkpoint file :

-  A directory containing the second most recent checkpoint file:

-  A file containing the topology of the system :

-  The original parameter file with which you ran the simulation :

-  A original parameter file with which you ran the simulation :

-  Files containing analysis output for example :

-  Force and energy tables :

-  VMD script for visualising the warmup :

-  VMD script for visualising the main trajectory :

To visualise your results using the vmd scripts you need to make sure
that you have vmd installed properly and that you have the special vmd
procedures used by the espresso team (i.e. support for the loadseries
command). Then you can visualise by typing:

vmd -e vmd\_animation.script

The script
----------

The file provided in the directory is a relatively complete script
written using mbtools. It is designed to run all of the examples
provided but no more. No doubt you will need to extend it for your own
purposes.

Variables used by 
~~~~~~~~~~~~~~~~~~

expects the user to set various parameters in a file (e.g. ). Some of
these parameters have defaults and generally don’t need to be worried
about except for specific cases. In the following list variables that
have no default and therefore must be set in the parameter file are
noted with an asterisk.

-  [] The type of thermostat to be used. Set to for a dpd thermostat.
   Any other value gives a langevin

-  Required if you set the thermo to

-  Required if you set the thermo to

-  [] The temperature at which the warmup is performed. The default
   behaviour is to use the system temperature

-  [100] Number of integrate steps per warmup cycle

-  [20] Number of calls to integrate over which the warmup occurs

-  [0] Warmup steps to be used for the warmup that occurs after
   particles are freed of any temporary constraints.

-  [0] Warmup times to be used for the warmup that occurs after
   particles are freed of any temporary constraints.

-  [] Whether to use the constant pressure barostat

-  The pressure you want to simulate at. Required if npt is set to

-  box mass. Required if npt is set to “on”

-  Required if npt is . Usually set to 1 as for langevin gamma

-  Required if npt is . Box friction

-  [] vmd mode

-  [8] The number of meshpoints per side for dividing the bilayer plane
   into a grid

-  [1000.0] Distance of the end tail bead from the bilayer midplane
   beyond which a lipid is deemed to have strayed from the membrane
   bulk.

-  The temperature of the simulation during the main run

-  Directory for output

-  Directory where forcetables are kept

-  a name for the simulation

-  the name of the file where the topology is written. Usually

-  A list of forcetable names to be used

-  Box dimensions

-  A complete list of the bonded interactions required

-  A complete list of the non-bonded interactions required

-  A list of system specifications (see documentation for the command in
   [mbtools::systemg])

-  A list of molecule types (see documentation in [mbtools::systemg])

-  timestep to be used during warmup integration

-  timestep for the main integration run

-  skin used for verlet nesting list criterion

-  langevin friction term

-  number of times to do main integration

-  number of steps in each main integration

-  How often to calculate the analysis

-  How often to print out configurations

-  a list of additional lines of commands to be written to the file

Analysis
--------

The analysis package is designed to help organise the many possible
analysis routines that can be performed during a simulation. This
documentation describes the basic user interface commands and then all
of the possible analysis routines. Instructions on how to add a new
analysis routine are given at the end of this section.

Basic commands
~~~~~~~~~~~~~~

The following commands comprise the user interface to the analysis
package.

At the start of a simulation all of the analysis that is to be performed
is specified using the command. Each time you want the analysis
performed a call to should be made. One can then call to write results
to file.

| ::mbtools::analysis::setup\_analysis : -outputdir.arg -suffix.arg
| -iotype.arg -g.arg -str.arg

-  [] A tcl list where each element of the list is a string specifying
   the name and complete argument list for a particular analysis to be
   carried out.

-  [] The directory where analysis output files will be created

-  [] Suffix that will be appended to standard file names for analysis
   output

-  [] The iotype that will be used when opening files for analysis. For
   an explanation of the different iotypes see the documentation for the
   standard tcl command open

-  [8] Number of grid points per side with which to divide the bilayer
   for height profile analyses

-  [4.0] Distance of a tail bead from bilayer midplane beyond which a
   lipid is deemed to be a stray lipid.

Sets up the analysis package for a simulation run or analysis run that
is about to be performed. This routine needs to be called before any
analysis can be performed.

::mbtools::analysis::do\_analysis :

Calls all of the routines corresponding to commands setup in . should be
called only after has already been called.

::mbtools::analysis::reset\_averages :

Calls all of the routines corresponding to commands setup in . These
routines vary from command to command but they typically reset the
storage and counter variables used for analysis results. is typically
only called internally by

::mbtools::analysis::print\_averages :

Calls all of the routines corresponding to commands setup in . These
routines typically print results to a file buffer. After printing the
routine is called internally. should be called only after has already
been called.

Available analysis routines
~~~~~~~~~~~~~~~~~~~~~~~~~~~

boxl : -verbose : output \|\| time\_vs\_boxl

Simply obtains the box dimensions from ESPResSo.

| clusters : -alipid.arg -verbose : output \|\| time\_vs\_clust,
| sizehisto.[format

-  alipid [1.29] Value for the area per lipid to be used in a making a
   rough calculation of the area of clusters

Calls the espresso command which groups molecules in the system into
aggregates. Output to is: maximum cluster size, minimum cluster size,
average size of clusters including those of size 2 or greater, standard
deviation of clusters including those of size 2 or greater, number of
clusters of size 2 or greater, total average cluster size, total cluster
size standard deviation, total number of clusters, length of the
interface between clusters, standard deviation of the interface length,
number of clusters for which length was calculate.

Additionally, at each call of the complete size histogram is printed to
a file with the formatted name ].

| density\_profile : -nbins.arg -hrange.arg -beadtypes.arg
| -colloidmoltypes.arg -r.arg -nogrid
| -verbose : output \|\| av\_zprof

-  [100] Number of slices into which the height range is divided for the
   purpose of calculating densities

-  [6] The maximum vertical distance from the bilayer midplane for which
   to calculate densities. Note that the complete vertical range is
   therefore 2\*varhrange

-  [0] A tcl list of the bead types for which to calculate a density
   profile

-  [] A tcl list of molecule types identifying the molecules which are
   colloids in the system. The default value is a null list

-  [0] A tcl list of sphere radii corresponding to the radii for each
   colloid type in the system. If this is non-zero the density profile
   will be calculated in spherical shells about the colloids in the
   system identified via colloidmoltypes or if colloidmoltypes is not
   set then the system center of mass is assumed for the colloid/vesicle
   center

-  If this is set a grid mesh will not be used to refine the density
   profile calculation by taking into account vertical differences
   between mesh points

Calculates the number density of each of the beadtypes given in
beadtypes as a function of the vertical distance from the bilayer
midplane. Lipids are also sorted according to their orientation and
assigned to upper or lower leaflets accordingly. Thus for a system with
3 beadtypes we would obtain 6 columns of output corresponding to 0
(lower) 1 (lower) 2 (lower) 2 (upper) 1 (upper) 0 (upper) where the
number refers to the bead type and upper or lower refers to the bilayer
leaflet.

energy : -verbose : output \|\| time\_vs\_energy

Obtains the internal energies of the system from the command of
ESPResSo.

flipflop : -verbose : output \|\| time\_vs\_flip

Makes a call to the command of ESPResSo and compares this with a
reference set of lipid orients obtained at the start of the simulation
with . Based on this comparison the number of lipids which have flipped
from their original positions is calculated

fluctuations : -verbose : output \|\| powav.dat

Routine for calculating the power spectrum of height and thickness
fluctuations for a flat bilayer sheet. Uses the routine in ESPResSo to
calculate the height and thickness functions and perform the fft. See
the documentation in the file for detail on what is calculated and how
to obtain a stiffness value from the resulting output. Note that this
routine causes a crash if it detects a large hole in the bilayer.

| localheights : -range.arg -nbins.arg -rcatch.arg -verbose :
| output \|\| av\_localh

-  [1.0] Range of local height deviations over which to bin

-  [100] Number of slices to divide up the height range into for the
   purposes of creating a profile

-  [1.9] The distance about a single lipid to use a starting value for
   finding the 6 closest neighbours

For each lipid we calculate its 6 nearest neighbours and then calculate
the height difference between the central lipid and these neighbours.
Taking these 6 values for each lipid we then create a histogram of
number densities as a function of the height difference.

localorients : -range.arg -nbins.arg -verbose : output \|\| av\_localo

-  range [1.0] Range of orientation deviations to consider

-  nbins [100] Number of bins to use for histogram

Calculates the projection of the lipid orientation vector onto the plane
for each lipid and then bins the absolute values of these vectors.

orient\_order : -verbose : output \|\| time\_vs\_oop

Calculates the orientational order parameter for each lipid through a
call to the espresso command .

stress\_tensor : -verbose : output \|\| time\_vs\_stress\_tensor

Calculates all 9 elements of the pressure tensor for the system through
a call to the espresso command

pressure : -verbose : output \|\| time\_vs\_pressure

Calculates the isotropic pressure through a call to . Results are
printed as a list of the various contributions in the following order: ,
, , , , . Where is the instantaneous pressure obtained directly from the
barostat.

stray : -verbose : output \|\| time\_vs\_stray

Calculates the number of stray lipids based on a call to .

Adding a new routine
~~~~~~~~~~~~~~~~~~~~

To add a new analysis routine you should create a new file called which
will contain all of your code. At the top of this file you should
declare a namespace for your analysis code and include all of the
internal variables inside that namespace as follows;

| namespace eval ::mbtools::analysis::myanalysis { variable av\_myresult
| variable av\_myresult\_i
| variable f\_tvsresult
| variable verbose
| namespace export setup\_myanalysis
| namespace export analyze\_myanalysis
| namespace export printav\_myanalysis
| namespace export resetav\_myanalisis
| }

Import your new file into the analysis package by adding a line like the
following to the file.

source [file join [file dirname [info script]] myanalysis.tcl]

You then need to implement the following essential functions within your
new namespace.

-  Typically you would use this function to initialise variables and
   open files.

   Called by . Arguments are allowed.

-  This function should print results to a file.

   Called by . Arguments are not allowed.

-  This function performs the actual analysis and should update the
   storage and averaging variables. Called by . Arguments are not
   allowed.

-  This function should update averages and reset variables accordingly
   depending on your requirements.

   Called by . Arguments are not allowed.

If any of these functions is not implemented the program will probably
crash.

System generation
-----------------

Package for setting up lipid membrane systems in a variety of
geometrical shapes.

Basic commands
~~~~~~~~~~~~~~

| ::mbtools::system\_generation::setup\_system : [system\_specs]
| [iboxl] [moltypes]

-  This is a list of structures called system specifications. Each such
   system specification in turn should be a list consisting of a
   geometry and a list detailing the number of each molecule type i.e.

   set system\_spec { geometry n\_molslist }

   The should be specified as a list with two elements. The first
   element should be a string “geometry” identifying this list as a
   geometry. The second element is a string containing the name of a
   geometry type followed by arguments to be passed to the routine .

   The should be specified as a list with two elements. The first
   element should be a string “n\_molslist” identifying this list as an
   n\_molslist. The second element is a list each element of which
   specifies a molecule type and the number of such molecules.

-  A list containing the lengths of each of the box side lengths.

-  A list, each element of which specifies a molecule type and type
   information. The exact format and requirements of this list are
   detailed for each molecule separately (see below for a list of
   molecule types and their requirements) however regardless of mol type
   the first two elements of the list must be a and a string specifying
   the moltype respectively.

Sets up the system including generating topologies and placing molecules
into specified geometries. Each geometry and list of molecules to be
placed into that geometry are grouped into a system spec.

Example:

The following code sets out the molecule types to be used in the
simulation by setting a list called . In this case two different lipid
types are setup and assigned to moltypeids 0 and 1 respectively. Moltype
0 will consist of three beads per lipid, the first of which is of
atomtype 0 and the second and third of which are of atomtype 1. Bonds in
the lipid will be of type 0 and 1 respectively.(see the function for
further details).

| set moltypes [list { 0 lipid { 0 1 1 } { 0 1 } }
| { 1 lipid { 0 2 2 2 } { 0 2 } } ]

We then construct system specs for a flat bilayer and a spherical
bilayer and group these into a list.

First the spherical

| set geometry { geometry “sphere -shuffle -c { 0.0 0.0 15.0 } ” }
| set n\_molslist { n\_molslist { { 0 1000 } } }
| lappend spherespec $geometry
| lappend spherespec $n\_molslist

The flat system\_spec

| set geometry { geometry “flat -fixz” }
| set n\_molslist { n\_molslist { { 1 3000 } } }
| lappend bilayerspec $geometry
| lappend bilayerspec $n\_molslist

Now group together the into a master list

| lappend system\_specs $spherespec
| lappend system\_specs $bilayerspec

Make the call to

| ::mbtools::system\_generation::setup\_system $system\_specs
| [setmd box\_l] $moltypes

::mbtools::system\_generation::get\_trappedmols :

returns the internal list variable which keeps track of all molecules
that have been trapped by their center of mass. This function should be
called after setup and would then typically be passed to the function .

::mbtools::system\_generation::get\_userfixedparts :

returns the internal list variable which keeps track of all particles
that have been fixed in position during the setup. This is useful for
later releasing particles after warmup routines have been completed.

::mbtools::system\_generation::get\_middlebead :

returns the internal variable .

Available geometries
~~~~~~~~~~~~~~~~~~~~

flat : -fixz -bondl.arg -crystal -half -pancake -shuffle

-  Fix the vertical positions of all particles. The ids of these
   particles are added to the list of which can later be obtained
   through a call to .

-  Sets lipids on a grid, instead of randomly.

-  Creates a halfbilayer (i.e. periodic only along one direction).
   Useful to measure a line tension.

-  Creates a spherical and flat bilayer. The diameter of the pancake
   cannot exceed the box\_l.

-  shuffle the topology prior to placing the lipids. This is required
   for a random lipid distribution because otherwise the lipids will be
   placed on the sphere in the order they appear in the topology

Creates a flat bilayer in the XY plane by random placement of lipids.

sphere : -c.arg -initarea.arg -bondl.arg -shuffle

-  [] The location of the center of the sphere relative to the center of
   the box

-  [1.29] An initial guess for the area per lipid. This guess is used to
   compute initial sphere dimensions based on the number of lipids. This
   initial guess is then iteratively refined until all lipids can be fit
   uniformly on the sphere.

-  shuffle the topology prior to placing the lipids. This is required
   for a random lipid distribution because otherwise the lipids will be
   placed on the sphere in the order they appear in the topology

Creates a spherical vesicle by placing molecules in an ordered manner at
uniform density on the surface of the sphere. Molecules are assumed to
have a uniform cross sectional area and closely matched (though not
identical) lengths. The radius of the vesicle will depend on the number
of lipids and the area per lipid.

sphere\_cap : -r.arg -half -c.arg -initarea.arg -bondl.arg -shuffle

-  [10.0] The radius of the whole sphere where the cap is shaped

-  Create a half of sphere with the amount of molecules available

-  [] The location of the center of the sphere relative to the center of
   the box

-  [1.29] An initial guess for the area per lipid. This guess is used to
   compute initial sphere dimensions based on the number of lipids. This
   initial guess is then iteratively refined until all lipids can be fit
   uniformly on the sphere.

-  shuffle the topology prior to placing the lipids. This is required
   for a random lipid distribution because otherwise the lipids will be
   placed on the sphere in the order they appear in the topology

Creates a spherical cap which is part of a vesicle of a radius , by
placing molecules in an ordered manner at uniform density on the surface
of the sphere. Molecules are assumed to have a uniform cross sectional
area and closely matched (though not identical) lengths. If the option
is defined, the radius of the vesicle will depend on the number of
lipids and the area per lipid.

torus : -c.arg -initarea.arg -ratio.arg -bondl.arg -shuffle

-  [] The location of the center of the torus relative to the center of
   the box.

-  [1.29] An initial guess for the area per lipid. This guess is used to
   compute initial radii based on the number of lipids. This initial
   guess is then iteratively refined until all lipids can be fit
   uniformly on the torus.

-  [1.4142] Ratio of major toroidal radius to minor toroidal radius.
   Default value is for the Clifford torus.

-  shuffle the topology prior to placing the lipids. This is required
   for a random lipid distribution because otherwise the lipids will be
   placed on the torus in the order they appear in the topology.

Creates a toroidal vesicle by placing molecules in an ordered manner at
uniform density on the surface of the torus. Molecules are assumed to
have a uniform cross sectional area and closely matched (though not
identical) lengths. The two radii of the torus will depend on the number
of lipids, the area per lipid and the ratio between radii.

cylinder : -c.arg -initarea.arg -bondl.arg -shuffle

-  [0.0 0.0 0.0]

-  [1.29]

-  shuffle the topology prior to placing the lipids.

Creates a cylinder which spans the box along one dimension by placing
molecules uniformly on its surface. Works in a similar way to the sphere
routine.

random : -exclude.arg -inside.arg -shuffle -bondl.arg

-  | [] an exclusion zone definition suitable for passing to
   | .

-  | [] an inclusion zone definition suitable for passing to
   | .

-  shuffle the topology prior to placing the lipids.

Places molecules randomly in space with a (sortof) random orientation
vector. If an exclusion zone is defined, then no molecules will be
placed such that their positions are within the zone. If an inclusion
zone if defined, then no molecules will be place outside this zone. For
instance,

| set geometry { geometry “random -exclude { sphere { 0.0 0.0 0.0 } 4.0
  }
  -inside { cuboid { 0.0 0.0 0.0 } { 15.0 15.0 15.0 } }” }

will randomly place molecules within the volume between a sphere with a
radius of :math:`4.0` and a cuboid with dimension
:math:`15.0 \times 15.0 \times 15.0` at the origin.

readfile : -ignore.arg -f.arg -t.arg

-  [] particle properties to be ignored during the file read.

-  [] The file containing the configuration to be used for setup. Must
   be an espresso blockfile with box length, particle and bonding
   information.

-  [] The topology file corresponding to the file to be read.

-  [0.000001] Tolerance for comparison of box dimensions.

Use particle positions contained in a file to initialise the locations
of particles for a particular geometry. The box dimensions in the file
and those set by the user are compared and an error is returned if they
are not the same to within a tolerance value of . Even though we read
from a file we also generate a topology from the and this topology is
compared with the topology that is read in to check if the number of
particles are the same.

| singlemol : -c.arg -o.arg -trapflag.arg -ctrap.arg
| -trapspring.arg -bondl.arg

-  [ 0.0 0.0 0.0 ] The molecule center. Exactly what this means depends
   on the molecule type.

-  [ 0.0 0.0 1.0 ] The orientation vector for the molecule. This is also
   molecule type dependent

-  [ 0 0 0 ] Set this optional argument to cause a molecule to be
   trapped by its center of mass. You should give three integers
   corresponding to each of the three coordinate axes. If a value of 1
   is given then motion in that axis is trapped.

-  [ “” ] Set this optional argument to the central point of the trap.
   This works much like an optical trap in that molecules will be
   attracted to this point via a simple harmonic spring force

-  [ 20 ] The spring constant for the trap potential (harmonic spring).

Simply place a single molecule at the desired position with the desired
orientation.

Adding a new geometry
~~~~~~~~~~~~~~~~~~~~~

To create a routine for setting up a system with a new type of geometry
. Start by creating a new file inside the directory. The new file should
declare a new namespace as a sub namespace of and export the proceedure
. Thus your file should begin with the lines

| namespace eval ::mbtools::system\_generation::mygeom {
| namespace export create\_mygeom
| }

Import your new file into the package by adding a line like the
following to the file

source [file join [file dirname [info script]] mygeom.tcl]

You then need to implement the proceedure within your new namespace as
follows

::mbtools::system\_generation::mygeom::create\_mygeom args

Available molecule types
~~~~~~~~~~~~~~~~~~~~~~~~

| lipid : typeinfo : { moltypeid “lipid” particletypelist
| bondtypelist }

-  A list of the particle types for each atom in the lipid. The
   particles are placed in the order in which they appear in this list.

-  A list of two s. The first id is used for bonds between consecutive
   beads in the lipid. The second defines the pseudo bending potential
   which is a two body bond acting across beads separated by exactly one
   bead.

Places atoms in a line to create a lipid molecule.

| hollowsphere : typeinfo : { moltypeid “hollowsphere”
| sphereparticlelist bondtype natomsfill }

-  A list of the particle types for each atom in the hollowsphere. The
   atoms that make up the outer shell must be listed first followed by
   the atoms that make up the inner filling.

-  The typeid for bonds linking atoms in the outer shell.

-  Number of filler atoms. The atom types for these will be obtained
   from the last in the .

Creates a sphere of beads arranged such that they have an approximate
spacing of and such that they optimally cover the sphere. The optimal
covering is obtained using the routines which are copyright R. H.
Hardin, N. J. A. Sloane and W. D. Smith, 1994, 2000. Thus the routine
will only work if you have installed icover and if you can successfully
run it from the command line in the directory that you started your
espresso job. These routines are serious overkill so if anybody can
think of a nice simple algorithm for generating a covering of the sphere
let us know.

| protein : typeinfo : { moltypeid “protein” particletypelist
| bondtypelist }

-  A list of the particle types for each atom in the protein.

-  A list of bondtypeids.

Create a protein molecule.

| spanlipid : typeinfo : { moltypeid “protein” particletypelist
| bondtypelist }

-  A list of the particle types for each atom in the lipid. Since this
   is a spanning lipid the first and last elements of this list would
   typically be head beads.

-  A list of two s with the same meaning as explained above for standard
   lipids.

Create a lipid which spans across the bilayer.

Adding a new molecule type
~~~~~~~~~~~~~~~~~~~~~~~~~~

To add a new molecule type you need to define a proceedure which
determines how the atoms that make up the molecule should be placed.
This proc will live directly in the namespace. Examples can be found in
.

In order to register your new molecule type to allow placement in any
geometry you need to add a call to it in the function . Make sure that
all arguments to your routine are included in this function call.

Utils
-----

Useful utilities routines for various types. Includes file management,
basic geometry and math procedures.

Setup commands
~~~~~~~~~~~~~~

| ::mbtools::utils::setup\_outputdir : [outputdir] -paramsfile.arg
| -tabdir.arg -tabnames.arg -startf.arg -ntabs.arg

-  Complete path of the directory to be setup. At least the parent of
   the directory must exist

-  [] Name of a file to be copied to the output directory

-  [] Full path name of the directory where forcetables are kept

-  [] Complete list of forcetables to be used in the simulation. These
   will be copied to the output directory

This routine is designed to setup a directory for simulation output. It
copies forcetables and the parameter file to the directory after
creating it if necessary.

::mbtools::utils::read\_startfile : [file]

-  Complete path of the file to be read. Should be an espresso
   blockfile.

Read in particle configuration from an existing file or simulation
snapshot

::mbtools::utils::read\_checkpoint : [dir]

-  Directory containing the checkpoint file which must be called .

Read in a checkpoint and check for success. Warn if the checkpoint does
not exist.

::mbtools::utils::read\_topology : [file]

-  Complete path of the file that contains the topology information.

Read in the topology from a file and then execute the command of
ESPResSo.

::mbtools::utils::set\_topology : [topo]

-  A valid topology.

Set the given topology and then execute the command of ESPResSo.

::mbtools::utils::set\_bonded\_interactions : [bonded\_parms]

-  A list of bonded interactions. Each element of this list should
   contain all the appropriate arguments in their correct order for a
   particular call to the espresso command. See the espresso command for
   a list of possible bonded interactions and correct syntax.

Set all the bonded interactions.

::mbtools::utils::set\_nb\_interactions : [nb\_parms]

-  A list of interactions. Each element of this list should contain all
   the appropriate arguments in their correct order for a particular
   call to the espresso command. See the espresso command for a list of
   possible non-bonded interactions and correct syntax.

Set all the bonded interactions.

::mbtools::utils::init\_random : [n\_procs]

-  The number of processors used in this job.

Initialize the random number generators on each processor based on the
current time with a fixed increment to the time seed used for each proc.

| ::mbtools::utils::initialize\_vmd : [flag] [outputdir]
| [ident] -extracommands.arg

-  Depending on the value of this parameter initialize vmd to one of its
   possible states:

   -  interactive : VMD is started and a connection to espresso
      established for immediate viewing of the current espresso process.
      With some luck this might even work sometimes! If VMD doesn’t get
      a proper connection to espresso then it will crash.

   -  offline : Just constructs the appropriate and files and writes
      them to the output directory so that files generated with writepdb
      can be viewed with .

   -  default : Any value other than those above for flag will just
      result in vmd not being initialized.

-  The directory where vmd output will be written.

-  A basename to be be given to vmd files.

-  [] A list of strings each of which will be written to the end of the
   . Use this to give additional commands to vmd.

Prepare for vmd output.

Warmup commands
~~~~~~~~~~~~~~~

| ::mbtools::utils::warmup : [steps] [times] -mindist.arg
| -cfgs.arg -outputdir.arg -vmdflag.arg -startcap.arg
| -capgoal.arg

-  number of integration steps used in each call to integrate.

-  number of times to call the integrate function during warmup.

-  [0] Terminate the warmup when the minimum particle distance is
   greater than this criterion. A value of 0 (default) results in this
   condition being ignored. If a condition is imposed this routine can
   become very very slow for large systems.

-  [-1] Write out a configuration file every cfgs calls to integrate.

-  [./] The directory for writing output.

-  [offline] If this flag is set to “offline” (default) pdb files will
   be generated for each configuration file generated.

-  [5] Starting value for the forcecap.

-  [1000] For the purposes of calculating a cap increment this value is
   used as a goal. The final forcecap will have this value.

Perform a series of integration steps while increasing forcecaps from an
initially small value.

Topology procs
~~~~~~~~~~~~~~

::mbtools::utils::maxpartid : [topo]

-  A valid topology.

Find the maximum particle id in a given topology.

::mbtools::utils::maxmoltypeid : [topo]

-  A valid topology.

Find the maximum molecule type id.

::mbtools::utils::listnmols : [topo]

-  A valid topology.

Construct a list with the number of molecules of each molecule type.

::mbtools::utils::minpartid : [topo]

-  A valid topology.

Minimum particle id for the given topology.

::mbtools::utils::minmoltype : [topo]

-  A valid topology/

Minimum molecule type id for this topology.

::mbtools::utils::listmoltypes : [topo]

-  A valid topology.

Make a list of all the molecule types in a topology. Makes a check for
duplication which would occur for an unsorted topology.

::mbtools::utils::listmollengths : [topo]

-  A valid topology.

Works out the length (number of atoms) of each molecule type and returns
a list of these lengths.

Math procs
~~~~~~~~~~

::mbtools::utils::dot\_product : A B

Returns A dot B

::mbtools::utils::matrix\_vec\_multiply : A B

Return the product of a matrix A with a vector B

::mbtools::utils::calc\_proportions : ilist

Calculate the number of times each integer occurs in the list ilist.

::mbtools::utils::average : data from to

-  A list of numbers to be averaged

-  Optional starting index in data

-  Optional ending index in data

Calculate the mean of a list of numbers starting from going up to .

::mbtools::utils::stdev : data from to

-  A list of numbers to find the std deviation of

-  Optional starting index in data

-  Optional ending index in data

Calculate the standard deviation of a list of numbers starting from
going up to .

::mbtools::utils::acorr : data

-  Data for which an autocorrelation is to be calculated

Calculate an autocorrelation function on a set of data.

::mbtools::utils::distance : pos1 pos2

-  A position vector

-  A position vector

Calculate the distance between two points whose position vectors are
given.

::mbtools::utils::distance\_min : pos1 pos2

-  A position vector

-  A position vector

Calculate the minimum image distance between two position vectors.

::mbtools::utils::min\_vec : pos1 pos2

-  A position vector

-  A position vector

Calculate the minimum image vector from position vector2 to postition 1,
*i.e.* pos1 - pos2.

::mbtools::utils::normalize : vec

-  The vector to be normalised

Normalize a vector

::mbtools::utils::scalevec : vec scale

-  The vector to be scaled

-  Scaling factor

Multiply all elements of a vector by a scaling factor

::mbtools::utils::uniquelist : original

-  A list possibly containing duplicate elements

Construct a list of all the unique elements in the original list
removing all duplication.

Miscellaneous procs
~~~~~~~~~~~~~~~~~~~

::mbtools::utils::trap\_mols : molstotrap

-  A list of trap values for molecules. This list would typically be
   obtained by calling immediately after the system has been setup.

Set the trap value for a list of molecules.

::mbtools::utils::isoutside : [pos] [zone]

-  The point whose status is to be determined

-  This will be a tcl list. The first element of the list must be a
   string with the name of the zone type and subsequent elements will be
   further information about the zone. Available zones are:

   -  : center radius

   -  : center {L W H}

Determines whether the point at is outside the zone. Parameter center
should be a tcl list. Returns 1 if it is and 0 if it is not.

::mbtools::utils::calc\_com : mol

-  The molecule

Calculate the center of mass of a molecule.

::mbtools::utils::centersofmass\_bymoltype : [moltypes]

-  A list of molecule type ids

Determine the center of mass of every molecule whose type matches an
item in the list moltypes. Returns a nested list where each element in
the list is itself a list of centers of mass for a given moltype.

mmsg
----

mmsg is designed to provide a more controlled way of printing messages
than the simple commands of Tcl. It has an ability to turn on or off
messages from particular namespaces.

Basic commands
~~~~~~~~~~~~~~

The following commands represent the standard interface for the package.
For consistency one should use these instead of a bare puts to standard
out. makes extensive use of these commands.

::mmsg::send : [namespace] [string] { [newline] }

-  A namespace. Typically this should be the current namespace which one
   can get via namespace current

-  The message you want printed

-  [yes] Set this to anything other than “yes” and no carriage return
   will be used after the message

The mmsg equivalent of puts. Designed for printing of simple status or
progress messages.

::mmsg::err : [namespace] [string] { [newline] }

-  A namespace. Typically this should be the current namespace which one
   can get via namespace current

-  The message you want printed

-  [yes] Set this to anything other than “yes” and no carriage return
   will be used after the message

Prints error messages and causes program to exit.

::mmsg::warn : [namespace] [string] { [newline] }

-  A namespace. Typically this should be the current namespace which one
   can get via namespace current

-  The message you want printed

-  [yes] Set this to anything other than “yes” and no carriage return
   will be used after the message

Prints warning messages.

::mmsg::debug : [namespace] [string] { [newline] }

-  A namespace. Typically this should be the current namespace which one
   can get via namespace current

-  The message you want printed

-  [yes] Set this to anything other than “yes” and no carriage return
   will be used after the message

Prints debug messages.

Control commands
~~~~~~~~~~~~~~~~

does several checks before it decides to print a message. For any given
message type it checks if that message type is allowed. It also checks
to see if the namespace given as an argument is in the allowable
namespaces list. The default behaviour is to print from the main mbtools
namespaces and the global namespace

{ :: ::mbtools::system\_generation ::mbtools::utils ::mbtools::analysis
}

Note that children of these namespaces must be explicitly enabled. All
message types except debug are also enabled by default. The following
commands allow this default behaviour to be changed.

::mmsg::setnamespaces : namespacelist

-  A list of all namespaces from which messages are to be printed

Allows control over which namespaces messages can be printed from.

::mmsg::enable : type

-  A string indicating a single message type to enable. Allowable values
   are “err”, “debug”, “send” and “warn”

Allows particular message types to be enabled: For example one could
enable debug output with

mmsg::enable “debug”

::mmsg::disable : type

-  A string indicating a single message type to disable. Allowable
   values are “err”, “debug”, “send” and “warn”

Allows particular message types to be disabled: For example one could
disable warning output with

mmsg::enable “warn”

.. [1]
   This documentation was written by Ira R. Cooke and published on his
   website. It has been transcripted by Tristan Bereau.
