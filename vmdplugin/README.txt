This directory contains the source files of the VMD molfile reader
plugin for the "VMD trajectory format":

       vtfplugin.c    - source code of the plugin
       vtfplugin.txt  - description of the format
       vtftest.vmd    - VMD script to test the plugin
       vtftest.c      - program to test the plugin
       Makefile	      - GNU Makefile to build the plugin
       colloid_pe.vtf - sample VTF file of the colloid-pe system
       colloid_pe.tcl - Espresso script of a colloid-pe system, used
		        to generate colloid_pe.vtf

COMPILING
---------
To compile the plugin, adapt the Makefile to your system. Usually,
it should be enough to modify the VMDDIR line to match the path where
you have installed VMD. Then you can compile the plugin via
"make". This will generate the plugin "vtfplugin.so" and the test
program "vtftest".

TESTING
-------
When the plugin has been compiled successfully, you can test the
plugin by running the VMD test script vtftest.vmd:

       vmd -e vtftest.vmd

INSTALLING
----------
If you want to install the plugin, simply copy the file vtftest.so to
an arbitrary directory. To activate the plugin in VMD, add the
following command to your ~/.vmdrc file:

       vmd_plugin_scandirectory <dir> *.so

where <dir> should be substituted with the directory where the plugin
was installed. Once the plugin has been activated in VMD, you can load
VTF, VSF and VCF files in VMD.

GENERATING VTF FILES
--------------------
Espresso contains two functions to generate VTF,VSF and VTF files:

writevsf <file> [short|verbose] [radius <radii>] [typedesc <typedesc>]

  Writes a structure block describing the system's structure to
  <file>. 
  The atom ids used in the file are identical to Espresso's particle
  ids. This makes it easy to write additional structure lines (as
  described in vtfplugin.txt) to the file, e.g. to add "resname"s to
  particle compounds, like chains.
  The output of this file can be used in a standalone VSF file, or at
  the beginning of a trajectory VTF file that contains a trajectory of
  a whole simulation.

  OPTIONS:
    short|verbose
      The keywords "short" and "verbose" specify, whether the output
      is in a human-readable, but somewhat longer format ("verbose"),
      or in a more compact form ("short"). The default is "verbose".

    radius <radii>
      <radii> is either "auto", or a Tcl-list describing the radii of
      the different particle types. When the keyword "auto" is used
      and a Lennard-Jones interaction between two particles of the
      given type is defined, the radius is set to be half the LJ sigma
      plus the LJ shift. Otherwise, the radius 0.5 is substituted. The
      default is "auto".
      Example: radius {{0 2.0} {1 auto} {2 1.0}}

   typedesc <typedesc>
     <typedesc> is a Tcl-list giving additional VTF atom-keywords (see
     vtfplugin.txt) to specify additional VMD characteristics of the
     atoms of the given type.
     Example: typedesc {{0 "segid colloid"} {1 "segid pe"}}

writevcf <file> [short|verbose] [folded|unfolded]

  Writes a coordinate (or timestep) block that contains all
  coordinates of the system's particles to <file>.
  
  OPTIONS:
    short|verbose
      The keywords "short" and "verbose" specify, whether the output
      is in a human-readable, but somewhat longer format ("verbose"),
      or in a more compact form ("short"). The default is "verbose".
    folded|absolute
      Specify whether the particle positions are written in absolute
      coordinates ("absolute") or folded into the central image of a
      periodic system ("folded"). The default is "absolute".
