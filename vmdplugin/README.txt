This directory contains the source files of the VMD molfile reader
plugin for the "VMD trajectory format" VTF:

       vtfplugin.c    - source code of the plugin
       vtfplugin.txt  - description of the format
       vtftest.vmd    - VMD script to test the plugin
       vtftest.c      - program to test the plugin
       Makefile	      - GNU Makefile to build the plugin
       colloid_pe.vtf - sample VTF file of the colloid-pe system
       colloid_pe.tcl - Espresso script of a colloid-pe system, used
		        to generate colloid_pe.vtf

WHY VTF?
--------
The VMD trajectory format has a number of advantages when used
together with Espresso instead of the PSF and PDB formats:

1. A whole trajectory can be written to a single ".vtf" file, instead
   of a PSF file and a series of PDB files that needs to be loaded into
   VMD via the "loadseries" command. However, it is also possible to
   write every timestep and the structure to separate files (".vsf"
   and ".vcf"). Of course, the VSF-files can also be used in
   conjunction with Espresso's IMD capabilities.

2. In the structure description, the VDW-radii of the atoms can be
   specified. In fact, all characteristics that a VMD atom may have
   can be specified (name, type, resid, resname, segid, ...)

3. The system length information is passed to VMD, so that periodic
   images and scripts like "pbcwrap" and "pbcbox" from the VMD script
   library can be used. 
   VMD script library: http://www.ks.uiuc.edu/Research/vmd/script_library/

4. If a system has fixed particles, it is possible to write only the
   coordinates of the particles that have been moved to the file. This
   can save some disk space for large systems.

5. The files are human readable.


COMPILATION
-----------
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

INSTALLATION
------------
If you want to install the plugin, simply copy the file vtftest.so to
an arbitrary directory. To activate the plugin in VMD, add the
following command to your ~/.vmdrc file:

       vmd_plugin_scandirectory <dir> *.so

where <dir> should be substituted with the directory where the plugin
was installed. Once the plugin has been activated in VMD, you can load
VTF, VSF and VCF files in VMD.

GENERATING VTF FILES
--------------------
The VTF format is described in the file "vtfplugin.txt". 
Espresso contains two functions to generate VTF,VSF and VTF files:
writevsf <file> <OPTIONS>
writevcf <file> <OPTIONS>

The synatx and usage of the commands is desribed in the Espresso
User's Guide.
