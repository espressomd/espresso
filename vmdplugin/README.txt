This directory contains the source files of the VMD molfile reader
plugin for the "VMD trajectory format":

       vtfplugin.c    - source code of the plugin
       vtfplugin.txt  - description of the format
       vtftest.c      - program to test the plugin
       Makefile	      - GNU Makefile to build the plugin
       colloid_pe.vtf - sample VTF file of the colloid-pe system
       colloid_pe.tcl - Espresso script of a colloid-pe system, used
		        to generate colloid_pe.vtf

To compile the plugin, adapt the Makefile to your system. Usually,
it should be enough to modify the VMDDIR line to match the path where
you have installed VMD. Then you can compile the plugin via
"make". This will generate the plugin "vtfplugin.so" and the test
program "vtftest".

When the plugin has been compiled successfully, you can activate the
plugin in VMD by adding the following command to the ~/.vmdrc file or
by executing this command in in the VMD console:

       vmd_plugin_scandirectory <dir> *.so

where <dir> should be substituted with the directory where the plugin
was compiled.

Once the plugin has been activated in VMD, you can load VTF, VSF and
VCF files in VMD.
