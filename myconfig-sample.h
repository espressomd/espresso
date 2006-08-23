/** \file myconfig-sample.h

    This is a sample for the file myconfig.h.

    Uncomment and add any of the following lines to myconfig.h to
    activate the corresponding feature of Espresso. It is recommended
    to turn only those features on that you actually need to optimize
    the performance of Espresso for your problem. 

    To access the information on the compilation status of the code
    you are working with in your Espresso Tcl-script, use the
    corresponding \ref tcl_features "Tcl-commands".

    If you add a new feature to Espresso, you also have to add the
    corresponding lines in the function \ref compilation_callback and
    to add documentation in <tt>doc/text/features.doc</tt>.
 
    <b>Responsible:</b> 
    <a href="mailto:arnolda@mpip-mainz.mpg.de">Axel</a>

*/
/* #define PARTIAL_PERIODIC */
/* #define ELECTROSTATICS */
/* #define ROTATION */
/* #define EXTERNAL_FORCES */
/* #define CONSTRAINTS */
/* #define MASS */
/* #define EXCLUSIONS */
/* #define COMFORCE */
/* #define COMFIXED */
/* #define MOLFORCES */
/* #define BOND_CONSTRAINT */

/* #define TABULATED */
/* #define LENNARD_JONES */
/* #define LJ_WARN_WHEN_CLOSE */
/* #define MORSE */
/* #define LJCOS */
/* #define LJCOS2 */
/* #define BUCKINGHAM */
/* #define SOFT_SPHERE */

/* Note: Activate ONLY ONE bonded angle potential out of the following! */
/* #define BOND_ANGLE_HARMONIC */
/* #define BOND_ANGLE_COSINE */
/* #define BOND_ANGLE_COSSQUARE */

/* #define NEMD */
/* #define NPT */ 
/* #define DPD */
/* #define LB */
