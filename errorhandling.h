#ifndef ERRORHANDLING_H
#define ERRORHANDLING_H

#include <tcl.h>

/** buffer for error messages during the integration process. */
extern char *error_msg;
extern int n_error_msg;


/* request space for leaving an error message to be passed to the master node.
   Also takes care of the error counter.
   @param errlen maximal length of the error message. If you use sprintf to create the error
   message, remember to use TCL_INTEGER/DOUBLE_SPACE as usual
   @return where to put the (null-terminated) string */
char *runtime_error(int errlen);

/** check for runtime errors on all nodes. This has to be called on all nodes synchronously.
    @return the number of characters in the error messages of all nodes together. */
int check_runtime_errors();

#endif
