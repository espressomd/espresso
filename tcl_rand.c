#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "config.h"
#include "global.h"
#include "debug.h"
#include "random.h"

/** \file tcl_rand A random generator for tcl.
 Usage: tcl_rand() for uniform double in ]0;1[
 tcl_rand(i <n>) for integer between 0 and n-1*/

int tcl_rand(ClientData data, Tcl_Interp *interp,
	  int argc, char **argv)
{
  int i_out;
  double d_out;
  char   buffer[TCL_DOUBLE_SPACE + 5];

  if (argc > 1) {
    switch(argv[1][0])
      {
      case 'i':
	if(argc < 3){
	  Tcl_AppendResult(interp, "wrong # args:  should be \"",
			   argv[0], " ?variable type? ?parameter?\"",
			   (char *) NULL);
	  return (TCL_ERROR);
	}else {
	  Tcl_GetInt(interp, argv[2], &i_out);
	  i_out = i_random(i_out);
	  sprintf(buffer, "%d", i_out);
	  Tcl_AppendResult(interp, buffer, (char *) NULL);
	  return (TCL_OK);
	} 
      case 'd':
	d_out = d_random();
	sprintf(buffer, "%f", d_out);
	Tcl_AppendResult(interp, buffer, (char *) NULL);
	return (TCL_OK);
      default:
	Tcl_AppendResult(interp, "wrong # args:  should be \"",
			 argv[0], " ?variable type? ?parameter?\"",
			 (char *) NULL);
	return (TCL_ERROR);
      }
  }else {
   d_out = d_random();
   sprintf(buffer, "%f", d_out);
   Tcl_AppendResult(interp, buffer, (char *) NULL); 
   return (TCL_OK);
  }
  printf("Error in tcl_rand(); this should never been shown\n");
  return(TCL_ERROR);
}
