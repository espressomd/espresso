/** \file blockfile_tcl.c
    Implementation of \ref blockfile_tcl.h "blockfile_tcl.h".
*/
#include <ctype.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "global.h"
#include "blockfile.h"
#include "blockfile_tcl.h"

int blockfile(ClientData data, Tcl_Interp *interp,
	      int argc, char **argv)
{
  char title[MAXBLOCKTITLE];
  char buffer[1024], *name;
  int tcl_file_mode;
  Tcl_Channel channel;
  Tcl_CmdInfo cmdInfo;
  int i, j, openbrackets, len, exists;

  if (argc < 4) {
    Tcl_AppendResult(interp, "wrong # args:  should be \"",
		     argv[0], " <file> read|write <what> <param>? ...\"",
		     (char *) NULL);
    return (TCL_ERROR);
  }

  if ((channel = Tcl_GetChannel(interp, argv[1], &tcl_file_mode)) == NULL)
    return (TCL_ERROR);

  /* the write commands */
  if (!strncmp(argv[2], "write", strlen(argv[2]))) {
    if (!(tcl_file_mode & TCL_WRITABLE)) {
      Tcl_AppendResult(interp, "\"", argv[1], "\" not writeable", (char *) NULL);
      return (TCL_ERROR);
    }

    /* write start tag */
    if (!strncmp(argv[3], "start", strlen(argv[3]))) {
      if (argc != 5 || strlen(argv[4]) > MAXBLOCKTITLE - 1) {
	Tcl_AppendResult(interp, "please give a short (< 63 bytes) title for your block and nothing more",
			 (char *) NULL);
	return (TCL_ERROR);
      }
      if (block_writestart(channel, argv[4]) != 0) {
	Tcl_AppendResult(interp, "\"", argv[1], "\" could not write start tag",
			 (char *) NULL);
	return (TCL_ERROR);
      }
      return (TCL_OK);
    }
    /* write end tag */
    else if (!strncmp(argv[3], "end", strlen(argv[3]))) {
      if (argc != 4) {
	Tcl_AppendResult(interp, "\"", argv[0],
			 " write end\" takes no parameters",
			 (char *) NULL);
	return (TCL_ERROR);
      }
      if (block_writeend(channel) != 0) {
	Tcl_AppendResult(interp, "\"", argv[1], "\" could not write end tag",
			 (char *) NULL);
	return (TCL_ERROR);
      }
      return (TCL_OK);
    }
    /* dump a global variable */
    else if (!strncmp(argv[3], "variable", strlen(argv[3]))) {
      if (argc != 5) {
	Tcl_AppendResult(interp, "variable name missing",
			 (char *) NULL);
	return (TCL_ERROR);
      }

      for (i = 0; fields[i].data != NULL; i++) {
	int len = strlen(argv[4]);
	if (len >= 1024) {
	  Tcl_AppendResult(interp, "\"", argv[4],
			   "\" is a brain-damaging long variable name and cannot be written to a file.",
			   (char *) NULL);
	  return (TCL_ERROR);
	}
	if (!strncmp(argv[4], fields[i].name, len)) {
	  len = strlen(fields[i].name);
	  if (block_writestart(channel, "variable") ||
	      (Tcl_Write(channel, (char *)fields[i].name, len) != len) ||
	      (Tcl_Write(channel, " ", 1) != 1)) {
	    Tcl_AppendResult(interp, "\"", argv[1], "\" could not write variable",
			     (char *) NULL);
	    return (TCL_ERROR);
	  }
	  block_write_data(channel, fields[i].type, fields[i].dimension, fields[i].data);

	  if (block_writeend(channel) != 0 ||
	      (Tcl_Write(channel, "\n", 1) != 1)) {
	    Tcl_AppendResult(interp, "\"", argv[1], "\" could not write data",
			     (char *) NULL);
	    return (TCL_ERROR);
	  }
	  return (TCL_OK);
	}
      }
      Tcl_AppendResult(interp, "unknown variable \"",
		       argv[4], "\"", (char *) NULL);
      return (TCL_ERROR);
    }
  }
  /* the read commands */
  else if (!strncmp(argv[2], "read", strlen(argv[2]))) {
    if (!(tcl_file_mode & TCL_READABLE)) {
      Tcl_AppendResult(interp, "\"", argv[1], "\" not readable", (char *) NULL);
      return (TCL_ERROR);
    }
    /* consume start tag only */
    if (!strncmp(argv[3], "start", strlen(argv[3]))) {
      if (argc != 4) {
	Tcl_AppendResult(interp, "\"", argv[0],
			 " read start\" takes no parameters",
			 (char *) NULL);
	return (TCL_ERROR);
      }
      if ((openbrackets = block_startread(channel, title)) == 1) {
	Tcl_AppendResult(interp, title, (char *)NULL);	
	return (TCL_OK);
      }
      else {
	Tcl_AppendResult(interp, "\"", argv[1], "\" could not read block start",
			 (char *) NULL);
	return (TCL_ERROR);
      }
    }
    /* consume next block and return data or variable contained */
    else if (!strncmp(argv[3], "auto", strlen(argv[3]))) {
      if (argc != 4) {
	Tcl_AppendResult(interp, "\"", argv[0],
			 " read start\" takes no parameters",
			 (char *) NULL);
	return (TCL_ERROR);
      }
      switch (openbrackets = block_startread(channel, title)) {
      case 0:	
      case 1:
	break;
      case -1:
	Tcl_AppendResult(interp, "eof", (char *)NULL);
	return (TCL_OK);
      default:
	Tcl_AppendResult(interp, "\"", argv[1], "\" could not read block start",
			 (char *) NULL);
	return (TCL_ERROR);
      }
      /* set global variable */
      if (!strcmp(title, "variable")) {
	if ((openbrackets = block_continueread(channel, openbrackets,
					       buffer, sizeof(buffer), ' ')) != 1) {
	  Tcl_AppendResult(interp, "\"", argv[1], "\" could not read variable name",
			   (char *) NULL);
	  return (TCL_ERROR);
	}
	for (i = 0; fields[i].data != NULL; i++)
	  if (!strcmp(fields[i].name, buffer)) {
	    if (block_read_data(channel, fields[i].type,
				fields[i].dimension, fields[i].data) != 0) {
	      Tcl_AppendResult(interp, "\"", argv[1], "\" could not read data",
			       (char *) NULL);
	      return (TCL_ERROR);
	    }
	    Tcl_AppendResult(interp, "setvar ", buffer, (char *) NULL);
	    if (block_continueread(channel, openbrackets, buffer, sizeof(buffer), 0) != 0) {
	      Tcl_ResetResult(interp);
	      Tcl_AppendResult(interp, "variable block for does not terminate correctly",
			       (char *) NULL);
	      return (TCL_ERROR);
	    }
	    for (j = 0; buffer[j] != 0; j++)
	      if (!isspace(j)) {
		Tcl_ResetResult(interp);
		Tcl_AppendResult(interp, "variable block for ", fields[i].name,
				 " contains garbage \"", buffer, "\"",
				 (char *) NULL);
		return (TCL_ERROR);
	      }
	    return (TCL_OK);
	  }
	/* not a tcl_md variable, so just dump variable / value pair */
	Tcl_AppendResult(interp, "uservar ", buffer, " ", (char *) NULL);
	while ((openbrackets =
		block_continueread(channel, openbrackets,
				   buffer, sizeof(buffer), 0)) > 0) {
	  Tcl_AppendResult(interp, buffer, (char *) NULL);
	}
	if (openbrackets < 0) {
	  Tcl_AppendResult(interp, "\"", argv[1], "\" could not read data",
			   (char *) NULL);
	  return (TCL_ERROR);
	}
	Tcl_AppendResult(interp, buffer, (char *) NULL);
	return (TCL_OK);
      }
      /* unknown field */
      else {
	len = 21; /* blockfile_read_auto_ + \0 */
	len += strlen(title);
	name = malloc(len);
	strcpy(name, "blockfile_read_auto_");
	strcat(name, title);
	exists = Tcl_GetCommandInfo(interp, name, &cmdInfo);
	free(name);
	if (exists) {
	  if (cmdInfo.proc(cmdInfo.clientData, interp,
			   argc, argv) != TCL_OK)
	    return (TCL_ERROR);
	  if (block_continueread(channel, openbrackets, buffer, sizeof(buffer), 0) != 0) {
	    Tcl_ResetResult(interp);
	    Tcl_AppendResult(interp, "variable block for does not terminate correctly",
			     (char *) NULL);
	    return (TCL_ERROR);
	  }
	  for (j = 0; buffer[j] != 0; j++)
	    if (!isspace(j)) {
	      Tcl_ResetResult(interp);
	      Tcl_AppendResult(interp, "variable block for ", fields[i].name,
			       " contains garbage \"", buffer, "\"",
			       (char *) NULL);
	      return (TCL_ERROR);
	    }
	  return (TCL_OK);
	}
	
	Tcl_AppendResult(interp, "usertag ", title, (char *) NULL);
	if (openbrackets != 0) {
	  Tcl_AppendResult(interp, " ", (char *) NULL);	  
	  while ((openbrackets =
		  block_continueread(channel, openbrackets,
				     buffer, sizeof(buffer), 0)) > 0) {
	    Tcl_AppendResult(interp, buffer, (char *) NULL);
	  }
	  if (openbrackets < 0) {
	    Tcl_AppendResult(interp, "\"", argv[1], "\" could not read data",
			     (char *) NULL);
	    return (TCL_ERROR);
	  }
	  Tcl_AppendResult(interp, buffer, (char *) NULL);
	}
	return (TCL_OK);	
      }
    }
    /* consume until end tag */
    else if (!strncmp(argv[3], "toend", strlen(argv[3]))) {
      if (argc != 4) {
	Tcl_AppendResult(interp, "\"", argv[0],
			 " read  toend\" takes no parameters",
			 (char *) NULL);
	return (TCL_ERROR);
      }
      openbrackets = 1;
      while ((openbrackets =
	      block_continueread(channel, openbrackets,
				 buffer, sizeof(buffer), 0)) > 0) {
	Tcl_AppendResult(interp, buffer, (char *) NULL);
      }
      if (openbrackets < 0) {
	Tcl_AppendResult(interp, "\"", argv[1], "\" could not read data",
			 (char *) NULL);
	return (TCL_ERROR);
      }
      Tcl_AppendResult(interp, buffer, (char *) NULL);
      return (TCL_OK);
    }
  }

  /* not a native action, try script support */
  len = 12; /* blockfile_ + _ + \0 */
  len += strlen(argv[2]) + strlen(argv[3]);
  name = malloc(len);
  strcpy(name, "blockfile_");
  strcat(name, argv[2]);
  strcat(name, "_");
  strcat(name, argv[3]);
  exists = Tcl_GetCommandInfo(interp, name, &cmdInfo);
  free(name);
  if (exists) {
    return cmdInfo.proc(cmdInfo.clientData, interp,
			argc, argv);
  }
  Tcl_AppendResult(interp, "unknown action \"", argv[2], " ", argv[3],
		   "\"", (char *)NULL);  
  return (TCL_ERROR);
}
