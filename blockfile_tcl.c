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
  char buffer[1024];
  int tcl_file_mode;
  Tcl_Channel channel;
  int i, j, openbrackets;

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
	  switch (fields[i].type) {
	  case TYPE_INT:
	    if (fields[i].dimension == 1)
	      block_write_int(channel, *(int *)fields[i].data);
	    else {
	      int dim[3] = { 0, 1, 1 };
	      dim[0] = fields[i].dimension;
	      block_write_int_array(channel, dim, (int *)fields[i].data);
	    }
	    break;
	  case TYPE_DOUBLE:
	    if (fields[i].dimension == 1)
	      block_write_double(channel, *(int *)fields[i].data);
	    else {
	      int dim[3] = { 0, 1, 1 };
	      dim[0] = fields[i].dimension;
	      block_write_double_array(channel, dim, (double *)fields[i].data);
	    }
	    break;
	  default: ;
	  }
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
    Tcl_AppendResult(interp, "unknown write action \"", argv[3],
		     "\"", (char *)NULL);
    return (TCL_ERROR);
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
	if (block_continueread(channel, 1, buffer, sizeof(buffer),' ') != 1) {
	  Tcl_AppendResult(interp, "\"", argv[1], "\" could not read variable name",
			   (char *) NULL);
	  return (TCL_ERROR);
	}
	for (i = 0; fields[i].data != NULL; i++)
	  if (!strcmp(fields[i].name, buffer)) {
	    int dim[3] = {1, 1, 1};
	    dim[0] =  fields[i].dimension;
	    if (block_read_any_data(channel, fields[i].type, dim, fields[i].data) != 0) {
	      Tcl_AppendResult(interp, "\"", argv[1], "\" could not read data",
			       (char *) NULL);
	      return (TCL_ERROR);
	    }
	    Tcl_AppendResult(interp, "setvar ", buffer, (char *) NULL);
	    if (block_continueread(channel, 1, buffer, sizeof(buffer), 0) != 0) {
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
	Tcl_AppendResult(interp, "unknownvar ", buffer, " ", (char *) NULL);
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
	return (TCL_OK);
      }
      /* unknown field */
      else {
	Tcl_AppendResult(interp, "unknowntag ", title, " ", (char *) NULL);
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

  Tcl_AppendResult(interp, "unknown action \"", argv[2],
		   "\"", (char *)NULL);  
  return (TCL_ERROR);
}
