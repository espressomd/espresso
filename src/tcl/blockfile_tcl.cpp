/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
/** \file blockfile_tcl.cpp
    Implements the blockfile command for writing Tcl-formatted data files.
*/

#include <cstdio>
#include <ctype.h>
#include <cstdlib>
#include <cstring>

#include "parser.hpp"
#include "communication.hpp"

/** \name Error return codes
    All functions return non-negative values on success. The following
    codes are issued on errors. */
/*@{*/

/** end of file. This will only be returned by \ref block_startread if the
    file is read completely to allow convenient parsing through all blocks.
    All other functions will treat EOF as an error.
*/
#define RETURN_CODE_EOF   -1
/** I/O error on the file or unexpected EOF. */
#define RETURN_CODE_ERROR -2
/** file format wrong */
#define RETURN_CODE_FILE_FORMAT -3

/*@}*/

/** The maximal size allowed for block titles. */
#define MAXBLOCKTITLE 64
/** Used in the write commands as buffer size for sprintf.
    Possible incompatability. The current value allows for
    up to 20 digits precision */
#define DOUBLE_SPACE 32
/* the format used for formatted double IO */
#define DOUBLE_FORMAT "%.10e"
/** Used in the write commands as buffer size for sprintf.
    Possible incompatability. The current value this allows
    up to 64 bit integers */
#define INT_SPACE 32

static int readchar(Tcl_Channel f)
{
  char c;
  if (Tcl_Read(f, &c, 1) != 1) {
    return (Tcl_Eof(f)) ? 0 : -1;
  }
  return c;
}

static int writestring(Tcl_Channel f, const char *s)
{
  int l = strlen(s);
  return (Tcl_Write(f, s, l) == l);
}

static int findNonWs(Tcl_Channel f)
{
  int c, cont = 1;
  while (cont) {
    switch ((c = readchar(f))) {
    case 0:
      c = RETURN_CODE_EOF;
      break;
    case -1:
      c = RETURN_CODE_ERROR;
      break;
    default: ;
    }
    if (!isspace(c))
      break;
  }
  return c;
}

static int readString(Tcl_Channel f, char *buffer, int size)
{
  int i = 0;
  char c;
  while(i < size - 1) {
    c = (i == 0) ? findNonWs(f) : readchar(f);

    switch (c) {
    case RETURN_CODE_EOF:
    case RETURN_CODE_ERROR:
      buffer[i] = 0;
      return c;
    default: ;
    }
    if (c == '}') {
      buffer[i] = 0;
      return 1;
    }
    else if (isspace(c)) {
      buffer[i] = 0;
      return 0;
    }
    buffer[i++] = c;
  }
  buffer[i] = 0;
  return ((i == size - 1) ? 1 : 0);
}

/** read the title of a block. 
 @param f the file
 @param index where to store the index if round
 @return the number of open braces on success (1 if data follows or 0
         if no data was contained in the block) or
         \ref RETURN_CODE_EOF or \ref RETURN_CODE_ERROR
         or \ref RETURN_CODE_FILE_FORMAT if no "{" is found. In that case
	 index contains the offending character.
*/
static int block_startread(Tcl_Channel f, char index[MAXBLOCKTITLE])
{
  int c;

  index[0] = 0;

  /* find the block start "{" */
  switch (c = findNonWs(f)) {
  case '{':
    break;
  case RETURN_CODE_EOF:  return RETURN_CODE_EOF;
  case RETURN_CODE_ERROR:  return RETURN_CODE_ERROR;
  default:
    index[0] = c;
    index[1] = 0;
    return RETURN_CODE_FILE_FORMAT;
  }

  /* since a block started, we consider from now on eof an error */

  switch (readString(f, index, MAXBLOCKTITLE)) {
  case RETURN_CODE_EOF:
  case RETURN_CODE_ERROR:
    return RETURN_CODE_ERROR;
  case 0:
    return 1;
  case 1:
    return 0;
  }
  return RETURN_CODE_ERROR;
}

/** read a portion of the data of a block.
 @param f the file
 @param brace_count the number of open braces. If open_braces is one, reading
        terminates at the matching closing bracket, else open_braces-1 additional
	braces will be consumed.
 @param data where to store the contained data
 @param size the size of the buffer pointed to by "data".
 @param spacer if this character is read and no brace is open, i. e. open_braces=1,
        reading also terminates. A spacer of 0 disables this feature. This can be used
	to read in for example space separated lists:
	\verbatim {demoblock 1 2 3}\endverbatim
 @return returns the number of open braces or \ref RETURN_CODE_EOF or \ref RETURN_CODE_ERROR.
         The number of open_braces will be non-zero if the data space was exhausted.
*/
static int block_continueread(Tcl_Channel f, int brace_count, char *data, int size,
			      char spacer)
{
  char c;
  int i;

  if (!data || !size)
    return RETURN_CODE_ERROR;

  data[0] = 0;

  if (brace_count == 0)
    return 0;

  /* scan block data until brace_count = 0 or space eaten up */
  i = 0;
  while (i < size - 1) {
    if ((c = readchar(f)) <= 0) {
      data[i] = 0;
      return RETURN_CODE_ERROR;
    }
 
    if (c == '{') brace_count++;
    if (c == '}') {
      if (--brace_count == 0) {
	/* read complete block, strip trailing whitespaces */
	while (i > 1 && isspace(data[i - 1]))
	  i--;
	data[i] = 0;
	return 0;
      }
    }
    if (c == spacer && brace_count == 1) {
      data[i] = 0;
      return brace_count;
    }
    data[i++] = c;
  }
  data[i] = 0;

  return brace_count;
}

/** write a start tag to the file f.
    @param f the file
    @param index the tag to write (if index is a string longer than \ref MAXBLOCKTITLE,
    an error is returned.
    @return 0 on success or RETURN_CODE_ERROR
*/
static int block_writestart(Tcl_Channel f, char index[MAXBLOCKTITLE])
{
  if (strlen(index) >= MAXBLOCKTITLE)
    return RETURN_CODE_ERROR;
  if (!writestring(f, "{") ||
      !writestring(f, index) ||
      !writestring(f, " "))
    return RETURN_CODE_ERROR;
  return 0;
}

/** write a end tag to the file f.
    @param f the file
    @return 0 on success or RETURN_CODE_ERROR
*/
static int block_writeend(Tcl_Channel f)
{
  if (!writestring(f, "} "))
    return RETURN_CODE_ERROR;

  return 0;
}

int tclcommand_blockfile(ClientData data, Tcl_Interp *interp,
			 int argc, char *argv[])
{
  char title[MAXBLOCKTITLE];
  char buffer[1024], *name;
  int tcl_file_mode;
  Tcl_Channel channel;
  Tcl_CmdInfo cmdInfo;
  int openbrackets, len, exists;

  if (argc < 4) {
    Tcl_AppendResult(interp, "wrong # args:  should be \"",
		     argv[0], " <file> read|write <what> <param>? ...\"",
		     (char *) NULL);
    return (TCL_ERROR);
  }

  if ((channel = Tcl_GetChannel(interp, argv[1], &tcl_file_mode)) == NULL)
    return (TCL_ERROR);

  /* ------------------ the write commands ------------------------------------ */
  if (!strncmp(argv[2], "write", strlen(argv[2]))) {
    if (!(tcl_file_mode & TCL_WRITABLE)) {
      Tcl_AppendResult(interp, "\"", argv[1], "\" not writeable", (char *) NULL);
      return (TCL_ERROR);
    }

    /* ------------------ write start tag ------------------------------------ */
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
    /* ------------------ write end tag ------------------------------------ */
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
  }
  /* ------------------- the read commands ------------------------------------ */
  else if (!strncmp(argv[2], "read", strlen(argv[2]))) {
    if (!(tcl_file_mode & TCL_READABLE)) {
      Tcl_AppendResult(interp, "\"", argv[1], "\" not readable", (char *) NULL);
      return (TCL_ERROR);
    }
    /* ----------------- read start tag only ------------------------------- */
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
    /* -------------- read next block and auto execute or return data --*/
    else if (!strncmp(argv[3], "auto", strlen(argv[3]))) {
      if (argc != 4) {
	Tcl_AppendResult(interp, "\"", argv[0],
			 " read auto\" takes no parameters",
			 (char *) NULL);
	return (TCL_ERROR);
      }
      switch (openbrackets = block_startread(channel, title)) {
      case 0:	
      case 1:
	break;
      case RETURN_CODE_EOF:
	Tcl_AppendResult(interp, "eof", (char *)NULL);
	return (TCL_OK);
      case RETURN_CODE_FILE_FORMAT:
	Tcl_AppendResult(interp, "illstring ", title, (char *)NULL);
	return (TCL_OK);
      default:
	Tcl_AppendResult(interp, "\"", argv[1], "\" could not read block start",
			 (char *) NULL);
	return (TCL_ERROR);
      }
      /* -------------- read unknown field ---------------------------- */
      len = 21; /* blockfile_read_auto_ + \0 */
      len += strlen(title);
      name = (char*)malloc(len);
      strcpy(name, "blockfile_read_auto_");
      strcat(name, title);
      exists = Tcl_GetCommandInfo(interp, name, &cmdInfo);
      free(name);

      if (exists) {
	int err = cmdInfo.proc(cmdInfo.clientData, interp,
			       argc, argv);
	return gather_runtime_errors(interp, err);
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
    /* ------------------- read until end tag ------------------------------ */
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
    /* -------- (probably) read any tag implemented in tcl ----------------- */
    else {
      len = 21; /* blockfile_read_auto_ + \0 */
      len += strlen(argv[3]);
      name = (char*)malloc(len);
      strcpy(name, "blockfile_read_auto_");
      strcat(name, argv[3]);
      exists = Tcl_GetCommandInfo(interp, name, &cmdInfo);
      free(name);
      if (exists) {
	int err;
	if (!((openbrackets = block_startread(channel, title)) == 1) ||
	    strncmp(title,argv[3], strlen(argv[3]))) {
	  Tcl_AppendResult(interp, "\"", argv[1], "\" did not contain block\"", argv[3],"\" you indicated",
			   (char *) NULL);
	  return (TCL_ERROR);
	}

        err = cmdInfo.proc(cmdInfo.clientData, interp, argc, argv);
	return gather_runtime_errors(interp, err);
      }
    }
  }

  /* not a native action, try script support */
  len = 12; /* blockfile_ + _ + \0 */
  len += strlen(argv[2]) + strlen(argv[3]);
  name = (char*)malloc(len);
  strcpy(name, "blockfile_");
  strcat(name, argv[2]);
  strcat(name, "_");
  strcat(name, argv[3]);
  exists = Tcl_GetCommandInfo(interp, name, &cmdInfo);
  free(name);
  if (exists) {
    int err = cmdInfo.proc(cmdInfo.clientData, interp, argc, argv);
    return gather_runtime_errors(interp, err);
  }
  Tcl_AppendResult(interp, "unknown action \"", argv[2], " ", argv[3],
		   "\"", (char *)NULL);  
  return (TCL_ERROR);
}
