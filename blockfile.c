// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
/** \file blockfile.c
    Implementation of \ref blockfile.h "blockfile.h". If you compile this
    file with the option -DTCL_FILE_IO instead of the C-like FILE *
    operations Tcl channels are used for IO. This is used in Espresso
    itself, while the FILE * operations can be used in subsidiary code.
    The FILE * version is also contained in the \ref libEspresso.
*/
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include "blockfile.h"
#include "utils.h"

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

#ifndef BLOCKFILE_STDIO
static int readchar(Tcl_Channel f)
{
  char c;
  if (Tcl_Read(f, &c, 1) != 1) {
    return (Tcl_Eof(f)) ? 0 : -1;
  }
  return c;
}

static int writestring(Tcl_Channel f, char *s)
{
  int l = strlen(s);
  return (Tcl_Write(f, s, l) == l);
}

#else
static int readchar(FILE *f)
{
  int c;
  c = getc(f);
  if (c == EOF)
    return (feof(f)) ? 0 : -1;
  return c;
}

static int writestring(FILE *f, char *s)
{
  return (fputs(s, f) != EOF);
}
#endif

static int findNonWs(FILETYPE f)
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

static int readString(FILETYPE f, char *buffer, int size)
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

int block_startread(FILETYPE f, char index[MAXBLOCKTITLE])
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

int block_continueread(FILETYPE f, int brace_count, char *data, int size,
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

int block_writestart(FILETYPE f, char index[MAXBLOCKTITLE])
{
  if (strlen(index) >= MAXBLOCKTITLE)
    return RETURN_CODE_ERROR;
  if (!writestring(f, "{") ||
      !writestring(f, index) ||
      !writestring(f, " "))
    return RETURN_CODE_ERROR;
  return 0;
}

int block_writeend(FILETYPE f)
{
  if (!writestring(f, "} "))
    return RETURN_CODE_ERROR;

  return 0;
}

int block_write_data(FILETYPE f, int type, int dim, void *data)
{
  int error, i;
  char buffer[INT_SPACE + DOUBLE_SPACE];
  
  switch (type) {
  case TYPE_DOUBLE:
    if ((error = block_writestart(f, "_dval_")) != 0)
      return error;
    break;
  case TYPE_INT:
    if ((error = block_writestart(f, "_ival_")) != 0)
      return error;
    break;
  default: return RETURN_CODE_ERROR;
  }


  for (i = 0; i < dim; i++) {
    switch (type) {
    case TYPE_DOUBLE:
      sprintf(buffer, DOUBLE_FORMAT " ", ((double *)data)[i]);
      break;
    case TYPE_INT:
      sprintf(buffer, "%d ", ((int *)data)[i]);
      break;
    default: /* keep the compilers happy */;
    }
    if (!writestring(f, buffer))
      return RETURN_CODE_ERROR;
  }
  return block_writeend(f);
}


int block_read_data(FILETYPE f, int type, int dim, void *data)
{
  int readType;
  int err, i;
  char index[MAXBLOCKTITLE];
  char valbuffer[DOUBLE_SPACE + INT_SPACE];

  if ((err = block_startread(f, index)) < 0)
    return err;
  if (!strcmp(index, "_ival_"))
    readType = TYPE_INT;
  else if (!strcmp(index, "_dval_"))
    readType = TYPE_DOUBLE;
  else
    return RETURN_CODE_FILE_FORMAT;

  if (readType != type)
    return RETURN_CODE_FILE_FORMAT;

  for (i = 0; i < dim; i++) {
    if ((err = readString(f, valbuffer, sizeof(valbuffer))) < 0)
      return RETURN_CODE_ERROR;
    /* premature end of data input */
    if (err == 1 && i != dim - 1)
      return RETURN_CODE_FILE_FORMAT;

    switch (type) {
    case TYPE_DOUBLE:
      ((double *)data)[i] = atof(valbuffer);
      break;
    case TYPE_INT:
      ((int *)data)[i] = atol(valbuffer);
      break;
    default: /* keep the compilers happy */;
    }
  }

  /* } already consumed by readToWsOrBracket */
  if (err == 1)
    return 0;

  /* consume } and detect excessive data */
  return (findNonWs(f) == '}') ? 0 : RETURN_CODE_FILE_FORMAT;
}
