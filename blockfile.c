/** \file binary_file.c
    Implementation of \ref blockfile.h "blockfile.h". If you compile this file with the option
    -DTCL_FILE_IO instead of the C-like FILE * operations Tcl channels are used for IO. This is
    used in tcl_md itself, while the FILE * operations can be used in subsidiary code.
*/
#include <stdio.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include "utils.h"
#include "blockfile.h"

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

#ifdef TCL_FILE_IO
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

static char findNonWs(FILETYPE f)
{
  char c;
  for(;;) {
    switch ((c = readchar(f))) {
    case 0:
      return RETURN_CODE_EOF;
    case -1:
      return RETURN_CODE_ERROR;
    default: ;
    }
    if (!isspace(c))
      return c;
  }
  /* never reached, Mr. cxx */
  return 0;
}

static int readToWsOrBracket(FILETYPE f, char *buffer, int size)
{
  int i = 0;
  char c;
  while(i < size - 1) {
    switch ((c = readchar(f))) {
    case 0:
      buffer[i] = 0;
      return RETURN_CODE_EOF;
    case -1:
      buffer[i] = 0;
      return RETURN_CODE_ERROR;
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
  char c;
  int i;

  index[0] = 0;

  /* find the block start "{" */
  switch (c = findNonWs(f)) {
  case '{':
    break;
  case RETURN_CODE_EOF:  return RETURN_CODE_EOF;
  default: return RETURN_CODE_ERROR;
  }

  /* since a block started, we consider from now on eof an error */

  /* find start of block title */
  if ((c = findNonWs(f)) <= 0)
    return RETURN_CODE_ERROR;
  index[0] = c;
  /* read block title */
  i = 1;
  for(;;) {
    if ((c = readchar(f)) <= 0) {
      index[i] = 0;
      return RETURN_CODE_ERROR;
    }

    if (isspace(c))
      break;

    if (i < MAXBLOCKTITLE - 1)
      index[i++] = c;
  }
  index[i] = 0;

  return 1;
}

int block_continueread(FILETYPE f, int brace_count, char *data, int size,
		       char spacer)
{
  char c;
  int i;

  if (!data || !size)
    return RETURN_CODE_ERROR ;

  data[0] = 0;

  /* scan block data until brace_count = 0 or space eaten up */
  i = 0;
  while (i < size - 1) {
    if (i == 0) {
      if ((c = findNonWs(f)) <= 0)
	return RETURN_CODE_ERROR;
    }
    else
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
    if (c == spacer) {
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

int block_write_double(FILETYPE f, double data)
{
  int error;
  char buffer[DOUBLE_SPACE];
  
  if ((error = block_writestart(f, "_dval_")) != 0)
    return error;

  sprintf(buffer, DOUBLE_FORMAT, data);
  if (!writestring(f, buffer))
    return RETURN_CODE_ERROR;

  return block_writeend(f);
}

int block_write_int(FILETYPE f, int data)
{
  int error;
  char buffer[INT_SPACE];

  if ((error = block_writestart(f, "_ival_")) != 0)
    return error;

  sprintf(buffer, "%d", data);
  if (!writestring(f, buffer))
    return RETURN_CODE_ERROR;

  return block_writeend(f);
}

int block_write_double_array(FILETYPE f, int dim[3], double *data)
{
  int error;
  int sx, sy, sz;
  char buffer[3 + INT_SPACE*3 + DOUBLE_SPACE];

  if ((error = block_writestart(f, "_darray_")) != 0)
    return error;

  sprintf(buffer, "%d %d %d ", dim[0], dim[1], dim[2]);
  if (!writestring(f, buffer))
    return RETURN_CODE_ERROR;
  for (sx = 0; sx < dim[0]; sx++)
    for (sy = 0; sy < dim[1]; sy++)
      for (sz = 0; sz < dim[2]; sz++) {
	sprintf(buffer, DOUBLE_FORMAT " ",
		data[get_linear_index(sx, sy, sz, dim)]);
	if (!writestring(f, buffer))
	  return RETURN_CODE_ERROR;
      }
  return block_writeend(f);
}

int block_write_int_array(FILETYPE f, int dim[3], int *data)
{
  int error;
  int sx, sy, sz;
  char buffer[3 + INT_SPACE*3];

  if ((error = block_writestart(f, "_iarray_")) != 0)
    return error;

  sprintf(buffer, "%d %d %d ", dim[0], dim[1], dim[2]);
  if (!writestring(f, buffer))
    return RETURN_CODE_ERROR;
  for (sx = 0; sx < dim[0]; sx++)
    for (sy = 0; sy < dim[1]; sy++)
      for (sz = 0; sz < dim[2]; sz++) {
	sprintf(buffer, "%d ",
		data[get_linear_index(sx, sy, sz, dim)]);
	if (!writestring(f, buffer))
	  return RETURN_CODE_ERROR;
      }
  return block_writeend(f);
}

int block_read_any_data(FILETYPE f, int type, int dim[3], void *data)
{
  int readType, array;
  int err, i, sx, sy, sz;
  char index[MAXBLOCKTITLE];
  char valbuffer[DOUBLE_SPACE + INT_SPACE];

  if ((err = block_startread(f, index)) < 0)
    return err;
  if (!strcmp(index, "_ival_")) {
    readType = TYPE_INT;
    array = 0;
  }
  else if (!strcmp(index, "_dval_")) {
    readType = TYPE_DOUBLE;
    array = 0;
  }
  else if (!strcmp(index, "_iarray_")) {
    readType = TYPE_INT;
    array = 1;
  }
  else if (!strcmp(index, "_darray_")) {
    readType = TYPE_DOUBLE;
    array = 1;
  }
  else
    readType = -1;
  if (readType != type)
    return RETURN_CODE_WDATA;

  if (array) {
    for (i = 0; i < 3; i++) {
      if (readToWsOrBracket(f, valbuffer, sizeof(valbuffer)) != 0)
	return RETURN_CODE_ERROR;
      if (dim[i] != atol(valbuffer))
	return RETURN_CODE_WDATA;
    }
  }
  else {
    for (i = 0; i < 3; i++) {
      if (dim[i] != 1)
	return RETURN_CODE_WDATA;
    }
  }

  for (sx = 0; sx < dim[0]; sx++)
    for (sy = 0; sy < dim[1]; sy++)
      for (sz = 0; sz < dim[2]; sz++) {
	if ((err = readToWsOrBracket(f, valbuffer, sizeof(valbuffer))) < 0)
	  return RETURN_CODE_ERROR;

	switch (type) {
	case TYPE_DOUBLE:
	  ((double *)data)[get_linear_index(sx, sy, sz, dim)] = atof(valbuffer);
	  break;
	case TYPE_INT:
	  ((int *)data)[get_linear_index(sx, sy, sz, dim)] = atol(valbuffer);
	  break;
	  break;
	default: ;
	}
      }

  /* } already consumed by readToWsOrBracket */
  if (err == 1)
    return 0;

  /* consume } */
  return (findNonWs(f) == '}') ? 0 : RETURN_CODE_ERROR;
}
