#ifndef BLOCKFILE_H
#define BLOCKFILE_H

#define MAXBLOCKTITLE 64

/* keep in sync with global.h ! */
#ifndef TYPE_INT
#define TYPE_INT    0
#endif
#ifndef TYPE_DOUBLE
#define TYPE_DOUBLE 1
#endif

#define RETURN_CODE_EOF   -1
#define RETURN_CODE_ERROR -2
#define RETURN_CODE_WDATA -3


#ifdef TCL_FILE_IO
#include <tcl.h>
#define FILETYPE Tcl_Channel
#else
#define FILETYPE FILE *
#endif

int block_startread(FILETYPE f, char index[MAXBLOCKTITLE]);
int block_continueread(FILETYPE f, int open_braces, char *data, int size,
		       char spacer);

int block_writestart(FILETYPE f, char index[MAXBLOCKTITLE]);
int block_writeend(FILETYPE f);
int block_write_double(FILETYPE f, double data);
int block_write_double_array(FILETYPE f, int dim[3], double *data);
int block_write_int(FILETYPE f, int data);
int block_write_int_array(FILETYPE f, int dim[3], int *data);
int block_read_any_data(FILETYPE f, int type, int dim[3], void *data);

#endif
