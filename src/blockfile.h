/*
  Copyright (C) 2010 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany
  
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
/** \file blockfile.h
    This contains routines to access files formatted as described in \ref blockformat.
    The file can be either a Tcl channel or a FILE *, which is called FILETYPE here.
    See \ref blockfile.c "blockfile.c" for more information.
*/
 
#ifndef BLOCKFILE_H
#define BLOCKFILE_H
 
/** The maximal size allowed for block titles. */
#define MAXBLOCKTITLE 64

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


/* keep in sync with global.h! These are only here for
   compilation without Tcl support. */
#ifndef TYPE_INT
/** \internal */
#define TYPE_INT    0
#endif
#ifndef TYPE_DOUBLE
/** \internal */
#define TYPE_DOUBLE 1
#endif

#ifdef BLOCKFILE_STDIO
#include <stdio.h>
/** \internal */
#define FILETYPE FILE *
#else
#include <tcl.h>
/** \internal */
#define FILETYPE Tcl_Channel
#endif

/** read the title of a block. 
 @param f the file
 @param index where to store the index if round
 @return the number of open braces on success (1 if data follows or 0
         if no data was contained in the block) or
         \ref RETURN_CODE_EOF or \ref RETURN_CODE_ERROR
         or \ref RETURN_CODE_FILE_FORMAT if no "{" is found. In that case
	 index contains the offending character.
*/
int block_startread(FILETYPE f, char index[MAXBLOCKTITLE]);
/** read a portion of the data of a block.
 @param f the file
 @param open_braces the number of open braces. If open_braces is one, reading
        terminates at the matching closing bracket, else open_braces-1 additional
	braces will be consumed.
 @param data where to store the contained data
 @param size the size of the buffer pointed to by \ref data.
 @param spacer if this character is read and no brace is open, i. e. open_braces=1,
        reading also terminates. A spacer of 0 disables this feature. This can be used
	to read in for example space separated lists:
	\verbatim {demoblock 1 2 3}\endverbatim
 @return returns the number of open braces or \ref RETURN_CODE_EOF or \ref RETURN_CODE_ERROR.
         The number of open_braces will be non-zero if the data space was exhausted.
*/
int block_continueread(FILETYPE f, int open_braces, char *data, int size,
		       char spacer);
/** write a start tag to the file f.
    @param f the file
    @param index the tag to write (if index is a string longer than \ref MAXBLOCKTITLE,
    an error is returned.
    @return 0 on success or RETURN_CODE_ERROR
*/
int block_writestart(FILETYPE f, char index[MAXBLOCKTITLE]);
/** write a end tag to the file f.
    @param f the file
    @return 0 on success or RETURN_CODE_ERROR
*/
int block_writeend(FILETYPE f);
/** write a data field to the file f. \ref type, \ref dim and \ref data are fully compatible
    to the fields of \ref Datafield.
    @param f the file
    @param type the field type, either \ref TYPE_INT or \ref TYPE_DOUBLE
    @param dim the fields number of elements
    @param data the buffer containing the data
    @return 0 on success or RETURN_CODE_ERROR
*/
int block_write_data(FILETYPE f, int type, int dim, void *data);
/** read a data field from file f. \ref type, \ref dim and \ref data are fully compatible
    to the fields of \ref Datafield. The type information is also stored in the file and
    counterchecked.
    @param f the file
    @param type the field type, either \ref TYPE_INT or \ref TYPE_DOUBLE
    @param dim the fields number of elements
    @param data the buffer where to store the data
    @return 0 on success or RETURN_CODE_ERROR or RETURN_CODE_FILE_FORMAT if the size of field
            or the type does not match the given parameters.
*/
int block_read_data(FILETYPE f, int type, int dim, void *data);

#endif
