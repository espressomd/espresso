#ifndef BINARYFILE_H
#define BINARYFILE_H
/**************************************************************
 * binary file format
 **************************************************************/

#define MDMAGIC "MD01"

/** first a header. Magic is used to verify data has correct format */
struct MDHeader {
  char magic[4]; /* must be MDMAGIC as above */
  int  n_rows;   /* number of rows contained */
};

/** next n_rows bytes describe the row contents */
#define POSX  0
#define POSY  1
#define POSZ  2
#define VX    3
#define VY    4
#define VZ    5
#define FX    6
#define FY    7
#define FZ    8
#define Q     9
#define TYPE  10

/** following is the particle data, always particle number (int) and n_rows ints/doubles
 *  The data types are: double for POS*,V*,F*,Q and int for type.
 *  DATA IS WRITTEN IN MACHINE FORMAT AND MAY NOT BE PORTABLE TO DIFFERENT ARCHITECTURES */

/** The end of the data is signalled by a single int -1 instead of a particle number */
#endif
