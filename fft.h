#ifndef FFT_H
#define FFT_H
/** \file fft.h
 *
 *  For more information about FFT usage, see \ref fft.c "fft.c".
 */

/************************************************
 * data types
 ************************************************/

/** Structure for performing a 1D FFT.  
 *
 *  This includes the information about the redistribution of the 3D
 *  FFT *grid before the actual FFT.  
*/
typedef struct {
  /** plan direction: 0 = Forward FFT, 1 = Backward FFT. */
  int dir;
  /** row direction of that FFT. */
  int row_dir;
  /** permutations from normal coordinate system. */
  int n_permute;
  /** number of 1D FFTs. */ 
  int n_ffts;
  /** plan for fft. */
  void *fft_plan;
  /** function for fft. */
  void (*fft_function)();

  /** size of local mesh before communication. */
  int old_mesh[3];
  /** size of local mesh after communication, also used for actual FFT. */
  int new_mesh[3];
  /** lower left point of local FFT mesh in global FFT mesh coordinates. */
  int start[3];
  /** size of new mesh (number of mesh points). */
  int new_size;

  /** number of nodes which have to communicate with each other. */ 
  int g_size;
  /** group of nodes which have to communicate with each other. */ 
  int *group;

  /** packing function for send blocks. */
  void (*pack_function)();
  /** Send block specification. 6 integers for each node: start[3], size[3]. */ 
  int *send_block;
  /** Send block communication sizes. */ 
  int *send_size;
  /** Recv block specification. 6 integers for each node: start[3], size[3]. */ 
  int *recv_block;
  /** Recv block communication sizes. */ 
  int *recv_size;
  /** size of send block elements. */
  int element;
} fft_forw_plan;

/** Additional information for backwards FFT.*/
typedef struct {
  /** plan direction. (e.g. fftw makro)*/
  int dir;
  /** plan for fft. */
  void *fft_plan;
  /** function for fft. */
  void (*fft_function)();

  /** packing function for send blocks. */
  void (*pack_function)(); 
} fft_back_plan;

extern fft_forw_plan fft_plan[4];

/************************************************
 * public functions
 ************************************************/


int fft_init(double *data, int *ca_mesh_dim, int *ca_mesh_margin);

void fft_perform_forw(double *data);

void fft_perform_back(double *data);

void fft_exit();

/** pack a block (size[3] starting at start[3]) of an input 3d-grid
 *  with dimension dim[3] into an output 3d-grid with dimension size[3].
 *
 *    The block with dimensions (size[0], size[1], size[2]) is stored
 *    in 'row-major-order' or 'C-order', that means the first index is
 *    changing slowest when running through the linear array. The
 *    element (i0 (slow), i1 (mid), i2 (fast)) has the linear index 
 *    li = i2 + size[2] * (i1 + (size[1]*i0)) 
 *
 *  @param in     pointer to input 3d-grid.
 *  @param out    pointer to outpu 3d-gird (block).
 *  @param start  start index of the block in the in-grid.
 *  @param size   size of the block (=dimension of the out-grid).
 *  @param dim    size of the in-grid.
 *  @param element size of a grid element (e.g. 1 for Real, 2 for Complex).
 */
void pack_block(double *in, double *out, int start[3], int size[3], 
		int dim[3], int element);

void pack_block_permute1(double *in, double *out, int start[3], int size[3], 
			 int dim[3], int element);

void pack_block_permute2(double *in, double *out, int start[3], int size[3], 
			 int dim[3],int element);
void unpack_block(double *in, double *out, int start[3], int size[3], 
		  int dim[3], int element);
#endif
