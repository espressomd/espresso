/** \file config.h 
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:arnolda@mpip-mainz.mpg.de">Axel</a>
*/


/** if defined, the code will be slower, but with the \ref periodic
    array you can choose which coordinates are bound to p.b.c and
    which are not. If not defined, all coordinates are bound to
    p.b.c. 

    Has effect on: \ref per_callback, \ref find_node, \ref fields, 
    \ref cells_init and \ref sort_particles_into_cells.
*/
#define PARTIAL_PERIODIC
