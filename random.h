#ifndef RANDOM_H
#define RANDOM_H

/** \file random.h a random generator
    <b>Responsible:</b>
    <a href="mailto:muehlbac@mpip-mainz.mpg.de">Frank</a>
    \todo Put any debug output in a DEBUG statement!!!
 */

/*----------------------------------------------------------*/

extern long   l_random(void);
extern int    i_random(int maxint);
extern double random(void);
extern void   init_random(void);

/*----------------------------------------------------------*/

#endif






