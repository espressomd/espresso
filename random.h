#ifndef RANDOM_H
#define RANDOM_H

/* provides a random generator, either a long integer or a uniform
   double between 0 and 1 is drawn. The random generator ran1
   from numerical recipies was choosen with some modifications.
   Be sure to run init_random before use!!! */

/*----------------------------------------------------------*/

extern long   l_random(void);
extern double random(void);
extern void   init_random(void);

/*----------------------------------------------------------*/

#endif






