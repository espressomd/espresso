/** \file debug.h
 This file controls debug facilities. 

 <b>Responsible:</b>
 <a href="mailto:arnolda@mpip-mainz.mpg.de">Axel</a>

 The implementation is found in
 \ref debug.c "debug.c".
 This file contains a couple of commented out defines. If they are commented in,
 they activate debug output for various aspects:
 <ul>
 <li> \verbatim #define COMM_DEBUG \endverbatim activate printing of communicated actions.
 <li> \verbatim #define INTEG_DEBUG \endverbatim activate integration debug output.
 <li> \verbatim #define CELL_DEBUG \endverbatim activate cell code debug output.
 <li> \verbatim #define GHOST_DEBUG \endverbatim activate ghost code debug output.
 <li> \verbatim #define GRID_DEBUG \endverbatim activate grid debug output.
 <li> \verbatim #define FORCE_DEBUG \endverbatim activate force debug output.
 <li> \verbatim #define VERLET_DEBUG \endverbatim activate verlet debug output.
 <li> \verbatim #define PARTICLE_DEBUG \endverbatim activate particle data related debug output.
 <li> \verbatim #define P3M_DEBUG \endverbatim activate P3M related debug output.
 <li> \verbatim #define FFT_DEBUG \endverbatim activate FFT related debug output.
 <li> \verbatim #define MPI_CORE \endverbatim generate a core dump when exiting abnormally due
 to MPI errors.
 <li> \verbatim #define FORCE_CORE \endverbatim generate a core dump even on regular termination.
 <li> \verbatim #define MALLOC_DEBUG \endverbatim replaces malloc, realloc and free by version that
 log their actions. May help debugging issues and shows memory leaks. BUT YOU GET A LOT OF OUTPUT!!
 <li> \verbatim #define ADDITIONAL_CHECKS \endverbatim adds additional checks which indicate failure
 of the code.
 </ul>

 For every define there exists a macro that can be used to encapsulate short lines (like printf("...",...);)
 of code that should be executed iff the respective *_DEBUG macro is defined.
*/

/* #define COMM_DEBUG */
/* #define INTEG_DEBUG */
/* #define CELL_DEBUG */
/* #define GHOST_DEBUG */
/* #define GRID_DEBUG */
/* #define FORCE_DEBUG */
/* #define VERLET_DEBUG */
/* #define PARTICLE_DEBUG */
/* #define P3M_DEBUG */
/* #define FFT_DEBUG */
/* #define RANDOM_DEBUG */
/* #define FENE_DEBUG */

#define MPI_CORE
#define FORCE_CORE

#define ADDITIONAL_CHECKS

/* #define MALLOC_DEBUG */

#if defined FORCE_CORE || defined MPI_CORE
/** this functions kills the task with SIGSEGV */
extern void core();

/** by setting this variable to 1, a regular exit is
    indicated. In that case, no core dump is generated.
*/
extern int regular_exit;
#endif

#ifdef MALLOC_DEBUG
extern void *_debug_malloc(int size);
#define malloc(x) _debug_malloc(x)
extern void *_debug_realloc(void *p, int size);
#define realloc(x,y) _debug_realloc(x,y)
extern void _debug_free(void *p);
#define free(x) _debug_free(x)
#endif

#ifdef COMM_DEBUG
#define COMM_TRACE(cmd) { cmd; }
#else
/** Equals { cmd } iff COMM_DEBUG is set. */
#define COMM_TRACE(cmd)
#endif

#ifdef PARTICLE_DEBUG
#define PART_TRACE(cmd) { cmd; }
#else
/** Equals { cmd } iff PARTICLE_DEBUG is set. */
#define PART_TRACE(cmd)
#endif

#ifdef INTEG_DEBUG
#define INTEG_TRACE(cmd) { cmd; }
#else
/** Equals { cmd } iff INTEG_DEBUG is set. */
#define INTEG_TRACE(cmd)
#endif

#ifdef CELL_DEBUG
#define CELL_TRACE(cmd) { cmd;  }
#else
/** Equals { cmd } iff CELL_DEBUG is set. */
#define CELL_TRACE(cmd)
#endif

#ifdef GHOST_DEBUG
#define GHOST_TRACE(cmd) { cmd;  }
#else
/** Equals { cmd } iff GHOST_DEBUG is set. */
#define GHOST_TRACE(cmd)
#endif

#ifdef GRID_DEBUG
#define GRID_TRACE(cmd) { cmd;  }
#else
/** Equals { cmd } iff GRID_DEBUG is set. */
#define GRID_TRACE(cmd)
#endif

#ifdef FORCE_DEBUG
#define FORCE_TRACE(cmd) { cmd;  }
#else
/** Equals { cmd } iff FORCE_DEBUG is set. */
#define FORCE_TRACE(cmd)
#endif

#ifdef VERLET_DEBUG
#define VERLET_TRACE(cmd) { cmd;  }
#else
/** Equals { cmd } iff VERLET_DEBUG is set. */
#define VERLET_TRACE(cmd)
#endif

#ifdef P3M_DEBUG
#define P3M_TRACE(cmd) { cmd;  }
#else
/** Equals { cmd } iff P3M_DEBUG is set. */
#define P3M_TRACE(cmd)
#endif

#ifdef FFT_DEBUG
#define FFT_TRACE(cmd) { cmd; }
#else
/** Equals { cmd } iff FFT_DEBUG is set. */
#define FFT_TRACE(cmd)
#endif

#ifdef RANDOM_DEBUG
#define RANDOM_TRACE(cmd) { cmd; }
#else
/** Equals { cmd } iff FFT_DEBUG is set. */
#define RANDOM_TRACE(cmd)
#endif

#ifdef FENE_DEBUG
#define FENE_TRACE(cmd) { cmd; }
#else
#define FENE_TRACE(cmd)
#endif
