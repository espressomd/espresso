/* #define COMM_DEBUG */
/* #define INTEG_DEBUG */
/* #define CELL_DEBUG */
/* #define GHOST_DEBUG */
/* #define GRID_DEBUG */
/* #define FORCE_DEBUG */
/* #define VERLET_DEBUG */
/* #define PARTICLE_DEBUG */

/* #define FORCE_CORE */

/* #define MALLOC_DEBUG */

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
#define COMM_TRACE(cmd)
#endif

#ifdef PARTICLE_DEBUG
#define PART_TRACE(cmd) { cmd; }
#else
#define PART_TRACE(cmd)
#endif

#ifdef INTEG_DEBUG
#define INTEG_TRACE(cmd) { if (this_node < 2) { cmd; } }
#else
#define INTEG_TRACE(cmd)
#endif

#ifdef CELL_DEBUG
#define CELL_TRACE(cmd) { if (this_node < 2) { cmd; } }
#else
#define CELL_TRACE(cmd)
#endif

#ifdef GHOST_DEBUG
#define GHOST_TRACE(cmd) { if (this_node < 2) { cmd; } }
#else
#define GHOST_TRACE(cmd)
#endif

#ifdef GRID_DEBUG
#define GRID_TRACE(cmd) { if (this_node < 2) { cmd; } }
#else
#define GRID_TRACE(cmd)
#endif

#ifdef FORCE_DEBUG
#define FORCE_TRACE(cmd) { if (this_node < 9) { cmd; } }
#else
#define FORCE_TRACE(cmd)
#endif

#ifdef VERLET_DEBUG
#define VERLET_TRACE(cmd) { if (this_node < 2) { cmd; } }
#else
#define VERLET_TRACE(cmd)
#endif
