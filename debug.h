/* #define COMM_DEBUG */
#define INTEG_DEBUG
#define CELL_DEBUG
#define GHOST_DEBUG
#define GRID_DEBUG
#define FORCE_DEBUG
#define VERLET_DEBUG
#define PARTICLE_DEBUG

/* #define FORCE_CORE */

/* #define MALLOC_DEBUG */

#ifdef MALLOC_DEBUG
extern inline void *___malloc(int size) {
  void *res = malloc(size);
  fprintf(stderr, "allocated %p+%d\n", res, size);
  return res;
}

extern inline void *___realloc(void *p, int size) {
  void *res = realloc(p, size);
  fprintf(stderr, "reallocated %p -> %p+%d\n", p, res, size);
  return res;
}

extern inline void ___free(void *p) {
  free(p);
  fprintf(stderr, "freed %p\n", p);
}

#define malloc(x) ___malloc(x)
#define realloc(x,y) ___realloc(x,y)
#define free(x) ___free(x)
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
