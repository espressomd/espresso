/** \file debug.c
    Implements the malloc replacements as described in \ref debug.h "debug.h". */

/* do NOT include debug.h !!!!!!! */
#include "stdio.h"
#include "stdlib.h"
#include "string.h"

/** memory from malloc will initially be filled with this value. */
static unsigned char fill_code = 0xaa;

/** Replacement for malloc that logs allocation size and returned address. */
extern void *_debug_malloc(int size);
/** Replacement for realloc that logs allocation size, original and returned address. */
extern void *_debug_realloc(void *p, int size);
/** Replacement for free that logs the freed address. */
extern void _debug_free(void *p);

void *_debug_malloc(int size)
{
  void *res = malloc(size);
  memset(res, fill_code, size);
  fprintf(stderr, "allocated %p+%d\n", res, size);
  return res;
}

void *_debug_realloc(void *p, int size)
{
  void *res = realloc(p, size);
  fprintf(stderr, "reallocated %p -> %p+%d\n", p, res, size);
  return res;
}

void _debug_free(void *p)
{
  free(p);
  fprintf(stderr, "freed %p\n", p);
}

