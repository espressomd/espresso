/* do NOT include debug.h !!!!!!! */
#include "stdio.h"
#include "stdlib.h"

static char fill_code = 0xaa;

extern void *_debug_malloc(int size);
extern void *_debug_realloc(void *p, int size);
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

