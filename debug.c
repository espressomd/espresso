/** \file debug.c
    Implements the malloc replacements as described in \ref debug.h "debug.h". */

/* do NOT include debug.h !!!!!!! */
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "communication.h"

/** memory from malloc will initially be filled with this value. */
static unsigned char fill_code = 0xaa;

/** Replacement for malloc that logs allocation size and returned address. */
extern void *_debug_malloc(int size);
/** Replacement for realloc that logs allocation size, original and returned address. */
extern void *_debug_realloc(void *p, int size);
/** Replacement for free that logs the freed address. */
extern void _debug_free(void *p);

int regular_exit = 0;
static int core_done = 0;

void core()
{
  if (!core_done && !regular_exit) {
    core_done = 1;
    fprintf(stderr, "%d: forcing core dump on irregular exit\n", this_node);
    kill(getpid(), SIGSEGV);
  }
}

void *_debug_malloc(int size)
{
  void *res = malloc(size);
  memset(res, fill_code, size);
  fprintf(stderr, "%d: allocated %p+%d\n", this_node, res, size);
  return res;
}

void *_debug_realloc(void *p, int size)
{
  void *res = realloc(p, size);
  fprintf(stderr, "%d: reallocated %p -> %p+%d\n", this_node, p, res, size);
  return res;
}

void _debug_free(void *p)
{
  free(p);
  fprintf(stderr, "%d: freed %p\n", this_node, p);
}

