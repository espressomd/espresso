/** \file debug.c
    Implements the malloc replacements as described in \ref debug.h "debug.h". */

#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "communication.h"
#include "debug.h"

#if defined FORCE_CORE || defined MPI_CORE
int regular_exit = 0;
#else
int regular_exit = 1;
#endif
static int core_done = 0;

void core()
{
  if (!core_done && !regular_exit) {
    core_done = 1;
    fprintf(stderr, "%d: forcing core dump on irregular exit\n", this_node);
    kill(getpid(), SIGSEGV);
  }
}
