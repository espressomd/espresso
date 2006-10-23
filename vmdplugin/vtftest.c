/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2003 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile$
 *      $Author$       $Locker$             $State$
 *      $Revision$       $Date$
 *
 ***************************************************************************/

/*
 * A general main for testing plugins.  
 * Compile using: gcc main.c plugin.c -I../../include -o plugintest
 * Replace plugin.c with the plugin file you want to test.
 * Usage: plugintest <filetype> <file> 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "molfile_plugin.h"

static molfile_plugin_t *plugin = 0;
static const char *filetype = NULL;

static int register_cb(void *v, vmdplugin_t *p) {
  if (!strcmp(p->type, MOLFILE_PLUGIN_TYPE) && !strcmp(p->name, filetype))
    plugin = (molfile_plugin_t *)p;

  return VMDPLUGIN_SUCCESS;
}

int main(int argc, char *argv[]) {
  const char *filename;
  int rc, natoms;
  molfile_timestep_t timestep;
  void *handle;

  if (argc < 3) {
    fprintf(stderr, "Usage: %s <filetype> <filename>\n", argv[0]);
    return 1;
  }
  filetype = argv[1];
  filename = argv[2];
 
  vmdplugin_init();
  vmdplugin_register(NULL, register_cb);
  if (!plugin) {
    fprintf(stderr, "No plugin for filetype %s was linked in!\n", filetype);
    return 1;
  }
  
  if (!plugin->open_file_read) {
    fprintf(stdout, "FAILED: No open_file_read found.\n");
    return 1;
  } 
  handle = plugin->open_file_read(filename, filetype, &natoms);
  if (!handle) {
    fprintf(stderr, "FAILED: open_file_read returned NULL\n");
    return 1;
  }
  printf("Opened file %s; found %d atoms\n", 
    filename, natoms);
  if (plugin->read_structure) {
    int optflags;
    molfile_atom_t *atoms;
    atoms = (molfile_atom_t *)malloc(natoms * sizeof(molfile_atom_t));
    rc = plugin->read_structure(handle, &optflags, atoms);
    if (rc) {
      fprintf(stderr, "FAILED: read_structure returned %d\n", rc);
      plugin->close_file_read(handle);
      return 1;
    } else {
      printf("Succesfully read atom structure information.\n");
    }
    if (plugin->read_bonds) {
      int nbonds, *from, *to;
      if ((rc = plugin->read_bonds(handle, &nbonds, &from, &to
#if vmdplugin_ABIVERSION >= 9
				   ,NULL
#endif
				   ))) {
        fprintf(stderr, "FAILED: read_bonds returned %d\n", rc);
      } else {
        printf("read_bonds read %d bonds\n", nbonds);
        free(from);
        free(to);
      }
    } else {
      printf("File contains no bond information\n");
    }
  }
  if (plugin->read_next_timestep) {
    int nsteps = 0;
    timestep.coords = (float *)malloc(3*natoms*sizeof(float));
    while (!(rc = plugin->read_next_timestep(handle, natoms, &timestep))) 
      nsteps++;
    free(timestep.coords);
    if (rc != -1) {
      fprintf(stderr, "FAILED: read_next_timestep returned %d\n", rc);
    } else {
      printf("successfully read %d timesteps\n", nsteps);
    }
  }
  plugin->close_file_read(handle); 

  vmdplugin_fini();
  printf("Tests finished.\n");
  return 0;
}

