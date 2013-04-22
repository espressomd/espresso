#include "cells.h" //I can't go in extern C
#include "communication.h" //I can't go in extern C
#include "cuda_common.h" //I can't go in extern C

#include "grid.h"
#include "config.h"
#include "random.h"
#include "particle_data.h"
#include "interaction_data.h"


  static void mpi_get_particles_slave_lb();
  static void mpi_send_forces_slave_lb();

  /*************** REQ_GETPARTS ************/
  void mpi_get_particles_lb(LB_particle_gpu *host_data)
  {
    int n_part;
    int g, pnode;
    Cell *cell;
    int c;
    MPI_Status status;

    int i;  
    int *sizes;
    sizes = (int*) malloc(sizeof(int)*n_nodes);

    n_part = cells_get_n_particles();

    /* first collect number of particles on each node */
    MPI_Gather(&n_part, 1, MPI_INT, sizes, 1, MPI_INT, 0, comm_cart);

    /* just check if the number of particles is correct */
    if(this_node > 0){
      /* call slave functions to provide the slave datas */
      mpi_get_particles_slave_lb();
    }
    else {
      /* master: fetch particle informations into 'result' */
      g = 0;
      for (pnode = 0; pnode < n_nodes; pnode++) {
        if (sizes[pnode] > 0) {
          if (pnode == 0) {
            for (c = 0; c < local_cells.n; c++) {
              Particle *part;
              int npart;  
              int dummy[3] = {0,0,0};
              double pos[3];
              cell = local_cells.cell[c];
              part = cell->part;
              npart = cell->n;
              for (i=0;i<npart;i++) {
                memcpy(pos, part[i].r.p, 3*sizeof(double));
                fold_position(pos, dummy);
                host_data[i+g].p[0] = (float)pos[0];
                host_data[i+g].p[1] = (float)pos[1];
                host_data[i+g].p[2] = (float)pos[2];
                
                host_data[i+g].v[0] = (float)part[i].m.v[0];
                host_data[i+g].v[1] = (float)part[i].m.v[1];
                host_data[i+g].v[2] = (float)part[i].m.v[2];
                
  #ifdef LB_ELECTROHYDRODYNAMICS
                host_data[i+g].mu_E[0] = (float)part[i].p.mu_E[0];
                host_data[i+g].mu_E[1] = (float)part[i].p.mu_E[1];
                host_data[i+g].mu_E[2] = (float)part[i].p.mu_E[2];
  #endif

  #ifdef ELECTROSTATICS
                if (coulomb.method == COULOMB_P3M_GPU) {
                  host_data[i+g].q = (float)part[i].p.q;
                }
  #endif
              }  
              g += npart;
            }  
          }
          else {
            MPI_Recv(&host_data[g], sizes[pnode]*sizeof(LB_particle_gpu), MPI_BYTE, pnode, REQ_GETPARTS,
            comm_cart, &status);
            g += sizes[pnode];
          }
        }
      }
    }
    COMM_TRACE(fprintf(stderr, "%d: finished get\n", this_node));
    free(sizes);
  }

  static void mpi_get_particles_slave_lb(){
   
    int n_part;
    int g;
    LB_particle_gpu *host_data_sl;
    Cell *cell;
    int c, i;

    n_part = cells_get_n_particles();

    COMM_TRACE(fprintf(stderr, "%d: get_particles_slave, %d particles\n", this_node, n_part));

    if (n_part > 0) {
      /* get (unsorted) particle informations as an array of type 'particle' */
      /* then get the particle information */
      host_data_sl = (LB_particle_gpu*) malloc(n_part*sizeof(LB_particle_gpu));
      
      g = 0;
      for (c = 0; c < local_cells.n; c++) {
        Particle *part;
        int npart;
        int dummy[3] = {0,0,0};
        double pos[3];
        cell = local_cells.cell[c];
        part = cell->part;
        npart = cell->n;

        for (i=0;i<npart;i++) {
          memcpy(pos, part[i].r.p, 3*sizeof(double));
          fold_position(pos, dummy);  
      
          host_data_sl[i+g].p[0] = (float)pos[0];
          host_data_sl[i+g].p[1] = (float)pos[1];
          host_data_sl[i+g].p[2] = (float)pos[2];

          host_data_sl[i+g].v[0] = (float)part[i].m.v[0];
          host_data_sl[i+g].v[1] = (float)part[i].m.v[1];
          host_data_sl[i+g].v[2] = (float)part[i].m.v[2];
          
  #ifdef LB_ELECTROHYDRODYNAMICS
          host_data_sl[i+g].mu_E[0] = (float)part[i].p.mu_E[0];
          host_data_sl[i+g].mu_E[1] = (float)part[i].p.mu_E[1];
          host_data_sl[i+g].mu_E[2] = (float)part[i].p.mu_E[2];
  #endif

  #ifdef ELECTROSTATICS
          if (coulomb.method == COULOMB_P3M_GPU) {
            host_data_sl[i+g].q = (float)part[i].p.q;
          }
  #endif
        }
        g+=npart;
      }
      /* and send it back to the master node */
      MPI_Send(host_data_sl, n_part*sizeof(LB_particle_gpu), MPI_BYTE, 0, REQ_GETPARTS, comm_cart);
      free(host_data_sl);
    }  
  }

  void mpi_send_forces_lb(LB_particle_force_gpu *host_forces){
  
    int n_part;
    int g, pnode;
    Cell *cell;
    int c;
    int i;  
    int *sizes;
    sizes = (int *) malloc(sizeof(int)*n_nodes);
    n_part = cells_get_n_particles();
    /* first collect number of particles on each node */
    MPI_Gather(&n_part, 1, MPI_INT, sizes, 1, MPI_INT, 0, comm_cart);

    /* call slave functions to provide the slave datas */
    if(this_node > 0) {
      mpi_send_forces_slave_lb();
    }
    else{
    /* fetch particle informations into 'result' */
    g = 0;
      for (pnode = 0; pnode < n_nodes; pnode++) {
        if (sizes[pnode] > 0) {
          if (pnode == 0) {
            for (c = 0; c < local_cells.n; c++) {
              int npart;  
              cell = local_cells.cell[c];
              npart = cell->n;
              for (i=0;i<npart;i++) {
                cell->part[i].f.f[0] += (double)host_forces[i+g].f[0];
                cell->part[i].f.f[1] += (double)host_forces[i+g].f[1];
                cell->part[i].f.f[2] += (double)host_forces[i+g].f[2];
              }
        g += npart;
            }
          }
          else {
          /* and send it back to the slave node */
          MPI_Send(&host_forces[g], sizes[pnode]*sizeof(LB_particle_force_gpu), MPI_BYTE, pnode, REQ_GETPARTS, comm_cart);      
          g += sizes[pnode];
          }
        }
      }
    }
    COMM_TRACE(fprintf(stderr, "%d: finished send\n", this_node));

    free(sizes);
  }

  static void mpi_send_forces_slave_lb(){

    int n_part;
    LB_particle_force_gpu *host_forces_sl;
    Cell *cell;
    int c, i;
    MPI_Status status;

    n_part = cells_get_n_particles();

    COMM_TRACE(fprintf(stderr, "%d: send_particles_slave, %d particles\n", this_node, n_part));


    if (n_part > 0) {
      int g = 0;
      /* get (unsorted) particle informations as an array of type 'particle' */
      /* then get the particle information */
      host_forces_sl = (LB_particle_force_gpu *) malloc(n_part*sizeof(LB_particle_force_gpu));
      MPI_Recv(host_forces_sl, n_part*sizeof(LB_particle_force_gpu), MPI_BYTE, 0, REQ_GETPARTS,
      comm_cart, &status);
      for (c = 0; c < local_cells.n; c++) {
        int npart;  
        cell = local_cells.cell[c];
        npart = cell->n;
        for (i=0;i<npart;i++) {
          cell->part[i].f.f[0] += (double)host_forces_sl[i+g].f[0];
          cell->part[i].f.f[1] += (double)host_forces_sl[i+g].f[1];
          cell->part[i].f.f[2] += (double)host_forces_sl[i+g].f[2];
        }
        g += npart;
      }
      free(host_forces_sl);
    } 
  }
  /*@}*/