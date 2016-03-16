/*
  Copyright (C) 2010,2011,2012,2016 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
/** \file lees_edwards_domain_decomposition.cpp
 *
 *  This file contains modifications to the domain decomposition functions 
 *  which are specific to Lees-Edwards boundary conditions.
 *  See also \ref lees_edwards_domain_decomposition.hpp ,
 *  \ref lees_edwards_comms_manager.hpp and  \ref domain_decomposition.hpp
 */



#include "domain_decomposition.hpp"
#include "lees_edwards_domain_decomposition.hpp"
#include "lees_edwards_comms_manager.hpp"
#include "errorhandling.hpp"
#include "initialize.hpp"
#include "lees_edwards.hpp"
#include "grid.hpp"


#ifdef LEES_EDWARDS

#define LE_TRACE(x) 
#ifdef  LE_DEBUG
FILE *comms_log = NULL;
int   comms_count = 0;
#endif

/** Init cell interactions for the Lees-Edwards cell system.
 * initializes the interacting neighbor cell list of a cell. The
 * created list of interacting neighbor cells is used by the verlet
 * algorithm (see verlet.cpp) to build the verlet lists.
 */
void le_dd_init_cell_interactions()
{
  int m,n,o,p,q,r,ind1,ind2,c_cnt=0,n_cnt=0;
  int extra_cells = 0;

  /* initialize cell neighbor structures */
  dd.cell_inter = (IA_Neighbor_List *) Utils::realloc(dd.cell_inter,local_cells.n*sizeof(IA_Neighbor_List));
  for(m=0; m<local_cells.n; m++) { 
    dd.cell_inter[m].nList = NULL; 
    dd.cell_inter[m].n_neighbors=0; 
  }

  /* loop over non-ghost cells */
  for(o=1; o<=dd.cell_grid[2]; o++) {
    for(n=1; n<=dd.cell_grid[1]; n++) {
      for(m=1; m<=dd.cell_grid[0]; m++) {

    /* plenty for most cases */
    dd.cell_inter[c_cnt].nList = (IA_Neighbor *) Utils::realloc(dd.cell_inter[c_cnt].nList, 14*sizeof(IA_Neighbor));
    
    n_cnt=0;
    ind1 = get_linear_index(m,n,o,dd.ghost_cell_grid);

    /* loop all 'conventional' neighbor cells */
    for(p=o-1; p<=o+1; p++) {      /*z-loop*/
      for(q=n-1; q<=n+1; q++) {    /*y-loop*/
        for(r=m-1; r<=m+2; r++) {  /*x-loop*/

            /* Extra neighbours in x only for some cases */
            if(    (q == 0                 && node_pos[1] == 0)
                || (q == dd.cell_grid[1]+1 && node_pos[1] == node_grid[1]-1) ){
                extra_cells++;
                dd.cell_inter[c_cnt].nList = (IA_Neighbor *) Utils::realloc(dd.cell_inter[c_cnt].nList, (extra_cells+14)*sizeof(IA_Neighbor));
            }else{
                if( r == m + 2 )
                    continue;
            }
                
            ind2 = get_linear_index(r,q,p,dd.ghost_cell_grid);
            
            if(ind2 >= ind1) {
                dd.cell_inter[c_cnt].nList[n_cnt].cell_ind = ind2;
                dd.cell_inter[c_cnt].nList[n_cnt].pList    = &cells[ind2];
                init_pairList(&dd.cell_inter[c_cnt].nList[n_cnt].vList);
#ifdef LE_DEBUG
    dd.cell_inter[c_cnt].nList[n_cnt].my_pos[0] = my_left[0] + r * dd.cell_size[0];
    dd.cell_inter[c_cnt].nList[n_cnt].my_pos[1] = my_left[1] + q * dd.cell_size[1];
    dd.cell_inter[c_cnt].nList[n_cnt].my_pos[2] = my_left[2] + p * dd.cell_size[2];
#endif
              n_cnt++;
            }
          }
        }
      }
        dd.cell_inter[c_cnt].n_neighbors = n_cnt;
        c_cnt++;
        }
      }
   }

#ifdef LE_DEBUG
  FILE *cells_fp;
  char cLogName[64];
  int  c,nn,this_n;
  double myPos[3];
  sprintf(cLogName, "cells_map%i.dat", this_node);
  cells_fp = fopen(cLogName,"w");

  /* print out line segments showing the vector from each cell to each neighbour cell*/
  for(c=0;c<c_cnt;c++){
     myPos[0] = my_left[0] + dd.cell_size[0] * ( 1 + c % dd.cell_grid[0] );  
     myPos[1] = my_left[1] + dd.cell_size[1] * ( 1 + (c / dd.cell_grid[0]) % dd.cell_grid[1]);  
     myPos[2] = my_left[2] + dd.cell_size[2] * ( 1 + (c / (dd.cell_grid[0] * dd.cell_grid[1])));  

     for(nn=0;nn<dd.cell_inter[c].n_neighbors;nn++){
        
        this_n = dd.cell_inter[c].nList[nn].cell_ind;
        fprintf(cells_fp,"%i %i %i %f %f %f %f %f %f\n",c,nn,this_n,
            myPos[0], myPos[1], myPos[2], 
            dd.cell_inter[c].nList[nn].my_pos[0], 
            dd.cell_inter[c].nList[nn].my_pos[1], 
            dd.cell_inter[c].nList[nn].my_pos[2]);
     }
  }  
  fclose(cells_fp);
#endif

}


/** update the 'shift' member of those GhostCommunicators, which use
    that value to speed up the folding process of its ghost members
    (see \ref dd_prepare_comm for the original), i.e. all which have
    GHOSTTRANS_POSSHFTD or'd into 'data_parts' upon execution of \ref
    dd_prepare_comm. */
void le_dd_update_communicators_w_boxl(le_dd_comms_manager *mgr)
{
  int cnt=0;
  int neighbor_index, neighbor_rank, dir, i;
  
  for( i = 0; i < my_neighbor_count; i++ ){
      
     neighbor_index = mgr->comms_order[i];

     /* find some geometry about this comm, comms order is y-x-z */
     dir = 1;
     if( neighbor_index / 2 == 0 )      { dir = 0;}
     else if ( neighbor_index / 2 == 2 ){ dir = 2;}

     neighbor_rank = node_neighbors[neighbor_index];

     if( neighbor_rank == this_node ) { 
        /* prepare folding of ghost positions */
        if(node_neighbor_wrap[neighbor_index] != 0) {
            cell_structure.exchange_ghosts_comm.comm[cnt].shift[dir]  = node_neighbor_wrap[neighbor_index]*box_l[dir];
            cell_structure.update_ghost_pos_comm.comm[cnt].shift[dir] = node_neighbor_wrap[neighbor_index]*box_l[dir];
            if( dir == 1 ){
                cell_structure.exchange_ghosts_comm.comm[cnt].shift[0]   = node_neighbor_wrap[neighbor_index]*lees_edwards_offset;
                cell_structure.update_ghost_pos_comm.comm[cnt].shift[0]  = node_neighbor_wrap[neighbor_index]*lees_edwards_offset;
            }
        }
        cnt++;
     }
     else {
    /* send/recv loop: easiest to just set the shift for both of them, only sends will actually be shifted. */
    for(int send_rec=0; send_rec<2; send_rec++) {

        /* prepare folding of ghost positions */
        if(node_neighbor_wrap[neighbor_index] != 0) {
            cell_structure.exchange_ghosts_comm.comm[cnt].shift[dir]  = node_neighbor_wrap[neighbor_index]*box_l[dir];
            cell_structure.update_ghost_pos_comm.comm[cnt].shift[dir] = node_neighbor_wrap[neighbor_index]*box_l[dir];
            if( dir == 1 ){
                cell_structure.exchange_ghosts_comm.comm[cnt].shift[0]   = node_neighbor_wrap[neighbor_index]*lees_edwards_offset;
                cell_structure.update_ghost_pos_comm.comm[cnt].shift[0]  = node_neighbor_wrap[neighbor_index]*lees_edwards_offset;
            }
        }
        cnt++;
        
      }
    }
  }

  
}

#ifdef LE_DEBUG
void printAllComms(FILE *ofile, GhostCommunicator *comm, int backwards){

  int i,j,cnt,ccnt, num;
  int neighbor_index;
  char ***comGrid;
  int  l1, l2;
  
  l1          = dd.cell_grid[0]+3;
  l2          = dd.cell_grid[1]+2;
  comGrid     = new char**[l1];
  
  for( i = 0; i < l1;i++){
    comGrid[i] = new char*[l2];
    for( j = 0; j < l2;j++){
        comGrid[i][j]=new char[8];
        sprintf(&comGrid[i][j][0],".......");
        comGrid[i][j][3]='\0';
        comGrid[i][j][7]='\0';
    }
  }
  
  /* prepare communicator... need 2 comms for node-node, 1 comm for intranode exchange */
  num = my_neighbor_count * 2;
  for(i = 0; i < my_neighbor_count; i++) {
        if( node_neighbors[i] == this_node ) num--;
  }

  fprintf(ofile,"printing comms for node %i dimensions %i-%i-%i\n",this_node,dd.cell_grid[0],dd.cell_grid[1],dd.cell_grid[2]);
  fprintf(ofile,"list of partners and sizes:\n");
  for(ccnt = 0; ccnt < num; ccnt++){
    if( backwards ) cnt = num - 1 - ccnt;
    else            cnt = ccnt;
    fprintf(ofile,"  %i %i\n",comm->comm[cnt].node,comm->comm[cnt].n_part_lists);
  }
  fprintf(ofile,"\n");
  
  /*loop over comms */
  for(ccnt = 0; ccnt < num; ccnt++) {
      
    int  neighbor_coords[3];
    
    if( backwards == LE_COMM_BACKWARDS ) cnt = num - 1 - ccnt;
    else                                 cnt = ccnt;

    neighbor_index   = comm->comm[cnt].node;
    map_node_array(neighbor_index, neighbor_coords);

    if( ( comm->comm[cnt].type == GHOST_SEND
       || comm->comm[cnt].type == GHOST_RECV )
       && comm->comm[cnt].node == this_node ){
        fprintf(ofile,"serious comms error\n");
        fprintf(stderr,"serious comms error\n");
        errexit();
    }
        
    /* clear the representation of this comm */
    for( i = 0; i < l1;i++){
        for( j = 0; j < l2;j++){
            sprintf(&comGrid[i][j][0],".......");
            comGrid[i][j][3]='\0';
            comGrid[i][j][7]='\0';
        }
    }
    
    if( comm->comm[cnt].type == GHOST_LOCL ){
        for(  i = 0; i < comm->comm[cnt].n_part_lists/2; i++){
            if( comm->comm[cnt].part_lists[i]->myIndex[2] > 1
             && comm->comm[cnt].part_lists[i]->myIndex[2] < 10){
                //fprintf(stderr,"%i: local send %i %i\n",this_node,comm->comm[cnt].part_lists[i]->myIndex[0],comm->comm[cnt].part_lists[i]->myIndex[1]);
                sprintf( comGrid[comm->comm[cnt].part_lists[i]->myIndex[0]]
                                [comm->comm[cnt].part_lists[i]->myIndex[1]],
                        "%is%i",this_node,comm->comm[cnt].part_lists[i]->myIndex[2]);
            }
        }
        for(  i = comm->comm[cnt].n_part_lists/2; i < comm->comm[cnt].n_part_lists; i++){
            if( comm->comm[cnt].part_lists[i]->myIndex[2] > 1
             && comm->comm[cnt].part_lists[i]->myIndex[2] < 10){
                //fprintf(stderr,"%i: local rec  %i %i\n",this_node,comm->comm[cnt].part_lists[i]->myIndex[0],comm->comm[cnt].part_lists[i]->myIndex[1]);
                sprintf(&comGrid[comm->comm[cnt].part_lists[i]->myIndex[0]]
                                [comm->comm[cnt].part_lists[i]->myIndex[1]][4],
                        "%ir%i",this_node,comm->comm[cnt].part_lists[i]->myIndex[2]);
            }
        }
    }else if (comm->comm[cnt].type == GHOST_SEND){
        for(  i = 0; i < comm->comm[cnt].n_part_lists; i++){
            if( comm->comm[cnt].part_lists[i]->myIndex[2] > 1
             && comm->comm[cnt].part_lists[i]->myIndex[2] < 10){
                //fprintf(stderr,"%i: ghost send %i %i\n",this_node,comm->comm[cnt].part_lists[i]->myIndex[0],comm->comm[cnt].part_lists[i]->myIndex[1]);
                sprintf( comGrid[comm->comm[cnt].part_lists[i]->myIndex[0]]
                                [comm->comm[cnt].part_lists[i]->myIndex[1]],
                        "%iS%i",neighbor_index,comm->comm[cnt].part_lists[i]->myIndex[2]);
            }
        }
    }else if (comm->comm[cnt].type == GHOST_RECV){
        for(  i = 0; i < comm->comm[cnt].n_part_lists; i++){
            if( comm->comm[cnt].part_lists[i]->myIndex[2] > 1
             && comm->comm[cnt].part_lists[i]->myIndex[2] < 10){
                //fprintf(stderr,"%i: ghost rec  %i %i\n",this_node,comm->comm[cnt].part_lists[i]->myIndex[0],comm->comm[cnt].part_lists[i]->myIndex[1]);
                sprintf(&comGrid[comm->comm[cnt].part_lists[i]->myIndex[0]]
                                [comm->comm[cnt].part_lists[i]->myIndex[1]][4],
                        "%iR%i",neighbor_index,comm->comm[cnt].part_lists[i]->myIndex[2]);
            }
        }
    }
    /* print a diagram of this comm, only if it is in a particular z-slice */

    if(   comm->comm[cnt].type == GHOST_LOCL
      ||  comm->comm[cnt].type == GHOST_SEND
      ||  comm->comm[cnt].type == GHOST_RECV ){
        fprintf(ofile,"COMM TYPE: %i partner %i, size %i             \n",
                comm->comm[cnt].type,comm->comm[cnt].node,comm->comm[cnt].n_part_lists);
        for( j = dd.cell_grid[1]+1; j >=0; j--){
            for( i = 0; i < dd.cell_grid[0]+3;i++){
                fprintf(ofile,"|%s%s",comGrid[i][j],&comGrid[i][j][4]);
            }   
            fprintf(ofile,"|\n");
        }
        fprintf(ofile,"\n*****************\n");

        if( comm->comm[cnt].type == GHOST_LOCL ){
        fprintf(ofile,"COMM TYPE: %i\n",comm->comm[cnt].type);
        for( j = dd.cell_grid[1]+1; j >=0; j--){
            fprintf(ofile,".\n");
        }
        fprintf(ofile,"\n*****************\n");
            
        }
    }
  }

  for( i = 0; i < dd.cell_grid[0]+3;i++){
    for( j = 0; j < dd.cell_grid[1]+2;j++){
        delete [] comGrid[i][j];
    }
    delete [] comGrid[i];
  }
  delete [] comGrid;
}
#endif

/** Create communicators for cell structure domain decomposition. (see \ref GhostCommunicator) */
void  le_dd_prepare_comm(le_dd_comms_manager *mgr, GhostCommunicator *comm, int data_parts)
{
  static int le_cells_state_physical = 1;
  int dir,lr,i,cnt, num, n_comm_cells[3], send_rec, thisCommCount;
  int lc[3],hc[3], neighbor_index;

#ifdef LE_DEBUG
  if( comms_log != NULL ){ fclose(comms_log);comms_log=NULL;}
  char vLogName[64];
  sprintf(vLogName, "%i_comms_%i.dat", comms_count++,this_node);
  comms_log = fopen(vLogName, "w");
#endif
  
  CELL_TRACE(fprintf(stderr,"%d: neighbours:", this_node));
  CELL_TRACE(for(i = 0; i < my_neighbor_count; i++)fprintf(stderr," %d",node_neighbors[i]);)
  CELL_TRACE(fprintf(stderr,"\n"));
  
  /* prepare communicator... need 2 comms for node-node, 1 comm for intranode exchange */
  num = my_neighbor_count * 2;
  for(i = 0; i < my_neighbor_count; i++) {
        if( node_neighbors[i] == this_node ) num--;
  }

  /* flag data_parts is global to all individual sub-comms prepared in the communicator.*/
  CELL_TRACE(fprintf(stderr,"%d: Create Communicator: Aprep_comm ghost flags %d num_neighbor_comms: %d\n",this_node,data_parts,num));
  prepare_comm(comm, data_parts, num);

  /* number of cells owned by this node to communicate in each direction:
     the sequence of comms decides the size of the surface to share, because ghost cells are acquired at each comm. */
  n_comm_cells[1] =       dd.cell_grid[0]       *  dd.cell_grid[2];       /* order of comms is y-z-x */
  n_comm_cells[2] =       dd.cell_grid[0]       * (dd.cell_grid[1] + 2);  /*x-dir has one extra ghost layer (not all used).*/
  n_comm_cells[0] =  2 * (dd.cell_grid[1]+2)    * (dd.cell_grid[2] + 2);  /*x-dir has one extra ghost layer (not all used).*/

  
  /* loop over neighbours */
  cnt=0;
  for(i = 0; i < my_neighbor_count; i++) {
       
    int neighbor_coords[3], neighbor_rank;
    neighbor_index = mgr->comms_order[i];

    /* must loop on order y-x-z for Lees Edwards: each imaging must carry ghost cells already genned in the
       plane perpendicular to it. */
    if( neighbor_index == 2 || neighbor_index == 3 || neighbor_index >= 6 ) {
        dir   = 1; /* sharing a y-plane of the cube */
        lc[0] = 1;
        lc[2] = 1;
        hc[0] = dd.cell_grid[0];
        hc[2] = dd.cell_grid[2];
    }
    else if( neighbor_index == 0 || neighbor_index == 1 ){
        dir   = 0; /* sharing an x-plane of the cube */
        /* set lc and hc at the specific comm... not here.*/
    }
    else { /*if(neighbor_index == 4 || neighbor_index == 5){*/
        dir   = 2; /* sharing a z-plane of the cube */
        lc[0] = 1;
        hc[0] = dd.cell_grid[0];
        lc[1] = 0;
        hc[1] = dd.cell_grid[1]+1;
    } 
    lr            = node_neighbor_lr[neighbor_index]; 

    /* find out where this neighbor is */
    neighbor_rank = node_neighbors[neighbor_index];
    map_node_array(neighbor_rank, neighbor_coords);
 
    /* transfers along a given axis have lc[axis] == hc[axis] == fixed: you are sending/receiving
     * a layer of cells defined by that value of x,y or z. 
     *
     * For example, a receive from the cell above you in y would the range:
     * (??, dd.cell_grid[1] + 1, ??) for both lc and hc
     * while a send to the cell above you would be the range:
     * (??, dd.cell_grid[1], ??)
     * because the arriving cells go into the "ghost" plane, while the sent ones are 'real'.
     *
     * For example, a (send to /receive from) the cell below in y would be:
     * send: (??, 1, ??)
     * recv: (??, 0, ??)
     */
    
    if( neighbor_rank == this_node ) { /* if copying cells on a single node, then comms are a bit simpler... */

        comm->comm[cnt].type          = GHOST_LOCL;
        comm->comm[cnt].node          = this_node;
        
        /* prepare folding of ghost positions */
        if( data_parts & GHOSTTRANS_POSSHFTD )
            comm->comm[cnt].shift[dir] = node_neighbor_wrap[neighbor_index]*box_l[dir];

        /*if top plane */
        if( lr == 1 ) {
                    lc[dir] = dd.cell_grid[dir];
                    hc[dir] = dd.cell_grid[dir];
        }
        /*bottom plane*/
        else {
                    lc[dir] = 1;
                    hc[dir] = 1;
        }
        
        /* Buffer has to contain both Send and Recv cells -> factor 2 */
        comm->comm[cnt].part_lists    = (ParticleList **)Utils::malloc(2*n_comm_cells[dir]*sizeof(ParticleList *));

        switch( dir ){
          //send the single-thickness x-ghost layer, and then parts of the double-thickness
          //layer only when needed.
          case 0:
            lc[1] = 0;
            hc[1] = dd.cell_grid[1] + 1;
            lc[2] = 0;
            hc[2] = dd.cell_grid[2] + 1;
        
            /* double-thickness send */
            if( lr == 0 ) hc[0] = 2;
            thisCommCount  = dd_fill_comm_cell_lists(comm->comm[cnt].part_lists,lc,hc);
            break;
            
          case 1:
            if( data_parts & GHOSTTRANS_POSSHFTD ){
                comm->comm[cnt].shift[0]  = node_neighbor_wrap[neighbor_index]*lees_edwards_offset;
            }
            /* Wrap on send */
            thisCommCount  = mgr->wrap_fill_comm_cell_lists(comm->comm[cnt].part_lists,
                                                            lc,
                                                            hc,
                                                            node_pos[0],
                                                            node_pos[0],
                                                            lr == 1,
                                                            1);
            break;
          case 2:
            thisCommCount  = dd_fill_comm_cell_lists(comm->comm[cnt].part_lists,lc,hc);
            break;
        }
        comm->comm[cnt].n_part_lists  = thisCommCount;

        //if( dir == 1 || dir == 0 )
        CELL_TRACE(fprintf(stderr,"%d: Yprep_comm %d send %d    dir %i grid-range (%d,%d,%d)-(%d,%d,%d) shift: %f %f %f\n",
                           this_node,cnt,thisCommCount,dir,
                           lc[0],lc[1],lc[2],hc[0],hc[1],hc[2],
                           comm->comm[cnt].shift[0],comm->comm[cnt].shift[1],comm->comm[cnt].shift[2]);)
        
        /* fill recv comm cells */
        if( lr == 1 ) {
                    lc[dir] = 0;
                    hc[dir] = 0;
        }
        else {
                    lc[dir] = dd.cell_grid[dir] + 1;
                    hc[dir] = dd.cell_grid[dir] + 1;
        }

        
        switch( dir ){
          case 0:
                lc[1] = 0;
                hc[1] = dd.cell_grid[1] + 1;
                lc[2] = 0;
                hc[2] = dd.cell_grid[2] + 1;

                if( lr == 0 ) hc[0] = dd.cell_grid[0] + 2;
                
                thisCommCount = dd_fill_comm_cell_lists(&comm->comm[cnt].part_lists[thisCommCount],lc,hc);
                break;

          case 1:
                thisCommCount = mgr->wrap_fill_comm_cell_lists(&comm->comm[cnt].part_lists[thisCommCount],
                                                            lc,
                                                            hc,
                                                            node_pos[0],
                                                            node_pos[0],
                                                            lr != 1,
                                                            0);
                break;
          case 2:
                thisCommCount = dd_fill_comm_cell_lists(&comm->comm[cnt].part_lists[thisCommCount],lc,hc);
                break;
        }
        comm->comm[cnt].n_part_lists += thisCommCount;


        LE_TRACE(fprintf(comms_log,"%d send/rcv %d lists within node %d grid (%d,%d,%d)-(%d,%d,%d) dir %d shift (%f,%f,%f)\n",
                          cnt,comm->comm[cnt].n_part_lists,comm->comm[cnt].node,lc[0],lc[1],lc[2],hc[0],hc[1],hc[2],
                                   dir,comm->comm[cnt].shift[0],comm->comm[cnt].shift[1],comm->comm[cnt].shift[2]);)

        cnt++;

    }else{
        /* send/recv loop: sends and receives must synchronise, or at worst be in the correct order to avoid deadlocks. */
        for(send_rec=0; send_rec<2; send_rec++) {  
        
            comm->comm[cnt].node          = neighbor_rank;
 
            if( (send_rec == 0  && neighbor_rank > this_node) 
              ||(send_rec == 1  && neighbor_rank < this_node)){
                comm->comm[cnt].type          = GHOST_SEND;
                /* prepare fold-on-send of ghost positions */
                if(node_neighbor_wrap[neighbor_index] != 0 ){
                   if( data_parts & GHOSTTRANS_POSSHFTD ){ 
                        comm->comm[cnt].shift[dir]    = node_neighbor_wrap[neighbor_index]*box_l[dir];
                        if( dir == 1 )
                            comm->comm[cnt].shift[0]  = node_neighbor_wrap[neighbor_index]*lees_edwards_offset;
                   }
                }

                /*choose the plane to send*/
                if( lr == 1 ) {
                        lc[dir] = dd.cell_grid[dir];
                        hc[dir] = dd.cell_grid[dir];
                }else {
                        lc[dir] = 1;
                        hc[dir] = 1;
                }
                comm->comm[cnt].part_lists    = (ParticleList **)Utils::malloc(n_comm_cells[dir]*sizeof(ParticleList *));


                switch( dir ){
                     case 0:
                            lc[1] = 0;
                            hc[1] = dd.cell_grid[1] + 1;
                            lc[2] = 0;
                            hc[2] = dd.cell_grid[2] + 1;
            
                            if( 2 > dd.cell_grid[0] ){
                                if( le_cells_state_physical ){
                                    fprintf(stderr,"%i: WARNING: Cells per node is currently too small (%d of %d) for Lees-Edwards Simulation.\n",
                                            this_node,hc[dir],dd.cell_grid[dir]);
                                    fprintf(stderr,"%i: WARNING: This may fix itself during setup, so will continue and notify if it does.\n",this_node);
                                    le_cells_state_physical = 0;
                                }
                            }else if( le_cells_state_physical == 0 ){
                                fprintf(stderr,"%i: UN-WARNING: Cell system is now large enough for LE simulation.\n",this_node);
                                le_cells_state_physical = 1;
                            }


                            /* just send and receive an entire extra layer */
                            if( lr == 0 ) hc[0] = 2;
                            thisCommCount = dd_fill_comm_cell_lists(comm->comm[cnt].part_lists,lc,hc);
                            break;
                            
                    case 1:
                        if( node_neighbor_wrap[neighbor_index] != 0 )
                                thisCommCount = mgr->wrap_fill_comm_cell_lists(comm->comm[cnt].part_lists,
                                                                    lc,
                                                                    hc,
                                                                    lr != 1 ? node_pos[0] : neighbor_coords[0],
                                                                    lr == 1 ? node_pos[0] : neighbor_coords[0],
                                                                    lr == 1,
                                                                    1);
                        else
                                thisCommCount = dd_fill_comm_cell_lists(comm->comm[cnt].part_lists,lc,hc);
                        break;                        
                    case 2:
                        thisCommCount = dd_fill_comm_cell_lists(comm->comm[cnt].part_lists,lc,hc);
                        break;
                }
                comm->comm[cnt].n_part_lists  = thisCommCount;
                
                LE_TRACE(fprintf(comms_log,"%d send %d lists to node %d grid (%d,%d,%d)-(%d,%d,%d) dir %d shift (%f,%f,%f)\n",
                          cnt,comm->comm[cnt].n_part_lists,comm->comm[cnt].node,lc[0],lc[1],lc[2],hc[0],hc[1],hc[2],
                                   dir,comm->comm[cnt].shift[0],comm->comm[cnt].shift[1],comm->comm[cnt].shift[2]);)
                
                cnt++;
            }else {
                comm->comm[cnt].type          = GHOST_RECV;
          
                if( lr == 1 ) {
                    /* if top plane */
                    lc[dir] = dd.cell_grid[dir] + 1;
                    hc[dir] = dd.cell_grid[dir] + 1;
                }
                else {
                /* bottom plane */
                    lc[dir] = 0;
                    hc[dir] = 0;
                }
                comm->comm[cnt].part_lists    = (ParticleList **)Utils::malloc(n_comm_cells[dir]*sizeof(ParticleList *));

                switch( dir ){
                    case 0:
                        lc[1] = 0;
                        lc[2] = 0;
                        hc[1] = dd.cell_grid[1] + 1;
                        hc[2] = dd.cell_grid[2] + 1;

                        if( 2 > dd.cell_grid[0] ){
                            if( le_cells_state_physical ){
                                fprintf(stderr,"%i: WARNING: Cells per node is currently too small (%d of %d) for Lees-Edwards Simulation.\n",
                                        this_node,hc[dir],dd.cell_grid[dir]);
                                fprintf(stderr,"%i: WARNING: This may fix itself during setup, so will continue and notify if it does.\n",this_node);
                                le_cells_state_physical = 0;
                            }
                        }else if( le_cells_state_physical == 0 ){
                            fprintf(stderr,"%i: UN-WARNING: Cell system is now large enough for LE simulation.\n",this_node);
                            le_cells_state_physical = 1;
                        }


                        /* just send and receive an entire extra layer */
                        if( lr == 1 ) hc[0] = dd.cell_grid[0] + 2;
                        thisCommCount = dd_fill_comm_cell_lists(comm->comm[cnt].part_lists,lc,hc);
                        break;
                     
                    case 1:
                        if( node_neighbor_wrap[neighbor_index] != 0)
                                    thisCommCount = mgr->wrap_fill_comm_cell_lists(comm->comm[cnt].part_lists,
                                                                    lc,
                                                                    hc,
                                                                    lr != 1 ? node_pos[0] : neighbor_coords[0],
                                                                    lr == 1 ? node_pos[0] : neighbor_coords[0],
                                                                    lr == 1,
                                                                    0);
                        else
                                thisCommCount = dd_fill_comm_cell_lists(comm->comm[cnt].part_lists,lc,hc);
                        break;
                    case 2:
                        thisCommCount = dd_fill_comm_cell_lists(comm->comm[cnt].part_lists,lc,hc);
                        break;
                }
                comm->comm[cnt].n_part_lists  = thisCommCount;
          
                LE_TRACE(fprintf(comms_log,"%i recv %d lists from node %d grid (%d,%d,%d)-(%d,%d,%d) dir %d shift (%f,%f,%f)\n",
                          cnt,comm->comm[cnt].n_part_lists,comm->comm[cnt].node,lc[0],lc[1],lc[2],hc[0],hc[1],hc[2],
                                   dir,comm->comm[cnt].shift[0],comm->comm[cnt].shift[1],comm->comm[cnt].shift[2]);)
                
                cnt++;
            }
        }
      }
    }

#ifdef LE_DEBUG
    printAllComms(comms_log, comm, 0);
    if( comms_log != NULL ){ fclose(comms_log);comms_log=NULL;}
#endif
    
}

/** Dynamic update of communicators for cell structure domain decomposition. (see \ref GhostCommunicator) */
void  le_dd_dynamic_update_comm(le_dd_comms_manager *mgr, GhostCommunicator *comm, int data_parts, int backwards)
{
  int lr,i,cnt, num, send_rec, thisCommCount;
  int lc[3],hc[3], neighbor_index;

#ifdef LE_DEBUG
  if( comms_log != NULL ){ fclose(comms_log);comms_log=NULL;}
  char vLogName[64];
  sprintf(vLogName, "%i_recomm_%i.dat", comms_count++, this_node);
  comms_log = fopen(vLogName, "w");
#endif

  LE_TRACE(fprintf(stderr,"%d: neighbours:", this_node);)
  LE_TRACE(for(i = 0; i < my_neighbor_count; i++)fprintf(stderr," %d",node_neighbors[i]);)
  LE_TRACE(fprintf(stderr,"\n");)

  /* calculate the total number of communicators */
  num = my_neighbor_count * 2;
  for(i = 0; i < my_neighbor_count; i++) {
        if( node_neighbors[i] == this_node ) num--;
  }

  /* loop over neighbours */
  cnt               = 0;
  if(backwards == LE_COMM_BACKWARDS) cnt = num - 1;
  
  for(i = 0; i < my_neighbor_count; i++) {

    int neighbor_coords[3], neighbor_rank;
    neighbor_index  = mgr->comms_order[i];
    neighbor_rank   = node_neighbors[neighbor_index];
        
    /* non-imaged comms or comms in x or z directions do not need to be updated */
    if( (neighbor_rank != this_node && node_neighbor_wrap[neighbor_index] == 0)
     ||   neighbor_index == 0
     ||   neighbor_index == 1
     ||   neighbor_index == 4
     ||   neighbor_index == 5 ) {
        /* skip one comm if this node is talking to itself, otherwise two comms (send+recv)*/
        if( neighbor_rank == this_node )
            backwards == LE_COMM_BACKWARDS ? cnt--  : cnt++;
        else
            backwards == LE_COMM_BACKWARDS ? cnt-=2 : cnt+=2;
        continue;
    }
    
    LE_TRACE(fprintf(stderr,"%d: Not Skipping comm: %d  %d dir 1 shift: %f %f %f\n\n",
                           this_node,cnt,thisCommCount,
                           comm->comm[cnt].shift[0],comm->comm[cnt].shift[1],comm->comm[cnt].shift[2]);)
        
    /* sharing a y-plane of the cube */
    lc[0] = 1;
    lc[2] = 1;
    hc[0] = dd.cell_grid[0];
    hc[2] = dd.cell_grid[2];
    
    lr            = node_neighbor_lr[neighbor_index];

    /* find out where this neighbor is */
    map_node_array(neighbor_rank, neighbor_coords);

    if( neighbor_rank == this_node ) { /* if copying cells on a single node, then comms are a bit simpler... */

        /*if top plane */
        if( lr == 1 ) {
                    lc[1] = dd.cell_grid[1];
                    hc[1] = dd.cell_grid[1];
        }
        /*bottom plane*/
        else {
                    lc[1] = 1;
                    hc[1] = 1;
        }

        /* sending or receiving coords or coords+vels through a y-wrap */
        if( data_parts & GHOSTTRANS_POSSHFTD ){
                comm->comm[cnt].shift[1]  = node_neighbor_wrap[neighbor_index]*box_l[1];
                comm->comm[cnt].shift[0]  = node_neighbor_wrap[neighbor_index]*lees_edwards_offset;
        }
        /* Wrap on send */
        LE_TRACE(fprintf(stderr,"%d: Comm to reinit: %d send %d    dir %i grid-range (%d,%d,%d)-(%d,%d,%d) shift: %f %f %f\n",
                           this_node,cnt,thisCommCount,dir,
                           lc[0],lc[1],lc[2],hc[0],hc[1],hc[2],
                           comm->comm[cnt].shift[0],comm->comm[cnt].shift[1],comm->comm[cnt].shift[2]);)
        
        thisCommCount  = mgr->wrap_fill_comm_cell_lists(comm->comm[cnt].part_lists,
                                                             lc,
                                                             hc,
                                                             node_pos[0],
                                                             node_pos[0],
                                                             lr == 1,
                                                             1);
        comm->comm[cnt].n_part_lists  = thisCommCount;

        LE_TRACE(fprintf(stderr,"%d: Reinitted comm %d send %d    dir %i grid-range (%d,%d,%d)-(%d,%d,%d) shift: %f %f %f\n\n",
                           this_node,cnt,thisCommCount,dir,
                           lc[0],lc[1],lc[2],hc[0],hc[1],hc[2],
                           comm->comm[cnt].shift[0],comm->comm[cnt].shift[1],comm->comm[cnt].shift[2]);)


        /* fill recv comm cells */
        if( lr == 1 ) {
                    lc[1] = 0;
                    hc[1] = 0;
        }
        else {
                    lc[1] = dd.cell_grid[1] + 1;
                    hc[1] = dd.cell_grid[1] + 1;
        }

        LE_TRACE(fprintf(stderr,"%d: Comm to reinit: %d receive %d grid-range (%d,%d,%d)-(%d,%d,%d) shift: %f %f %f\n",
                           this_node,cnt,thisCommCount,
                           lc[0],lc[1],lc[2],hc[0],hc[1],hc[2],
                           comm->comm[cnt].shift[0],comm->comm[cnt].shift[1],comm->comm[cnt].shift[2]);)

        thisCommCount = mgr->wrap_fill_comm_cell_lists(&comm->comm[cnt].part_lists[thisCommCount],
                                                            lc,
                                                            hc,
                                                            node_pos[0],
                                                            node_pos[0],
                                                            lr != 1,
                                                            0);
        comm->comm[cnt].n_part_lists += thisCommCount;



        LE_TRACE(fprintf(comms_log,"%i send/rec %d lists within node %d grid (%d,%d,%d)-(%d,%d,%d) dir %d shift (%f,%f,%f)\n",
                          cnt,comm->comm[cnt].n_part_lists,comm->comm[cnt].node,lc[0],lc[1],lc[2],hc[0],hc[1],hc[2],
                                   dir,comm->comm[cnt].shift[0],comm->comm[cnt].shift[1],comm->comm[cnt].shift[2]);)
                
        if( backwards == LE_COMM_BACKWARDS ){
            /* exchange the send and receive parts of the comm */
            int nlist2=comm->comm[cnt].n_part_lists/2;
            for(int j=0;j<nlist2;j++) {
                ParticleList* tmplist = comm->comm[cnt].part_lists[j];
                comm->comm[cnt].part_lists[j]        = comm->comm[cnt].part_lists[j+nlist2];
                comm->comm[cnt].part_lists[j+nlist2] = tmplist;
            }
            cnt--;
        }else{
            cnt++;
        }
    }else{
        /* send/recv loop: sends and receives must synchronise, or at worst be in the correct order to avoid deadlocks. */
        for(send_rec=0; send_rec<2; send_rec++) {

            /* prepare fold-on-send/fold-on-rec of ghost positions (or pos+vels)*/
            if(node_neighbor_wrap[neighbor_index] != 0 ){
                if( data_parts & GHOSTTRANS_POSSHFTD ){
                        comm->comm[cnt].shift[1]  = node_neighbor_wrap[neighbor_index]*box_l[1];
                        comm->comm[cnt].shift[0]  = node_neighbor_wrap[neighbor_index]*lees_edwards_offset;
                }
            }
                
            if( (send_rec == 0  && neighbor_rank > this_node)
              ||(send_rec == 1  && neighbor_rank < this_node)){


                /*choose the plane to send*/
                if( lr == 1 ) {
                        lc[1] = dd.cell_grid[1];
                        hc[1] = dd.cell_grid[1];
                }else {
                        lc[1] = 1;
                        hc[1] = 1;
                }

                thisCommCount = mgr->wrap_fill_comm_cell_lists(comm->comm[cnt].part_lists,
                                                                    lc,
                                                                    hc,
                                                                    lr != 1 ? node_pos[0] : neighbor_coords[0],
                                                                    lr == 1 ? node_pos[0] : neighbor_coords[0],
                                                                    lr == 1,
                                                                    1);
                comm->comm[cnt].n_part_lists  = thisCommCount;


                LE_TRACE(fprintf(comms_log,"%i send %d lists to node %d grid (%d,%d,%d)-(%d,%d,%d) dir %d shift (%f,%f,%f)\n",
                          cnt,comm->comm[cnt].n_part_lists,comm->comm[cnt].node,lc[0],lc[1],lc[2],hc[0],hc[1],hc[2],
                                   dir,comm->comm[cnt].shift[0],comm->comm[cnt].shift[1],comm->comm[cnt].shift[2]);)
                

                backwards == LE_COMM_BACKWARDS ? cnt-- : cnt++;
            }else {
                if( lr == 1 ) {
                    /* if top plane */
                    lc[1] = dd.cell_grid[1] + 1;
                    hc[1] = dd.cell_grid[1] + 1;
                }
                else {
                /* bottom plane */
                    lc[1] = 0;
                    hc[1] = 0;
                }
                
                thisCommCount = mgr->wrap_fill_comm_cell_lists(comm->comm[cnt].part_lists,
                                                                    lc,
                                                                    hc,
                                                                    lr != 1 ? node_pos[0] : neighbor_coords[0],
                                                                    lr == 1 ? node_pos[0] : neighbor_coords[0],
                                                                    lr == 1,
                                                                    0);
                comm->comm[cnt].n_part_lists  = thisCommCount;

                LE_TRACE(fprintf(comms_log,"%i recv %d lists from node %d grid (%d,%d,%d)-(%d,%d,%d) dir %d shift (%f,%f,%f)\n",
                          cnt,comm->comm[cnt].n_part_lists,comm->comm[cnt].node,lc[0],lc[1],lc[2],hc[0],hc[1],hc[2],
                                   dir,comm->comm[cnt].shift[0],comm->comm[cnt].shift[1],comm->comm[cnt].shift[2]);)

                backwards == LE_COMM_BACKWARDS ? cnt-- : cnt++;
            }
        }
      }
    }
#ifdef LE_DEBUG
    printAllComms(comms_log, comm, backwards);
    if( comms_log != NULL ){ fclose(comms_log);comms_log=NULL;}
#endif
}



#endif //LEES_EDWARDS
