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
/** \file lees_edwards_comms_manager.cpp
 *
 *  This file defined methods for a class to hold some
 *  state for Lees-Edwards Comms updates.  Comms structures
 *  need to be partly rebuilt whenever the Lees-Edwards offset
 *  changes.
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

/** Allocate and fill a short int array with the order of comms.
    The neighbor index is into \ref node_neighbors is stored,
    not the actual node id (a comms system with one node only
    will have 6 entries in \ref node_neighbours, all pointing to itself).*/
void le_dd_comms_manager :: init(int neighbor_count){

  int                        i, cnt;

  comms_order = new unsigned short int [neighbor_count];

  i           = 0;
  comms_order[i++] = 2;  /* regular y neighbors come first */
  comms_order[i++] = 3;

  /* extra LE y-neighbors */
  for( cnt = 6; cnt < neighbor_count; cnt++ )
      comms_order[i++] = cnt;

  /* z and then x neighbors */
  comms_order[i++] = 4;
  comms_order[i++] = 5;
  comms_order[i++] = 0;
  comms_order[i++] = 1;

  le_offset_cells_up   = -1;
  le_offset_cells_down = -1;

  update_on_le_offset_change();

}

/** Change some state based on a new lees-edwards offset,
    variables saved are \ref total_xcells,
    \ref le_offset_cells_down, \ref le_offset_cells_up*/
int le_dd_comms_manager :: update_on_le_offset_change(){

        double le_dx;
        int    old_off_up, old_off_down;

        old_off_up   = le_offset_cells_up;
        old_off_down = le_offset_cells_down;

        total_xcells = dd.cell_grid[0] * node_grid[0];

        le_dx        = lees_edwards_offset;
        while( le_dx < 0.0 )      {le_dx += box_l[0]; }
        while( le_dx >= box_l[0] ){le_dx -= box_l[0]; }

        /* round to the right (up) on transfer from lower to upper
           because the extra neighbor cell is provided to the right. */
        le_offset_cells_up = int(le_dx * dd.inv_cell_size[0]);
        if( le_dx * dd.inv_cell_size[0] > (float)le_offset_cells_up )
            le_offset_cells_up++;

        /* round to the right (down) on transfer from upper to lower as well
           because the extra neighbor cell is also provided to the right. */
        le_offset_cells_down = -1 * int(le_dx * dd.inv_cell_size[0]);
        while(le_offset_cells_down < 0) le_offset_cells_down += total_xcells;

        if( le_offset_cells_down == old_off_down
         && le_offset_cells_up   == old_off_up ){
             cells_shifted = LE_CELLS_SAME;
             return( LE_CELLS_SAME );
        }else{
             cells_shifted = LE_CELLS_SHIFTED;
             return( LE_CELLS_SHIFTED );
        }
}


/**Sending or receiving cells in Y direction with LE wrap means looking at the
 * existing cell list with an offset.  This function fills the comm list
 * using cells with an offset or no cells at all, depending on what is required.*/
int le_dd_comms_manager :: wrap_fill_comm_cell_lists(Cell        **part_lists,
                                                      int           lc[3],
                                                      int           hc[3],
                                                      int    lower_xIndex,
                                                      int    upper_xIndex,
                                                      int    amUpper,
                                                      int    amSender)
{

  int    i,m,n,o,c=0;
  int    imaged_o, lower_base_o, lower_top_o;
  int    upper_base_o, upper_top_o;
  int    doIt;

  lower_base_o    = lower_xIndex    * dd.cell_grid[0];
  lower_top_o     = lower_base_o    + dd.cell_grid[0];
  upper_base_o    = upper_xIndex    * dd.cell_grid[0];
  upper_top_o     = upper_base_o    + dd.cell_grid[0];

  /* actually copy a different cell to o (or no cells),
   * based on Lees-Edwards imaging*/
  doIt = 0;
  for(o=lc[0]; o<=hc[0]; o++){

      if( amUpper == 1 ){

        /* find cell pos in the global list of cells */
        if( amSender )
            imaged_o = (upper_base_o + (o-1) + le_offset_cells_down) % total_xcells;
        else
            imaged_o = (lower_base_o + (o-1) + le_offset_cells_up) % total_xcells;


        /* test if the received cell is in the range of the receiving node */
        if( ( amSender && ( imaged_o < lower_base_o || imaged_o >= lower_top_o ) )
          ||(!amSender && ( imaged_o < upper_base_o || imaged_o >= upper_top_o ) ) ){
            doIt = 0;
        }
        else {
            doIt = 1;
            if( !amSender ) {
                imaged_o = imaged_o - upper_base_o + 1;
            }
        }
      }else{
        /* find cell pos in the global list of cells */
        if( amSender )
            imaged_o = (lower_base_o + (o-1) + le_offset_cells_up)   % total_xcells;
        else
            imaged_o = (upper_base_o + (o-1) + le_offset_cells_down) % total_xcells;

        /* test if the received cell is in the range of the receiving node */
        if( ( amSender && ( imaged_o < upper_base_o || imaged_o >= upper_top_o ) )
          ||(!amSender && ( imaged_o < lower_base_o || imaged_o >= lower_top_o ) ) ){
            doIt = 0;
        }
        else{
            doIt = 1;
            if( !amSender ) {
                imaged_o = imaged_o - lower_base_o + 1;
            }
        }
      }


    CELL_TRACE(fprintf(stderr,"%i: linking cell with range %f -- %f %+f to cell with x range: %f -- %f\n",
                    this_node, (o-1)*dd.cell_size[0],(o)*dd.cell_size[0],lees_edwards_offset,
                    (imaged_o-1)*dd.cell_size[0],(imaged_o)*dd.cell_size[0]);)

    for(n=lc[1]; n<=hc[1]; n++) {
      for(m=lc[2]; m<=hc[2]; m++) {

        if( amSender )
            i = get_linear_index(o,n,m,dd.ghost_cell_grid);
        else
            i = get_linear_index(imaged_o,n,m,dd.ghost_cell_grid);

        if( doIt ){
            part_lists[c] = &cells[i];
            c++;
        }
      }
    }
  }
  return c;
}
#endif