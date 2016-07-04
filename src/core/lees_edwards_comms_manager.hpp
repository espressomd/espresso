/*
  Copyright (C) 2010,2012,2016 The ESPResSo project
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
/** \file lees_edwards_comms_manager.hpp
 Class definition and a few macros for Lees Edwards comms rebuild manager.
*/
#ifndef LEES_EDWARDS_CM_H
#define LEES_EDWARDS_CM_H

#include "utils.hpp"
#include "cells.hpp"
#include "integrate.hpp"
#include "verlet.hpp"
#include "thermostat.hpp"

/** \name Switches and state flags for Lees Edwards */
/************************************************************/
/*@{*/
typedef enum{
    LE_COMM_FORWARDS = 0,
    LE_COMM_BACKWARDS
} LE_COMM_REVERSED_T;

typedef enum{
    LE_CELLS_SAME = 0,
    LE_CELLS_SHIFTED
} LE_CELLS_SHIFTED_T;

/*@}*/


/************************************************************/
/** \name Class Definitions */
/************************************************************/
/*@{*/
/** Class to manage rebuilds
 *  of the communicators after a Lees-Edwards offset
 *  change
 */
class le_dd_comms_manager {
    public:
    /** Setup the comms order.*/
        void                init(int neighbor_count);

   /** change some state based on a new lees-edwards offset */
        int                 update_on_le_offset_change();
        
    /** Update a y-direction comm after a Lees-Edwards offset change.*/
        int                 wrap_fill_comm_cell_lists(Cell **part_lists,
                                                     int    lc[3],
                                                     int    hc[3],
                                                     int    lower_xIndex,
                                                     int    upper_xIndex,
                                                     int    amUpper,
                                                     int    amSender);
        
    /** Array of neighbour indexes in order to carry out the comms.*/
        unsigned short int *comms_order;

    /** destructor for the comms manager */
       ~le_dd_comms_manager(){
            delete [] comms_order;
        }
        
    private:
    /** Handy to have the number of cells in x for the system.*/
        int                 total_xcells;
    /** Positive integer number of cells which the current LE offset implies.*/
        int                 le_offset_cells_up, le_offset_cells_down;
    /** Has there been a shift in the cell wrap lately?*/
        int                 cells_shifted;
};
/*@}*/
#endif //LEES_EDWARDS_CM_H
