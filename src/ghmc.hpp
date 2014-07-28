/*
  Copyright (C) 2010,2012,2013 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany
  
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
#ifndef GHMC_H
#define GHMC_H
/** \file ghmc.hpp

    This file contains the implementation of the GHMC (Generalized 
    Hybrid Monte Carlo) thermostat.
    
 */

#define GHMC_MFLIP_OFF       0
#define GHMC_MFLIP_ON        1
#define GHMC_MFLIP_RAND      2

#define GHMC_TSCALE_OFF      0
#define GHMC_TSCALE_ON       1

#define GHMC_MOVE_REJECT     0
#define GHMC_MOVE_ACCEPT     1

#ifdef GHMC

/************************************************/
/** \name Data Types */
/************************************************/
/*@{*/
/** Data structure describing a slab and the velocities occuring their in. */
typedef struct {

	/** MC  statistics variables */
  int att, acc;
	double hmlt_old, hmlt_new;
	
} Ghmc;
/*@}*/

/************************************************************/
/** \name Exported Variables */
/************************************************************/
/*@{*/

/** Structure containing the ghmc relevant information */
extern Ghmc ghmcdata;
/*@}*/

#endif

extern int ghmc_mflip;
extern int ghmc_tscale;
extern int ghmc_att, ghmc_acc;
extern int ghmc_mc_res;

/************************************************************/
/** \name Exported Functions */
/************************************************************/
/*@{*/

#ifdef GHMC 

/**  initilize global thermostat parameters*/
void thermo_init_ghmc();

/**  initilize ghmc during integration*/
void ghmc_init();

/**  momentum update step of ghmc */
void ghmc_momentum_update();

/**  ghmc MC step*/
void ghmc_mc();

/**  ghmc close proc*/
void ghmc_close();

/**  save particles state */
void save_last_state();

/**  load particles last saved state */
void load_last_state();

/*@}*/

/* endif GHMC */
#endif

/* endif GHMC_H */
#endif
