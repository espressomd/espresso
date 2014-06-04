/*
  Copyright (C) 2010,2012,2013 The ESPResSo project
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
/** \file thermostat.cpp
    Implementation of \ref thermostat.hpp "thermostat.h"
 */
#include "thermostat.hpp"
#include "communication.hpp"
#include "lattice.hpp"
#include "npt.hpp"
#include "ghmc.hpp"

#if (!defined(FLATNOISE) && !defined(GAUSSRANDOMCUT) && !defined(GAUSSRANDOM))
#define FLATNOISE
#endif


/* thermostat switch */
int thermo_switch = THERMO_OFF;
/** Temperature */
double temperature = -1.0;

/* LANGEVIN THERMOSTAT */
/* Langevin friction coefficient gamma. */
double langevin_gamma = 0.0;
/* Friction coefficient gamma for rotation */
double langevin_gamma_rotation;

/* NPT ISOTROPIC THERMOSTAT */
// INSERT COMMENT
double nptiso_gamma0 = 0.0;
// INSERT COMMENT
double nptiso_gammav = 0.0;

/* GHMC THERMOSTAT */
// Number of NVE-MD steps in each GHMC cycle
int ghmc_nmd = 1;
// phi parameter for partial momentum update step in GHMC
double ghmc_phi = 0;

double langevin_pref1, langevin_pref2, langevin_pref2_rotation;
/** buffers for the work around for the correlated random values which cool the system,
    and require a magical heat up whenever reentering the integrator. */
static double langevin_pref2_buffer, langevin_pref2_rotation_buffer;

#ifdef NPT
double nptiso_pref1;
double nptiso_pref2;
double nptiso_pref3;
double nptiso_pref4;
#endif


void thermo_init_langevin() 
{
  langevin_pref1 = -langevin_gamma/time_step;
#if defined (FLATNOISE)
  langevin_pref2 = sqrt(24.0*temperature*langevin_gamma/time_step);
#elif defined (GAUSSRANDOMCUT) || defined (GAUSSRANDOM)
  langevin_pref2 = sqrt(2.0*temperature*langevin_gamma/time_step);
#else
#error No Noise defined
#endif
  
#ifdef ROTATION 
  langevin_gamma_rotation = langevin_gamma/3;
#if defined (FLATNOISE)
  langevin_pref2_rotation = sqrt(24.0*temperature*langevin_gamma_rotation/time_step);
#elif defined (GAUSSRANDOMCUT) || defined (GAUSSRANDOM)
  langevin_pref2_rotation = sqrt(2.0*temperature*langevin_gamma_rotation/time_step);
#else
#error No Noise defined
#endif
  THERMO_TRACE(fprintf(stderr,"%d: thermo_init_langevin: langevin_gamma_rotation=%f, langevin_pref2_rotation=%f",this_node, langevin_gamma_rotation,langevin_pref2_rotation));
#endif
  THERMO_TRACE(fprintf(stderr,"%d: thermo_init_langevin: langevin_pref1=%f, langevin_pref2=%f",this_node,langevin_pref1,langevin_pref2));  
}

#ifdef NPT
void thermo_init_npt_isotropic()
{
  if (nptiso.piston != 0.0) {
#if defined (FLATNOISE)
    nptiso_pref1 = -nptiso_gamma0*0.5 * time_step;
    nptiso_pref2 = sqrt(12.0*temperature*nptiso_gamma0*time_step) * time_step;
    nptiso_pref3 = -nptiso_gammav*(1.0/nptiso.piston)*0.5*time_step;
    nptiso_pref4 = sqrt(12.0*temperature*nptiso_gammav*time_step);
#elif defined (GAUSSRANDOMCUT) || defined (GAUSSRANDOM)
    nptiso_pref1 = -nptiso_gamma0*0.5 * time_step;
    nptiso_pref2 = sqrt(1.0*temperature*nptiso_gamma0*time_step) * time_step;
    nptiso_pref3 = -nptiso_gammav*(1.0/nptiso.piston)*0.5*time_step;
    nptiso_pref4 = sqrt(1.0*temperature*nptiso_gammav*time_step);
#else
#error No Noise defined
#endif
    THERMO_TRACE(fprintf(stderr,"%d: thermo_init_npt_isotropic: nptiso_pref1=%f, nptiso_pref2=%f, nptiso_pref3=%f, nptiso_pref4=%f \n",this_node,nptiso_pref1,nptiso_pref2,nptiso_pref3,nptiso_pref4));
  }
  else {
    thermo_switch = ( thermo_switch ^ THERMO_NPT_ISO );
    THERMO_TRACE(fprintf(stderr,"%d: thermo_init_npt_isotropic: switched off nptiso (piston=%f; thermo_switch=%d) \n",this_node,nptiso.piston,thermo_switch));
  }
}
#endif

void thermo_init()
{
  if(thermo_switch == THERMO_OFF){
    return;
  }
#ifdef INTER_DPD
  if(thermo_switch & THERMO_INTER_DPD)  inter_dpd_init();
#endif
  if(thermo_switch & THERMO_LANGEVIN ) thermo_init_langevin();
#ifdef DPD
  if(thermo_switch & THERMO_DPD) thermo_init_dpd();
#endif
#ifdef NPT
  if(thermo_switch & THERMO_NPT_ISO)   thermo_init_npt_isotropic();
#endif
#ifdef GHMC
  if(thermo_switch & THERMO_GHMC) thermo_init_ghmc();
#endif
}

void thermo_heat_up()
{
  if(thermo_switch & THERMO_LANGEVIN) {
    langevin_pref2_buffer          = langevin_pref2;
    langevin_pref2_rotation_buffer = langevin_pref2_rotation;
    langevin_pref2 *= sqrt(3);
    langevin_pref2_rotation *= sqrt(3);
  }
#ifdef DPD
  else if (thermo_switch & THERMO_DPD){dpd_heat_up();}
#endif
#ifdef INTER_DPD
  else if (thermo_switch & THERMO_INTER_DPD) {inter_dpd_heat_up();}
#endif
}

void thermo_cool_down()
{
  if(thermo_switch & THERMO_LANGEVIN) {
    langevin_pref2          = langevin_pref2_buffer;
    langevin_pref2_rotation = langevin_pref2_rotation_buffer;
  }
#ifdef DPD
  else if (thermo_switch & THERMO_DPD){dpd_cool_down();}
#endif
#ifdef INTER_DPD
  else if (thermo_switch & THERMO_INTER_DPD) inter_dpd_cool_down();
#endif
}

