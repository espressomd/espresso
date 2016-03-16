/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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
#ifndef DPD_H
#define DPD_H
/** \file dpd.hpp
 *  Routines to use dpd as thermostat or pair force
 *  T. Soddemann, B. Duenweg and K. Kremer, Phys. Rev. E 68, 046702 (2003)
 *  \ref forces.cpp
*/

#include "utils.hpp"
#include "thermostat.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "virtual_sites.hpp"

/** Flag to decide wether to allow for fixed particles with DPD */
extern int dpd_ignore_fixed_particles;

/** DPD Friction coefficient gamma. */
extern double dpd_gamma;
/** DPD thermostat cutoff */
extern double dpd_r_cut;
/** DPD thermostat weight function */
extern int dpd_wf;

/** DPD transversal Friction coefficient gamma. */
extern double dpd_tgamma;
/** trans DPD thermostat cutoff */
extern double dpd_tr_cut;
/** trans DPD thermostat weight function */
extern int dpd_twf;

#ifdef DPD
extern double dpd_r_cut_inv;
extern double dpd_pref1;
extern double dpd_pref2;

#ifdef TRANS_DPD 
extern double dpd_tr_cut_inv;
extern double dpd_pref3;
extern double dpd_pref4;
#endif

void dpd_switch_off(void);
void thermo_init_dpd();
void dpd_heat_up();
void dpd_cool_down();

/** Calculate Random Force and Friction Force acting between particle
    p1 and p2 and add them to their forces. */
void add_dpd_thermo_pair_force(Particle * p1, Particle * p2,
                               double d[3], double dist, double dist2);
#endif

#ifdef INTER_DPD
void inter_dpd_heat_up();
void inter_dpd_cool_down();
void inter_dpd_switch_off();
int inter_dpd_set_params(int part_type_a, int part_type_b,
			 double gamma, double r_c, int wf,
			 double tgamma, double tr_c,
			 int twf);
void inter_dpd_init();
void inter_dpd_update_params(double pref2_scale);

void add_inter_dpd_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
                              double d[3], double dist, double dist2);
#endif

#endif

