/*
  Copyright (C) 2016,2017 The ESPResSo project

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

#ifndef OBSERVABLES_OBSERVABLE_HPP
#define OBSERVABLES_OBSERVABLE_HPP

#include <fstream>
#include <vector>
#include <string>

#include "core/PartCfg.hpp"

namespace Observables {
class Observable {
public:
  Observable();
  virtual ~Observable() = default;
  int calculate();

  /* IO functions for observables */
  void set_filename(std::string const &filename, bool binary);
  bool writable() const;
  void write();

  virtual int n_values() const { return 0; }
  std::vector<double> last_value;

private:
  virtual int actual_calculate(PartCfg & partCfg) = 0;

  int n;
  double last_update;

  virtual void do_write();
  std::ofstream m_ofile;
  std::string m_filename;
  bool m_binary;
};

//\//\typedef struct {
//\  Observable* reference_observable;
//\  unsigned int n_sweeps;
//\} observable_average_container;
//\
//\
//\typedef struct {
//\// FIXME finish the implementation of scattering length
//\  IntList* id_list;
//\  DoubleList *scattering_length; // Scattering lengths of particles
//\  int order;
//\  int dim_sf; // number of q vectors
//\  int *q_vals; // values of q vectors
//\  double *q_density; // number of q vectors per bin
//\  // entries for fast version
//\  //int num_k_vecs;
//\  int k_density;
//\} observable_sf_params;
//\
//\/** See if particles from idList1 interact with any of the particles in idList2 
//\input parameters are passed via struct iw_params
//\*/
//\
//\typedef struct {
//\  double cutoff;
//\  IntList *ids1;
//\  IntList *ids2;
//\} iw_params;
//\
//\
//\void mpi_observable_lb_radial_velocity_profile_slave_implementation();
//\
//\
//\typedef struct { 
//\    IntList *id_list;
//\    int type;
//\    double minr;
//\    double maxr;
//\    int rbins;
//\    int start_point_id;
//\    int end_point_id;
//\    // id_flag == 0 : actual positions given, otherwise two particle ids for the start- and 
//\    // end-point are given
//\    int id_flag;
//\    double start_point[3];
//\    double end_point[3];
//\} radial_density_data;
//\
 //\
 //\typedef struct { 
//\    IntList *id_list;
//\    int npoly;
//\    int cut_off;
//\} spatial_polym_data;
//\
 //\// uses the same data as spatial_polymer_properties
//\
 //\typedef struct {
//\    IntList *id_list;
//\    int poly_len;
//\    int npoly;
//\    int k;
//\    int n_bins;
//\    double r_min;
//\    double r_max;
//\} k_dist_data;
//\
 //\typedef struct {
//\  int *p1_types;
//\  int n_p1;
//\  int *p2_types;
//\  int n_p2;
//\  double r_min;
//\  double r_max;
//\  int r_bins;
//\} rdf_profile_data;

} // Namespace Observables
#endif
