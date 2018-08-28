/*
Copyright (C) 2010-2018 The ESPResSo project

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
int ObservableRdf::actual_calculate(PartCfg &partCfg) {
  if (!sortPartCfg()) {
    return 1;
  }
  double *last = last_value;
  rdf_profile_data *rdf_data = (rdf_profile_data *)container;
  calc_rdf(rdf_data->p1_types, rdf_data->n_p1, rdf_data->p2_types,
           rdf_data->n_p2, rdf_data->r_min, rdf_data->r_max, rdf_data->r_bins,
           last);
  return 0;
}
