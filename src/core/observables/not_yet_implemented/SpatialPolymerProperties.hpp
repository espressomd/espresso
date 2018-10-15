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
int ObbservableSpatialPolymerProperties::actual_calculate(PartCfg &partCfg) {
  if (!sortPartCfg()) {
    ostringstream errtxt; //= runtime_error(128);
    errtxt << "{094 could not sort partCfg}";
    runtimeError(errtxt);
    return -1;
  }
  double *A = last_value;
  int poly_len = n;
  for (int i = 0; i < poly_len; i++)
    A[i] = 0.0;

#ifdef ELECTROSTATICS
  spatial_polym_data *p_data = (spatial_polym_data *)container;
  IntList *ids = p_data->id_list;
  for (int i = 0; i < ids->n; i++) {
    A[i % poly_len] += partCfg[ids->e[i]].e->p.q / ((double)p_data->npoly);
  }
#endif
  return 0;
}
