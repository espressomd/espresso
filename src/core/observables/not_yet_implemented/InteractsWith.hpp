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
int ObservableInteractsWith::actual_calculate(PartCfg &partCfg) {
  double *A = last_value;
  iw_params *params = (iw_params *)container;
  IntList *ids1;
  IntList *ids2;
  int i, j;
  //  double dx,dy,dz;
  double dist2;
  double cutoff2 = params->cutoff * params->cutoff;
  double pos1[3], pos2[3], dist[3];
  ids1 = params->ids1;
  ids2 = params->ids2;
  if (!sortPartCfg()) {
    runtimeErrorMsg() << "could not sort partCfg";
    return -1;
  }
  for (i = 0; i < ids1->n; i++) {
    if (ids1->e[i] >= n_part)
      return 1;
    pos1[0] = partCfg[ids1->e[i]].r.p[0];
    pos1[1] = partCfg[ids1->e[i]].r.p[1];
    pos1[2] = partCfg[ids1->e[i]].r.p[2];
    for (j = 0; j < ids2->n; j++) {
      if (ids2->e[j] >= n_part)
        return 1;
      if (ids2->e[j] == ids1->e[i]) // do not count self-interaction :-)
        continue;
      A[i] = 0;
      pos2[0] = partCfg[ids2->e[j]].r.p[0];
      pos2[1] = partCfg[ids2->e[j]].r.p[1];
      pos2[2] = partCfg[ids2->e[j]].r.p[2];
      get_mi_vector(dist, pos1, pos2);
      dist2 = dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2];
      if (dist2 < cutoff2) {
        A[i] = 1;
        break;
        // interaction found for i, go for next
      }
    }
  }
  return 0;
}
