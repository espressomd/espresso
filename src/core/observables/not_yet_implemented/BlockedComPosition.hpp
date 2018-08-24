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
int ObservableBlockedComPosition::actual_calculate(PartCfg &partCfg) {
  double *A = last_value;
  unsigned int i;
  unsigned int block;
  unsigned int n_blocks;
  unsigned int blocksize;
  unsigned int id;
  double total_mass = 0;
  IntList *ids;
  if (!sortPartCfg()) {
    runtimeErrorMsg() << "could not sort partCfg";
    return -1;
  }
  ids = (IntList *)container;
  n_blocks = n / 3;
  blocksize = ids->n / n_blocks;
  for (block = 0; block < n_blocks; block++) {
    total_mass = 0;
    for (i = 0; i < blocksize; i++) {
      id = ids->e[block * blocksize + i];
      if (ids->e[i] >= n_part)
        return 1;
      A[3 * block + 0] += (partCfg[id]).p.mass * partCfg[id].r.p[0];
      A[3 * block + 1] += (partCfg[id]).p.mass * partCfg[id].r.p[1];
      A[3 * block + 2] += (partCfg[id]).p.mass * partCfg[id].r.p[2];
      total_mass += (partCfg[ids->e[i]]).p.mass;
    }
    A[3 * block + 0] /= total_mass;
    A[3 * block + 1] /= total_mass;
    A[3 * block + 2] /= total_mass;
  }
  return 0;
}
