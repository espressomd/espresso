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
int ObservableRadialDensityProfile::actual_calculate(PartCfg &partCfg) {
  double *A = last_value;
  int binr, binphi, binz;
  double ppos[3];
  double r, phi, z;
  int img[3];
  double bin_volume;
  IntList *ids;

  if (!sortPartCfg()) {
    runtimeErrorMsg() << "could not sort partCfg";
    return -1;
  }
  radial_profile_data *pdata;
  pdata = (radial_profile_data *)container;
  ids = pdata->id_list;
  double rbinsize = (pdata->maxr - pdata->minr) / pdata->rbins;
  double phibinsize = (pdata->maxphi - pdata->minphi) / pdata->phibins;
  double zbinsize = (pdata->maxz - pdata->minz) / pdata->zbins;

  for (int i = 0; i < n; i++) {
    A[i] = 0;
  }
  for (int i = 0; i < ids->n; i++) {
    if (ids->e[i] >= n_part)
      return 1;
    /* We use folded coordinates here */
    memmove(ppos, partCfg[ids->e[i]].r.p, 3 * sizeof(double));
    memmove(img, partCfg[ids->e[i]].l.i, 3 * sizeof(int));
    fold_position(ppos, img);
    transform_to_cylinder_coordinates(ppos[0] - pdata->center[0],
                                      ppos[1] - pdata->center[1],
                                      ppos[2] - pdata->center[2], &r, &phi, &z);
    // printf("%f %f %f %f %f %f\n", ppos[0], ppos[1], ppos[2],
    // r*cos(phi)+pdata->center[0], r*sin(phi)+pdata->center[1],
    // z+pdata->center[2]);
    binr = (int)floor((r - pdata->minr) / rbinsize);
    binphi = (int)floor((phi - pdata->minphi) / phibinsize);
    binz = (int)floor((z - pdata->minz) / zbinsize);

    if (binr >= 0 && binr < pdata->rbins && binphi >= 0 &&
        binphi < pdata->phibins && binz >= 0 && binz < pdata->zbins) {
      bin_volume = PI * ((pdata->minr + (binr + 1) * rbinsize) *
                             (pdata->minr + (binr + 1) * rbinsize) -
                         (pdata->minr + (binr)*rbinsize) *
                             (pdata->minr + (binr)*rbinsize)) *
                   zbinsize * phibinsize / 2 / PI;
      A[binr * pdata->phibins * pdata->zbins + binphi * pdata->zbins + binz] +=
          1. / bin_volume;
    }
  }
  return 0;
}
