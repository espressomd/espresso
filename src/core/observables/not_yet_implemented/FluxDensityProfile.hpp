int ObservableFluxDensityProfile::actual_calculate(PartCfg & partCfg) {
  double* A = last_value;
  int binx, biny, binz;
  double ppos[3];
  double x, y, z;
  int img[3];
  double bin_volume;
  IntList* ids;

  if (!sortPartCfg()) {
      runtimeErrorMsg() <<"could not sort partCfg";
    return -1;
  }
  profile_data* pdata;
  pdata=(profile_data*) container;
  ids=pdata->id_list;
  double xbinsize=(pdata->maxx - pdata->minx)/pdata->xbins;
  double ybinsize=(pdata->maxy - pdata->miny)/pdata->ybins;
  double zbinsize=(pdata->maxz - pdata->minz)/pdata->zbins;
  double v[3];
  double v_x, v_y, v_z;
    
  for (int i = 0; i< n; i++ ) {
    A[i]=0;
  }

  for (int i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_part)
      return 1;
    /* We use folded coordinates here */
    v[0]=partCfg[ids->e[i]].m.v[0]/time_step;
    v[1]=partCfg[ids->e[i]].m.v[1]/time_step;
    v[2]=partCfg[ids->e[i]].m.v[2]/time_step;
    memmove(ppos, partCfg[ids->e[i]].r.p, 3*sizeof(double));
    memmove(img, partCfg[ids->e[i]].l.i, 3*sizeof(int));
    fold_position(ppos, img);
    // The position of the particle is by definition the middle of old and new position
  
    x=ppos[0];
    y=ppos[1];
    z=ppos[2];

    binx  =(int)floor((x-pdata->minx)/xbinsize);
    biny  =(int)floor((y-pdata->miny)/ybinsize);
    binz  =(int)floor((z-pdata->minz)/zbinsize);


    if (binx>=0 && binx < pdata->xbins && biny>=0 && biny < pdata->ybins && binz>=0 && binz < pdata->zbins) {
      bin_volume=xbinsize*ybinsize*zbinsize;
      v_x=v[0];
      v_y=v[1];
      v_z=v[2];
      A[3*(binx*pdata->ybins*pdata->zbins + biny*pdata->zbins + binz) + 0] += v_x/bin_volume;
      A[3*(binx*pdata->ybins*pdata->zbins + biny*pdata->zbins + binz) + 1] += v_y/bin_volume;
      A[3*(binx*pdata->ybins*pdata->zbins + biny*pdata->zbins + binz) + 2] += v_z/bin_volume;
    } 
  }
  return 0;
}
