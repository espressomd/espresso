int ObservableRadialFluxDensityProfile::actual_calculate(PartCfg & partCfg) {
  double* A = last_value;
  int binr, binphi, binz;
  double ppos[3];
  double unfolded_ppos[3];
  double r, phi, z;
  int img[3];
  double bin_volume;
  IntList* ids;

  if (!sortPartCfg()) {
      runtimeErrorMsg() <<"could not sort partCfg";
    return -1;
  }
  radial_profile_data* pdata;
  pdata=(radial_profile_data*) container;
  ids=pdata->id_list;
  double rbinsize=(pdata->maxr - pdata->minr)/pdata->rbins;
  double phibinsize=(pdata->maxphi - pdata->minphi)/pdata->phibins;
  double zbinsize=(pdata->maxz - pdata->minz)/pdata->zbins;
  double v[3];
  double v_r, v_phi, v_z;

  if (last_update==sim_time) {
    return ES_ERROR;
  }
    
  for (int i = 0; i< n; i++ ) {
    A[i]=0;
  }
  double* old_positions=(double*) pdata->container;
  if (old_positions[0] == CONST_UNITITIALIZED) {
    for (int i = 0; i<ids->n; i++ ) {
      memmove(unfolded_ppos, partCfg[ids->e[i]].r.p, 3*sizeof(double));
      memmove(img, partCfg[ids->e[i]].l.i, 3*sizeof(int));
      unfold_position(unfolded_ppos, img);
      old_positions[3*i+0]=unfolded_ppos[0];
      old_positions[3*i+1]=unfolded_ppos[1];
      old_positions[3*i+2]=unfolded_ppos[2];
    }
    return 0;
  }
  for (int i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_part)
      return 1;
/* We use folded coordinates here */
    memmove(unfolded_ppos, partCfg[ids->e[i]].r.p, 3*sizeof(double));
    memmove(img, partCfg[ids->e[i]].l.i, 3*sizeof(int));
    unfold_position(unfolded_ppos, img);
    v[0]=(unfolded_ppos[0] - old_positions[3*i+0]);
    v[1]=(unfolded_ppos[1] - old_positions[3*i+1]);
    v[2]=(unfolded_ppos[2] - old_positions[3*i+2]);
    memmove(ppos, partCfg[ids->e[i]].r.p, 3*sizeof(double));
    memmove(img, partCfg[ids->e[i]].l.i, 3*sizeof(int));
    fold_position(ppos, img);
    // The position of the particle is by definition the middle of old and new position
    ppos[0]+=0.5*v[0]; ppos[1]+=0.5*v[1]; ppos[2]+=0.5*v[2];
    fold_position(ppos, img);
    v[0]/=(sim_time - last_update);
    v[1]/=(sim_time - last_update);
    v[2]/=(sim_time - last_update);
    if (i==0) {
//      printf("(%3.4f) %f %f %f\n", sim_time-last_update, v[2], partCfg[ids->e[i]].m.v[2]/time_step,v[2]* partCfg[ids->e[i]].m.v[2]/time_step/time_step);
//      printf("(%3.3f) %f %f", sim_time, old_positions[3*i+2], unfolded_ppos[2]);
    }
    old_positions[3*i+0]=unfolded_ppos[0];
    old_positions[3*i+1]=unfolded_ppos[1];
    old_positions[3*i+2]=unfolded_ppos[2];
    transform_to_cylinder_coordinates(ppos[0]-pdata->center[0], ppos[1]-pdata->center[1], ppos[2]-pdata->center[2], &r, &phi, &z);
    binr  =(int)floor((r-pdata->minr)/rbinsize);
    binphi=(int)floor((phi-pdata->minphi)/phibinsize);
    binz  =(int)floor((z-pdata->minz)/zbinsize);

    if (binr>=0 && binr < pdata->rbins && binphi>=0 && binphi < pdata->phibins && binz>=0 && binz < pdata->zbins) {
      bin_volume=PI*((pdata->minr+(binr+1)*rbinsize)*(pdata->minr+(binr+1)*rbinsize) - (pdata->minr+(binr)*rbinsize)*(pdata->minr+(binr)*rbinsize)) *zbinsize * phibinsize/2/PI;
      v_r = 1/r*((ppos[0]-pdata->center[0])*v[0] + (ppos[1]-pdata->center[1])*v[1]); 
      v_phi = 1/r/r*((ppos[0]-pdata->center[0])*v[1]-(ppos[1]-pdata->center[1])*v[0]);
      v_z = v[2];
      A[3*(binr*pdata->phibins*pdata->zbins + binphi*pdata->zbins + binz) + 0] += v_r/bin_volume;
      A[3*(binr*pdata->phibins*pdata->zbins + binphi*pdata->zbins + binz) + 1] += v_phi/bin_volume;
      A[3*(binr*pdata->phibins*pdata->zbins + binphi*pdata->zbins + binz) + 2] += v_z/bin_volume;
    } 
  }
  return 0;
}
