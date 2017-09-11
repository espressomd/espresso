int ObservableComVelocity::actual_calculate(PartCfg & partCfg) {
  double* A = last_value;
  double v_com[3] = { 0. , 0., 0. } ;
  double total_mass = 0;
  IntList* ids;
  if (!sortPartCfg()) {
      runtimeErrorMsg() <<"could not sort partCfg";
    return -1;
  }
  ids=(IntList*) container;
  for (int i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_part)
      return 1;
    v_com[0] += (partCfg[ids->e[i]]).p.mass*partCfg[ids->e[i]].m.v[0]/time_step;
    v_com[1] += (partCfg[ids->e[i]]).p.mass*partCfg[ids->e[i]].m.v[1]/time_step;
    v_com[2] += (partCfg[ids->e[i]]).p.mass*partCfg[ids->e[i]].m.v[2]/time_step;
    total_mass += (partCfg[ids->e[i]]).p.mass;
  }
  A[0]=v_com[0]/total_mass;
  A[1]=v_com[1]/total_mass;
  A[2]=v_com[2]/total_mass;
  return 0;
}
