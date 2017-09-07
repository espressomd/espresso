int ObservableComForce::actual_calculate(PartCfg & partCfg) {
  double* A = last_value;
  double f_com[3] = { 0. , 0., 0. } ;
  IntList* ids;
  if (!sortPartCfg()) {
      runtimeErrorMsg() <<"could not sort partCfg";
    return -1;
  }
  ids=(IntList*) container;
  for (int i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_part)
      return 1;
    f_com[0] += partCfg[ids->e[i]].f.f[0]/time_step/time_step*2;
    f_com[1] += partCfg[ids->e[i]].f.f[1]/time_step/time_step*2;
    f_com[2] += partCfg[ids->e[i]].f.f[2]/time_step/time_step*2;
  }
  A[0]=f_com[0];
  A[1]=f_com[1];
  A[2]=f_com[2];
  return 0;
}
