#ifdef DIPOLES
int ObservableComDipoleMoment::actual_calculate(PartCfg & partCfg) {
  double* A = last_value;
  double d[3] = {0. , 0., 0. } ;
  IntList* ids;
  if (!sortPartCfg()) {
      runtimeErrorMsg() <<"could not sort partCfg";
    return -1;
  }
  ids=(IntList*) container;
  for (int i = 0; i<ids->n; i++ ) {
    if (ids->e[i] > n_part)
      return 1;
    d[0] += partCfg[ids->e[i]].r.dip[0];
    d[1] += partCfg[ids->e[i]].r.dip[1];
    d[2] += partCfg[ids->e[i]].r.dip[2];
  }
  A[0]=d[0];
  A[1]=d[1];
  A[2]=d[2];
  return 0;
}
#endif
