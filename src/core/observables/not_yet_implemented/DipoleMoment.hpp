int ObservableDipoleMoment::actual_calculate(PartCfg & partCfg) {
  double* A = last_value;
  double charge;
  double j[3] = {0. , 0., 0. } ;
  IntList* ids;
  if (!sortPartCfg()) {
      runtimeErrorMsg() <<"could not sort partCfg";
    return -1;
  }
  ids=(IntList*) container;
  for (int i = 0; i<ids->n; i++ ) {
    if (ids->e[i] > n_part)
      return 1;
    charge = partCfg[ids->e[i]].p.q;
    j[0] += charge * partCfg[ids->e[i]].r.p[0];
    j[1] += charge * partCfg[ids->e[i]].r.p[1];
    j[2] += charge * partCfg[ids->e[i]].r.p[2];
  }
  A[0]=j[0];
  A[1]=j[1];
  A[2]=j[2];
  return 0;
}
