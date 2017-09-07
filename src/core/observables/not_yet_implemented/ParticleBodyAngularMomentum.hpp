int ObservableParticleBodyAngularMomentum::actual_calculate(PartCfg & partCfg) {
  double* A = last_value;
  IntList* ids;
  if (!sortPartCfg()) {
      runtimeErrorMsg() <<"could not sort partCfg";
    return -1;
  }
  ids=(IntList*) container;
  for ( int i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_part)
      return 1;

#ifdef ROTATION

    A[3*i + 0] = partCfg[ids->e[i]].m.omega[0];
    A[3*i + 1] = partCfg[ids->e[i]].m.omega[1];
    A[3*i + 2] = partCfg[ids->e[i]].m.omega[2];

#else

    A[3*i + 0] = 0.0;
    A[3*i + 1] = 0.0;
    A[3*i + 2] = 0.0;

#endif

  }
  return 0;
}
