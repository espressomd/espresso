int ObservableParticleAngularMomentum::actual_calculate() {
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

    double RMat[9];
    double omega[3];
    define_rotation_matrix(&partCfg[ids->e[i]], RMat);
    omega[0] = RMat[0 + 3*0]*partCfg[ids->e[i]].m.omega[0] + RMat[1 + 3*0]*partCfg[ids->e[i]].m.omega[1] + RMat[2 + 3*0]*partCfg[ids->e[i]].m.omega[2];
    omega[1] = RMat[0 + 3*1]*partCfg[ids->e[i]].m.omega[0] + RMat[1 + 3*1]*partCfg[ids->e[i]].m.omega[1] + RMat[2 + 3*1]*partCfg[ids->e[i]].m.omega[2];
    omega[2] = RMat[0 + 3*2]*partCfg[ids->e[i]].m.omega[0] + RMat[1 + 3*2]*partCfg[ids->e[i]].m.omega[1] + RMat[2 + 3*2]*partCfg[ids->e[i]].m.omega[2];

    A[3*i + 0] = omega[0];
    A[3*i + 1] = omega[1];
    A[3*i + 2] = omega[2];

#else

    A[3*i + 0] = 0.0;
    A[3*i + 1] = 0.0;
    A[3*i + 2] = 0.0;

#endif

  }
  return 0;
}
