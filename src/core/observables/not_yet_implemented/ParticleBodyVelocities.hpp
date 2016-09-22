int ObservableParticleBodyVelocities::actual_calculate() {
  double* A = last_value;
  IntList* ids;
  if (!sortPartCfg()) {
      runtimeErrorMsg() <<"could not sort partCfg";
    return -1;
  }
  ids=(IntList*) container;
  for (int i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_part)
      return 1;


#ifdef ROTATION

    double RMat[9];
    double vel_lab[3];
    double vel_body[3];

    vel_lab[0] = partCfg[ids->e[i]].m.v[0]/time_step;
    vel_lab[1] = partCfg[ids->e[i]].m.v[1]/time_step;
    vel_lab[2] = partCfg[ids->e[i]].m.v[2]/time_step;
    define_rotation_matrix(&partCfg[ids->e[i]], RMat);

    vel_body[0] = RMat[0 + 3*0]*vel_lab[0] + RMat[0 + 3*1]*vel_lab[1] + RMat[0 + 3*2]*vel_lab[2];
    vel_body[1] = RMat[1 + 3*0]*vel_lab[0] + RMat[1 + 3*1]*vel_lab[1] + RMat[1 + 3*2]*vel_lab[2];
    vel_body[2] = RMat[2 + 3*0]*vel_lab[0] + RMat[2 + 3*1]*vel_lab[1] + RMat[2 + 3*2]*vel_lab[2];

    A[3*i + 0] = vel_body[0];
    A[3*i + 1] = vel_body[1];
    A[3*i + 2] = vel_body[2];

#else

    A[3*i + 0] = 0.0;
    A[3*i + 1] = 0.0;
    A[3*i + 2] = 0.0;

#endif

  }
  return 0;
}
