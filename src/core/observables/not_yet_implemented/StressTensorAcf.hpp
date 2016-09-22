int ObservableStressTensorAcfObs::actual_calculate() {
  double* A = last_value;
  double stress_tensor[9];
  if (!sortPartCfg()) {
      runtimeErrorMsg() <<"could not sort partCfg";
    return -1;
  }
  observable_compute_stress_tensor(1,stress_tensor,9);
  A[0]=stress_tensor[1];
  A[1]=stress_tensor[5];
  A[2]=stress_tensor[6];
  A[3]=stress_tensor[0]-stress_tensor[4];
  A[4]=stress_tensor[0]-stress_tensor[8];
  A[5]=stress_tensor[4]-stress_tensor[8];
  return 0;
}
