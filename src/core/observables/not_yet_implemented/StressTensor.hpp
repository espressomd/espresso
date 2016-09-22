int ObservableStressTensor::actual_calculate() {
  if (!sortPartCfg()) {
      runtimeErrorMsg() <<"could not sort partCfg";
    return -1;
  }
  observable_compute_stress_tensor(1,last_value,n);
  return 0;
}


