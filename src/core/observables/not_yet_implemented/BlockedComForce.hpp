int ObservableBlockedComForce::actual_calculate(PartCfg & partCfg) {
  double* A = last_value;
  unsigned int i;
  unsigned int block;
  unsigned int n_blocks;
  unsigned int blocksize;
  unsigned int id;
  IntList* ids;
  if (!sortPartCfg()) {
      runtimeErrorMsg() <<"could not sort partCfg";
    return -1;
  }
  ids=(IntList*) container;
  n_blocks=n/3; 
  blocksize=ids->n/n_blocks;
  for ( block = 0; block < n_blocks; block++ ) {
    for ( i = 0; i < blocksize; i++ ) {
      id = ids->e[block*blocksize+i];
      if (ids->e[i] >= n_part)
        return 1;
      A[3*block+0] +=  partCfg[id].f.f[0]/time_step/time_step*2;
      A[3*block+1] +=  partCfg[id].f.f[1]/time_step/time_step*2;
      A[3*block+2] +=  partCfg[id].f.f[2]/time_step/time_step*2;
    }
  }
  return 0;
}
