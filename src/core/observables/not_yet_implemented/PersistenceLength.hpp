int ObservablePersistenceLength::actual_calculate(PartCfg & partCfg) {
    if (!sortPartCfg()) {
      ostringstream errtxt; // = runtimeError(128);
      errtxt << "{094 could not sort partCfg}";
	  runtimeError(errtxt);
      return -1;
    }
	double* A = last_value;
	spatial_polym_data *p_data = (spatial_polym_data *) container;
	IntList *ids=p_data->id_list;

	double v1[3];
	double v2[3];
	double abs1, abs2;
	int num_parts = ids->n;
	int cut_off = p_data->cut_off;
	int n_A = n;

	for (int i = 0; i<n_A; i++ )
		A[i]=0.0;

	for (int i = 0; i<n_A; i++) {
		// i : distance between polymer segments
		for ( int j=cut_off; j<num_parts - i -cut_off; j++ ) {
			vecsub(partCfg[ids->e[j]].r.p, partCfg[ids->e[j+1]].r.p, v1);
			abs1 = normr(v1);
			vecsub(partCfg[ids->e[i+j]].r.p, partCfg[ids->e[i+j+1]].r.p, v2);
			abs2 = normr(v2);
			A[i] += scalar(v1, v2)/(abs1*abs2 * (double) (num_parts - i - 2*cut_off) );
		}
	}
	return 0;
}	
