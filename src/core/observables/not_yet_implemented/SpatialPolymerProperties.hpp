int ObbservableSpatialPolymerProperties::actual_calculate(PartCfg & partCfg) {
    if (!sortPartCfg()) {
      ostringstream errtxt ; //= runtime_error(128);
	  errtxt << "{094 could not sort partCfg}";
	  runtimeError(errtxt);
      return -1;
    }
	double* A = last_value;
	int poly_len = n;
	for (int i = 0; i<poly_len; i++ )
		A[i]=0.0;

#ifdef ELECTROSTATICS
	spatial_polym_data *p_data = (spatial_polym_data *) container;
	IntList *ids=p_data->id_list;
	for (int i = 0; i<ids->n; i++){
		A[i%poly_len] += partCfg[ids->e[i]].p.q/( (double) p_data->npoly);
	}
#endif
	return 0;
}
