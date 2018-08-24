int ObservablePolymerKDistribution::actual_calculate(PartCfg & partCfg) {
	if (!sortPartCfg()) {
		ostringstream errtxt;
		errtxt << "{094 could not sort partCfg}";
		runtimeError(errtxt);
		return -1;
	}
	double* A = last_value;
	for (int i = 0; i < n; i++) {
		A[i] = 0.0;
	}
	k_dist_data *data = (k_dist_data *) container;
	IntList *ids = data->id_list;
	int npoly = data->npoly;
	int poly_len = data-> poly_len;
	int k = data->k;
	double dist_vec[3];
	double dist;
	double r_min = data->r_min;
	double r_max = data->r_max;
	int n_bins   = data->n_bins;
	int bin_id =0;
	double bin_size = (r_max - r_min) / n_bins;
	int number_of_pairs =(int) ( floor(poly_len / (k +1 ) ) + ( (poly_len - 1)%(k+1) == 0 ? 1 : 0 ) - 1);
	for (int i = 0; i < npoly; i++) 
		for (int j = 0; j < poly_len - k; j+= k) {
			get_mi_vector(dist_vec, partCfg[ids->e[i*poly_len + j]].r.p, partCfg[ids->e[i*poly_len + j + k]].r.p);
			dist = normr(dist_vec);
			bin_id = (int) floor( (dist - r_min)/bin_size );
			if (bin_id < n_bins && bin_id >= 0) {
				A[bin_id] += 1.0/(npoly * number_of_pairs);
			}
		}
	return 0;
}	
