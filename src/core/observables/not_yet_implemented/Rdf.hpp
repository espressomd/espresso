int ObservableRdf::actual_calculate(PartCfg & partCfg) {
  if (!sortPartCfg()) {
    return 1;
  }
  double * last = last_value;
  rdf_profile_data * rdf_data = (rdf_profile_data *) container;
  calc_rdf(rdf_data->p1_types, rdf_data->n_p1,
	   rdf_data->p2_types, rdf_data->n_p2,
	   rdf_data->r_min,    rdf_data->r_max,
	   rdf_data->r_bins,   last);
  return 0;
}
