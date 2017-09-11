int ObservableRadialDensityDistribution::actual_calculate(PartCfg & partCfg) {
  if (!sortPartCfg()) {
    ostringstream errtxt; // = runtime_error(128);
    errtxt << "{094 could not sort partCfg} ";
	runtimeError(errtxt);
    return -1;
  }

  radial_density_data *r_data = (radial_density_data *) container;
  IntList *ids;  
  if ( GC_init && Type_array_init ) {
	  ids = (IntList *) Utils::malloc(sizeof(IntList));

	  //using the grandcanonical scheme, always update the particle id list
	  ids->e = (int *) Utils::malloc(sizeof(int)*type_array[Index.type[r_data->type]].max_entry);
	  memmove(ids->e, type_array[Index.type[r_data->type]].id_list, type_array[Index.type[r_data->type]].max_entry*sizeof(int));
	  ids->n = type_array[Index.type[r_data->type]].max_entry;
	  ids->max = type_array[Index.type[r_data->type]].cur_size;
  } else { 
	  ids = r_data->id_list;
  }
  //r_data->id_list = ids;
  //ids = r_data->id_list;

  double* A = last_value;
  int n_A   = n;
  double start_point[3];
  double end_point[3];
  int image_box[3];
  if ( r_data->id_flag ) {
	  // Using particle_ids to specify the start and endpoints
	  memmove(start_point, partCfg[r_data->start_point_id].r.p, 3*sizeof(double));
	  memmove(image_box, partCfg[r_data->start_point_id].l.i, 3*sizeof(int));
	  unfold_position(start_point, image_box);
	  memmove(end_point, partCfg[r_data->end_point_id].r.p, 3*sizeof(double));
	  memmove(image_box, partCfg[r_data->end_point_id].l.i, 3*sizeof(int));
	  unfold_position(end_point, image_box);
  } else {
	  memmove(start_point, r_data->start_point, 3*sizeof(double));
	  memmove(end_point, r_data->end_point, 3*sizeof(double));
  }

  double *bin_volume = (double *) Utils::malloc(sizeof(double)*r_data->rbins);
 
  double part_pos[3];
  double AB[3];		// normalized normal vector pointing to start point
  double BA[3];		// ...									end point
  double Apart[3];	// vector difference start_point - part_pos
  double norm_Apart;
  double normAB;
  double dist;
  double angle;
  double r_dist;
  double tmp[3];
  double tmp2[3];
  char periodic_points = 0;
  int frac = 0;
  double epsilon = 1e-9;

  double rbin_size = (r_data->maxr - r_data->minr)/(r_data->rbins);
  int rbin_id;

  if ( r_data->id_flag ) {
	  // in the case of particles use minimal image convention for axis definition
	  get_mi_vector(AB, start_point, end_point);
	  get_mi_vector(BA, end_point, start_point);
	  if ( normr(AB) < epsilon )  {
		  //start and end point are periodic images of one another
		  periodic_points = 1;
	  }
	  end_point[0] = start_point[0] - AB[0];
	  end_point[1] = start_point[1] - AB[1];
	  end_point[2] = start_point[2] - AB[2];
  } else {
	  // otherwise use the given absolute positions but check whether those points are too close
	  get_mi_vector(AB, start_point, end_point);
	  if ( normr(AB) < epsilon )  {
		  //start and end point are periodic images of one another
		  periodic_points = 1;
	  }
	  vecsub(start_point, end_point, AB);
	  vecsub(end_point, start_point, BA);
  }

  //printf("id: %d, lenAB: %1.3e, lenBA: %1.3e, %f %f %f, %f %f %f\n", r_data->id_flag, normr(AB), normr(BA), start_point[0], start_point[1], start_point[2], end_point[0], end_point[1], end_point[2]);
  normAB = normr(AB); 
  for (int i=0; i<r_data->rbins; i++) {
	  bin_volume[i] = normAB*( pow(r_data->minr + (i+1)*rbin_size, 2) - pow(r_data->minr + i*rbin_size, 2) ) * PI;
  }

  unit_vector(AB, AB);
  unit_vector(BA, BA);

  for (int i=0; i<n_A; i++)
	  A[i] = 0.0;

  for (int i=0; i < ids->n; i++){
	  part_pos[0]=partCfg[ids->e[i]].r.p[0];
	  part_pos[1]=partCfg[ids->e[i]].r.p[1];
	  part_pos[2]=partCfg[ids->e[i]].r.p[2];
	  // that might seem weird, but this ensures that the used particle position is the
	  // closest image to one of the two points that define the axis
	  get_mi_vector(tmp, part_pos, start_point);
	  get_mi_vector(tmp2, part_pos, end_point);

	  if ( periodic_points ) { 
		  // the axis spans the whole box, so any particle is in 
		  // scope
		  dist = -1;
	  } else if ( normr(tmp) < normr(tmp2) ) {
		  part_pos[0] = start_point[0] + tmp[0];
		  part_pos[1] = start_point[1] + tmp[1];
		  part_pos[2] = start_point[2] + tmp[2];
		  dist = scalar(AB, tmp);
	  } else {
		  part_pos[0] = end_point[0] + tmp2[0];
		  part_pos[1] = end_point[1] + tmp2[1];
		  part_pos[2] = end_point[2] + tmp2[2];
		  dist =  scalar(BA, tmp2);
	  }

	  if (dist > 0 ) { 

		  continue;
	  }

	  // particle in scope
	  // now calculate the distance from the given axis using simple geometry
	  vecsub(start_point, part_pos, Apart);

	  norm_Apart = normr(Apart);
	  angle = acos(scalar(Apart, AB)/norm_Apart);
	  r_dist= sin(angle)*norm_Apart;

	  rbin_id = (int) floor( (r_dist - r_data->minr)/rbin_size );
	  if ( rbin_id >= r_data->rbins || rbin_id < 0 )
		  continue;

	  A[rbin_id]+=1.0/bin_volume[rbin_id];
	  frac++;
  }
//  printf("fraction of parts: %d %d\n", frac, ids->n);
  free(bin_volume);
  if ( GC_init && Type_array_init ) {
	  free(ids->e);
	  free(ids);
	}
  return 0;
}
