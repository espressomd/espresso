/* this is included in statistics.c */

// list of the currently specified box boundaries
DoubleList boundaries = { NULL, 0 };
// number of boxes
int n_part_in_bin = 0;
// the boxes with the particle identities
IntList *part_in_bin = NULL;

static void wall_sort_particles()
{
  // 1. reallocate the boxes for the particle identities
  // free old boxes
  for (int i = 0; i < n_part_in_bin; ++i) {
    realloc_intlist(&part_in_bin[i], part_in_bin[i].n = 0);
  }
  part_in_bin = realloc(part_in_bin, (boundaries.n - 1)*sizeof(IntList));
  // initialize new ones
  for (int i = n_part_in_bin; i < boundaries.n-1; ++i) {
    init_intlist(&part_in_bin[i]);
  }
  n_part_in_bin = boundaries.n-1;

  // 2. for each particle, find the box and put its
  // identity there
  for(int i=0; i<n_total_particles; i++) {
    double x = partCfg[i].r.p[0];
    // ignore particles outside the boundaries
    if (x < boundaries.e[0] || x > boundaries.e[boundaries.n-1])
      continue;
    // simple bisection on the particle's x-coordinate
    int s = 0;
    int e = boundaries.n - 1;
    while (e - s > 1) {
      int c = (e + s)/2;
      if (x >= boundaries.e[c]) s = c; else e = c;
    }
    // and add the particle to the resulting list
    realloc_grained_intlist(&part_in_bin[s], part_in_bin[s].n + 1, 8);
    part_in_bin[s].e[part_in_bin[s].n++] = i;
  }
}

static void calc_wallmsdyz(double *g, int bin)
{
  // loop over all stored configurations
  // and calculate the MSD with respect to the current configuration
  // MSD of configuration with itself is always 0
  g[0] = 0;
  for(int k = 1; k < n_configs; k++) {
    g[k] = 0.0;
    // loop over all particles in the specified bin and add up MSD
    for (int i = 0; i < part_in_bin[bin].n; ++i) {
      int p = part_in_bin[bin].e[i];
      g[k] +=
	+ SQR(configs[n_configs-1][3*p + 1]-configs[n_configs-1-k][3*p + 1])
	+ SQR(configs[n_configs-1][3*p + 2]-configs[n_configs-1-k][3*p + 2]);
    }
    // normalize
    g[k] /= part_in_bin[bin].n;
  }
}

static void calc_wallmsdx(double *g, int bin)
{
  // see calc_wallmsdyz, just for x
  g[0] = 0;
  for(int k = 1; k < n_configs; k++) {
    g[k] = 0.0;
    for (int i = 0; i < part_in_bin[bin].n; ++i) {
      int p = part_in_bin[bin].e[i];
      g[k] += SQR(configs[n_configs-1][3*p]-configs[n_configs-1-k][3*p]);
    }
    g[k] /= part_in_bin[bin].n;
  }
}

static void calc_wallbondyz(double *g, int bin, double rclocal, double rmax, int rbins)
{
  // char buffer[TCL_INTEGER_SPACE + TCL_DOUBLE_SPACE + 2];
  double rclocal2=rclocal*rclocal, neigh[part_in_bin[bin].n];
  double sum_sai_r=0.0, sum_sai_m=0.0, sai_r[part_in_bin[bin].n], sai_m[part_in_bin[bin].n];
  int grdf[rbins];

  for (int i = 0; i < part_in_bin[bin].n; ++i) {
    int nb=0;
    for (int j = 0; j < part_in_bin[bin].n; ++j) {
      if (j==i) continue;
      int p = part_in_bin[bin].e[i];
      int q = part_in_bin[bin].e[j];
      //printf("%d   %d\n",p,q);
      // minimum image vector between the two particles
      double diff[3];
      get_mi_vector(diff, partCfg[p].r.p, partCfg[q].r.p);

      double dist = SQR(diff[1]) + SQR(diff[2]);
      if (dist <rclocal2) {
	neigh[nb] = atan2(diff[2], diff[1]);
	nb++;  // counting neighbours
      }

      // end of inner particle 
    }
    printf ("part %d neibours %d\n",i, nb);
   
    //counting angle and sai
    if (nb>3) {
      double sinus=0.0;
      double cosinus=0.0;
      for (int l=0; l<nb; l++){
	
	double teta= 6*neigh[l];
	cosinus = cosinus+cos (teta);
	sinus =sinus+ sin (teta);
      }
      sai_r[i] = cosinus/nb;
      sai_m[i] = sinus/nb;
      //printf("bin %d   nb %d  i %d  sin %e cos %e\n", bin,nb, i, sinus, cosinus);
    } else sai_r[i]= sai_m[i]=0.0;
    
    sum_sai_r = sum_sai_r + sai_r[i];
    sum_sai_m = sum_sai_m + sai_m[i];
    //printf("bin %d  part %d  sai %e  saim %e\n", bin,i, sum_sai_r, sum_sai_m);
  
    // end of outer particle
  }

  double tot_sai_m= sqrt((sum_sai_r * sum_sai_r) + (sum_sai_m * sum_sai_m))/part_in_bin[bin].n;
  g[rbins]=tot_sai_m;
  //printf("bin %d    sai %e\n", bin, tot_sai_m);

  double bin_width = rmax / (double)rbins;
  double inv_bin_width = 1.0 / bin_width;

  // initially set all counts to 0
  for(int i = 0; i < rbins; ++i) {
    g[i] = 0.0;
    grdf[i]=0;
  }

  // number of samples gathered
  //int cnt = 0;

  // loop over all distinct particle pairs in the bin
  //printf("here1");
  for (int i = 0; i < part_in_bin[bin].n; ++i) {
    for (int j = i + 1; j < part_in_bin[bin].n; ++j) {
      int p = part_in_bin[bin].e[i];
      int q = part_in_bin[bin].e[j];

      // minimum image vector between the two particles
      double diff[3];
      get_mi_vector(diff, partCfg[p].r.p, partCfg[q].r.p);
      double dist = sqrt(SQR(diff[1]) + SQR(diff[2]));

      // binning
      int ind = (int) (dist*inv_bin_width);
      //printf("index %d\n", ind);
      if (ind >=0 && ind < rbins) {
	double sai_j1j2 = ((sai_r[i] * sai_r[j]) + (sai_m[i] * sai_m[j]));
	g[ind] += sai_j1j2;
	grdf[ind]++;
	// printf("%e\n",g[ind]);
      }
      //cnt++;
    }
  }

  // normalization
  // the width of the slice we look at
  // double width = boundaries.e[bin+1] - boundaries.e[bin];
  // average density of particles
  //double av_density = cnt/(width*box_l[1]*box_l[2]);

  for(int i = 0; i < rbins; ++i) {
    // normalize
    if (grdf[i]==0) {
      g[i]=0;
      continue;
    }
    g[i] /= grdf[i];
    //printf(" %e  rdf %d\n",g[i], grdf[i]);
    // Tcl_AppendResult(interp, buffer, (char *)NULL); 
  }
}


static void calc_wallrdfyz(double *g, int bin, double rmin, double rmax, int rbins)
{
  double bin_width = (rmax-rmin) / (double)rbins;
  double inv_bin_width = 1.0 / bin_width;

  // initially set all counts to 0
  for(int i = 0; i < rbins; ++i) g[i] = 0.0;

  // number of samples gathered
  int cnt = 0;

  // loop over all distinct particle pairs in the bin
  for (int i = 0; i < part_in_bin[bin].n; ++i) {
    for (int j = i + 1; j < part_in_bin[bin].n; ++j) {
      int p = part_in_bin[bin].e[i];
      int q = part_in_bin[bin].e[j];

      // minimum image vector between the two particles
      double diff[3];
      get_mi_vector(diff, partCfg[p].r.p, partCfg[q].r.p);
      double dist = sqrt(SQR(diff[1]) + SQR(diff[2]));

      // binning
      int ind = (int) ((dist - rmin)*inv_bin_width);
      if (ind >=0 && ind < rbins) {
        g[ind]++;
      }
      cnt++;
    }
  }

  // normalization
  // the width of the slice we look at
  double width = boundaries.e[bin+1] - boundaries.e[bin];
  // average density of particles
  double av_density = cnt/(width*box_l[1]*box_l[2]);

  for(int i = 0; i < rbins; ++i) {
    double r_in  =     i*bin_width + rmin; 
    double r_out = (i+1)*bin_width + rmin;
    double bin_volume = PI*(r_out*r_out - r_in*r_in)*width;
    // normalize first to density, and then to relative density
    g[i] = (g[i]/bin_volume) / av_density;
  }
}

int parse_wallstuff(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze wallmsd -xy|-z <min> <max>' */
  /******************************************************************************/
  char buffer[TCL_INTEGER_SPACE + TCL_DOUBLE_SPACE + 2];
  DoubleList g;
  int job, bin;
  double rmin, rmax,rclocal;
  int rbins;
  enum { BINS, MX, MYZ, RDFYZ,BONDYZ, PRINT };

  if (argc < 2) {
    Tcl_AppendResult(interp, "expected: analyze wallstuff -bins <binboundaries> | -myz <bin> |-mx <bin> | -rdfyz <bin> <rmin> <rmax> <rdfbins>",
		     (char *)NULL);
    return TCL_ERROR;
  }

  // 1. what do we do?
  if (ARG0_IS_S("-bins")) {
    job = BINS;
  }
  else if (ARG0_IS_S("-mx") && argc == 2) {
    job = MX;
  }
  else if (ARG0_IS_S("-myz") && argc == 2) {
    job = MYZ;
  }
  else if (ARG0_IS_S("-rdfyz") && argc == 5) {
    job = RDFYZ;
  }
  else if (ARG0_IS_S("-bondyz") && argc == 5) {
    job = BONDYZ;
  }
  else if (ARG0_IS_S("-print") && argc == 2) {
    job = PRINT;
  }
  else {
    Tcl_AppendResult(interp, ": analyze wallstuff -bins|-myz|-mx|-rdfyz|-bondyz ...", (char *)NULL);
    return TCL_ERROR;
  }
  
  // 2. parameters
  // 2. a) 1. parameter, bin or boundaries
  switch (job) {
  case BINS:
    realloc_doublelist(&boundaries, boundaries.n = 0);
    if (!ARG_IS_DOUBLELIST(1, boundaries)) {
      return TCL_ERROR;
    }
    if (boundaries.n < 2) {
      return (TCL_ERROR);
    }
    break;
  case MX:
  case MYZ:
  case RDFYZ:
  case BONDYZ:
  case PRINT:
    if (!ARG_IS_I(1, bin)) {
      return (TCL_ERROR);
    }
    if (bin < 0 || bin >= boundaries.n-1) {
      return (TCL_ERROR);
    }
    break;
  }

  // 2. other parameters, only for rdf
  switch (job) {
  case RDFYZ:
    if (!ARG_IS_D(2, rmin)) {
      return (TCL_ERROR);
    }
    if (!ARG_IS_D(3, rmax)) {
      return (TCL_ERROR);
    }
    if (!ARG_IS_I(4, rbins)) {
      return (TCL_ERROR);
    }
    break;
  case BONDYZ:
    if (!ARG_IS_D(2, rclocal)) {
      return (TCL_ERROR);
    }
    if (!ARG_IS_D(3, rmax)) {
      return (TCL_ERROR);
    }
    if (!ARG_IS_I(4, rbins)) {
      return (TCL_ERROR);
    }
    break;
  
  }

  // result double list
  init_doublelist(&g);

  // check that data is there
  switch (job) {
  case BINS:
  case RDFYZ:
    // these cases use partCfg
    updatePartCfg(WITHOUT_BONDS);
    break;
  case BONDYZ:
    // these cases use partCfg
    updatePartCfg(WITHOUT_BONDS);
    break;
  case MX:
  case MYZ:
    // these cases use the positions array
    if (n_configs == 0) {
      Tcl_AppendResult(interp, "no configurations found! Use 'analyze append' to save some!",
		       (char *)NULL);
      return TCL_ERROR;
    }
    break;
  }

  // finally, do what is necessary
  switch (job) {
  case BINS:
    wall_sort_particles();
    break;
  case MX:
    realloc_doublelist(&g, g.n = n_configs);
    calc_wallmsdx(g.e, bin);    
    break;
  case MYZ:
    realloc_doublelist(&g, g.n = n_configs);
    calc_wallmsdyz(g.e, bin);    
    break;
  case RDFYZ:
    realloc_doublelist(&g, g.n = rbins);
    calc_wallrdfyz(g.e, bin, rmin, rmax, rbins);
    break;
  case BONDYZ:
    realloc_doublelist(&g, g.n = rbins+1);
    calc_wallbondyz(g.e, bin, rclocal, rmax, rbins);
    break;
  case PRINT:
    // just write out what wall_sort_particles has put into
    // this bin
    for (int i = 1; i < part_in_bin[bin].n; i++) { 
      sprintf(buffer," %d",part_in_bin[bin].e[i]);
      Tcl_AppendResult(interp, buffer, (char *)NULL); 
    }
    break;
  }

  // print out double results, if any
  if (g.n) {
    sprintf(buffer,"%f",g.e[0]);
    Tcl_AppendResult(interp, buffer, (char *)NULL); 

    for (int i = 1; i < g.n; i++) { 
      sprintf(buffer," %f",g.e[i]);
      Tcl_AppendResult(interp, buffer, (char *)NULL); 
    }
    realloc_doublelist(&g, g.n = 0);
  }

  return (TCL_OK);
}

double histogram[3][200];
double count[3];

int parse_sqr(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze dipol' */
  /******************************************************************************/
  int img[3];
  char buffer[3*(TCL_DOUBLE_SPACE + 1) + 4];

  if (argc > 0 && ARG0_IS_S("-init")) {

    for (int j = 0; j < 3; ++j) {
      count[j] = 0;
      for (int i = 0; i < 200; ++i)
        histogram[j][i] = 0;
    }
    return TCL_OK;
  }

  double volume = 125*box_l[0]*box_l[1]*box_l[2];
  if (argc > 0 && ARG0_IS_S("-retrieve")) {
    for (int i = 0; i < 200; ++i) {
      double r_in = (i*sqrt(3)*box_l[0])/200;
      double r_out = ((i+1)*sqrt(3)*box_l[0])/200;
      double bin_volume = (4.0/3.0) * PI * ((r_out*r_out*r_out) - (r_in*r_in*r_in));
      sprintf(buffer,"%f", 0.5*(r_out + r_in));
      Tcl_AppendResult(interp, "{ ", buffer, (char *)NULL); 
      for (int j = 0; j < 3; j++) {
        double rdf = histogram[j][i] * volume / (bin_volume * count[j]);
        sprintf(buffer," %f", rdf);
        Tcl_AppendResult(interp, buffer, (char *)NULL);
      }
      Tcl_AppendResult(interp, " } ", (char *)NULL); 
    }
    return TCL_OK;
  }

  updatePartCfg(WITHOUT_BONDS);

  for(int i=0; i<n_total_particles; i++) {
    int t1 = partCfg[i].p.type;

    double fpos1[3];
    memcpy(fpos1, partCfg[i].r.p, 3*sizeof(double));
    fold_position(fpos1, img);

    for(int j=0; j < n_total_particles; j++) {
      int sort;
      int t2 =  partCfg[j].p.type;
      double fpos2[3];

      if (t1 == t2) {
        if (t1 == 1)
          sort = 1;
        else
          sort = 0;
      }
      else
        sort = 2;
      
      memcpy(fpos2, partCfg[j].r.p, 3*sizeof(double));
      fold_position(fpos2, img);

      for(double imx = -2.0; imx <= 2.0; imx += 1.0)
        for(double imy = -2.0; imy <= 2.0; imy += 1.0)
          for(double imz = -2.0; imz <= 2.0; imz += 1.0) {
            double rdst2 =
              sqrt(SQR(imx*box_l[0] + fpos2[0]-fpos1[0]) +
                   SQR(imy*box_l[1] + fpos2[1]-fpos1[1]) +
                   SQR(imz*box_l[2] + fpos2[2]-fpos1[2]));

            int bin = (int)(floor(200*rdst2/(sqrt(3)*box_l[0])));
            if (bin < 0) bin = 0;
            if (bin > 200-1) bin = 200-1;

            if (i!=j) histogram[sort][bin]++;
            count[sort]++;
          }
    }
  }

  return TCL_OK;
}

int parse_dipol(Tcl_Interp *interp, int argc, char **argv)
{
#ifdef ELECTROSTATICS
  /* 'analyze dipol' */
  /**************************************************************************/
  char buffer[3*(TCL_DOUBLE_SPACE + 1) + 4];

  updatePartCfg(WITHOUT_BONDS);
  
  double dipole[3], dipfold[3], diff[3];
  
  for (int c = 0; c < 3; ++c) {
    dipole[c] = dipfold[c] = diff[c] = 0;
  }

  for(int i=0; i<n_total_particles; i++) {
    double q = partCfg[i].p.q;
    double pos[3], fpos[3]; int img[3];
    memcpy(pos, partCfg[i].r.p, 3*sizeof(double));
    memcpy(fpos, pos , 3*sizeof(double));
    fold_position(fpos, img);

    for (int c = 0; c < 3; ++c) {
      dipole[c]  += q*pos[c];
      dipfold[c] += q*fpos[c];
    }
#if 0
    for(int j=i+1; j < n_total_particles; j++) {
      double q2 = partCfg[j].p.q;
      double fpos2[3];
      
      memcpy(fpos2, partCfg[j].r.p, 3*sizeof(double));
      fold_position(fpos2, img);

      for (int c = 0; c < 3; ++c) {
	diff[c] += q2*q*fabs(fpos2[c] - fpos[c]);
      }
    }
#endif
  }

  sprintf(buffer,"%f %f %f",
	  dipole[0], dipole[1], dipole[2]);
  
#if 0
  (SQR(dipole[0]) + SQR(dipole[1]) + SQR(dipole[2]),
   SQR(dipfold[0]) + SQR(dipfold[1]) + SQR(dipfold[2]),
   - box_l[0]*diff[0] - box_l[1]*diff[1] - box_l[2]*diff[2]
   - SQR(dipfold[0]) - SQR(dipfold[1]) - SQR(dipfold[2])
   );
#endif

  Tcl_AppendResult(interp, buffer, (char *)NULL); 
#endif  
  return TCL_OK;
}

#define BINS 100
#define RMAX 5.0
static double force[BINS];
static int cnt[BINS];

int parse_pairpot(Tcl_Interp *interp, int argc, char **argv)
{
#ifdef ELECTROSTATICS
  /* 'analyze pairpot' */
  /**************************************************************************/
  char buffer[4*(TCL_DOUBLE_SPACE + 1) + 1];

  updatePartCfg(WITHOUT_BONDS);
  
  if (argc > 0 && ARG0_IS_S("-init")) {
    for (int i = 0; i < BINS; ++i) {
      force[i] = 0.0;
      cnt[i] = 0;
    }
    return TCL_OK;
  }
  else if (argc > 0 && ARG0_IS_S("-get")) {
    for (int i = 0; i < BINS; ++i) {
      if (cnt[i]) {
	sprintf(buffer,"%.8e %.8e ", (0.5 + i)*(RMAX/BINS),
		force[i]/cnt[i]/(0.5*time_step*time_step));
      } else {
	sprintf(buffer,"%.8e * ", (0.5 + i)*(RMAX/BINS));
      }
      Tcl_AppendResult(interp, buffer, (char *)NULL); 
    }
    return TCL_OK;
  }

  for(int i=0; i<n_total_particles; i++) {
    double q = partCfg[i].p.q;

    for(int j=i+1; j < n_total_particles; j++) {
      double q2 = partCfg[j].p.q;
      double diff[3];
      get_mi_vector(diff, partCfg[i].r.p, partCfg[j].r.p);
      double r = 0.0, s = 0.0;
      for (int k = 0; k < 3; ++k) {
	double f = partCfg[i].f.f[k] - partCfg[j].f.f[k];
	s += f*diff[k];
	r += SQR(diff[k]);
      }
      r = sqrt(r);
      int bin = (int)(floor(r*(BINS/RMAX)));
      
      if (bin >= 0 && bin < BINS) {
	force[bin] += s/r/(q*q2);
	cnt[bin]++;
      }
    }
  }
#endif
  return TCL_OK;
}

int parse_current(Tcl_Interp *interp, int argc, char **argv)
{
#ifdef ELECTROSTATICS
  /* 'analyze current' */
  /***************************************************************************/
  char buffer[3*(TCL_DOUBLE_SPACE + 1) + 4];

  updatePartCfg(WITHOUT_BONDS);
  
  double current[3];
  
  for (int c = 0; c < 3; ++c) {
    current[c] = 0;
  }

  for(int i=0; i<n_total_particles; i++) {
    double q = partCfg[i].p.q/time_step;

    for (int c = 0; c < 3; ++c) {
      current[c] += q*partCfg[i].m.v[c];
    }
  }

  sprintf(buffer,"%f %f %f",
	  current[0], current[1], current[2]);

  Tcl_AppendResult(interp, buffer, (char *)NULL); 
#endif  
  return TCL_OK;
}
