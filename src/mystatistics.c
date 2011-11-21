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
	neigh[nb] = atan2(diff[2], diff[1]); // saving the angle
	//if (neigh[nb]> 1.) neigh[nb]=1.; 
        //if (neigh[nb]<-1.) neigh[nb]=-1.; 
	nb++;  // counting neighbours
      }

      // end of inner particle 
    }
    //printf ("part %d neibours %d\n",i, nb);
   
    //counting angle and sai
    if (nb>3) {
      double sinus=0.0;
      double cosinus=0.0;
      for (int l=0; l<nb; l++){
	
	double teta= 6.*neigh[l];
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
  //printf("tot_sai %e parts   %d\n",  tot_sai_m, part_in_bin[bin].n);
  g[rbins]=tot_sai_m;
  g[rbins+1]=tot_sai_m*tot_sai_m;
  //printf("bin %d    sai %e sqr %e\n", bin, tot_sai_m,  g[rbins+1]);
  
  double bin_width = rmax / (double)rbins;
  double inv_bin_width = 1.0 / bin_width;

  // initially set all counts to 0
  for(int i = 0; i < rbins; ++i) {
    g[i] = 0.0;
    grdf[i]=0;
  }

  // number of samples gathered
  // int cnt = 0;

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

static void calc_scaling (double *g, int bin, int boxes, double rclocal)
{
  // char buffer[TCL_INTEGER_SPACE + TCL_DOUBLE_SPACE + 2];
  double rclocal2=rclocal*rclocal, neigh[part_in_bin[bin].n];
  // 12  was set up in initialization as max neighbours in 2d
  int limit=pow(2,boxes);
  double sum_sai_r[limit][limit], sum_sai_m[limit][limit], sai_r, sai_m;
  int division[limit][limit] ;
  double ystep=box_l[1]/(double)limit, zstep=box_l[2]/(double)limit;
  //printf("boxyz %e %e  limit %d  ystep %e  zstep %e check %d\n", box_l[1],box_l[2], limit, ystep,zstep, ((int)99.999)/100);
  
  //setting sums and counters to zero
  for (int i = 0; i < limit; ++i) {
    for (int j = 0; j <limit; ++j) {
      sum_sai_r[i][j]=0.0;
      sum_sai_m[i][j]=0.0;
      division[i][j]=0;
    }
  }
 
  for (int i = 0; i < part_in_bin[bin].n; ++i) {
    int nb=0;
    int p = part_in_bin[bin].e[i];
    for (int j = 0; j < part_in_bin[bin].n; ++j) {
      if (j==i) continue;
      int q = part_in_bin[bin].e[j];
      //printf("%d   %d\n",p,q);
      // minimum image vector between the two particles
      double diff[3];
      get_mi_vector(diff, partCfg[p].r.p, partCfg[q].r.p);
      
      double dist = SQR(diff[1]) + SQR(diff[2]);
      if (dist <rclocal2) {
	neigh[nb] = atan2(diff[2], diff[1]); // saving the angle
	nb++;  // counting neighbours
      }
      
      // end of inner particle 
    }
    
    //printf ("part %d neibours %d\n",i, nb);
    
    //counting angle and sai
    if (nb > 0) {
      double sinus=0.0;
      double cosinus=0.0;
      for (int l=0; l<nb; l++){
	double teta= 6.*neigh[l];
	cosinus = cosinus+cos (teta);
	sinus = sinus+ sin (teta);
      }
      sai_r = cosinus/nb;
      sai_m = sinus/nb;
    }
    else {
      sai_r = sai_m = 0.0;
    }
    //printf("bin %d   nb %d  i %d  sin %e cos %e\n", bin,nb, i, sinus, cosinus);
    
    fold_position(partCfg[p].r.p, partCfg[p].l.i);
    double y=partCfg[p].r.p[1];
    double z=partCfg[p].r.p[2];
    int line=(int)(y/ystep);
    int col=(int)(z/zstep);
    //if (line!=0||col!=0){
    //printf("id %d  line %d col %d y %e  z %e  yc %e  zc %e \n",p,line,col,y,z,ycheck,zcheck);
    // break;
    //}
    sum_sai_r[line][col] = sum_sai_r[line][col] + sai_r;
    sum_sai_m[line][col] = sum_sai_m[line][col] + sai_m;
    division[line][col]++;
    //printf("bin %d  part %d  sai %e  saim %e\n", bin,i, sum_sai_r, sum_sai_m);
    
    // end of outer particle
  }
  //printf("tot_sai %e parts   %d\n",  tot_sai_m, part_in_bin[bin].n);
  int cnt=0;
  for (int j = 0; j < limit; ++j) {
    for (int i = 0; i <limit; ++i) {
      g[cnt]=sum_sai_r[i][j];
      g[cnt+1]=sum_sai_m[i][j];
      g[cnt+2]=division[i][j];
      //printf("cnt %d  bin%d  sqr %e\n", cnt,bin,  g[cnt]);
      cnt=cnt+3;
    }
  }
  
  // g[rbins]=tot_sai_m;
  //g[rbins+1]=tot_sai_m*tot_sai_m;
  //printf("bin %d    sai %e\n", bin, tot_sai_m);
  
}

static void calc_scaling2 (double *g, int bin, int boxes, double rclocal)
{
  // char buffer[TCL_INTEGER_SPACE + TCL_DOUBLE_SPACE + 2];
  double rclocal2=rclocal*rclocal, neigh[part_in_bin[bin].n];
 
  int limit=pow(2,boxes);
  double sum_sai_r[limit][limit], sum_sai_m[limit][limit], sai_r, sai_m;
  // 12  was set up in initialization as max neighbours in 2d
  int division[limit][limit],neighbour[limit][limit][13] ;
  double ystep=box_l[1]/(double)limit, zstep=box_l[2]/(double)limit;
  //printf("boxyz %e %e  limit %d  ystep %e  zstep %e check %d\n", box_l[1],box_l[2], limit, ystep,zstep, ((int)99.999)/100);
  
  //setting sums and counters to zero
  int cnt=0;
  for (int i = 0; i < limit; ++i) {
    for (int j = 0; j <limit; ++j) {
      sum_sai_r[i][j]=0.0;
      sum_sai_m[i][j]=0.0;
      division[i][j]=0; 
      g[cnt]=0.0;
      cnt++;
      for (int l = 0; l < 13; ++l) {
	neighbour[i][j][l]=0;
	g[cnt+l]=0.0;
      }
      cnt=cnt+13;
    }
  }
      
  
  for (int i = 0; i < part_in_bin[bin].n; ++i) {
    int nb=0;
    int p = part_in_bin[bin].e[i];
    for (int j = 0; j < part_in_bin[bin].n; ++j) {
      if (j==i) continue;
      int q = part_in_bin[bin].e[j];
      //printf("%d   %d\n",p,q);
      // minimum image vector between the two particles
      double diff[3];
      get_mi_vector(diff, partCfg[p].r.p, partCfg[q].r.p);
      
      double dist = SQR(diff[1]) + SQR(diff[2]);
      if (dist <rclocal2) {
	neigh[nb] = atan2(diff[2], diff[1]); // saving the angle
	nb++;  // counting neighbours
      }
      
      // end of inner particle 
    }
    
    //printf ("part %d neibours %d\n",i, nb);
    
    //counting angle and sai
    if (nb > 0) {
      double sinus=0.0;
      double cosinus=0.0;
      for (int l=0; l<nb; l++){
	double teta= 6.*neigh[l];
	cosinus = cosinus+cos (teta);
	sinus = sinus+ sin (teta);
      }
      sai_r = cosinus/nb;
      sai_m = sinus/nb;
    }
    else {
      sai_r = sai_m = 0.0;
    }
    //printf("bin %d   nb %d  i %d  sin %e cos %e\n", bin,nb, i, sinus, cosinus);
    
    fold_position(partCfg[p].r.p, partCfg[p].l.i);
    double y=partCfg[p].r.p[1];
    double z=partCfg[p].r.p[2];
    int line=(int)(y/ystep);
    int col=(int)(z/zstep);
    //if (line!=0||col!=0){
    //printf("id %d  line %d col %d y %e  z %e  yc %e  zc %e \n",p,line,col,y,z,ycheck,zcheck);
    // break;
    //}
    sum_sai_r[line][col] = sum_sai_r[line][col] + sai_r;
    sum_sai_m[line][col] = sum_sai_m[line][col] + sai_m;
    division[line][col]++;
    if (nb<13) {neighbour[line][col][nb]++;} 
    //printf("bin %d  part %d  sai %e  saim %e  nbr%d\n", bin,i, sum_sai_r[line][col], sum_sai_m[line][col], neighbour[line][col][nb]);
    
    // end of outer particle
  }
  //printf("tot_sai %e parts   %d\n",  tot_sai_m, part_in_bin[bin].n);
  cnt=0;
  for (int j = 0; j < limit; ++j) {
    for (int i = 0; i <limit; ++i) { 
      if (division[i][j]!=0) {
	g[cnt]=(sum_sai_r[i][j]*sum_sai_r[i][j]+sum_sai_m[i][j]*sum_sai_m[i][j])/(double)(division[i][j]*division[i][j]);
	//printf("cnt %d  bin %d  sqr %e %e %d\n", cnt,bin, g[cnt], sum_sai_r[i][j]*sum_sai_r[i][j],division[i][j] );
	cnt++;
	for (int nbr = 0; nbr <13; ++nbr) {
	  g[cnt+nbr]=(double)neighbour[i][j][nbr]/division[i][j];
	  // printf("counter %d  bin %d  sqr %e\n", nbr,bin, g[cnt+nbr]);
	}
	//printf("cnt %d  bin %d  sqr %e\n", cnt,bin, g[cnt-1]);
	cnt=cnt+13;
      } else { cnt=cnt+14;}
    }
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
  int rbins, boxes;
  enum { BINS, MX, MYZ, RDFYZ,BONDYZ,SCALE,SCALE2, PRINT };

  if (argc < 2) {
    Tcl_AppendResult(interp, "expected: analyze wallstuff -bins <binboundaries> | -myz <bin> |-mx <bin> | -rdfyz <bin> <rmin> <rmax> <rdfbins> |-bondyz| -scale| -scale2 ",
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
  else if (ARG0_IS_S("-scale") && argc == 4) {
    job = SCALE;
  }
  else if (ARG0_IS_S("-scale2") && argc == 4) {
    job = SCALE2;
  }
  else if (ARG0_IS_S("-print") && argc == 2) {
    job = PRINT;
  }
  else {
    Tcl_AppendResult(interp, ": analyze wallstuff -bins|-myz|-mx|-rdfyz|-bondyz|-scale|-scale2...", (char *)NULL);
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
  case SCALE:
  case SCALE2:
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
  case SCALE:
    if (!ARG_IS_I(2, boxes)) {
      return (TCL_ERROR);
    }
  case SCALE2:
    if (!ARG_IS_I(2, boxes)) {
      return (TCL_ERROR);
    }
    if (!ARG_IS_D(3, rclocal)) {
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
  case SCALE:
    // these cases use partCfg
    updatePartCfg(WITHOUT_BONDS);
    break;
  case SCALE2:
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
    realloc_doublelist(&g, g.n = rbins+2);
    calc_wallbondyz(g.e, bin, rclocal, rmax, rbins);
    break;
  case SCALE:
    realloc_doublelist(&g, g.n = 3*pow(4,boxes));
    calc_scaling (g.e,bin, boxes, rclocal);
    break;
  case SCALE2:
    realloc_doublelist(&g, g.n = 14*pow(4,boxes));
    calc_scaling2 (g.e,bin, boxes, rclocal);
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


// order parameter, no types!!!!!!!!!!!!!!!!!!!!!!!!!
// dcut - min q6 for assosiating for the cluster 
static void calc_q6(double *g, double rclocal, double dcut)
{
  // char buffer[TCL_INTEGER_SPACE + TCL_DOUBLE_SPACE + 2];
  double rclocal2=rclocal*rclocal, phi, theta;
  double   q6square[n_total_particles];
  // spherical harmonics im and re parts separetly
  // only positive values, since for positive m it's just minus in front (odd) and for psi for negative m
  // careful when summing up!
  double y0[n_total_particles], y1s[n_total_particles], y1c[n_total_particles], y2s[n_total_particles], y2c[n_total_particles];
  double y3s[n_total_particles], y3c[n_total_particles], y4s[n_total_particles], y4c[n_total_particles];
  double y5s[n_total_particles], y5c[n_total_particles], y6s[n_total_particles], y6c[n_total_particles];
  int neighb[n_total_particles];
  double pi=4*atan(1.);
  int maxneighb=35; // maximum of neighbours i can think of
  int cluster[n_total_particles*maxneighb], connect[n_total_particles];
 

  printf("starting\n");
  for (int i = 0; i <  n_total_particles; ++i) {
    //check for (int i = 0; i <  1; ++i) {
    int nb=0; //number of neighb for current

    //intializing for summing up
    y0[i]=0.0; y1s[i]=0.0; y1c[i]=0.0; y2s[i]=0.0; y2c[i]=0.0;
    y3s[i]=0.0; y3c[i]=0.0; y4s[i]=0.0; y4c[i]=0.0; y5s[i]=0.0; y5c[i]=0.0; y6s[i]=0.0; y6c[i]=0.0;

    for (int j = 0; j < n_total_particles; ++j) {
      if (j==i) continue;
      //printf("%d   %d\n",p,q);
      //check for (int j = 0; j < 1; ++j) {
      // minimum image vector between the two particles
      double diff[3];
      get_mi_vector(diff, partCfg[i].r.p, partCfg[j].r.p);
      double dist2 = SQR(diff[1]) + SQR(diff[2])+SQR(diff[0]);

      //neighbours withing cutting off raduis
      if (dist2 <rclocal2) {
	phi = atan2(diff[1], diff[0]); // polar angle
	theta = acos(diff[2]/sqrt(dist2)); // azimuth angle
	//check phi=pi*3./4.;
	//check theta=pi-atan(sqrt(2.));
	printf("i %d j %d   phi %e theta %e \n",i, j, phi/pi, theta/pi);
	//summing up separetly spherical harmonics
	//no coeff in front of harmonics, would be multiplied with division over neighb
	double cos1=cos(theta);
	double cos2=cos1*cos1;
	double cos3=cos2*cos1;
	double sin1=sin(theta);
	double sin2=sin1*sin1;
	double sin3=sin2*sin1;
	//check y0[i]+=231.*cos3*cos3-315.*cos2*cos2+105.*cos2-5.;
	//check y1c[i]+=sin1*(33.*cos2*cos3-30.*cos3+5.*cos1);
	//check y2c[i]+=sin2*(33.*cos2*cos2-18.*cos2+1.);	
	//check y3c[i]+=sin3*(11.*cos3-3.*cos1);
	//check y4c[i]+=sin2*sin2*(11.*cos2-1.);	
	//check y5c[i]+=sin2*sin3*cos1;
	//check y6c[i]+=sin3*sin3;

	y0[i]+=231.*cos3*cos3-315.*cos2*cos2+105.*cos2-5.;
	y1c[i]+=cos(phi)*sin1*(33.*cos2*cos3-30.*cos3+5.*cos1);
	y1s[i]+=sin(phi)*sin1*(33.*cos2*cos3-30.*cos3+5.*cos1);
	y2c[i]+=cos(2.*phi)*sin2*(33.*cos2*cos2-18.*cos2+1.);
	y2s[i]+=sin(2.*phi)*sin2*(33.*cos2*cos2-18.*cos2+1.);
	y3c[i]+=cos(3.*phi)*sin3*(11.*cos3-3.*cos1);
	y3s[i]+=sin(3.*phi)*sin3*(11.*cos3-3.*cos1);
	y4c[i]+=cos(4.*phi)*sin2*sin2*(11.*cos2-1.);
	y4s[i]+=sin(4.*phi)*sin2*sin2*(11.*cos2-1.);
	y5c[i]+=cos(5.*phi)*sin2*sin3*cos1;
	y5s[i]+=sin(5.*phi)*sin2*sin3*cos1;
	y6c[i]+=cos(6.*phi)*sin3*sin3;
	y6s[i]+=sin(6.*phi)*sin3*sin3;
	//if (neigh[nb]> 1.) neigh[nb]=1.; 
        //if (neigh[nb]<-1.) neigh[nb]=-1.; 
	nb++;  // counting neighbours
      }	
      // end of the inner particle 
    }

    //check	y0[i] = y0[i]/16.;
    //check y1c[i] = y1c[i]/16.;
    
    //check  y2c[i] = y2c[i]/64./2.;
    
    //check y3c[i] = y3c[i]/32./12.;
    //check  y4c[i] = y4c[i]/32./120.;
    //check y5c[i] = y5c[i]/32./120.;
    //check y6c[i] = y6c[i]/64./30./24.;
    //check  printf("y0 %e\n",y0[i]);
    //check  printf("y1 %e\n",y1c[i]);
    //check  printf("y2 %e\n",y2c[i]);
    //check printf("y3 %e\n",y3c[i]);
    //check  printf("y4 %e\n",y4c[i]);
    //check  printf("y5 %e\n",y5c[i]);
    //check  printf("y6 %e\n",y6c[i]);
    
    y0[i] = y0[i]*sqrt(13./pi)/nb/32.;
    y1c[i] = y1c[i]*sqrt(273.*0.5/pi)/nb/16.;
    y1s[i] = y1s[i]*sqrt(273.*0.5/pi)/nb/16.;
    y2c[i] = y2c[i]*sqrt(1365./pi)/nb/64.;
    y2s[i] = y2s[i]*sqrt(1365./pi)/nb/64.;
    y3c[i] = y3c[i]*sqrt(1365./pi)/nb/32.;
    y3s[i] = y3s[i]*sqrt(1365./pi)/nb/32.;
    y4c[i] = y4c[i]*sqrt(91.*0.5/pi)/nb/32.*3.;
    y4s[i] = y4s[i]*sqrt(91.*0.5/pi)/nb/32.*3.;
    y5c[i] = y5c[i]*sqrt(1001./pi)/nb/32.*3.;
    y5s[i] = y5s[i]*sqrt(1001./pi)/nb/32.*3.;
    y6c[i] = y6c[i]*sqrt(3003./pi)/nb/64.;
    y6s[i] = y6s[i]*sqrt(3003./pi)/nb/64.;
    //printf("y0 %e\n",y0[i]);
    // printf("y1c %e y1s %e\n",y1c[i],y1s[i] );
    //printf("y2c %e y2s %e\n",y2c[i],y2s[i] );
    //printf("y3c %e y3s %e\n",y3c[i],y3s[i] );
    //printf("y4c %e y4s %e\n",y4c[i],y4s[i] );
    //printf("y5c %e y5s %e\n",y5c[i],y5s[i] );
    //printf("y6c %e y6s %e\n",y6c[i],y6s[i] );

    q6square[i] = y0[i]*y0[i]+2.*SQR(y1c[i])+2.*SQR(y1s[i])+2.*SQR(y2c[i])+2.*SQR(y2s[i])+2.*SQR(y3c[i])+2.*SQR(y3s[i]);
    q6square[i] += 2.*(SQR(y4c[i])+SQR(y4s[i])+SQR(y5c[i])+SQR(y5s[i])+SQR(y6c[i])+SQR(y6s[i]));
    //!!q6square[i] *=4.*pi/13.;
    neighb[i] = nb;
    printf("neibours %d parameter %e \n",nb, q6square[i]);
    //q6square[i]=0.;
    //!!!!!!!!!!!!!!!! check only!
    g[i*3] = sqrt(q6square[i]*4.*pi/13.);
    g[i*3+1] = (double)neighb[i];
    //!!! nothing done about neighbours    
    //printf ("part %d neibours %d\n",i, nb)
    //printf("bin %d   nb %d  i %d  sin %e cos %e\n", bin,nb, i, sinus, cosinus);
    //printf("bin %d  part %d  sai %e  saim %e\n", bin,i, sum_sai_r, sum_sai_m);
    // end of outer particle
  }

  return;
  
  
  // loop over all distinct particle pairs 
  //printf("here1");
  double y0ij, y1ij,y2ij,y3ij,y4ij,y5ij,y6ij;
  
  //initialization, cant be 0, since would be used for clustering
  for(int i = 0; i < (n_total_particles*maxneighb); ++i) {
    cluster[i]=n_total_particles;
  }
  for (int i = 0; i < n_total_particles; ++i) {
    int nb=0;
    for (int j = 0; j < n_total_particles; ++j) {
      if (j==i) continue;
      // minimum image vector between the two particles
      double diff[3];
      get_mi_vector(diff, partCfg[i].r.p, partCfg[j].r.p);
      double dist2 = SQR(diff[1]) + SQR(diff[2])+SQR(diff[0]);
      if (dist2 < rclocal2) {
	// here are _all_ harmonics - and + m
	y0ij=y0[i]*y0[j];
	y1ij=y1c[i]*y1c[j]+y1s[i]*y1s[j];
	y2ij=y2c[i]*y2c[j]+y2s[i]*y2s[j];
	y3ij=y3c[i]*y3c[j]+y3s[i]*y3s[j];
	y4ij=y4c[i]*y4c[j]+y4s[i]*y4s[j];
	y5ij=y5c[i]*y5c[j]+y5s[i]*y5s[j];
	y6ij=y6c[i]*y6c[j]+y6s[i]*y6s[j];

	double q6ij= y0ij+2.*y1ij+2.*y2ij+2.*y3ij+2.*y4ij+2.*y5ij+2.*y6ij;
	//printf("i %d j %d   q6ij  %e  q6i2  %e  q6j2  %e\n",i, j, q6ij,q6square[i],q6square[j] );
	q6ij/=sqrt(q6square[i]*q6square[j]);
	//printf("i %d j %d   q6ij normalized  %e dcut %e \n",i, j, q6ij,dcut);
	if (q6ij<dcut) printf("Attention!! i %d j %d   q6ij normalized  %2.15e dcut %e \n",i, j, q6ij,dcut);
	if (q6ij>=dcut){
	  cluster[i*maxneighb+nb]=j;
	  //printf("i %d j %d   q6ij  %e dcut %e\n",i, j, q6ij,dcut);
	  nb++;
	}
	// printf("%e\n",g[ind]);
      }
    }
    connect[i]=nb;// total number of connections [solid] for particle i
    g[3*i+2]=(double)nb;
  }
  
  // printf(" %e  rdf %d\n",g[i], grdf[i]);
  // Tcl_AppendResult(interp, buffer, (char *)NULL); 
}




int parse_order(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze wallmsd -xy|-z <min> <max>' */
  /******************************************************************************/
  char buffer[TCL_INTEGER_SPACE + TCL_DOUBLE_SPACE + 2];
  DoubleList g;
  int job;
  double rclocal,dcut;
  
  enum {Q6, Q4};

  if (argc < 3) {
    Tcl_AppendResult(interp, "expected: analyze order -q4 <rclocal> <dcut>  |-q6 <rclocal> <dcut> ",
		     (char *)NULL);
    return TCL_ERROR;
  }

  // 1. what do we do?
  if (ARG0_IS_S("-q6") && argc == 3) {
    job = Q6;
  }
  else if (ARG0_IS_S("-q4") && argc == 3) {
    job = Q4;
  }
  else {
    Tcl_AppendResult(interp, ": analyze order -q4|-q6 ...", (char *)NULL);
    return TCL_ERROR;
  }
  
  // 2. parameters

  switch (job) {
  case Q4:
    if (!ARG_IS_D(1, rclocal)) {
      return (TCL_ERROR);
    }
    if (!ARG_IS_D(2, dcut)) {
      return (TCL_ERROR);
    }
    break;
  case Q6:
    if (!ARG_IS_D(1, rclocal)) {
      return (TCL_ERROR);
    }
    if (!ARG_IS_D(2, dcut)) {
      return (TCL_ERROR);
    }
    break;
  
  }
  
  // result double list
  init_doublelist(&g);
  
  // check that data is there
  switch (job) {
  case Q4:
    // these cases use partCfg
    updatePartCfg(WITHOUT_BONDS);
    break;
  case Q6:
    // these cases use partCfg
    updatePartCfg(WITHOUT_BONDS);
    break;
  }

  // finally, do what is necessary
  switch (job) {
  case Q4:
    //realloc_doublelist(&g, g.n = rbins);
    //calc_q4(g.e, rclocal, dcut);
    //break;
    Tcl_AppendResult(interp, "analyze order -q4 <rclocal>  is not written, sorry for inconvenience",
		     (char *)NULL);
    return TCL_ERROR;
  case Q6:
    realloc_doublelist(&g, g.n = (3*n_total_particles));
    calc_q6(g.e, rclocal,dcut);
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
    double pos[3], fpos[3]; int img[3] = {0, 0, 0};
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
