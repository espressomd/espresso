/*
  Copyright (C) 2013 The ESPResSo project
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

#include "statistics_wallstuff.hpp"

#include "utils.hpp"
#include "grid.hpp"

// list of the currently specified box boundaries
DoubleList wallstuff_boundaries = { NULL, 0 };
// the boxes with the particle identities
IntList *wallstuff_part_in_bin = NULL;

void wall_sort_particles()
{
  // 1. reallocate the boxes for the particle identities
  // save current number of bins
  int n_part_in_bin = wallstuff_boundaries.n-1;
  // free old boxes
  for (int i = 0; i < n_part_in_bin; ++i) {
    realloc_intlist(&wallstuff_part_in_bin[i], wallstuff_part_in_bin[i].n = 0);
  }
  wallstuff_part_in_bin = (IntList*)realloc(wallstuff_part_in_bin, 
                                            (wallstuff_boundaries.n-1)*sizeof(IntList));
  // initialize new ones
  for (int i = n_part_in_bin; i < wallstuff_boundaries.n-1; ++i) {
    init_intlist(&wallstuff_part_in_bin[i]);
  }

  // 2. for each particle, find the box and put its
  // identity there
  for(int i=0; i<n_part; i++) {
    double x = partCfg[i].r.p[0];
    // ignore particles outside the wallstuff_boundaries
    if (x < wallstuff_boundaries.e[0] || x > wallstuff_boundaries.e[wallstuff_boundaries.n-1])
      continue;
    // simple bisection on the particle's x-coordinate
    int s = 0;
    int e = wallstuff_boundaries.n - 1;
    while (e - s > 1) {
      int c = (e + s)/2;
      if (x >= wallstuff_boundaries.e[c]) s = c; else e = c;
    }
    // and add the particle to the resulting list
    realloc_grained_intlist(&wallstuff_part_in_bin[s], wallstuff_part_in_bin[s].n + 1, 8);
    wallstuff_part_in_bin[s].e[wallstuff_part_in_bin[s].n++] = i;
  }
}

void calc_wallmsdyz(double *g, int bin)
{
  // loop over all stored configurations
  // and calculate the MSD with respect to the current configuration
  // MSD of configuration with itself is always 0
  g[0] = 0;
  for(int k = 1; k < n_configs; k++) {
    g[k] = 0.0;
    // loop over all particles in the specified bin and add up MSD
    for (int i = 0; i < wallstuff_part_in_bin[bin].n; ++i) {
      int p = wallstuff_part_in_bin[bin].e[i];
      g[k] +=
	+ SQR(configs[n_configs-1][3*p + 1]-configs[n_configs-1-k][3*p + 1])
	+ SQR(configs[n_configs-1][3*p + 2]-configs[n_configs-1-k][3*p + 2]);
    }
    // normalize
    g[k] /= wallstuff_part_in_bin[bin].n;
  }
}

void calc_wallmsdx(double *g, int bin)
{
  // see calc_wallmsdyz, just for x
  g[0] = 0;
  for(int k = 1; k < n_configs; k++) {
    g[k] = 0.0;
    for (int i = 0; i < wallstuff_part_in_bin[bin].n; ++i) {
      int p = wallstuff_part_in_bin[bin].e[i];
      g[k] += SQR(configs[n_configs-1][3*p]-configs[n_configs-1-k][3*p]);
    }
    g[k] /= wallstuff_part_in_bin[bin].n;
  }
}

void calc_wallbondyz(double *g, int bin, double rclocal, double rmax, int rbins)
{
  // char buffer[TCL_INTEGER_SPACE + TCL_DOUBLE_SPACE + 2];
  double rclocal2=rclocal*rclocal, neigh[wallstuff_part_in_bin[bin].n];
  double sum_sai_r=0.0, sum_sai_m=0.0, sai_r[wallstuff_part_in_bin[bin].n], sai_m[wallstuff_part_in_bin[bin].n];
  int grdf[rbins];

  for (int i = 0; i < wallstuff_part_in_bin[bin].n; ++i) {
    int nb=0;
    for (int j = 0; j < wallstuff_part_in_bin[bin].n; ++j) {
      if (j==i) continue;
      int p = wallstuff_part_in_bin[bin].e[i];
      int q = wallstuff_part_in_bin[bin].e[j];
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

  double tot_sai_m= sqrt((sum_sai_r * sum_sai_r) + (sum_sai_m * sum_sai_m))/wallstuff_part_in_bin[bin].n;
  //printf("tot_sai %e parts   %d\n",  tot_sai_m, wallstuff_part_in_bin[bin].n);
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
  for (int i = 0; i < wallstuff_part_in_bin[bin].n; ++i) {
    for (int j = i + 1; j < wallstuff_part_in_bin[bin].n; ++j) {
      int p = wallstuff_part_in_bin[bin].e[i];
      int q = wallstuff_part_in_bin[bin].e[j];

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
  // double width = wallstuff_boundaries.e[bin+1] - wallstuff_boundaries.e[bin];
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

void calc_scaling(double *g, int bin, int boxes, double rclocal)
{
  // char buffer[TCL_INTEGER_SPACE + TCL_DOUBLE_SPACE + 2];
  double rclocal2=rclocal*rclocal, neigh[wallstuff_part_in_bin[bin].n];
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
 
  for (int i = 0; i < wallstuff_part_in_bin[bin].n; ++i) {
    int nb=0;
    int p = wallstuff_part_in_bin[bin].e[i];
    for (int j = 0; j < wallstuff_part_in_bin[bin].n; ++j) {
      if (j==i) continue;
      int q = wallstuff_part_in_bin[bin].e[j];
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
  //printf("tot_sai %e parts   %d\n",  tot_sai_m, wallstuff_part_in_bin[bin].n);
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

void calc_scaling2 (double *g, int bin, int boxes, double rclocal)
{
  // char buffer[TCL_INTEGER_SPACE + TCL_DOUBLE_SPACE + 2];
  double rclocal2=rclocal*rclocal, neigh[wallstuff_part_in_bin[bin].n];
 
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
      
  
  for (int i = 0; i < wallstuff_part_in_bin[bin].n; ++i) {
    int nb=0;
    int p = wallstuff_part_in_bin[bin].e[i];
    for (int j = 0; j < wallstuff_part_in_bin[bin].n; ++j) {
      if (j==i) continue;
      int q = wallstuff_part_in_bin[bin].e[j];
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
  //printf("tot_sai %e parts   %d\n",  tot_sai_m, wallstuff_part_in_bin[bin].n);
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

void calc_wallrdfyz(double *g, int bin, double rmin, double rmax, int rbins)
{
  double bin_width = (rmax-rmin) / (double)rbins;
  double inv_bin_width = 1.0 / bin_width;

  // initially set all counts to 0
  for(int i = 0; i < rbins; ++i) g[i] = 0.0;

  // number of samples gathered
  int cnt = 0;

  // loop over all distinct particle pairs in the bin
  for (int i = 0; i < wallstuff_part_in_bin[bin].n; ++i) {
    for (int j = i + 1; j < wallstuff_part_in_bin[bin].n; ++j) {
      int p = wallstuff_part_in_bin[bin].e[i];
      int q = wallstuff_part_in_bin[bin].e[j];

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
  double width = wallstuff_boundaries.e[bin+1] - wallstuff_boundaries.e[bin];
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
