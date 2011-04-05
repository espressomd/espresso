#include "statistics_observable.h"
#include "particle_data.h"
#include "parser.h"
#include "integrate.h"

observable** observables = 0;
int n_observables = 0; 


static int convert_types_to_ids(IntList * type_list, IntList * id_list); 
  
int tclcommand_parse_profile(Tcl_Interp* interp, int argc, char** argv, int* change, int* dim_A, profile_data** pdata_);
int tclcommand_parse_radial_profile(Tcl_Interp* interp, int argc, char** argv, int* change, int* dim_A, radial_profile_data** pdata);
int sf_print_usage(Tcl_Interp* interp);



int tclcommand_observable_print(Tcl_Interp* interp, int argc, char** argv, int* change, observable* obs) {
  char buffer[TCL_DOUBLE_SPACE];
  double* values=malloc(obs->n*sizeof(double));
  (*obs->fun)(obs->args, values, obs->n);
  for (int i = 0; i<obs->n; i++) {
    Tcl_PrintDouble(interp, values[i], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *)NULL );
  }
  return TCL_OK;
}


int tclcommand_observable_particle_velocities(Tcl_Interp* interp, int argc, char** argv, int* change, observable* obs) {
  IntList* ids;
  int temp;
  if (parse_id_list(interp, argc-1, argv+1, &temp, &ids) != TCL_OK ) 
    return TCL_ERROR;
  obs->fun=&observable_particle_velocities;
  obs->args=ids;
  obs->n=3*ids->n;
  *change=1+temp;
  return TCL_OK;
}

int tclcommand_observable_com_velocity(Tcl_Interp* interp, int argc, char** argv, int* change, observable* obs) {
  IntList* ids;
  int temp;
  if (parse_id_list(interp, argc-1, argv+1, &temp, &ids) != TCL_OK ) 
    return TCL_ERROR;
  obs->fun=&observable_com_velocity;
  obs->args=ids;
  obs->n=3;
  *change=1+temp;
  return TCL_OK;
}

int tclcommand_observable_particle_positions(Tcl_Interp* interp, int argc, char** argv, int* change, observable* obs) {
  IntList* ids;
  int temp;
  printf("parsing particle positions\n");
  if (parse_id_list(interp, argc-1, argv+1, &temp, &ids) != TCL_OK ) 
     return TCL_ERROR;
  obs->fun = &observable_particle_positions;
  obs->args=(void*)ids;
  obs->n=3*ids->n;
  *change=1+temp;
  return TCL_OK;
}


int tclcommand_observable_density_profile(Tcl_Interp* interp, int argc, char** argv, int* change, observable* obs){
  int temp;
  profile_data* pdata;
  obs->fun = &observable_density_profile;
  if (! tclcommand_parse_profile(interp, argc-1, argv+1, &temp, &obs->n, &pdata) == TCL_OK ) 
    return TCL_ERROR;
  obs->args=(void*)pdata;
  obs->n=pdata->nbins;
  *change=1+temp;
  return TCL_OK;
}

int tclcommand_observable_radial_density_profile(Tcl_Interp* interp, int argc, char** argv, int* change, observable* obs){
  int temp;
  radial_profile_data* rpdata;
  obs->fun = &observable_radial_density_profile;
  if (! tclcommand_parse_radial_profile(interp, argc-1, argv+1, &temp, &obs->n, &rpdata) == TCL_OK ) 
     return TCL_ERROR;
  obs->args=(void*)rpdata;
  obs->n=rpdata->nbins;
  *change=1+temp;
  return TCL_OK;
}

int tclcommand_observable_particle_currents(Tcl_Interp* interp, int argc, char** argv, int* change, observable* obs){
#ifdef ELECTROSTATICS
  int temp;
  IntList* ids;
  obs->fun = &observable_particle_currents;
  if (! parse_id_list(interp, argc-1, argv+1, &temp, &ids) == TCL_OK ) 
    return TCL_ERROR;
  obs->args=(void*)ids;
  obs->n=3*ids->n;
  *change=1+temp;
  return TCL_OK;
#else
  Tcl_AppendResult(interp, "Feature ELECTROSTATICS needed for observable particle_currents", (char *)NULL);
  return TCL_ERROR;
#endif
  
}

int tclcommand_observable_currents(Tcl_Interp* interp, int argc, char** argv, int* change, observable* obs) {
#ifdef ELECTROSTATICS
  int temp;
  IntList* ids;
  obs->fun = &observable_currents;
  if (! parse_id_list(interp, argc-1, argv+1, &temp, &ids) == TCL_OK ) 
    return TCL_ERROR;
  obs->args=(void*)ids;
  obs->n=3;
  *change=1+temp;
  return TCL_OK;
#else
  Tcl_AppendResult(interp, "Feature ELECTROSTATICS needed for observable currents", (char *)NULL);
  return TCL_ERROR;
#endif
} 

int tclcommand_observable_structure_factor(Tcl_Interp* interp, int argc, char** argv, int* change, observable* obs) {
  int order, order2;
  int* order_p;

  Tcl_AppendResult(interp, "Structure Factor not available yet!!", (char *)NULL);
  return TCL_ERROR;
  if (argc > 1 && ARG1_IS_I(order)) {
    obs->fun = &observable_structure_factor;
    order_p=malloc(sizeof(int));
    *order_p=order;
    obs->args=(void*) order_p;
    int order2,i,j,k,l,n ; 
    order2=order*order;
    l=0;
    // lets counter the number of entries for the DSF
    for(i=-order; i<=order; i++) 
      for(j=-order; j<=order; j++) 
        for(k=-order; k<=order; k++) {
          n = i*i + j*j + k*k;
          if ((n<=order2) && (n>=1)) 
            l=l+2;
  }
    obs->n=l;
    *change=2;
    return TCL_OK;
  } else { 
    sf_print_usage(interp);
    return TCL_ERROR; 
  }
}

int tclcommand_observable_interacts_with(Tcl_Interp* interp, int argc, char** argv, int* change, observable* obs) {
  IntList *ids, *ids1;
  int temp;
  double cutoff;
  obs->fun = &observable_interacts_with;
  ids=(IntList*)malloc(sizeof(IntList));
  ids1=(IntList*)malloc(sizeof(IntList));
  iw_params* iw_params_p;
  if (! parse_id_list(interp, argc-1, argv+1, &temp, &ids) == TCL_OK ) {
    free(ids);
    free(ids1);
    free(iw_params_p);
    return TCL_ERROR;
  }
  iw_params_p=(iw_params*)malloc(sizeof(iw_params));
  iw_params_p->ids1=ids;
  *change=1+temp;
  if (! parse_id_list(interp, argc-3, argv+3, &temp, &ids1) == TCL_OK ) {
    free(ids);
    free(ids1);
    free(iw_params_p);
    return TCL_ERROR;
  }
  *change+=temp;
  iw_params_p->ids2=ids1;
  if ( argc < 7 || !ARG_IS_D(6,cutoff)) {
    Tcl_AppendResult(interp, "Usage: analyze correlation ... interacts_with id_list1 id_list2 cutoff", (char *)NULL);
    free(ids);
    free(ids1);
    free(iw_params_p);
  return TCL_ERROR;
  }
  *change+=1;
  iw_params_p->cutoff=cutoff;
  obs->args=(void*)iw_params_p;
  obs->n=ids->n; // number of ids from the 1st argument
  return TCL_OK;
}

//  if (ARG0_IS_S("textfile")) {
//    // We still can only handle full files
//    if ( argc>1 ) {
//      fds= malloc(sizeof(file_data_source));
//      error = file_data_source_init(fds, argv[1], 0);
//      *change=2;
//      if (!error) {
//        *A_args=(void*) fds;
//        *dim_A = fds->n_columns;
//        *A_fun = (void*)&file_data_source_readline;
//        return TCL_OK;
//      } else {
//        Tcl_AppendResult(interp, "Error reading file ", argv[1] ,"\n", (char *)NULL);
//        Tcl_AppendResult(interp, file_data_source_init_errors[error] ,"\n", (char *)NULL);
//        return TCL_ERROR;
//      }
//    } else {
//      Tcl_AppendResult(interp, "Error in parse_observable textfile: no filename given" , (char *)NULL);
//      return TCL_ERROR;
//    }
//    return TCL_OK;
//  } 
//  if (ARG0_IS_S("tclinput")) {
//    if (argc>1 && ARG1_IS_I(temp)) {
//      *dim_A = temp;
//      *A_fun = &tcl_input;
//      *A_args=malloc(sizeof(tcl_input_data));
//      *change =2;
//      return TCL_OK;
//    } else {
//      Tcl_AppendResult(interp, "\nError in parse_observable tclinfo. You must pass the dimension of the observable." , (char *)NULL);
//      return TCL_ERROR;
//    }
//  }else {
//    return observable_usage(interp);
//  }
//  return 0 ;















#define REGISTER_OBSERVABLE(name,parser) \
  if (ARG1_IS_S(#name)) { \
    observables[n_observables]=malloc(sizeof(observable)); \
    if (parser(interp, argc-1, argv+1, &temp, observables[n_observables]) ==TCL_OK) { \
      n_observables++; \
      argc-=1+temp; \
      argv+=1+temp; \
      return TCL_OK; \
    } else { \
      free(observables[n_observables]);\
      Tcl_AppendResult(interp, "\nError parsing observable ", #name, "\n", (char *)NULL); \
      return TCL_ERROR; \
    } \
  } 



int tclcommand_observable(ClientData data, Tcl_Interp *interp, int argc, char **argv){
//  file_data_source* fds;
  char buffer[TCL_INTEGER_SPACE];
  int n;
  int error=0;
  int temp;
  int order=0;
  double cutoff;
  int *order_p;
  IntList* ids;
  IntList* ids1;
  iw_params* iw_params_p;
  
  profile_data* pdata=0;
  radial_profile_data* rpdata=0;

  observables=(observable**) realloc(observables, (n_observables+1)*sizeof(observable*)); 
  if (argc<2) {
    Tcl_AppendResult(interp, "Usage!!!\n", (char *)NULL);
    return TCL_ERROR;
  }

  if (argc > 1 && ARG_IS_S(1, "n_observables")) {
	  sprintf(buffer, "%d", n_observables);
    Tcl_AppendResult(interp, buffer, (char *)NULL );
    return TCL_OK;
  }

  REGISTER_OBSERVABLE(particle_velocities, tclcommand_observable_particle_velocities);
  REGISTER_OBSERVABLE(com_velocity, tclcommand_observable_com_velocity);
  REGISTER_OBSERVABLE(particle_positions, tclcommand_observable_particle_positions);
  REGISTER_OBSERVABLE(particle_currents, tclcommand_observable_particle_currents);
  REGISTER_OBSERVABLE(currents, tclcommand_observable_currents);
  REGISTER_OBSERVABLE(structure_factor, tclcommand_observable_structure_factor);
  REGISTER_OBSERVABLE(interacts_with, tclcommand_observable_interacts_with);
//  REGISTER_OBSERVABLE(obs_nothing, tclcommand_observable_obs_nothing);
//  REGISTER_OBSERVABLE(flux_profile, tclcommand_observable_flux_profile);
  REGISTER_OBSERVABLE(density_profile, tclcommand_observable_density_profile);
//  REGISTER_OBSERVABLE(lb_velocity_profile, tclcommand_observable_lb_velocity_profile);
  REGISTER_OBSERVABLE(radial_density_profile, tclcommand_observable_radial_density_profile);
  
  if (ARG1_IS_I(n)) {
    if (n>=n_observables) {
      Tcl_AppendResult(interp, "Observable %d does not yet exist.", n ,"\n", (char *)NULL);
      return TCL_ERROR;
    }
    if (argc > 2 && ARG_IS_S(2,"print")) {
      return tclcommand_observable_print(interp, argc-2, argv+2, &temp, observables[n]);
    }
  }
  Tcl_AppendResult(interp, "Unknown observable ", argv[1] ,"\n", (char *)NULL);
  return TCL_ERROR;
}


int parse_id_list(Tcl_Interp* interp, int argc, char** argv, int* change, IntList** ids ) {
  int i,ret;
//  char** temp_argv; int temp_argc;
//  int temp;
  IntList* input=malloc(sizeof(IntList));
  IntList* output=malloc(sizeof(IntList));
  init_intlist(input);
  alloc_intlist(input,1);
  init_intlist(output);
  alloc_intlist(output,1);
  if (argc<1) {
    Tcl_AppendResult(interp, "Error parsing id list\n", (char *)NULL);
    return TCL_ERROR;
  }


  if (ARG0_IS_S("id")) {
    if (!parse_int_list(interp, argv[1],input)) {
      Tcl_AppendResult(interp, "Error parsing id list\n", (char *)NULL);
      return TCL_ERROR;
    } 
    *ids=input;
    for (i=0; i<input->n; i++) {
      if (input->e[i] >= n_total_particles) {
        Tcl_AppendResult(interp, "Error parsing ID list. Given particle ID exceeds the number of existing particles\n", (char *)NULL);
        return TCL_ERROR;
      }
    }
    *change=2;
    return TCL_OK;

  } else if ( ARG0_IS_S("types") ) {
    if (!parse_int_list(interp, argv[1],input)) {
      Tcl_AppendResult(interp, "Error parsing types list\n", (char *)NULL);
      return TCL_ERROR;
    } 
    if( (ret=convert_types_to_ids(input,output))<=0){
        Tcl_AppendResult(interp, "Error parsing types list. No particle with given types found.\n", (char *)NULL);
        return TCL_ERROR;
    } else { 
      *ids=output;
    }
    *change=2;
    return TCL_OK;
  } else if ( ARG0_IS_S("all") ) {
    if( (ret=convert_types_to_ids(NULL,output))<=0){
        Tcl_AppendResult(interp, "Error parsing keyword all. No particle found.\n", (char *)NULL);
        return TCL_ERROR;
    } else { 
      *ids=output;
    }
    *change=1;
    return TCL_OK;
  }

  Tcl_AppendResult(interp, "unknown keyword given to observable: ", argv[0] , (char *)NULL);
  return TCL_ERROR;
}


static int convert_types_to_ids(IntList * type_list, IntList * id_list){ 
      int i,j,n_ids=0,flag;
      sortPartCfg();
      for ( i = 0; i<n_total_particles; i++ ) {
         if(type_list==NULL) { 
		/* in this case we select all particles */
               flag=1 ;
         } else {  
                flag=0;
                for ( j = 0; j<type_list->n ; j++ ) {
                    if(partCfg[i].p.type == type_list->e[j])  flag=1;
	        }
         }
	 if(flag==1){
              realloc_intlist(id_list, id_list->n=n_ids+1);
	      id_list->e[n_ids] = i;
	      n_ids++;
	 }
      }
      return n_ids;
}

int observable_particle_velocities(void* idlist, double* A, unsigned int n_A) {
  unsigned int i;
  IntList* ids;
  sortPartCfg();
  ids=(IntList*) idlist;
  for ( i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_total_particles)
      return 1;
    A[3*i + 0] = partCfg[ids->e[i]].m.v[0]/time_step;
    A[3*i + 1] = partCfg[ids->e[i]].m.v[1]/time_step;
    A[3*i + 2] = partCfg[ids->e[i]].m.v[2]/time_step;
  }
  return 0;
}

#ifdef ELECTROSTATICS
int observable_particle_currents(void* idlist, double* A, unsigned int n_A) {
  unsigned int i;
  double charge;
  IntList* ids;
  sortPartCfg();
  ids=(IntList*) idlist;
  for ( i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_total_particles)
      return 1;
    charge = partCfg[ids->e[i]].p.q;
    A[3*i + 0] = charge * partCfg[ids->e[i]].m.v[0]/time_step;
    A[3*i + 1] = charge * partCfg[ids->e[i]].m.v[1]/time_step;
    A[3*i + 2] = charge * partCfg[ids->e[i]].m.v[2]/time_step;
  }
  return 0;
}
int observable_currents(void* idlist, double* A, unsigned int n_A) {
  unsigned int i;
  double charge;
  double j[3] = {0. , 0., 0. } ;
  IntList* ids;
  sortPartCfg();
  ids=(IntList*) idlist;
  for ( i = 0; i<ids->n; i++ ) {
    if (ids->e[i] > n_total_particles)
      return 1;
    charge = partCfg[ids->e[i]].p.q;
    j[0] += charge * partCfg[ids->e[i]].m.v[0]/time_step;
    j[1] += charge * partCfg[ids->e[i]].m.v[1]/time_step;
    j[2] += charge * partCfg[ids->e[i]].m.v[2]/time_step;
  }
  A[0]=j[0];
  A[1]=j[1];
  A[2]=j[2];
  return 0;
}
#endif

int observable_com_velocity(void* idlist, double* A, unsigned int n_A) {
/* TODO: this does not work with MASS ... */
  unsigned int i;
  double v_com[3] = { 0. , 0., 0. } ;
  IntList* ids;
  sortPartCfg();
  ids=(IntList*) idlist;
  for ( i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_total_particles)
      return 1;
    v_com[0] += partCfg[ids->e[i]].m.v[0]/time_step;
    v_com[1] += partCfg[ids->e[i]].m.v[1]/time_step;
    v_com[2] += partCfg[ids->e[i]].m.v[2]/time_step;
  }
  A[0]=v_com[0]/ids->n;
  A[1]=v_com[1]/ids->n;
  A[2]=v_com[2]/ids->n;
  printf("v_com %f %f %f\n", A[0], A[1], A[2]);
  return 0;
}

int observable_density_profile(void* pdata_, double* A, unsigned int n_A) {
  unsigned int i;
  int bin;
  double ppos[3];
  int img[3];
  IntList* ids;
  profile_data* pdata;
  sortPartCfg();
  pdata=(profile_data*) pdata_;
  ids=pdata->id_list;
    
  for ( i = 0; i<pdata->nbins; i++ ) {
    A[i]=0;
  }
  for ( i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_total_particles)
      return 1;
/* We use folded coordinates here */
    memcpy(ppos, partCfg[ids->e[i]].r.p, 3*sizeof(double));
    memcpy(img, partCfg[ids->e[i]].l.i, 3*sizeof(int));
    fold_position(ppos, img);
    bin= (int) floor( pdata->nbins*  (ppos[2]-pdata->minz)/(pdata->maxz-pdata->minz));
/* uncomment this line for unfolded coordinates */
//    bin= (int) floor( pdata->nbins*  (partCfg[ids->e[i]].r.p[2]-pdata->minz)/(pdata->maxz-pdata->minz));
    if (bin>=0 && bin < pdata->nbins) {
      A[bin] += 1./box_l[0]/box_l[1]/(pdata->maxz - pdata->minz)*pdata->nbins;
    }
  }
  return 0;
}

int observable_radial_density_profile(void* pdata_, double* A, unsigned int n_A) {
  unsigned int i;
  int bin;
  double ppos[3];
  double r;
  int img[3];
  double bin_volume;
  IntList* ids;
  radial_profile_data* pdata;
  sortPartCfg();
  pdata=(radial_profile_data*) pdata_;
  ids=pdata->id_list;
    
  for ( i = 0; i<pdata->nbins; i++ ) {
    A[i]=0;
  }
  for ( i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_total_particles)
      return 1;
/* We use folded coordinates here */
    memcpy(ppos, partCfg[ids->e[i]].r.p, 3*sizeof(double));
    memcpy(img, partCfg[ids->e[i]].l.i, 3*sizeof(int));
    fold_position(ppos, img);
    r=sqrt( (ppos[0]-pdata->center[0])*(ppos[0]-pdata->center[0])+(ppos[1]-pdata->center[1])*(ppos[1]-pdata->center[1]));
    bin= (int) floor( pdata->nbins*  r/pdata->maxr);

    if (bin>=0 && bin < pdata->nbins) {
      bin_volume=PI*((bin+1)*(bin+1)-bin*bin)*(pdata->maxr*pdata->maxr)/pdata->nbins/pdata->nbins*box_l[2];
      A[bin] += 1./bin_volume;
    }
  }
  return 0;
}

int observable_particle_positions(void* idlist, double* A, unsigned int n_A) {
  unsigned int i;
  IntList* ids;
  sortPartCfg();
  ids=(IntList*) idlist;
  for ( i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_total_particles)
      return 1;
      A[3*i + 0] = partCfg[ids->e[i]].r.p[0];
      A[3*i + 1] = partCfg[ids->e[i]].r.p[1];
      A[3*i + 2] = partCfg[ids->e[i]].r.p[2];
  }
  return 0;
}

int observable_structure_factor(void* params_p, double* A, unsigned int n_A) {
  int i,j,k,l,p;
  int order, order2, n;
  double twoPI_L, C_sum, S_sum, qr;
  observable_sf_params* params;
  params = (observable_sf_params*)params_p;
  order = params->order;
  order2=order*order;
  twoPI_L = 2*PI/box_l[0];
  
  sortPartCfg();

    for(p=0; p<n_A; p++) {
       A[p]   = 0.0;
    }

    l=0;
    //printf("n_A: %d, dim_sf: %d\n",n_A, params.dim_sf); fflush(stdout);
    for(i=-order; i<=order; i++) {
      for(j=-order; j<=order; j++) {
        for(k=-order; k<=order; k++) {
	  n = i*i + j*j + k*k;
	  if ((n<=order2) && (n>=1)) {
	    C_sum = S_sum = 0.0;
            //printf("l: %d, n: %d %d %d\n",l,i,j,k); fflush(stdout);
	    for(p=0; p<n_total_particles; p++) {
	      qr = twoPI_L * ( i*partCfg[p].r.p[0] + j*partCfg[p].r.p[1] + k*partCfg[p].r.p[2] );
	      C_sum+= partCfg[p].p.scattering_length * cos(qr);
	      S_sum-= partCfg[p].p.scattering_length * sin(qr);
	    }
            A[l]   =C_sum;
            A[l+1] =S_sum;
            l=l+2;
	  }
	}
      }
    }
    //printf("finished calculating sf\n"); fflush(stdout);
    return 0;
}

int observable_interacts_with (void* params_p, double* A, unsigned int n_A) {
  iw_params *params=(iw_params*)params_p;
  IntList* ids1;
  IntList* ids2;
  int i,j;
//  double dx,dy,dz;
  double dist2;
  double cutoff2=params->cutoff*params->cutoff;
  double pos1[3], pos2[3], dist[3];
  ids1=params->ids1;
  ids2=params->ids2;
  sortPartCfg();
  for ( i = 0; i<ids1->n; i++ ) {
    if (ids1->e[i] >= n_total_particles)
      return 1;
    pos1[0]=partCfg[ids1->e[i]].r.p[0];
    pos1[1]=partCfg[ids1->e[i]].r.p[1];
    pos1[2]=partCfg[ids1->e[i]].r.p[2];
    for ( j = 0; j<ids2->n; j++ ) {
      if (ids2->e[j] >= n_total_particles)
        return 1;
      A[i] = 0;
      pos2[0]=partCfg[ids2->e[j]].r.p[0];
      pos2[1]=partCfg[ids2->e[j]].r.p[1];
      pos2[2]=partCfg[ids2->e[j]].r.p[2];
      get_mi_vector(dist,pos1,pos2);
      dist2= dist[0]*dist[0] + dist[1]*dist[1] + dist[2]*dist[2];
      if(dist2<cutoff2) {
        A[i] = 1;
	break;
	// interaction found for i, go for next
      }
    }
  }
  return 0;
}



//int file_data_source_init(file_data_source* self, char* filename, IntList* columns) {
//  int counter=1;
//  char* token;
//  if (filename==0)
//    return 1;
//  self->f = fopen(filename, "r");
//  if (! self->f )
//    return 2;
//  fgets(self->last_line, MAXLINELENGTH, self->f);
//  while (self->last_line && self->last_line[0] == 35) {
//    fgets(self->last_line, MAXLINELENGTH, self->f);
//  }
//  if (!self->last_line)
//    return 3;
//// Now lets count the tokens in the first line
//  token=strtok(self->last_line, " \t\n");
//  while (token) {
////    printf("reading token **%s**\n", token);
//    token=strtok(NULL, " \t\n");
//    counter++;
//  }
//  self->n_columns = counter;
//  rewind(self->f);
//  self->data_left=1;
////  printf("I found out that your file has %d columns\n", self->n_columns);
//  if (columns !=0)
//    /// Here we would like to check if we can handle the desired columns, but this has to be implemented!
//    return -1;
//  return 0;
//}

int file_data_source_readline(void* xargs, double* A, int dim_A) {
//  file_data_source* self = xargs;
//  int counter=0;
//  char* token;
//  char* temp;
//
//  temp=fgets(self->last_line, MAXLINELENGTH, self->f);
//  while (temp!= NULL && self->last_line && self->last_line[0] == 35) {
//    temp=fgets(self->last_line, MAXLINELENGTH, self->f);
//  }
//  if (!self->last_line || temp==NULL) {
////    printf("nothing left\n");
//    self->data_left=0;
//    return 3;
//  }
//  token=strtok(self->last_line, " \t\n");
//  while (token) {
////    printf("reading token: ");
//    A[counter]=atof(token);
////    printf("%f ", A[counter]);
//    token=strtok(NULL, " \t\n");
//    counter++;
//    if (counter >= dim_A) {
////      printf("urgs\n");
//      return 4;
//    }
//  }
////  printf("\n");
//  return 0;
//}
//
//int tcl_input(void* data, double* A, unsigned int n_A) {
//  tcl_input_data* input_data = (tcl_input_data*) data;
//  int i, tmp_argc;
//  const char  **tmp_argv;
//  Tcl_SplitList(input_data->interp, input_data->argv[0], &tmp_argc, &tmp_argv);
//  // function prototype from man page:
//  // int Tcl_SplitList(interp, list, argcPtr, argvPtr)
//  if (tmp_argc < n_A) {
//    Tcl_AppendResult(input_data->interp, "Not enough arguments passed to analyze correlation update", (char *)NULL);
//    return 1;
//  }
//  for (i = 0; i < n_A; i++) {
//    if (Tcl_GetDouble(input_data->interp, tmp_argv[i], &A[i]) != TCL_OK) {
//      Tcl_AppendResult(input_data->interp, "error parsing argument ", input_data->argv[i],"\n", (char *)NULL);
//      return 1;
//    }
//  }
//  return 0;
}


int tclcommand_parse_profile(Tcl_Interp* interp, int argc, char** argv, int* change, int* dim_A, profile_data** pdata_) {
  int temp;
  *change=0;
  profile_data* pdata=(profile_data*)malloc(sizeof(profile_data));
  *pdata_ = pdata;
  pdata->id_list=0;
  pdata->minz=0;
  pdata->maxz=box_l[2];
  pdata->nbins=0;
  printf("\n");
  if (ARG0_IS_S("id") || ARG0_IS_S("type") || ARG0_IS_S("all")) {
    if (!parse_id_list(interp, argc, argv, &temp, &pdata->id_list )==TCL_OK) {
      Tcl_AppendResult(interp, "Error reading profile: Error parsing particle id information\n" , (char *)NULL);
      return TCL_ERROR;
    } else {
      *change+=temp;
      argc-=temp;
      argv+=temp;
    }
  } 
  if ( ARG0_IS_S("minz")){
    if (argc>1 && ARG1_IS_D(pdata->minz)) {
      argc-=2;
      argv+=2;
      *change+=2;
    } else {
      Tcl_AppendResult(interp, "Error in profile: could not read minz\n" , (char *)NULL);
      return TCL_ERROR;
    } 
  } 
  if ( ARG0_IS_S("maxz") ) {
    if (argc>1 && ARG1_IS_D(pdata->maxz)) {
      argc-=2;
      argv+=2;
      *change+=2;
    } else {
      Tcl_AppendResult(interp, "Error in profile: could not read maxz\n" , (char *)NULL);
      return TCL_ERROR;
    } 
  } 
  if (ARG0_IS_S("nbins")) {
    if (argc>1 && ARG1_IS_I(pdata->nbins)) {
      argc-=2;
      argv+=2;
      *change+=2;
    } else {
      Tcl_AppendResult(interp, "Error in profile: could not read nbins\n" , (char *)NULL);
      return TCL_ERROR;
    } 
  }
  
  temp=0;
  if (pdata->id_list==0) {
    Tcl_AppendResult(interp, "Error in profile: particle ids/types not specified\n" , (char *)NULL);
    temp=1;
  }
  if (pdata->nbins<1) {
    Tcl_AppendResult(interp, "Error in profile: nbins not specified\n" , (char *)NULL);
    temp=1;
  }
  if (temp)
    return TCL_ERROR;
  else
    return TCL_OK;
}

int tclcommand_parse_radial_profile(Tcl_Interp* interp, int argc, char** argv, int* change, int* dim_A, radial_profile_data** pdata_) {
  int temp;
  *change=0;
  radial_profile_data* pdata=(radial_profile_data*)malloc(sizeof(radial_profile_data));
  *pdata_ = pdata;
  pdata->id_list=0;
  pdata->maxr=1e100;
  pdata->minz=0;
  pdata->maxz=box_l[2];
  pdata->center[0]=1e100;pdata->center[1]=1e100;pdata->center[2]=1e100;
  pdata->nbins=0;
  if (argc < 1) {
    Tcl_AppendResult(interp, "UsageL radial_profile id $ids center $x $y $z maxr $r_max nbins $n\n" , (char *)NULL);
    return TCL_ERROR;
  }
  if (ARG0_IS_S("id") || ARG0_IS_S("type") || ARG0_IS_S("all")) {
    if (!parse_id_list(interp, argc, argv, &temp, &pdata->id_list )==TCL_OK) {
      Tcl_AppendResult(interp, "Error reading profile: Error parsing particle id information\n" , (char *)NULL);
      return TCL_ERROR;
    } else {
      *change+=temp;
      argc-=temp;
      argv+=temp;
    }
  } 
  if ( ARG0_IS_S("center")){
    if (argc>3 && ARG1_IS_D(pdata->center[0]) && ARG_IS_D(2,pdata->center[1]) && ARG_IS_D(3,pdata->center[2])) {
      argc-=4;
      argv+=4;
      *change+=4;
    } else {
      Tcl_AppendResult(interp, "Error in radial_profile: could not read center\n" , (char *)NULL);
      return TCL_ERROR;
    } 
  } 
  if ( ARG0_IS_S("maxr") ) {
    if (argc>1 && ARG1_IS_D(pdata->maxr)) {
      argc-=2;
      argv+=2;
      *change+=2;
    } else {
      Tcl_AppendResult(interp, "Error in radial_profile: could not read maxr\n" , (char *)NULL);
      return TCL_ERROR;
    } 
  } 
  if ( ARG0_IS_S("minz")){
    if (argc>1 && ARG1_IS_D(pdata->minz)) {
      argc-=2;
      argv+=2;
      *change+=2;
    } else {
      Tcl_AppendResult(interp, "Error in profile: could not read minz\n" , (char *)NULL);
      return TCL_ERROR;
    } 
  } 
  if ( ARG0_IS_S("maxz") ) {
    if (argc>1 && ARG1_IS_D(pdata->maxz)) {
      argc-=2;
      argv+=2;
      *change+=2;
    } else {
      Tcl_AppendResult(interp, "Error in profile: could not read maxz\n" , (char *)NULL);
      return TCL_ERROR;
    } 
  } 
  if (ARG0_IS_S("nbins")) {
    if (argc>1 && ARG1_IS_I(pdata->nbins)) {
      argc-=2;
      argv+=2;
      *change+=2;
    } else {
      Tcl_AppendResult(interp, "Error in radial_profile: could not read nbins\n" , (char *)NULL);
      return TCL_ERROR;
    } 
  }
  
  temp=0;
  if (pdata->id_list==0) {
    Tcl_AppendResult(interp, "Error in radial_profile: particle ids/types not specified\n" , (char *)NULL);
    temp=1;
  }
  if (pdata->center[0]>1e90) {
    Tcl_AppendResult(interp, "Error in radial_profile: center not specified\n" , (char *)NULL);
    temp=1;
  }
  if (pdata->maxr>1e90) {
    Tcl_AppendResult(interp, "Error in radial_profile: maxr not specified\n" , (char *)NULL);
    temp=1;
  }
  if (pdata->nbins<1) {
    Tcl_AppendResult(interp, "Error in radial_profile: nbins not specified\n" , (char *)NULL);
    temp=1;
  }
  if (temp)
    return TCL_ERROR;
  else
    return TCL_OK;
}


int sf_print_usage(Tcl_Interp* interp) {
  Tcl_AppendResult(interp, "\nusage: structure_factor order delta_t tau_max  tau_lin", (char *)NULL);
  return TCL_ERROR;
}
