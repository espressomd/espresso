
#include "statistics_average.h"


average* averages;
int n_averages;

int averages_add_empty_average() {
  n_averages++;
  realloc(averages, n_averages*sizeof(Average*));
  averages[n_averages-1]=(Average*) malloc(sizeof(Average));
  return 0;
}

int average_init(Average* self, observable* obs, int with_variance) {
  self->obs=obs;
  self->n_obs=obs->n;
  self->with_variance=with_variance;
  self->n_sweeps=0;
  self->average=(double*) malloc(self->n_obs*sizeof(double));
  if (self->with_variance)
    self->variance=(double*) malloc(self->n_obs*sizeof(double)); 
  self->ready=self->obs->ready;
}

int 


