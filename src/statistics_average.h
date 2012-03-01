
#ifndef STATISTICS_AVERAGE_H
#define

#include "statistics_observable.h"



typedef struct {
  observable* obs;
  double* average;
  double* variance;
  int n_obs;
  int n_sweeps; 
  int ready;
  int with_variance;
} average;

/* The actual containers! */
extern average* averages;
extern int n_averages;

/* The Interface */
int averages_add_empty_average();
int average_init(average* self, observable* obs, int with_variance=1);
int average_update(average* self);
int average_print(average* self);
int average_print_formatted(average* self, Tcl_Interp* interp);


#endif

