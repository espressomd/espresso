
int ObservableAverage::actual_update() {
    observable_average_container* data = (observable_average_container*) container;
    data->n_sweeps++;
    int error = data->reference_observable->calculate();
    if ( error != 0)
      return 1;
    double factor = 1 / (double) data->n_sweeps;
    for (int i =0; i<n; i++) {
      last_value[i] = (1-factor)*last_value[i] + factor*data->reference_observable->last_value[i];
    }
    return 0;
}

int ObservableAverage::reset() {
    observable_average_container* data = (observable_average_container*) container;
    data->n_sweeps=0;
    int error = data->reference_observable->calculate();
    for (int i =0; i<n; i++) {
      last_value[i] = 0;
    }
    return error;
}


