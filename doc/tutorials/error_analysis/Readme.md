# Tutorial: error analysis

Statistical analysis of simulation results.

## Learning objectives

### Part 1: introduction and binning analysis

* Give a brief overview of common measures of dispersion (standard deviation,
  confidence interval, standard error of the mean)
* Teach binning analysis and apply it on well-behaved data and on data where
  it doesn't converge (synthetic data is generated using the AR1 process)

### Part 2: autocorrelation analysis

* Teach autocorrelation analysis
* Integrate the ACF to determine the standard error of the mean
* Extract the autocorrelation time and use that information to decrease the
  observable sampling frequency (and thus reduce the amount of data and
  improve performance) and increase the simulation time for better statistics
