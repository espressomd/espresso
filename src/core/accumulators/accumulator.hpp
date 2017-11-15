#ifndef _ACCUMULATORS_ACCUMULATOR_H
#define _ACCUMULATORS_ACCUMULATOR_H


#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>

#include "observables/Observable.hpp"


namespace Accumulators {

class Accumulator
{
public:
  // The observable we want to accumulate.
  std::shared_pointer <Observables::Observable> m_obs;
  void add_sample();
  double get_mean();
  double get_variance();
private:
  using namespace boost::accumulators;
  accumulator_set<std::vector<double>, stats<tag::mean, tag::variance>> m_acc;

};

}

#endif