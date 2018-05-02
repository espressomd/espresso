#ifndef CORE_ACCUMULATORS_ACCUMULATORBASE
#define CORE_ACCUMULATORS_ACCUMULATORBASE

namespace Accumulators {

class AccumulatorBase {
public:
  virtual void update() = 0;
};
}


#endif
