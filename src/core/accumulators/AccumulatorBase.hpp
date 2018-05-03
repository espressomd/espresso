#ifndef CORE_ACCUMULATORS_ACCUMULATORBASE
#define CORE_ACCUMULATORS_ACCUMULATORBASE

namespace Accumulators {

class AccumulatorBase {
public:
  AccumulatorBase(int delta_N)
      : delta_N(delta_N){};
  void auto_update();
  virtual void update() = 0;
  // Number of timesteps between automatic updates.
  int delta_N;
  int counter = 0;
};

inline void AccumulatorBase::auto_update() {
  if (counter % delta_N == 0) {
    update();
  }
  ++counter;
}
}

#endif
