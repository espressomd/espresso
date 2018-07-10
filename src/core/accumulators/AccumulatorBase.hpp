#ifndef CORE_ACCUMULATORS_ACCUMULATORBASE
#define CORE_ACCUMULATORS_ACCUMULATORBASE

namespace Accumulators {

class AccumulatorBase {
public:
  explicit AccumulatorBase(int delta_N=1)
      : m_delta_N(delta_N){};
  void auto_update();
  int &delta_N() {return m_delta_N;};
private:
  virtual void update() = 0;
  // Number of timesteps between automatic updates.
  int m_delta_N;
  int m_counter = 0;
};

inline void AccumulatorBase::auto_update() {
  if (m_counter % m_delta_N == 0) {
    update();
  }
  ++m_counter;
}
}

#endif
