#include "OneParticleForce.h"

template<class Adapter>
class HarmonicPotential : public OneParticleForce<Adapter>
{
 public:
  void init(double _x, double _k) {
    x0 = _x;
    k = _k;
  };

  
  void run(void);
  void computeForces();
  void writeForces();

 private:
  double x0;
  double k;
};

template<class Adapter>
void HarmonicPotential<Adapter>::writeForces() {
}

template<class Adapter>
void HarmonicPotential<Adapter>::computeForces() {
  typename Adapter::RType::const_iterator it;
  for(it = this->a.getR().begin(); it != this->a.getR().end(); ++it) {
    
  }
}

template<class Adapter>
void HarmonicPotential<Adapter>::run(void) {
  this->computeForces();
  this->writeForces();
}
