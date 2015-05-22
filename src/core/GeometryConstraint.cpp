#include "GeometryConstraint.hpp"

namespace ConstraintClass {

  void GeometryConstraint::reflect_particle(Particle *p, const double *distance_vector, const double *folded_pos) {
    double vec[3];
    double norm; 

    memcpy(vec, distance_vector, 3*sizeof(double));

    norm=sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
    p->r.p[0] = p->r.p[0]-2*vec[0];
    p->r.p[1] = p->r.p[1]-2*vec[1];
    p->r.p[2] = p->r.p[2]-2*vec[2];

    /* vec seams to be the vector that points from the wall to the particle*/
    /* now normalize it */ 
    switch(reflection_type) {
    case REFLECTION_NORMAL:
      vec[0] /= norm;
      vec[1] /= norm;
      vec[2] /= norm;
      /* calculating scalar product - reusing var norm */
      norm = vec[0] *  p->m.v[0] + vec[1] * p->m.v[1] + vec[2] * p->m.v[2];
      /* now add twice the normal component to the velcity */
      p->m.v[0] = p->m.v[0]-2*vec[0]*norm; /* norm is still the scalar product! */
      p->m.v[1] = p->m.v[1]-2*vec[1]*norm;
      p->m.v[2] = p->m.v[2]-2*vec[2]*norm;
      break;
    case REFLECTION_NORMAL_TANGENTIAL:
      /* if bounce back, invert velocity */
      p->m.v[0] =-p->m.v[0]; /* norm is still the scalar product! */
      p->m.v[1] =-p->m.v[1];
      p->m.v[2] =-p->m.v[2];
      break;  
    case REFLECTION_NONE:
      break;
    }
  }
};
