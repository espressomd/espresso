#include "InteractionConstraint.hpp"
#include "forces_inline.hpp"

namespace ConstraintClass {

  void InteractionConstraint::add_force(Particle *p, const double *folded_pos) {

    int n, j;
    double dist, vec[3], force[3], torque1[3], torque2[3];

    IA_parameters *ia_params;

    ia_params=get_ia_param(p->p.type, part_rep.p.type);
    dist=0.;
    for (j = 0; j < 3; j++) {
      force[j] = 0;
#ifdef ROTATION
      torque1[j] = torque2[j] = 0;
#endif
    }

    if(checkIfInteraction(ia_params)) {
      m_shape.calculate_dist(folded_pos, &dist, vec);
      if ( dist > 0 ) {
	calc_non_bonded_pair_force(p, &part_rep,
				   ia_params,vec,dist,dist*dist, force,
				   torque1, torque2);
	// #ifdef TUNABLE_SLIP
	// 	  if (constraints[n].c.wal.tunable_slip == 1 ) {
	// 		  add_tunable_slip_pair_force(p1, &constraints[n].part_rep,ia_params,vec,dist,force);
	// 	  }
	// #endif
      }
      else if (penetrable && (dist <= 0)) {
	if ( (!only_positive) && ( dist < 0 ) ) {
          calc_non_bonded_pair_force(p, &part_rep,
                                     ia_params,vec,-1.0*dist,dist*dist, force,
                                     torque1, torque2);
        }
      }
      else {
	if(reflection_type != REFLECTION_NONE){
          reflect_particle(p, vec, folded_pos);
        } else {
          ostringstream msg;
          msg <<"wall constraint "<< n<<" violated by particle "<<p->p.identity;
          runtimeError(msg);
        }
      }
    }
  }
};
