#include "TabulatedBondLength.hpp"
#include "debug.hpp"

/** Calculate a tabulated bond length force with number type_num (see
    \ref Bonded_ia_parameters) between particles p1 and p2 and add it
    to the particle forces. The force acts in the direction of the
    connecting vector between the particles. For distances smaller
    than the tabulated range it uses a linear extrapolation based on
    the first two tabulated force values.
    Needs feature TABULATED compiled in (see \ref config.hpp). */
int Bond::TabulatedBondLength::calc_bonded_pair_force(Particle *p1, Particle *p2, double dx[3], 
			  double force[3]) const 
{

  int i;
  double fac, dist = sqrt(sqrlen(dx));

  if(dist > m_maxval)
    return 1;

  fac = bonded_tab_force_lookup(dist);
  
  for(i=0;i<3;i++)
    force[i] = fac*dx[i];

  ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: TAB BOND f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist,fac));
  ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: TAB BOND f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist,fac));

  return 0;

}



/** Calculate and return a tabulated bond length energy with number
    type_num (see \ref Bonded_ia_parameters) between particles p1 and
    p2. For distances smaller than the tabulated range it uses a
    quadratic extrapolation based on the first two tabulated force
    values and the first tabulated energy value. 
    Needs feature TABULATED compiled in (see \ref config.hpp). */
int Bond::TabulatedBondLength::calc_bonded_pair_energy(Particle *p1, Particle *p2, 
			   double dx[3], double *_energy) const
{

  double dist = sqrt(sqrlen(dx));

  if(dist > m_maxval)
    return 1;

  /* For distances smaller than the tabulated minimim take quadratic
     extrapolation from first two values of the force table! 
     This corresponds to the linear extrapolation of the force at that point.
     This sould not occur too often, since it is quite expensive!
  */
  if( dist <  m_minval) {
    double x0, b;
    b = (m_f[1]-m_f[0])*m_invstepsize;
    x0 = m_minval - m_f[0]/b;
    *_energy = ( (m_e[0] + 0.5*b*SQR(m_minval-x0)) -
		 0.5*b*SQR(dist-x0) );
  }
  else
    *_energy = bonded_tab_energy_lookup(dist);

  return 0;

}
