#include "OverlapBondLength.hpp"
#include "debug.hpp"
#include "interaction_data.hpp"

/** Computes the two body overlapped bonded force.
    Adds this force to the particle forces in forces.hpp (see \ref
   tclcommand_inter).
    @param p1        Pointer to first particle.
    @param p2        Pointer to second/middle particle.
    @param iaparams  bond type number of the angle interaction (see \ref
   tclcommand_inter).
    @param dx        particle distance vector
    @param force     returns force of particle 1
    @return 0.
    Needs feature OVERLAPPED compiled in (see \ref config.hpp).
*/
int Bond::OverlapBondLength::calc_bonded_pair_force(Particle *p1, Particle *p2, double dx[3], 
						    double force[3]) const 
{

  int i;
  double fac;

  int ig;
  double ev;
  double prob = 0.;
  double Apart = 0.;

  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);

  /* compute fac = - 1/r * dU/dr */
  /* prob = sum_(i=1,N) [ a_i * exp(-(r -b_i)^2/c_i^2))] */
  /* Apart = sum_(i=1,N) [ a_i * exp(-(r -b_i)^2/c_i^2)) * (x-b_i) / c_i^2] */
  /* fac = - 2.0 * ( 1/r + Apart / prob ) */
  for (ig = 0; ig < m_noverlaps; ig++) {
    ev = (dist - m_para_b[ig]) /
         m_para_c[ig];
    prob = prob + m_para_a[ig] * exp(-1.0 * ev * ev);
    Apart =
        Apart +
        m_para_a[ig] * exp(-1.0 * ev * ev) *
            (dist - m_para_b[ig]) /
            (m_para_c[ig] * m_para_c[ig]);
  }
  fac = -2. * (1 / dist + Apart / prob);
  fac /= dist;

  /* compute force */
  for (i = 0; i < 3; i++)
    force[i] = fac * dx[i];

  ONEPART_TRACE(if (p1->p.identity == check_id)
                    fprintf(stderr, "%d: OPT: HARMONIC f = (%.3e,%.3e,%.3e) "
                                    "with part id=%d at dist %f fac %.3e\n",
                            this_node, p1->f.f[0], p1->f.f[1], p1->f.f[2],
                            p2->p.identity, dist2, fac));
  ONEPART_TRACE(if (p2->p.identity == check_id)
                    fprintf(stderr, "%d: OPT: HARMONIC f = (%.3e,%.3e,%.3e) "
                                    "with part id=%d at dist %f fac %.3e\n",
                            this_node, p2->f.f[0], p2->f.f[1], p2->f.f[2],
                            p1->p.identity, dist2, fac));

  return 0;

}


/** Computes the two body overlapped angle interaction energy (see \ref
   tclcommand_inter, \ref tclcommand_analyze).
    @param p1        Pointer to first particle.
    @param p2        Pointer to second/middle particle.
    @param iaparams  bond type number of the angle interaction (see \ref
   tclcommand_inter).
    @param dx        particle distance vector
    @param _energy   returns energy of this interaction
    @return 0.
    Needs feature OVERLAPPED compiled in (see \ref config.hpp).
*/
int Bond::OverlapBondLength::calc_bonded_pair_energy(Particle *p1, Particle *p2, 
						    double dx[3], double *_energy) const
{

  int ig;
  double ev;
  double prob = 0.;

  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);

  /*compute overlapped bond energy */
  /*prob = sum_(i=1,N) [ a_i * exp(-(r -b_i)^2/c_i^2))] */
  /*pot  = -1.0 * Log { prob / r^2 } */
  for (ig = 0; ig < m_noverlaps; ig++) {
    ev = (dist - m_para_b[ig]) /
         m_para_c[ig];
    prob = prob + m_para_a[ig] * exp(-1.0 * ev * ev);
  }
  *_energy = -log(prob / dist2);

  return 0;

}
