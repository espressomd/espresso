#ifndef CONSTRAINT_H
#define CONSTRAINT_H

/** \file constraint.h
 *  Routines for handling of constraints. Implemented are walls, cylinders and spheres.
 *  Only activ if you define CONSTRAINTS in \ref config.h.
 *  see also \ref interaction_data.h
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:limbach@mpip-mainz.mpg.de">Hanjo</a>
*/


#ifdef CONSTRAINTS
MDINLINE void add_constraints_forces(Particle *p1)
{
  int i,n;
  double dist, c_dist, vec[3], frac2, frac6, fac;
  
  for(n=0;n<n_constraints;n++) {
    switch(constraints[n].type) {
      
    case CONSTRAINT_WAL: 
      /* calc distance */
      dist = -constraints[n].c.wal.d;
      for(i=0;i<3;i++) dist += p1->r.p[i]*constraints[n].c.wal.n[i];
      if(dist<constraints[n].LJ_cut) {
	/* force factor */
	frac2 = SQR( constraints[n].LJ_sig / dist);
	frac6 = frac2*frac2*frac2;
	fac   = 48.0 * constraints[n].LJ_eps * frac6*(frac6 - 0.5)*frac2 * dist;
	/* apply force */
	for(i=0;i<3;i++) p1->f[i] += fac * constraints[n].c.wal.n[i];
	/* fprintf(stderr,"CONSTRAINT: part %d at pos (%.2f,%.2f,%.2f) dist %.2f force (%.2e,%.2e,%.2e)\n",
	p1->r.identity,p1->r.p[0],p1->r.p[1],p1->r.p[2],dist,fac*constraints[n].c.wal.n[0],fac*constraints[n].c.wal.n[1],fac*constraints[n].c.wal.n[2]); */
      }
      break;
      
    case CONSTRAINT_SPH:
      /* calc distance */
      dist = constraints[n].c.sph.rad;
      c_dist=0.0;
      fac=0.0;
      for(i=0;i<3;i++) {
	vec[i] = constraints[n].c.sph.pos[i] - p1->r.p[i];
	c_dist += SQR(vec[i]);
      }
      c_dist = sqrt(c_dist);
      dist -= c_dist;
      for(i=0;i<3;i++) vec[i] /= c_dist;
      if(dist<constraints[n].LJ_cut) {
	/* force factor */
	frac2 = SQR(constraints[n].LJ_sig/dist);
	frac6 = frac2*frac2*frac2;
	fac   = 48.0 * constraints[n].LJ_eps * frac6*(frac6 - 0.5)*frac2 * dist;
	/* apply force */
	for(i=0;i<3;i++) p1->f[i] += fac * vec[i];
	
	/* fprintf(stderr,"CONSTRAINT: part %d at pos (%.2f,%.2f,%.2f) dist %.2f force (%.2e,%.2e,%.2e)\n",
	   p1->r.identity,p1->r.p[0],p1->r.p[1],p1->r.p[2],dist,fac*vec[0],fac*vec[1],fac*vec[2]); */
      }
      break;
    case CONSTRAINT_CYL:
      
      /* NOT IMPLEMENTED */

      break;
    }
  }  
}
#endif


#endif
