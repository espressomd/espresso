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
MDINLINE void add_wall_force(Particle *p1, Particle *c_p, Constraint_wall *c, int type)
{
  int i;
  double dist, vec[3];
  IA_parameters *ia_params;

  ia_params=get_ia_param(p1->r.type, type);

  if (ia_params->LJ_eps > 0. ) {
    dist = -c->d;
    for(i=0;i<3;i++) dist += p1->r.p[i]*c->n[i];
  
    for(i=0;i<3;i++) vec[i] = c->n[i]*dist;
  
    if (dist > 0) {
        add_lj_pair_force(p1, c_p, ia_params, vec, dist);
    }
    else if ((dist * c->d)< 0) {
        fprintf(stderr,"CONSTRAINT: ERROR! part %d at (%.2e,%.2e,%.2e) out of constraint!\n",
		  p1->r.identity,p1->r.p[0],p1->r.p[1],p1->r.p[2]);
          fprintf(stderr,"%f %f %f %f %f\n", vec[0],vec[1],vec[2],dist,c->d);
	    errexit();    
    }
  }
}


MDINLINE void add_sphere_force(Particle *p1, Particle *c_p, Constraint_sphere *c, int type)
{
  int i;
  double dist, vec[3], c_dist;
  IA_parameters *ia_params;
 
  ia_params=get_ia_param(p1->r.type, type);

  if (ia_params->LJ_eps > 0. ) {
    dist = c->rad;
    c_dist=0.0;
  
    for(i=0;i<3;i++) {
        vec[i] =  c->pos[i] - p1->r.p[i] ;
        c_dist += SQR(vec[i]);
    }

    c_dist = sqrt(c_dist);
    dist -= c_dist;
  
    /* Assume that if there is a particle at the center, there is no interaction with the constraint */
    if (dist > 0) {
        for(i=0;i<3;i++) vec[i] *= (dist/c_dist);
        add_lj_pair_force(p1, c_p, ia_params, vec, dist);
    }
    else if (dist < 0) {
        fprintf(stderr,"CONSTRAINT: ERROR! part %d at (%.2e,%.2e,%.2e) out of constraint!\n",
		  p1->r.identity,p1->r.p[0],p1->r.p[1],p1->r.p[2]);
	    errexit();
    }
  }
}

MDINLINE void add_cylinder_force(Particle *p1, Particle *c_p, Constraint_cylinder *c, int type)
{
  int i;
  double dist, vec[3], c_dist;
  IA_parameters *ia_params;
  
  ia_params=get_ia_param(p1->r.type, type);

  if (ia_params->LJ_eps > 0. ) {
    /* put an infinite cylinder along z axis */
    dist = c->rad;
    c_dist = 0.0;
  
    for(i=0;i<2;i++) {
        vec[i] = c->pos[i] - p1->r.p[i];
        c_dist += SQR(vec[i]);
    }
    vec[2] = 0.;

    c_dist = sqrt(c_dist);
    dist -= c_dist;
    /* Assume that if there is a particle at the center, there is no interaction with the constraint */
    if ( dist > 0 ) {
        for(i=0;i<3;i++) vec[i] *= (dist/c_dist);
        add_lj_pair_force(p1, c_p, ia_params, vec, dist);
    }
    else if (dist < 0) {
        fprintf(stderr,"CONSTRAINT: ERROR! part %d at (%.2e,%.2e,%.2e) out of constraint!\n",
		  p1->r.identity,p1->r.p[0],p1->r.p[1],p1->r.p[2]);
	    errexit();
    }
  }
}

MDINLINE void add_constraints_forces(Particle *p1)
{
  int n;
  /* dummy particle for constraint. Later we will use this to
     get the force information on the constraint */
  Particle c_p;

  for(n=0;n<n_constraints;n++) {
    c_p.r.identity = -n;
    switch(constraints[n].type) {
    case CONSTRAINT_WAL: add_wall_force(p1, &c_p, &constraints[n].c.wal, constraints[n].particle_type); break;
    case CONSTRAINT_SPH: add_sphere_force(p1, &c_p, &constraints[n].c.sph, constraints[n].particle_type); break;
    case CONSTRAINT_CYL: add_cylinder_force(p1, &c_p, &constraints[n].c.cyl, constraints[n].particle_type); break;
    }
  }  
}
#endif


#endif
