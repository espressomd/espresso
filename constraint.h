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
MDINLINE void add_wall_force(Particle *p1, Particle *c_p, Constraint_wall *c)
{
  int i;
  double dist, vec[3];
  IA_parameters *ia_params;

  ia_params=get_ia_param(p1->r.type, c_p->r.type);

  if (ia_params->LJ_eps > 0. ) {
    dist = -c->d;
    for(i=0;i<3;i++) dist += p1->r.p[i]*c->n[i];
  
    for(i=0;i<3;i++) vec[i] = c->n[i]*dist;
  
    if (dist > 0)
      add_lj_pair_force(p1, c_p, ia_params, vec, dist);
    else {
      fprintf(stderr,"CONSTRAINT WALL : ERROR! part %d at (%.2e,%.2e,%.2e) out of constraint!\n",
	      p1->r.identity,p1->r.p[0],p1->r.p[1],p1->r.p[2]);
      errexit();
    }
  }
}


MDINLINE void add_sphere_force(Particle *p1, Particle *c_p, Constraint_sphere *c)
{
  int i;
  double fac, dist, vec[3], c_dist;

  IA_parameters *ia_params;
 
  ia_params=get_ia_param(p1->r.type, c_p->r.type);
 
  if (ia_params->LJ_eps > 0. ) {
    c_dist=0.0;
    for(i=0;i<3;i++) {
      vec[i] = c->pos[i] - p1->r.p[i];
      c_dist += SQR(vec[i]);
    }

    c_dist = sqrt(c_dist);
    dist = c->rad - c_dist;

    if (dist > 0) {
      fac = dist / c_dist;
      for(i=0;i<3;i++) vec[i] *= fac;  
      add_lj_pair_force(p1, c_p, ia_params, vec, dist);
    }
    else {
      fprintf(stderr,"CONSTRAINT SPHERE: ERROR! part %d at (%.2e,%.2e,%.2e) out of constraint!\n",
	      p1->r.identity,p1->r.p[0],p1->r.p[1],p1->r.p[2]);
      errexit();
    }
  }
}

MDINLINE void add_cylinder_force(Particle *p1, Particle *c_p, Constraint_cylinder *c)
{
  int i;
  double fac, dist, vec[3], c_dist;
  IA_parameters *ia_params;
  
  ia_params=get_ia_param(p1->r.type, c_p->r.type);

  if (ia_params->LJ_eps > 0. ) {
    /* put an infinite cylinder along z axis */
    c_dist = 0.0;
    
    for(i=0;i<2;i++) {
      vec[i] = c->pos[i] - p1->r.p[i];
      c_dist += SQR(vec[i]);
    }
    vec[2] = 0.;

    c_dist = sqrt(c_dist);
    dist = c->rad - c_dist;

    if (dist > 0) {
      fac = dist / c_dist;
      for(i=0;i<2;i++) vec[i] *= fac;
      add_lj_pair_force(p1, c_p, ia_params, vec, dist);
    }
    else {
      fprintf(stderr,"CONSTRAINT CYLINDER: ERROR! part %d at (%.2e,%.2e,%.2e) out of constraint!\n",
	      p1->r.identity,p1->r.p[0],p1->r.p[1],p1->r.p[2]);
      errexit();
    }
  }
}

MDINLINE void add_rod_force(Particle *p1, Particle *c_p, Constraint_rod *c)
{
  int i;
  double fac, dist, vec[3], c_dist_2, c_dist;
  IA_parameters *ia_params;
  
  ia_params=get_ia_param(p1->r.type, c_p->r.type);

  /* also needed for coulomb */
  c_dist_2 = 0.0;
  for(i=0;i<2;i++) {
    vec[i] = p1->r.p[i] - c->pos[i];
    c_dist_2 += SQR(vec[i]);
  }
  vec[2] = 0.;

  /* charge stuff. This happens even if the particle does not feel the constraint. The electrostatic
     formulas for pbc highly dislike partial interactions anyways.
     THIS HAS TO BE DONE FIRST SINCE LJ CHANGES vec!!!
  */
  /*  fprintf(stderr, "%d: bj %f q %f l %f\n", this_node, coulomb.bjerrum, p1->r.q, c->lambda);*/
  if (coulomb.bjerrum > 0.0 && p1->r.q != 0.0 && c->lambda != 0.0) {
    fac = 2*coulomb.bjerrum*c->lambda*p1->r.q/c_dist_2;
    if (temperature > 0)
      fac *= temperature;
    p1->f[0]  += fac*vec[0];
    p1->f[1]  += fac*vec[1];
    c_p->f[0] -= fac*vec[0];
    c_p->f[1] -= fac*vec[1];
    /* fprintf(stderr, "%d: vec %f %f -> f %f %f\n", this_node, vec[0], vec[1],
       fac*vec[0], fac*vec[1]); */
  }

  if (ia_params->LJ_eps > 0. ) {
    /* put an infinite cylinder along z axis */
    c_dist = sqrt(c_dist_2);
    dist = c_dist - c->rad;

    if (dist > 0) {
      fac = dist / c_dist;
      for(i=0;i<2;i++) vec[i] *= fac;
      add_lj_pair_force(p1, c_p, ia_params, vec, dist);
    }
    else {
      fprintf(stderr,"CONSTRAINT ROD: ERROR! part %d at (%.2e,%.2e,%.2e) within rod!\n",
	      p1->r.identity,p1->r.p[0],p1->r.p[1],p1->r.p[2]);
      errexit();
    }
  }
}

MDINLINE void add_constraints_forces(Particle *p1)
{
  int n;

  for(n=0;n<n_constraints;n++) {
    switch(constraints[n].type) {
    case CONSTRAINT_WAL: add_wall_force(p1, &constraints[n].part_rep, &constraints[n].c.wal); break;
    case CONSTRAINT_SPH: add_sphere_force(p1, &constraints[n].part_rep, &constraints[n].c.sph); break;
    case CONSTRAINT_CYL: add_cylinder_force(p1, &constraints[n].part_rep, &constraints[n].c.cyl); break;
    case CONSTRAINT_ROD: add_rod_force(p1, &constraints[n].part_rep, &constraints[n].c.rod); break;
    }
  }
}
#endif

MDINLINE void init_constraint_forces()
{
  int n, i;
  
  for (n = 0; n < n_constraints; n++)
    for (i = 0; i < 3; i++)
      constraints[n].part_rep.f[i] = 0;
}

#endif
