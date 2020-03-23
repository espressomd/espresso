#include "particle_index.hpp"

std::vector<Particle *> local_particles;

void update_local_particles(const ParticleList *pl) {
  Particle *p = pl->part;
  int n = pl->n, i;
  for (i = 0; i < n; i++)
    set_local_particle_data(p[i].p.identity, &p[i]);
}

Particle *append_indexed_particle(ParticleList *l, Particle &&part) {
  auto const re = l->resize(l->n + 1);
  auto p = new (&(l->part[l->n - 1])) Particle(std::move(part));

  if (re)
    update_local_particles(l);
  else
    set_local_particle_data(p->p.identity, p);
  return p;
}
