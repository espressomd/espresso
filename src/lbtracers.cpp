#include "lbtracers.hpp"
#include "lb.hpp"
#include "integrate.hpp"
#include "initialize.hpp"
#include "communication.hpp"

#ifdef LBTRACERS

// DEBUG: velocity fluctuations
/*double v2Avg = 0;
  double vAvg = 0;
  unsigned long int countAvg = 0;*/
// End DEBUG

//double time_step;
//extern double skin (?)

//Update Position ~ Euler
void update_mol_pos_particle(Particle *p) {
  //Do Euler for particle p; assume velocity has already been calculated 
  // & is stored in particle data
  int j;
  double skin2 = SQR(0.5 * skin);
		
  for(j=0;j<3;j++) 
    {
      // Euler
      p->r.p[j] = p->r.p[j] + p->m.v[j]*time_step;
    }
	    
  // Check if a particle might have crossed a box border (Verlet criterium); 
  //if possible resort_particles = 1
  //const double dist = distance2(p[j].r.p,p[j].l.p_old);
  const double dist2 = distance2( p->r.p, p->l.p_old);	
  if ( dist2 > skin2 ) { resort_particles = 1; }		
	
  //printf("Resort = %d\n");
	
	
  // DEBUG
  // Average fluct
  /*vAvg = ( vAvg*countAvg + p->m.v[0]) / (double) ( countAvg + 1);
    const double v2 = p->m.v[0]*p->m.v[0];
    v2Avg = ( v2Avg*countAvg + v2) / (double) ( countAvg + 1);
    countAvg++;
	
    if ( p->p.identity == 0) printf("vAvg = %e, v2Avg = %e\n", vAvg, v2Avg);*/
	
  // End DEBUG
}

//Update Velocity ~ Get interpolated velocity of LB
void update_mol_vel_particle(Particle *p) 
{
  int j;
  double p_temp[3];
	
  for(j=0;j<3;j++) {
    p_temp[j] = p->r.p[j];
  }
	
  if ( !(lattice_switch & LATTICE_LB_GPU) )
    {
      // Need to interpolate velocity here only for CPU
      // For GPU it is already stored
      double v_int[3] = {0,0,0};
		
      lb_lbfluid_get_interpolated_velocity_lbtrace(p_temp,v_int, p->p.identity);
      for ( j = 0; j < 3; j++)	
	p->m.v[j] = v_int[j];
    }
}

//Distribute forces
void distribute_mol_force() {
  //All forces these sites are subject to are influencing the LB fluid, not other
  //particles, therefore no forces need to be distributed here. => Do nothing
}


#endif
