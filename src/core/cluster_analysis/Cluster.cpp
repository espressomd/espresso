#include "particle_data.hpp"
#include "grid.hpp"
#ifdef GSL
  #include "gsl/gsl_fit.h"
#endif
#include "Vector.hpp"
#include <vector>

#include "Cluster.hpp" 

namespace ClusterAnalysis {

//Center of mass of an aggregate
Vector3d Cluster::center_of_mass() 
{
 Vector3d com;
  
  // The distances between the particles "folded", such that all distances
  // are smaller than box_l/2 in a periodic system. The 1st particle
  // of the cluster is arbitrarily chosen as reference.
  
  Vector3d reference_position;
  Vector3d dist_to_reference;
  for (int i=0; i<3; i++) {
    reference_position[i] = local_particles[particles[0]]->r.p[i];
  }
  
  double total_mass=0.;
  for (int pid : particles)  //iterate over all particle ids within a cluster
  {
    get_mi_vector(dist_to_reference, local_particles[pid]->r.p, reference_position); //add current particle positions
    com = com + dist_to_reference *local_particles[pid]->p.mass;
    total_mass += local_particles[pid]->p.mass;
  }

  // Normalize by numer of particles
  com = com * 1./total_mass;

  // Re-add reference position
  com = com +reference_position;

  // Fold into simulation box
  
  for (int i = 0; i < 3; i ++) {
    com[i] = fmod(com[i],box_l[i]);
  }
  return com;
}



double Cluster::longest_distance() {
  double ld=0.;
  for (auto a=particles.begin();a!=particles.end();a++) { 
    for (auto b=a;++b!=particles.end();) {
      double dist[3];
      get_mi_vector(dist, local_particles[*a]->r.p,local_particles[*b]->r.p);
       
      // Larger than previous largest distance?
      if (ld < sqrt(sqrlen(dist))) {
        ld=sqrt(sqrlen(dist)); //save bigger value as longest distance - ld
      }
    }
  }
  return ld;
}


//Radius of gyration
double Cluster::radius_of_gyration() {
  // Center of mass
  Vector3d com=center_of_mass();
  double sum_sq_dist=0.;
  for (auto const pid : particles) {
    double distance[3];
    get_mi_vector(distance, com.begin(), local_particles[pid]->r.p);
// calculate square length of this distance  
    sum_sq_dist += sqrlen(distance);
  }   
 
  return sqrt(sum_sq_dist/particles.size());
}



double Cluster::fractal_dimension(double dr, double& mean_sq_residual) {
#ifdef GSL
  Vector3d com = center_of_mass();  
// calculate Df using linear regression on the logarithms of diameters [__std::vector<double> diameters__] and num of particles [__std::vector<int> pcounts__] within the diameters
  
  std::vector<double> distances;

  for (auto const& it : particles) {
    double dist[3];
    get_mi_vector(dist, com.begin(), local_particles[it]->r.p); 
    distances.push_back(sqrt(sqrlen(dist))); //add distance from the current particle to the com in the distances vectors
  }
  
  // Increment radius in steps of dr and count particles within a sphere of that radius
  // Stop, when the sphere contains all particles of the cluster
  double rad = 0.0;
  int particles_in_sphere=0;
  const int cluster_size=particles.size();
  std::vector<double> log_diameters;
  std::vector<double> log_pcounts;
  while (particles_in_sphere < cluster_size) 
  { 

    // Count particles in the sphere
    particles_in_sphere=0;
    for (double d: distances) {
     if (d <= rad) particles_in_sphere+=1;
    }
    if (particles_in_sphere > 0) 
    {
      log_pcounts.push_back(log(particles_in_sphere)); 
      log_diameters.push_back(log(rad*2.0));
    }
    rad += dr;  //increase the radius
  }

// usage: Function: int gsl_fit_linear (const double * x, const size_t xstride, const double * y, const size_t ystride, size_t n, double * c0, double * c1, double * cov00, double * cov01, double * cov11, double * sumsq) 
  const int n=log_pcounts.size();
  double c0, c1, cov00, cov01, cov11, sumsq;
  gsl_fit_linear (&log_diameters.front(), 1, &log_pcounts.front(), 1, n, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);  
  mean_sq_residual =sumsq/log_diameters.size();
  return c1;
#else
  runtimeErrorMsg()<< "GSL (gnu scientific library) is required for fractal dimension calculation.";
#endif
}


}
