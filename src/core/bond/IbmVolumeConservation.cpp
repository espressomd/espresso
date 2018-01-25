#include "immersed_boundary/ibm_volume_conservation.hpp" //VolumesCurrent Variable
#include "IbmVolumeConservation.hpp"
#include "grid.hpp" //get_mi_vector

//gets bl_id of ibm triel 
int Bond::IbmVolumeConservation::add_bonded_force(Particle *p1, int bl_id)
{

  // Our particle is the leading particle of a triel
  // Get second and third particle of the triangle
  Particle *p2, *p3 = NULL;
  if(auto bond_partners =  get_n_bond_partners<2>(p1, bl_id)){
    p2 = (*bond_partners)[0];
    p3 = (*bond_partners)[1];

    // Unfold position of first node
    // this is to get a continuous trajectory with no jumps when box boundaries are crossed
    double x1[3] = { p1->r.p[0], p1->r.p[1], p1->r.p[2] };
    int img[3] = { p1->l.i[0], p1->l.i[1], p1->l.i[2] };
    unfold_position(x1,img);
            
    // Unfolding seems to work only for the first particle of a triel
    // so get the others from relative vectors considering PBC
    double a12[3];
    get_mi_vector(a12, p2->r.p, x1);
    double a13[3];
    get_mi_vector(a13, p3->r.p, x1);
            

            
    // Now we have the true and good coordinates
    // Compute force according to eq. C.46 Krüger thesis
    // It is the same as deriving Achim's equation w.r.t x
    /*                        const double fact = m_kappaV * 1/6. * (IBMVolumesCurrent[m_softID] - m_volRef) / IBMVolumesCurrent[m_softID];
             
			      double x2[3];
			      double x3[3];
            
			      for (int i=0; i < 3; i++)
			      {
			      x2[i] = x1[i] + a12[i];
			      x3[i] = x1[i] + a13[i];
			      }
             
			      double n[3];
			      vector_product(x3, x2, n);
			      for (int k=0; k < 3; k++) p1->f.f[k] += fact*n[k];
			      vector_product(x1, x3, n);
			      for (int k=0; k < 3; k++) p2->f.f[k] += fact*n[k];
			      vector_product(x2, x1, n);
			      for (int k=0; k < 3; k++) p3->f.f[k] += fact*n[k];*/
            

    // This is Dupin 2008. I guess the result will be very similar as the code above
    double n[3];
    vector_product(a12, a13, n);
    const double ln = sqrt( n[0]*n[0] + n[1]*n[1] + n[2]*n[2] );
    const double A = 0.5 * ln;
    const double fact = m_kappaV * (VolumesCurrent[m_softID] - m_volRef) / VolumesCurrent[m_softID];
    double nHat[3];
    nHat[0] = n[0] / ln;
    nHat[1] = n[1] / ln;
    nHat[2] = n[2] / ln;
            
    double force[3];
    force[0] = -fact * A * nHat[0];
    force[1] = -fact * A * nHat[1];
    force[2] = -fact * A * nHat[2];
            
    // Add forces
    for (int k=0; k < 3; k++)
      {
	p1->f.f[k] += force[k];
	p2->f.f[k] += force[k];
	p3->f.f[k] += force[k];
      };

  return 0;

  }
  else{
    return 2;
  };

}

// no energy contribution
int Bond::IbmVolumeConservation::add_bonded_energy(Particle *p1, int bl_id)
{
  return 0;
}

int Bond::IbmVolumeConservation::get_soft_ID(Particle *p1, int bl_id, int *softID, int *bond_map_id)
{

#ifdef IMMERSED_BOUNDARY
  if( p1->p.isVirtual){
    *softID = m_softID;
    *bond_map_id = bl_id;
    return 0;
  }
  else { 
    printf("Error. Encountered non-virtual particle with VOLUME_CONSERVATION_IBM\n");
    return 2; // see BondContainer loop -> it return immediately and we can call exit(1)
  };
#endif
#ifndef IMMERSED_BOUNDARY
  return 0;
#endif
  
}

int Bond::IbmVolumeConservation::calc_volumes(Particle *p1, int bl_id, double *tempVol)
{

  // Our particle is the leading particle of a triel
  // Get second and third particle of the triangle
  Particle *p2, *p3 = NULL;
  if(auto bond_partners =  get_n_bond_partners<2>(p1, bl_id)){            
    p2 = (*bond_partners)[0];
    p3 = (*bond_partners)[1];
    // Unfold position of first node
    // this is to get a continuous trajectory with no jumps when box boundaries are crossed
    double x1[3] = { p1->r.p[0], p1->r.p[1], p1->r.p[2] };
    int img[3] = { p1->l.i[0], p1->l.i[1], p1->l.i[2] };
    unfold_position(x1,img);
            
    // Unfolding seems to work only for the first particle of a triel
    // so get the others from relative vectors considering PBC
    double a12[3];
    get_mi_vector(a12, p2->r.p, x1);
    double a13[3];
    get_mi_vector(a13, p3->r.p, x1);
            
    double x2[3];
    double x3[3];
            
    for (int i=0; i < 3; i++)
      {
	x2[i] = x1[i] + a12[i];
	x3[i] = x1[i] + a13[i];
      }
            
    // Volume of this tetrahedron
    // See Cha Zhang et.al. 2001, doi:10.1109/ICIP.2001.958278
    // http://research.microsoft.com/en-us/um/people/chazhang/publications/icip01_ChaZhang.pdf
    // The volume can be negative, but it is not necessarily the "signed volume" in the above paper
    // (the sign of the real "signed volume" must be calculated using the normal vector; the result of the calculation here
    // is simply a term in the sum required to calculate the volume of a particle). Again, see the paper.
    // This should be equivalent to the formulation using vector identities in Krüger thesis
            
    const double v321 = x3[0] * x2[1] * x1[2];
    const double v231 = x2[0] * x3[1] * x1[2];
    const double v312 = x3[0] * x1[1] * x2[2];
    const double v132 = x1[0] * x3[1] * x2[2];
    const double v213 = x2[0] * x1[1] * x3[2];
    const double v123 = x1[0] * x2[1] * x3[2];
            
    tempVol[m_softID] += 1.0/6.0 * (-v321 + v231 + v312 - v132 - v213 + v123);
    return 0;
  }
  else{
    return 2;
  };

}
