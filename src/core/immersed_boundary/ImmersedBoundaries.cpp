
#include "config.hpp"

#ifdef IMMERSED_BOUNDARY

#include "particle_data.hpp"
#include "interaction_data.hpp"
#include "grid.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "ImmersedBoundaries.hpp" 





/************
  IBM_VolumeConservation
Calculate (1) volumes, (2) volume force and (3) add it to each virtual particle
This function is called from integrate_vv
 **************/

void ImmersedBoundaries::volume_conservation()
{
  // Calculate volumes
  calc_volumes();
  
  calc_volume_force();
  
  // Center-of-mass output
  //if ( numWriteCOM > 0 )
//    IBM_CalcCentroids(frameNum, simTime);
}

/************
  IBM_InitVolumeConservation
 *************/

void ImmersedBoundaries::init_volume_conservation()
{
  
  // Check since this function is called at the start of every integrate loop
  // Also check if volume has been set due to reading of a checkpoint
  if ( !VolumeInitDone )
  {
    
    // Calculate volumes
    calc_volumes();
    
    //numWriteCOM = 0;
    
    // Loop through all bonded interactions and check if we need to set the reference volume
    for (int i=0; i < bonded_ia_params.size(); i++)
    {
      if ( bonded_ia_params[i].type == BONDED_IA_IBM_VOLUME_CONSERVATION )
      {
        // This check is important because InitVolumeConservation may be called accidentally
        // during the integration. Then we must not reset the reference
        if ( bonded_ia_params[i].p.ibmVolConsParameters.volRef == 0 )
        {
          const int softID =bonded_ia_params[i].p.ibmVolConsParameters.softID;
          bonded_ia_params[i].p.ibmVolConsParameters.volRef = VolumesCurrent[softID];
        }
      }
      
    }
    
  }
  
  VolumeInitDone = true;
  
}


/****************
  IBM_VolumeConservation_ResetParams
 *****************/

int ImmersedBoundaries::volume_conservation_reset_params(const int bond_type, const double volRef)
{
  
  // Check if bond exists and is of correct type
  if ( bond_type >= bonded_ia_params.size() ) return ES_ERROR;
  if ( bonded_ia_params[bond_type].type != BONDED_IA_IBM_VOLUME_CONSERVATION ) return ES_ERROR;
  
  // Specific stuff
  // We need to set this here, since it is not re-calculated at the restarting of a sim as, e.g., triel
  bonded_ia_params[bond_type].p.ibmVolConsParameters.volRef = volRef;
  
  //Communicate this to whoever is interested
  mpi_bcast_ia_params(bond_type, -1);

  return ES_OK;
}

/***********
   IBM_VolumeConservation_SetParams
************/

int ImmersedBoundaries::volume_conservation_set_params(const int bond_type, const int softID, const double kappaV)
{
  // Create bond
  make_bond_type_exist(bond_type);
  
  // General bond parameters
  bonded_ia_params[bond_type].type = BONDED_IA_IBM_VOLUME_CONSERVATION;
  bonded_ia_params[bond_type].num = 0;        // This means that Espresso requires one bond partner. Here we simply ignore it, but Espresso cannot handle 0.
  
  // Specific stuff
  if ( softID > MaxNumIBM) { printf("Error: softID (%d) is larger than MaxNumIBM (%d)\n", softID, MaxNumIBM); return ES_ERROR; }
  if ( softID < 0) { printf("Error: softID (%d) must be non-negative\n", softID); return ES_ERROR; }
  
  bonded_ia_params[bond_type].p.ibmVolConsParameters.softID = softID;
  bonded_ia_params[bond_type].p.ibmVolConsParameters.kappaV = kappaV;
  bonded_ia_params[bond_type].p.ibmVolConsParameters.volRef = 0;
  // NOTE: We cannot compute the reference volume here because not all interactions are setup
  // and thus we do not know which triangles belong to this softID
  // Calculate it later in the init function
  
  //Communicate this to whoever is interested
  mpi_bcast_ia_params(bond_type, -1);
  
  return ES_OK;
}

/****************
   calc_volumes
Calculate partial volumes on all compute nodes
and call MPI to sum up
****************/

void ImmersedBoundaries::calc_volumes()
{
  
  // Partial volumes for each soft particle, to be summed up
  std::vector<double> tempVol(MaxNumIBM);
  
  // Loop over all particles on local node
  for (int c = 0; c < local_cells.n; c++)
  {
    const Cell *const cell = local_cells.cell[c];
    
    for (int i = 0; i < cell->n; i++)
    {
      Particle &p1 = cell->part[i];
      
      // Check if particle has a BONDED_IA_IBM_TRIEL and a BONDED_IA_IBM_VOLUME_CONSERVATION
      // Basically this loops over all triangles, not all particles
      // First round to check for volume conservation and virtual
      // Loop over all bonds of this particle
      // Actually j loops over the bond-list, i.e. the bond partners (see particle_data.hpp)
      int softID = -1;
      int j = 0;
      while ( j < p1.bl.n )
      {
        const int type_num = p1.bl.e[j];
        const Bonded_ia_parameters &iaparams = bonded_ia_params[type_num];
        const int type = iaparams.type;
        if ( type == BONDED_IA_IBM_VOLUME_CONSERVATION )
        {
          if ( p1.p.is_virtual) softID = iaparams.p.ibmVolConsParameters.softID;
          else { printf("Error. Encountered non-virtual particle with VOLUME_CONSERVATION_IBM\n"); exit(1); }
        }
        // Iterate, increase by the number of partners of this bond + 1 for bond type
        j += iaparams.num+1;
      }
      
      
      // Second round for triel
      if ( softID > -1 )
      {
        j = 0;
        while ( j < p1.bl.n)
        {
          const int type_num = p1.bl.e[j];
          const Bonded_ia_parameters &iaparams = bonded_ia_params[type_num];
          const int type = iaparams.type;
          
          if ( type == BONDED_IA_IBM_TRIEL )
          {
            // Our particle is the leading particle of a triel
            // Get second and third particle of the triangle
            Particle *p2 = local_particles[p1.bl.e[j+1]];
            if (!p2)
            {
              runtimeErrorMsg() << "{IBM_calc_volumes: 078 bond broken between particles "
                << p1.p.identity << " and " << p1.bl.e[j+1] << " (particles not stored on the same node)} ";
              return;
            }
            Particle *p3 = local_particles[p1.bl.e[j+2]];
            if (!p3)
            {
              runtimeErrorMsg() << "{IBM_calc_volumes: 078 bond broken between particles "
                << p1.p.identity << " and " << p1.bl.e[j+2] << " (particles not stored on the same node)} ";
              return;
            }
            
            // Unfold position of first node
            // this is to get a continuous trajectory with no jumps when box boundaries are crossed
            double x1[3] = { p1.r.p[0], p1.r.p[1], p1.r.p[2] };
            int img[3] = { p1.l.i[0], p1.l.i[1], p1.l.i[2] };
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
            
            tempVol[softID] += 1.0/6.0 * (-v321 + v231 + v312 - v132 - v213 + v123);
          }
          // Iterate, increase by the number of partners of this bond + 1 for bond type
          j += iaparams.num+1;
        }
      }
    }
  }
  
  for (int i = 0; i < MaxNumIBM; i++) VolumesCurrent[i] = 0;
  
  // Sum up and communicate
  MPI_Allreduce(&(tempVol.front()), &(VolumesCurrent.front()), MaxNumIBM, MPI_DOUBLE, MPI_SUM, comm_cart);
  
}

/*****************
  calc_volume_force
Calculate and add the volume force to each node
*******************/

void ImmersedBoundaries::calc_volume_force()
{
  // Loop over all particles on local node
  for (int c = 0; c < local_cells.n; c++)
  {
    const Cell *const cell = local_cells.cell[c];
    
    for (int i = 0; i < cell->n; i++)
    {
      Particle &p1 = cell->part[i];
      
      // Check if particle has a BONDED_IA_IBM_TRIEL and a BONDED_IA_IBM_VOLUME_CONSERVATION
      // Basically this loops over all triangles, not all particles
      // First round to check for volume conservation and virtual
      // Loop over all bonds of this particle
      // Actually j loops over the bond-list, i.e. the bond partners (see particle_data.hpp)
      int softID = -1;
      double volRef = 0.;
      double kappaV = 0.;
      int j = 0;
      while ( j < p1.bl.n )
      {
        const int type_num = p1.bl.e[j];
        const Bonded_ia_parameters &iaparams = bonded_ia_params[type_num];
        const int type = iaparams.type;
        if ( type == BONDED_IA_IBM_VOLUME_CONSERVATION )
        {
          if ( !p1.p.is_virtual) { printf("Error. Encountered non-virtual particle with VOLUME_CONSERVATION_IBM\n"); exit(1); }
          softID = iaparams.p.ibmVolConsParameters.softID;
          volRef = iaparams.p.ibmVolConsParameters.volRef;
          kappaV = iaparams.p.ibmVolConsParameters.kappaV;
        }
        // Iterate, increase by the number of partners of this bond + 1 for bond type
        j += iaparams.num+1;
      }
      
      
      // Second round for triel
      if ( softID > -1 )
      {
        j = 0;
        while ( j < p1.bl.n)
        {
          const int type_num = p1.bl.e[j];
          const Bonded_ia_parameters &iaparams = bonded_ia_params[type_num];
          const int type = iaparams.type;
          
          if ( type == BONDED_IA_IBM_TRIEL )
          {
            // Our particle is the leading particle of a triel
            // Get second and third particle of the triangle
            Particle *p2 = local_particles[p1.bl.e[j+1]];
            Particle *p3 = local_particles[p1.bl.e[j+2]];
            
            // Unfold position of first node
            // this is to get a continuous trajectory with no jumps when box boundaries are crossed
            double x1[3] = { p1.r.p[0], p1.r.p[1], p1.r.p[2] };
            int img[3] = { p1.l.i[0], p1.l.i[1], p1.l.i[2] };
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
            /*                        const double fact = kappaV * 1/6. * (IBMVolumesCurrent[softID] - volRef) / IBMVolumesCurrent[softID];
             
             double x2[3];
            double x3[3];
            
            for (int i=0; i < 3; i++)
            {
              x2[i] = x1[i] + a12[i];
              x3[i] = x1[i] + a13[i];
            }
             
             double n[3];
             vector_product(x3, x2, n);
             for (int k=0; k < 3; k++) p1.f.f[k] += fact*n[k];
             vector_product(x1, x3, n);
             for (int k=0; k < 3; k++) p2->f.f[k] += fact*n[k];
             vector_product(x2, x1, n);
             for (int k=0; k < 3; k++) p3->f.f[k] += fact*n[k];*/
            
            
            // This is Dupin 2008. I guess the result will be very similar as the code above
            double n[3];
            vector_product(a12, a13, n);
            const double ln = sqrt( n[0]*n[0] + n[1]*n[1] + n[2]*n[2] );
            const double A = 0.5 * ln;
            const double fact = kappaV * (VolumesCurrent[softID] - volRef) / VolumesCurrent[softID];
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
              p1.f.f[k] += force[k];
              p2->f.f[k] += force[k];
              p3->f.f[k] += force[k];
            }
            
          }
          // Iterate, increase by the number of partners of this bond + 1 for bond type
          j += iaparams.num+1;
        }
      }
    }
  }
}


#endif
