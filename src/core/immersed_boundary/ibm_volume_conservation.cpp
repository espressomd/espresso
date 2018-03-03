
#include "config.hpp"

const int MaxNumIBM = 1000;
double VolumesCurrent[MaxNumIBM] = {0};

#ifdef IMMERSED_BOUNDARY

#include "particle_data.hpp"
#include "interaction_data.hpp"
#include "grid.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "immersed_boundary/ibm_volume_conservation.hpp"

#include "bond/Bond.hpp" // for Bond::Bond
#include "bond/BondType.hpp"// for Bond::BondType
#include "bond/IbmVolumeConservation.hpp" // for Bond::IbmVolumeConservation

#include "utils/make_unique.hpp" //for creating a unique ptr to a bond class object

// ****** Internal variables & functions ********
bool VolumeInitDone = false;

void CalcVolumes();
void CalcVolumeForce();

using std::ostringstream;

/************
  IBM_VolumeConservation
Calculate (1) volumes, (2) volume force and (3) add it to each virtual particle
This function is called from integrate_vv
 **************/

void IBM_VolumeConservation()
{
  // Calculate volumes
  CalcVolumes();
  
  CalcVolumeForce();
  
  // Center-of-mass output
  //if ( numWriteCOM > 0 )
//    IBM_CalcCentroids(frameNum, simTime);
}

void IBM_InitVolumeConservation()
{
  
  // Check since this function is called at the start of every integrate loop
  // Also check if volume has been set due to reading of a checkpoint
  if ( !VolumeInitDone /*&& VolumesCurrent[0] == 0 */)
  {
    
    // Calculate volumes
    CalcVolumes();
    bond_container.init_Vol_Con();
    
  };
  
  VolumeInitDone = true;
  
}

/****************
  IBM_VolumeConservation_ResetParams
 *****************/

int IBM_VolumeConservation_ResetParams(const int bond_type, const double volRef)
{
  
  // Check if bond exists and is of correct type
  if ( bond_type >= n_bonded_ia ) return ES_ERROR;

  auto vol_con_bond = dynamic_cast<Bond::IbmVolumeConservation*>
    (bond_container.get_Bond(bond_type));
  
  if(vol_con_bond){
    // Specific stuff
  // We need to set this here, since it is not re-calculated at the restarting of a sim as, e.g., triel
    vol_con_bond->ResetParams(VolRef);
    
  }
  else{
    printf("Dynamic Cast of Bond Failed: IBM Volume Conservation!\n");
    return ES_ERROR;
  }
  return ES_OK;
}

/***********
   IBM_VolumeConservation_SetParams
************/

int IBM_VolumeConservation_SetParams(const int bond_type, const int softID, const double kappaV)
{

  // Specific stuff
  if ( softID > MaxNumIBM) {
    printf("Error: softID (%d) is larger than MaxNumIBM (%d)\n", softID, MaxNumIBM);
    return ES_ERROR;
  };
  if ( softID < 0) {
    printf("Error: softID (%d) must be non-negative\n", softID);
    return ES_ERROR;
  };
  
  // NOTE: We cannot compute the reference volume here because not all interactions are setup
  // and thus we do not know which triangles belong to this softID
  // Calculate it later in the init function

  bond_container.set_bond_by_type(bond_type, Utils::make_unique<Bond::IbmVolumeConservation>
				  (softId, 0, kappaV));
  
  return ES_OK;
}


void CalcVolumes()
{
  
  // Partial volumes for each soft particle, to be summed up
  double tempVol[MaxNumIBM] = {0};
  
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
      int ibm_vol_con_bl_id = -1;
      //loop over ibmvolcon and get softid
      if(bond_container.ibm_vol_con_softID_loop(&p1, &softID, &ibm_vol_con_bl_id) == 2){
	exit(1);
      };
      // Second round for triel
      if ( softID > -1 && ibm_vol_con_bl_id > -1)
      {
	int j = 0;
        while ( j < p1.bl.n)
        {
          int bond_map_id = p1.bl.e[j];
	  Bond::Bond* current_bond = bond_container.get_Bond(bond_map_id);
	  //if bond cannot be found -> exit
	  if(current_bond == NULL){
	    runtimeErrorMsg() << "IBMVolConservation: bond type of atom "
			      << p1.p.identity << " unknown\n";
	    exit(1);
	  }
	  else{
	    int n_partners = current_bond->get_number_of_bond_partners();
	    Bond::BondType type = current_bond->get_Bond_Type();
	    if (type == Bond::BondType::BONDED_IA_IBM_TRIEL){
	      Bond::IbmVolumeConservation* IbmVolConBond = bond_container.get_IBM_Vol_Con_Bond(ibm_vol_con_bl_id);
	      //if bond cannot be found -> exit
	      if(IbmVolConBond == NULL){
		runtimeErrorMsg() << "IBMVolConservation: IBMVolCon_bond type of atom "
				  << p1.p.identity << " unknown\n";
		exit(1);
	      }
	      else{
		//if bond partners don't exist -> exit
		// give now bond map id of ibm triel bond for getting bond partners of ibm triel
		if(IbmVolConBond->calc_volumes(&p1, bond_map_id, tempVol)==2){
		  exit(1);
		};
	      };
	    };
	    // Iterate, increase by the number of partners of this bond + 1 for bond type
	    j += n_partners+1;
	  };//else
	};//while j < p1
      };//if softID
    };//for cell
  
  for (int i = 0; i < MaxNumIBM; i++) VolumesCurrent[i] = 0;
  
  // Sum up and communicate
  MPI_Allreduce(tempVol, VolumesCurrent, MaxNumIBM, MPI_DOUBLE, MPI_SUM, comm_cart);
  
  };//for
}

void CalcVolumeForce()
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
      int ibm_vol_con_bl_id = -1;
      //loop over ibmvolcon and get softid
      if(bond_container.ibm_vol_con_softID_loop(&p1, &softID, &ibm_vol_con_bl_id) == 2){
	exit(1);
      };      
      // Second round for triel
      if ( softID > -1 && ibm_vol_con_bl_id > -1)
      {
	int j = 0;
        while ( j < p1.bl.n)
        {
          const int bond_map_id = p1.bl.e[j];
	  Bond::Bond* current_bond = bond_container.get_Bond(bond_map_id);
	  //if bond cannot be found -> exit
	  if(current_bond == NULL){
	    runtimeErrorMsg() << "IBMVolConservation: bond type of atom "
			      << p1.p.identity << " unknown\n";
	    exit(1);
	  };
	  int n_partners = current_bond->get_number_of_bond_partners();
	  Bond::BondType type = current_bond->get_Bond_Type();
          if ( type == Bond::BondType::BONDED_IA_IBM_TRIEL )
          {
            //take bond_map_id of ibm triel and not ibmvolcon for bond partners
	    if(current_bond->add_bonded_force(&p1, bond_map_id)==2){
	      exit(1);
	    };
            
          }
          // Iterate, increase by the number of partners of this bond + 1 for bond type
          j += n_partners+1;
        }//while
      }//if soft id
    }//for cell
  }// for cells
}


#endif
