// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
#ifndef GB_H
#define GB_H

/** \file gb.h
 *  Routines to calculate the Gay-Berne energy and force 
 *  for a pair of particles.
 *  \ref forces.c
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:antypov@mpip-mainz.mpg.de">Dmytro</a>
*/

#ifdef ROTATION

MDINLINE int gay_berne_set_params(int part_type_a, int part_type_b,
				  double eps, double sig, double cut,
				  double k1, double k2,
				  double mu, double nu)
{
  IA_parameters *data, *data_sym;

  make_particle_type_exist(part_type_a);
  make_particle_type_exist(part_type_b);
    
  data     = get_ia_param(part_type_a, part_type_b);
  data_sym = get_ia_param(part_type_b, part_type_a);

  if (!data || !data_sym) {
    return TCL_ERROR;
  }

  /* GB should be symmetrically */
  data->GB_eps    = data_sym->GB_eps    = eps;
  data->GB_sig    = data_sym->GB_sig    = sig;
  data->GB_cut    = data_sym->GB_cut    = cut;
  data->GB_k1     = data_sym->GB_k1     = k1;
  data->GB_k2     = data_sym->GB_k2     = k2;
  data->GB_mu     = data_sym->GB_mu     = mu;
  data->GB_nu     = data_sym->GB_nu     = nu;
 
   /* Calculate dependent parameters */

  data->GB_chi1 = data_sym->GB_chi1 = ((data->GB_k1*data->GB_k1) - 1) / ((data->GB_k1*data->GB_k1) + 1);
  data->GB_chi2 = data_sym->GB_chi2 = (pow(data->GB_k2,(1/data->GB_mu))-1)/(pow(data->GB_k2,(1/data->GB_mu))+1);

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);
  mpi_bcast_ia_params(part_type_b, part_type_a);

  return TCL_OK;
}

MDINLINE int printgbIAToResult(Tcl_Interp *interp, int i, int j)
{
  char buffer[TCL_DOUBLE_SPACE];
  IA_parameters *data = get_ia_param(i, j);

  Tcl_PrintDouble(interp, data->GB_eps, buffer);
  Tcl_AppendResult(interp, "gay-berne ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->GB_sig, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->GB_cut, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->GB_k1, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->GB_k2, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->GB_mu, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->GB_nu, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->GB_chi1, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->GB_chi2, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);

  return TCL_OK;
}

MDINLINE int gb_parser(Tcl_Interp * interp,
		       int part_type_a, int part_type_b,
		       int argc, char ** argv)
{
  double tmp;
  double eps, sig, cut;
  double k1, k2, mu, nu;
  int change;

  /* there are 9 parameters for gay-berne, but you read in only 7 of them.
     The rest is calculated in gay_berne_set_params.
  */

  if (argc < 8) {
    Tcl_AppendResult(interp, "gay-berne needs 7 parameters: "
		     "<gb_eps> <gb_sig> <gb_cut> <gb_k1> <gb_k2> <gb_mu> <gb_nu>",
		     (char *) NULL);
    return 0;
  }

  /* copy gay-berne parameters */
  if ((! ARG_IS_D(1, eps))   ||
      (! ARG_IS_D(2, sig))   ||
      (! ARG_IS_D(3, cut))   ||
      (! ARG_IS_D(4, k1 ))   ||
      (! ARG_IS_D(5, k2 ))   ||
      (! ARG_IS_D(6, mu ))   ||	
      (! ARG_IS_D(7, nu )    )) {
    Tcl_AppendResult(interp, "gay-berne needs 7 DOUBLE parameters: "
		     "<gb_eps> <gb_sig> <gb_cut> <gb_k1> <gb_k2> <gb_mu> <gb_nu>",
		     (char *) NULL);
    return 0;
  }
  change = 8;

  if (argc >= 10 && ARG_IS_D(8, tmp) && ARG_IS_D(9, tmp))
    change += 2;
  else
    Tcl_ResetResult(interp);

  if (gay_berne_set_params(part_type_a, part_type_b, eps, sig, cut, k1, k2, mu, nu) == TCL_ERROR) {
    Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
    return TCL_ERROR;
  }
  return change;
}

MDINLINE void add_gb_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
                                 double d[3], double dist)

{double a,b,c, X, Xcut,
	Brack,BrackCut,
	Bra12,Bra12Cut,
	u1x, u1y, u1z,
	u2x, u2y, u2z,	  
	E,E1,E2, Sigma,
	Brhi1,Brhi2,
	Plus1,Minus1,
	Plus2,Minus2,
	Koef1,Koef2,			/*  mu/E2  and  Sigma^3/2  */
	dU_dr, dU_da, dU_db, dU_dc,     /*  all derivatives        */
	FikX,FikY,FikZ,			/*  help for forces        */
	Gx,Gy,Gz;			/*  help for torques       */
 
  if (dist < ia_params->GB_cut) {  
     
  u1x = 2*(p1->r.quat[1]*p1->r.quat[3] + p1->r.quat[0]*p1->r.quat[2]);	     
  u1y = 2*(p1->r.quat[2]*p1->r.quat[3] - p1->r.quat[0]*p1->r.quat[1]);
  u1z =    p1->r.quat[0]*p1->r.quat[0] - p1->r.quat[1]*p1->r.quat[1] -
           p1->r.quat[2]*p1->r.quat[2] + p1->r.quat[3]*p1->r.quat[3] ;
	   
  u2x = 2*(p2->r.quat[1]*p2->r.quat[3] + p2->r.quat[0]*p2->r.quat[2]);	     
  u2y = 2*(p2->r.quat[2]*p2->r.quat[3] - p2->r.quat[0]*p2->r.quat[1]);
  u2z =    p2->r.quat[0]*p2->r.quat[0] - p2->r.quat[1]*p2->r.quat[1] -
           p2->r.quat[2]*p2->r.quat[2] + p2->r.quat[3]*p2->r.quat[3] ;   
    
        a = d[0]*u1x + d[1]*u1y + d[2]*u1z;
	b = d[0]*u2x + d[1]*u2y + d[2]*u2z;
	c =  u1x*u2x +  u1y*u2y +  u1z*u2z;	
       E1 = 1/sqrt(1-ia_params->GB_chi1*ia_params->GB_chi1*c*c);
    Plus1 = (a+b)/(1+ia_params->GB_chi1*c);
    Plus2 = (a+b)/(1+ia_params->GB_chi2*c);
   Minus1 = (a-b)/(1-ia_params->GB_chi1*c);
   Minus2 = (a-b)/(1-ia_params->GB_chi2*c);
    Brhi2 = (ia_params->GB_chi2/dist/dist)*(Plus2*(a+b) + Minus2*(a-b));
       E2 = 1-0.5*Brhi2;
        E = 4*ia_params->GB_eps*pow(E1,ia_params->GB_nu)*pow(E2,ia_params->GB_mu);
    Brhi1 = (ia_params->GB_chi1/dist/dist)*(Plus1*(a+b) + Minus1*(a-b));
    Sigma = 1/sqrt(1-0.5*Brhi1);
    Koef1 = ia_params->GB_mu/E2;
    Koef2 = Sigma*Sigma*Sigma*0.5;
    
        X = 1/(dist - Sigma + ia_params->GB_sig);
     Xcut = 1/(ia_params->GB_cut - Sigma + ia_params->GB_sig);
     
if (X < 1.25) { /* 1.25 corresponds to the interparticle penetration of 0.2 units of length.
                   If they are not that close, the GB forces and torques are calculated */
      
    Brack = X*X*X;
 BrackCut = Xcut*Xcut*Xcut;
    Brack = Brack*Brack;
 BrackCut = BrackCut*BrackCut;
 
    Bra12 = 6*Brack*X*(2*Brack-1);
 Bra12Cut = 6*BrackCut*Xcut*(2*BrackCut-1);
    Brack = Brack*(Brack-1);
 BrackCut = BrackCut*(BrackCut-1);
     
/*-------- Here we calculate derivatives -----------------------------*/ 

    dU_dr = E*(Koef1*Brhi2*(Brack-BrackCut)-Koef2*Brhi1*(Bra12-Bra12Cut)-Bra12*dist)/dist/dist;
    Koef1 = Koef1*ia_params->GB_chi2/dist/dist; 
    Koef2 = Koef2*ia_params->GB_chi1/dist/dist;
    dU_da = E*(Koef1*(Minus2+Plus2)*(BrackCut-Brack)+Koef2*(Plus1+Minus1)*(Bra12-Bra12Cut));
    dU_db = E*(Koef1*(Minus2-Plus2)*(Brack-BrackCut)+Koef2*(Plus1-Minus1)*(Bra12-Bra12Cut));
    dU_dc = E*((Brack-BrackCut)*(ia_params->GB_nu*E1*E1*ia_params->GB_chi1*ia_params->GB_chi1*c+
            0.5*Koef1*ia_params->GB_chi2*(Plus2*Plus2-Minus2*Minus2))-
            (Bra12-Bra12Cut)*0.5*Koef2*ia_params->GB_chi1*(Plus1*Plus1-Minus1*Minus1));
	   
/*--------------------------------------------------------------------*/
 
 FikX = -dU_dr*d[0] - dU_da*u1x - dU_db*u2x;
 FikY = -dU_dr*d[1] - dU_da*u1y - dU_db*u2y;
 FikZ = -dU_dr*d[2] - dU_da*u1z - dU_db*u2z;
 
    p1->f.f[0] += FikX;
    p1->f.f[1] += FikY;
    p1->f.f[2] += FikZ;
   
    p2->f.f[0] -= FikX;
    p2->f.f[1] -= FikY;
    p2->f.f[2] -= FikZ;
    
/* calculate torque:  torque = u_1 x G   */

  Gx = -dU_da*d[0] - dU_dc*u2x;
  Gy = -dU_da*d[1] - dU_dc*u2y;
  Gz = -dU_da*d[2] - dU_dc*u2z;
  
p1->f.torque[0]+= u1y*Gz - u1z*Gy;
p1->f.torque[1]+= u1z*Gx - u1x*Gz;
p1->f.torque[2]+= u1x*Gy - u1y*Gx;

/* calculate torque:  torque = u_2 x G     */

  Gx = -dU_db*d[0] - dU_dc*u1x;
  Gy = -dU_db*d[1] - dU_dc*u1y;
  Gz = -dU_db*d[2] - dU_dc*u1z;

p2->f.torque[0]+= u2y*Gz - u2z*Gy;
p2->f.torque[1]+= u2z*Gx - u2x*Gz;
p2->f.torque[2]+= u2x*Gy - u2y*Gx;
        }
else { /* the particles are too close to each other */
       Koef1  = 100;
       
    p1->f.f[0] += Koef1 * d[0];
    p1->f.f[1] += Koef1 * d[1];
    p1->f.f[2] += Koef1 * d[2];
   
    p2->f.f[0] -= Koef1 * d[0];
    p2->f.f[1] -= Koef1 * d[1];
    p2->f.f[2] -= Koef1 * d[2];
      }
  }
}

MDINLINE double gb_pair_energy(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist)
{double a,b,c, X, Xcut,
	Brack,BrackCut,
	u1x, u1y, u1z,
	u2x, u2y, u2z,	  
	E,E1,E2, Sigma,
	Plus1, Minus1,
	Plus2, Minus2;
	
  if (dist < ia_params->GB_cut) {  
    
  u1x = 2*(p1->r.quat[1]*p1->r.quat[3] + p1->r.quat[0]*p1->r.quat[2]);	     
  u1y = 2*(p1->r.quat[2]*p1->r.quat[3] - p1->r.quat[0]*p1->r.quat[1]);
  u1z =    p1->r.quat[0]*p1->r.quat[0] - p1->r.quat[1]*p1->r.quat[1] -
           p1->r.quat[2]*p1->r.quat[2] + p1->r.quat[3]*p1->r.quat[3] ;
	   
  u2x = 2*(p2->r.quat[1]*p2->r.quat[3] + p2->r.quat[0]*p2->r.quat[2]);	     
  u2y = 2*(p2->r.quat[2]*p2->r.quat[3] - p2->r.quat[0]*p2->r.quat[1]);
  u2z =    p2->r.quat[0]*p2->r.quat[0] - p2->r.quat[1]*p2->r.quat[1] -
           p2->r.quat[2]*p2->r.quat[2] + p2->r.quat[3]*p2->r.quat[3] ;   
    
        a = d[0]*u1x + d[1]*u1y + d[2]*u1z;
	b = d[0]*u2x + d[1]*u2y + d[2]*u2z;
	c =  u1x*u2x +  u1y*u2y +  u1z*u2z;
		
    Plus1 = (a+b)/(1+ia_params->GB_chi1*c);
    Plus2 = (a+b)/(1+ia_params->GB_chi2*c);
   Minus1 = (a-b)/(1-ia_params->GB_chi1*c);
   Minus2 = (a-b)/(1-ia_params->GB_chi2*c);
       E1 = 1/sqrt(1-ia_params->GB_chi1*ia_params->GB_chi1*c*c);
       E2 = 1-0.5*(ia_params->GB_chi2/dist/dist)*(Plus2*(a+b) + Minus2*(a-b));
        E = 4*ia_params->GB_eps*pow(E1,ia_params->GB_nu)*pow(E2,ia_params->GB_mu);  
    Sigma = 1/sqrt(1-0.5*(ia_params->GB_chi1/dist/dist)*(Plus1*(a+b) + Minus1*(a-b)));
        
        X = 1/(dist - Sigma + ia_params->GB_sig);
     Xcut = 1/(ia_params->GB_cut - Sigma + ia_params->GB_sig);
     
    Brack = X*X*X;
 BrackCut = Xcut*Xcut*Xcut;
    Brack = Brack*Brack;
 BrackCut = BrackCut*BrackCut;
    Brack = Brack*(Brack-1);
 BrackCut = BrackCut*(BrackCut-1);
   return E*(Brack-BrackCut);
  }  
return 0.0;
}

#endif
#endif
