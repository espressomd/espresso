#ifndef TRIEL_H
#define TRIEL_H

#include "utils.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "grid.hpp"
#include "lb.hpp"
#include "integrate.hpp"
#include <iostream>
#include <fstream>
using namespace std;

#ifdef TRIELASTIC

int triel_set_params(int bond_type, int ind1, int ind2, int ind3, double max, double ks, double ka);
int triel_reset_params(int bond_type, double lo, double lpo, double cospo, double sinpo, double Area0, double max, double ks, double ka);

static int time_step_p = -1;

//Use knowledge that the x-axis in rotates system is parallel to r(p1->p2) in original system;
//To find the corresponding unit vector to y in the rotated system, construct vector perpendicular to r(p1->p2);
inline void RotateForces(double f1_rot[2], double f2_rot[2], double f1[3], double f2[3], double r1[3], double r2[3]) {
	double xu[3];
	double y[3], yu[3];
	double sca;
	int i;
	unit_vector(r1,xu);
	sca = scalar(r2,xu);
	for(i=0; i<3; i++) {
		y[i] = r2[i] - sca * xu[i];
	}
	
	unit_vector(y,yu);

	f1[0] = f1_rot[0] * xu[0] + f1_rot[1] * yu[0]; f1[1] = f1_rot[0] * xu[1] + f1_rot[1] * yu[1]; f1[2] = f1_rot[0] * xu[2] + f1_rot[1] * yu[2]; 
    f2[0] = f2_rot[0] * xu[0] + f2_rot[1] * yu[0]; f2[1] = f2_rot[0] * xu[1] + f2_rot[1] * yu[1]; f2[2] = f2_rot[0] * xu[2] + f2_rot[1] * yu[2];
    
}

inline int calc_triel_force(Particle *p_ind1, Particle *p_ind2, Particle *p_ind3,
			      Bonded_ia_parameters *iaparams, double force1[3], double force2[3]) 
{
  double dxy, dxx, dyy, dyx; //Displacment gradient tensor matrix elements
  double gxy, gyx, gxx, gyy; //matrix calculations to get I1 and I2
    double e1, e2, i1, i2, i11, i12, i21, i22, i23, i24; // eigen values 1, 2 and D matrix elements  
    double A0, a1, a2, a3, b1, b2, b3;
    double l, lp, sinp, cosp;  // l, l' and cross and dot product between them 
    double vec1[3], vec2[3], vecpro[3];
    double f1_rot[2], f2_rot[2];
    
    //Calculate the current shape of the triangle (l,lp,cos(phi),sin(phi));
    //l = length between 1 and 3
    get_mi_vector(vec2, p_ind3->r.p, p_ind1->r.p);
    //vecsub(p_ind3->r.p,p_ind1->r.p,vec2);
    l = sqrt (sqrlen(vec2));
    //lp = lenght between 1 and 2
    get_mi_vector(vec1, p_ind2->r.p, p_ind1->r.p);
    //vecsub(p_ind2->r.p,p_ind1->r.p,vec1);
    lp = sqrt (sqrlen(vec1));
    //cosp / sinp angle functions between these vectors; calculated directly via the producs
    cosp = scalar(vec1,vec2)/(lp*l);
    vector_product(vec1, vec2, vecpro);
    sinp = sqrt(sqrlen(vecpro))/(l*lp);
    
    if( (lp-iaparams->p.triel.lpo > iaparams->p.triel.maxdist) ||  (l-iaparams->p.triel.lo > iaparams->p.triel.maxdist)) {
        return 1;
    }
    
    
    //Calculate forces in common plane (after assumed triangle rotation/translation in xy-plane);
    //Note that certain geometries and parameters (e.g. a3=0) can be used to speed the code up.
    //For now it will be played safe and done in detail. 
	dxx = lp/iaparams->p.triel.lpo;
	dxy = ((l*cosp/iaparams->p.triel.lo) - (lp*iaparams->p.triel.cospo/iaparams->p.triel.lpo)) / iaparams->p.triel.sinpo;
	dyx = 0.0;
	dyy = (l*sinp)/(iaparams->p.triel.lo * iaparams->p.triel.sinpo);
	gxx = SQR(dxx)+SQR(dyx);
	gxy = dxx*dxy + dyx*dyy;
	gyx = dxx*dxy + dyy*dyx;
	gyy = SQR(dxy) + SQR(dyy);
	i1 = (gxx + gyy) - 2;
	i2 = ((gxx * gyy) - (gxy * gyx)) - 1;
	i11 = 1.0; i12 = 1.0;
	i21 = gyy; i22 = -gyx; i23 = i22; i24 = gxx;
	
#ifdef TRIELNEOHOOKEAN	
	e1 = iaparams->p.triel.ks/6.0;
	e2 = (-1)*iaparams->p.triel.ks/(6.0*(i2+1.0)*(i2+1.0));
#else
	//Skalak-Law = Standard
	e1 = iaparams->p.triel.ks*(i1+1)/6.0;
	e2 = (-1)*iaparams->p.triel.ks/6.0 + iaparams->p.triel.ka*i2/6.0;
#endif	
    
    //For sake of better readability shorten the call for the triangle's constants:
	A0 = iaparams->p.triel.Area0; 
	a1 = iaparams->p.triel.a1; a2 = iaparams->p.triel.a2; a3 = iaparams->p.triel.a3;
	b1 = iaparams->p.triel.b1; b2 = iaparams->p.triel.b2; b3 = iaparams->p.triel.b3;
	
    f1_rot[0] = A0*((-1)*e1*((i11*2*a1*dxx)+(i12*2*b1*dxy))+ (-1)*e2*((i21*2*a1*dxx)+(i22*(a1*dxy+b1*dxx))+(i23*(a1*dxy+b1*dxx))+(i24*2*b1*dxy)));
    f1_rot[1] = A0*((-1)*e1*((i11*0.0)+(i12*2*b1*dyy))+ (-1)*e2*((i21*0.0)+(i22*a1*dyy)+(i23*a1*dyy)+(i24*2*b1*dyy)));
    f2_rot[0] = A0*((-1)*e1*((i11*2*a2*dxx)+(i12*2*b2*dxy))+ (-1)*e2*((i21*2*a2*dxx)+(i22*(a2*dxy+b2*dxx))+(i23*(a2*dxy+b2*dxx))+(i24*2*b2*dxy)));
    f2_rot[1] = A0*((-1)*e1*((i11*0.0)+(i12*2*b2*dyy))+ (-1)*e2*((i21*0.0)+(i22*a2*dyy)+(i23*a2*dyy)+(i24*2*b2*dyy)));
    //((int) ((double*)(fields[FIELD_TIMESTEP].data)));
      
    double vo[3], vpo[3];
    get_mi_vector(vo, p_ind3->m.v, p_ind1->m.v);
    get_mi_vector(vpo, p_ind2->m.v, p_ind1->m.v);
    double damping_coeff = 1; // kg/s
    double damping_force_o = (damping_coeff * (sqrt(sqrlen(vo))));
    double damping_force_po = (damping_coeff * (sqrt(sqrlen(vpo))));

     ofstream myfile ("c.dat", std::ofstream::out | std::ofstream::app);
    double v_int[3] = {0,0,0};
    double p_temp[3] = {p_ind1->r.p[0], p_ind1->r.p[1], p_ind1->r.p[0]};
    int time_step_c = step_counter;
    if( time_step_c > time_step_p ) { 
      
      // myfile << "Time Step: " << time_step_c << "\n";
      // myfile << "Particle's Velocity: " << p_ind1->m.v[0] << ", " <<  p_ind1->m.v[1] << ", " << p_ind1->m.v[2] << "\n";
      // myfile << "Particle's Relative Velocity: " << (p_ind2->m.v[0] - p_ind1->m.v[0]) << "\n";
      // lb_lbfluid_get_interpolated_velocity_lbtrace(p_temp,v_int, p_ind1->p.identity);
      // myfile << "Fluid Velocity around Particle: " << v_int[0] << ", " << v_int[1] << "," << v_int[2] << "\n";
      // myfile << "Force Rot0 : " << f1_rot[0] << ", " << f1_rot[1] << "\n";
      // myfile << "Force Rot1 : " << f2_rot[0] << ", " << f2_rot[1] << "\n";
      // myfile << "DampForce in X: " << damping_force_po << "\n";
      // myfile << "DampForce in Y: " << damping_force_o << "\n";
     myfile << "Particle's Old Velocity: " << p_ind1->m.v_old[0] << ", " <<  p_ind1->m.v_old[1] << ", " << p_ind1->m.v_old[2] << "\n";
     myfile << "Particle's Velocity: " << p_ind1->m.v[0] << ", " <<  p_ind1->m.v[1] << ", " << p_ind1->m.v[2] << "\n"; 
     
    }



    // f1_rot[0] += f1_rot[0] > 0 ? - damping_force_po : f1_rot[0] != 0 ? damping_force_po : 0;
    // f1_rot[1] += f1_rot[1] > 0 ? - damping_force_o : f1_rot[1] != 0 ? damping_force_o : 0;
    // f2_rot[0] += f2_rot[0] > 0 ? - damping_force_po : f2_rot[0] != 0 ? damping_force_po : 0;
    // f2_rot[1] += f2_rot[1] > 0 ? - damping_force_o : f2_rot[1] != 0 ? damping_force_o : 0;

    // f1_rot[0] += 1000;
    // f1_rot[1] += 1000;
    // f2_rot[0] += -1000;
    // f2_rot[1] += -1000;

    // f1_rot[0] *= -damping_force_po;
    // f1_rot[1] *= -damping_force_o;
    // f2_rot[0] *= damping_force_po;
    // f2_rot[1] *= damping_force_o;


    
    // printf("vx2, vx1: %d, %d", p_ind2->m.v[0], p_ind1->m.v[0]);
    // printf("vpo = %d \n", (sqrt(sqrlen(vpo)))); 

    //Rotate forces back into original position of triangle
    RotateForces(f1_rot,f2_rot,force1,force2, vec1, vec2); 

    // if( time_step_c > time_step_p  ) { 
    //   myfile << "Force Rot0 : " << f1_rot[0] << ", " << f1_rot[1] << "\n";
    //   myfile << "Force Rot1 : " << f2_rot[0] << ", " << f2_rot[1] << "\n";
    //   myfile << "Force0: " << force1[0] << ", " << force1[1] << ", " << force1[2] << "\n";
    //   myfile << "Force1: " << force2[0] << ", " << force2[1] << ", " << force2[2] << "\n\n";

    // }

    time_step_p = time_step_c;
   
   
    
    
    /*printf("Adding force to particles \n");
    force1[0] = 1;
    force2[0] = -1;*/
    
    return 0;
}

#endif
#endif
