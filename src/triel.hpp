#ifndef TRIEL_H
#define TRIEL_H

#include "utils.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "grid.hpp"

#ifdef TRIELASTIC

int triel_set_params(int bond_type, int ind1, int ind2, int ind3, double max, double ks, double ka);
int triel_reset_params(int bond_type, double lo, double lpo, double cospo, double sinpo, double Area0, double max, double ks, double ka);

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
    double dxy, dxx, dyy, dyx;
    double gxy, gyx, gxx, gyy;
    double e1, e2, i1, i2, i11, i12, i21, i22, i23, i24;
    double a1, a2, a3, b1, b2, b3;
    double l, lp, sinp, cosp;
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
	i1 = gxx + gyy - 2;
	i2 = gxx * gyy - gxy * gyx - 1;
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
	a1 = iaparams->p.triel.a1; a2 = iaparams->p.triel.a2; a3 = iaparams->p.triel.a3;
	b1 = iaparams->p.triel.b1; b2 = iaparams->p.triel.b2; b3 = iaparams->p.triel.b3;
	
    f1_rot[0] = iaparams->p.triel.Area0*((-1)*e1*((i11*2*a1*dxx)+(i12*2*b1*dxy))+ (-1)*e2*((i21*2*a1*dxx)+(i22*(a1*dxy+b1*dxx))+(i23*(a1*dxy+b1*dxx))+(i24*2*b1*dxy)));
    f1_rot[1] = iaparams->p.triel.Area0*((-1)*e1*((i11*0.0)+(i12*2*b1*dyy))+ (-1)*e2*((i21*0.0)+(i22*a1*dyy)+(i23*a1*dyy)+(i24*2*b1*dyy)));
    f2_rot[0] = iaparams->p.triel.Area0*((-1)*e1*((i11*2*a2*dxx)+(i12*2*b2*dxy))+ (-1)*e2*((i21*2*a2*dxx)+(i22*(a2*dxy+b2*dxx))+(i23*(a2*dxy+b2*dxx))+(i24*2*b2*dxy)));
    f2_rot[1] = iaparams->p.triel.Area0*((-1)*e1*((i11*0.0)+(i12*2*b2*dyy))+ (-1)*e2*((i21*0.0)+(i22*a2*dyy)+(i23*a2*dyy)+(i24*2*b2*dyy)));
    
    //Rotate forces back into original position of triangle
    RotateForces(f1_rot,f2_rot,force1,force2, vec1, vec2); 
    
    
    /*printf("Adding force to particles \n");
    force1[0] = 1;
    force2[0] = -1;*/
    
    return 0;
}

#endif
#endif
