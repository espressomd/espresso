#include "OifLocalForces.hpp"
#include "grid.hpp" // get_mi_vector
#include "utils.hpp" // unfold position
#include "interaction_data.hpp"

inline double KS(double lambda){ // Defined by (19) from Dupin2007
    double res;
    res = (pow(lambda,0.5) + pow(lambda,-2.5))/(lambda + pow(lambda,-3.));
    return res;
}

//force calculation
/** Computes the local forces (Dupin2007) and adds them
    to the particle forces (see \ref tclcommand_inter). 
    @param p1,p2,p3     Pointers to particles of triangle 1.
    @param p2,p3,p4     Pointers to particles of triangle 2.
    (triangles have particles p2 and p3 in common)
    @return 0
*/
int Bond::OifLocalForces::calc_bonded_four_particle_force(Particle *p1, Particle *p2, Particle *p3, 
				   Particle *p4, double force[3], double force2[3], 
							 double force3[3], double force4[3]) const {
	int i, img[3];
	double fp1[3],fp2[3],fp3[3],fp4[3];
	double AA[3],BB[3],CC[3];
    double n1[3],n2[3],dn1,dn2,phi,aa;
    double dx[3],fac,dr,len2,len,lambda;
    double A,h[3],rh[3],hn;

	// first find out which particle out of p1, p2 (possibly p3, p4) is not a ghost particle. In almost all cases it is p2, however, it might be other one. we call this particle reference particle.
	if (p2->l.ghost != 1) {
		//unfold non-ghost particle using image, because for physical particles, the structure p->l.i is correctly set
		memmove(fp2, p2->r.p, 3*sizeof(double));
		memmove(img, p2->l.i, 3*sizeof(int));
		unfold_position(fp2,img);
		// other coordinates are obtained from its relative positions to the reference particle
		get_mi_vector(AA, p1->r.p, fp2);
		get_mi_vector(BB, p3->r.p, fp2);
		get_mi_vector(CC, p4->r.p, fp2);
		for (i=0; i<3; i++) {
            fp1[i] = fp2[i] + AA[i];
            fp3[i] = fp2[i] + BB[i];
            fp4[i] = fp2[i] + CC[i];
        }
	} else {
		// in case  particle p2 is a ghost particle
		if (p1->l.ghost != 1) {
			memmove(fp1, p1->r.p, 3*sizeof(double));
			memmove(img, p1->l.i, 3*sizeof(int));
			unfold_position(fp1,img);
			get_mi_vector(AA, p2->r.p, fp1);
			get_mi_vector(BB, p3->r.p, fp1);
			get_mi_vector(CC, p4->r.p, fp1);
			for (i=0; i<3; i++) {
                fp2[i] = fp1[i] + AA[i];
                fp3[i] = fp1[i] + BB[i];
                fp4[i] = fp1[i] + CC[i];
            }
		} else {
			// in case the first and the second particle are ghost particles
			if (p3->l.ghost != 1) {
				memmove(fp3, p3->r.p, 3*sizeof(double));
				memmove(img, p3->l.i, 3*sizeof(int));
				unfold_position(fp3,img);
				get_mi_vector(AA, p1->r.p, fp3);
				get_mi_vector(BB, p2->r.p, fp3);
				get_mi_vector(CC, p4->r.p, fp3);
				for (i=0; i<3; i++) {
                    fp1[i] = fp3[i] + AA[i];
                    fp2[i] = fp3[i] + BB[i];
                    fp4[i] = fp3[i] + CC[i];
                }
			} else {
				// in case the first and the second particle are ghost particles
				if (p4->l.ghost != 1) {
					memmove(fp4, p4->r.p, 3*sizeof(double));
					memmove(img, p4->l.i, 3*sizeof(int));
					unfold_position(fp4,img);
					get_mi_vector(AA, p1->r.p, fp4);
					get_mi_vector(BB, p2->r.p, fp4);
					get_mi_vector(CC, p3->r.p, fp4);
					for (i=0; i<3; i++) {
                        fp1[i] = fp4[i] + AA[i];
                        fp2[i] = fp4[i] + BB[i];
                        fp3[i] = fp4[i] + CC[i];
                    }
				} else {
					printf("Something wrong in oif_local_forces.hpp: All particles in a bond are ghost particles, impossible to unfold the positions...");
					return 0;
				}
			}
		}
	}
	
    for(i=0; i<3; i++) {
        force[i] = 0;
        force2[i] = 0;
        force3[i] = 0;
        force4[i] = 0;
    }
    
    // non-linear stretching
    if (m_ks > TINY_OIF_ELASTICITY_COEFFICIENT) {
        vecsub(fp2,fp3,dx);
        len2 = sqrlen(dx);
        len = sqrt(len2);
        dr = len - m_r0;
        lambda = 1.0*len/m_r0;
        fac = -m_ks * KS(lambda) * dr; // no normalization
        for(i=0; i<3; i++) {
            force2[i] += fac*dx[i]/len;
            force3[i] += -fac*dx[i]/len;
        }
    }
    
    // linear stretching
    if (m_kslin > TINY_OIF_ELASTICITY_COEFFICIENT) {
        vecsub(fp2,fp3,dx);
        len2 = sqrlen(dx);
        len = sqrt(len2);
        dr = len - m_r0;
        fac = -m_kslin * dr; // no normalization
        for(i=0; i<3; i++) {
            force2[i] += fac*dx[i]/len;
            force3[i] += -fac*dx[i]/len;
        }
    }
    
    /* bending
       forceT1 is restoring force for triangle p1,p2,p3 and force2T restoring force for triangle p2,p3,p4
       p1 += forceT1; p2 -= 0.5*forceT1+0.5*forceT2; p3 -= 0.5*forceT1+0.5*forceT2; p4 += forceT2; */
    if (m_kb > TINY_OIF_ELASTICITY_COEFFICIENT) {
        if (m_phi0 < 0.001 || m_phi0 > 2*M_PI - 0.001)
            printf("oif_local_forces.hpp, calc_oif_local: Resting angle is close to zero!\n");
        get_n_triangle(fp2,fp1,fp3,n1);
        dn1=normr(n1);
        get_n_triangle(fp2,fp3,fp4,n2);
        dn2=normr(n2);
        phi = angle_btw_triangles(fp1,fp2,fp3,fp4);
        if (phi < 0.001 || phi > 2*M_PI - 0.001) printf("oif_local_forces.hpp, calc_oif_local: Angle approaches 0 or 2*Pi\n");
        aa = (phi - m_phi0)/m_phi0;
        fac = m_kb * aa;
        for(i=0; i<3; i++) {
            force[i] += fac * n1[i]/dn1;
            force2[i] -= (0.5 * fac * n1[i]/dn1 + 0.5 * fac * n2[i]/dn2);
            force3[i] -= (0.5 * fac * n1[i]/dn1 + 0.5 * fac * n2[i]/dn2);
            force4[i] += fac * n2[i]/dn2;
        }
    }

    /* local area
       for both triangles
       only 1/3 of calculated forces are added, because each triangle will enter this calculation 3 times (one time per edge)
     test
       normalisation by square root of the triangle's area, see "R.Tothova: Comparison of Different Formulas for Local Area Conservation Modulus in Spring Network Models"
    */
    if (m_kal > TINY_OIF_ELASTICITY_COEFFICIENT) {
        
        // triangle p1,p2,p3
        for(i=0; i<3; i++){            // centroid of triangle p1,p2,p3
            h[i]=1.0/3.0 *(fp1[i]+fp2[i]+fp3[i]);
        }
        A=area_triangle(fp1,fp2,fp3);
        aa=(A - m_A01)/sqrt(A);
        vecsub(h,fp1,rh);
        hn=normr(rh);
        for(i=0; i<3; i++) {          // local area force for p1
            force[i] += m_kal * aa * rh[i]/(3*hn);
        }
        vecsub(h,fp2,rh);
        hn=normr(rh);
        for(i=0; i<3; i++) {          // local area force for p2
            force2[i] += m_kal * aa * rh[i]/(3*hn);
        }
        vecsub(h,fp3,rh);
        hn=normr(rh);
        for(i=0; i<3; i++) {          // local area force for p3
            force3[i] += m_kal * aa * rh[i]/(3*hn);
        }
        
        // triangle p2,p3,p4
        for(i=0; i<3; i++) {         // centroid of triangle p2,p3,p4
            h[i]=1.0/3.0 *(fp2[i]+fp3[i]+fp4[i]);
        }
        A=area_triangle(fp2,fp3,fp4);
        aa=(A - m_A02)/sqrt(A);
        vecsub(h,fp4,rh);
        hn=normr(rh);
        for(i=0; i<3; i++) {          // local area force for p4
            force4[i] += m_kal * aa * rh[i]/(3*hn);
        }
        vecsub(h,fp2,rh);
        hn=normr(rh);
        for(i=0; i<3; i++) {          // local area force for p2
            force2[i] += m_kal * aa * rh[i]/(3*hn);
        }
        vecsub(h,fp3,rh);
        hn=normr(rh);
        for(i=0; i<3; i++) {          // local area force for p3
            force3[i] += m_kal * aa * rh[i]/(3*hn);
        }
    }
  return 0;

}
//energy calculation
int Bond::OifLocalForces::calc_bonded_four_particle_energy(Particle *p1, Particle *p2, Particle *p3, 
							  Particle *p4, double *_energy) const {
  *_energy = 0.0;
  return 0;

}

void Bond::OifLocalForces::write_force_to_particle(Particle *p1, Particle *p2, Particle *p3, 
			     Particle *p4, double force[3], double force2[3],
			     double force3[3], double force4[3]) const
{
  for (int j = 0; j < 3; j++) {
    p1->f.f[j] += force2[j];
    p2->f.f[j] += force[j];
    p3->f.f[j] += force3[j];
    p4->f.f[j] += force4[j];
  }
}

boost::any Bond::OifLocalForces::get_bond_parameters_from_bond() const
{

  Oif_local_forces_bond_parameters params = {m_r0, m_ks, m_kslin, m_phi0, m_kb, m_A01, m_A02, m_kal};
  return boost::any(params);
  
}
