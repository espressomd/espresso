#include"MembraneCollision.hpp"
#include"grid.hpp" //get_mi_vector unfold_position
#include"utils.hpp" //get_n_triangle

int Bond::MembraneCollision::add_bonded_four_particle_force(Particle *p1, Particle *p2, 
							    Particle *p3, Particle *p4, 
							    double force[3], double force2[3], 
							    double force3[3], double force4[3]) 
  const{

    double n[3],dn;
	int j;
	double fp1[3],fp2[3],fp3[3],fp4[3];
	int img[3];
    double AA[3],BB[3],CC[3];

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
        for (int i=0; i < 3; i++) { fp1[i] = fp2[i] + AA[i]; fp3[i] = fp2[i] + BB[i]; fp4[i] = fp2[i] + CC[i]; }
    } else {
        // in case  particle p2 is a ghost particle
        if (p1->l.ghost != 1) {
            memmove(fp1, p1->r.p, 3*sizeof(double));
            memmove(img, p1->l.i, 3*sizeof(int));
            unfold_position(fp1,img);
            get_mi_vector(AA, p2->r.p, fp1);
            get_mi_vector(BB, p3->r.p, fp1);
            get_mi_vector(CC, p4->r.p, fp1);
            for (int i=0; i < 3; i++) { fp2[i] = fp1[i] + AA[i]; fp3[i] = fp1[i] + BB[i];  fp4[i] = fp1[i] + CC[i];}
        } else {
            // in case the first and the second particle are ghost particles
            if (p3->l.ghost != 1) {
                memmove(fp3, p3->r.p, 3*sizeof(double));
                memmove(img, p3->l.i, 3*sizeof(int));
                unfold_position(fp3,img);
                get_mi_vector(AA, p1->r.p, fp3);
                get_mi_vector(BB, p2->r.p, fp3);
                get_mi_vector(CC, p4->r.p, fp3);
                for (int i=0; i < 3; i++) { fp1[i] = fp3[i] + AA[i]; fp2[i] = fp3[i] + BB[i]; fp4[i] = fp3[i] + CC[i]; }
            } else {
                // in case the first and the second particle are ghost particles
                if (p4->l.ghost != 1) {
                    memmove(fp4, p4->r.p, 3*sizeof(double));
                    memmove(img, p4->l.i, 3*sizeof(int));
                    unfold_position(fp4,img);
                    get_mi_vector(AA, p1->r.p, fp4);
                    get_mi_vector(BB, p2->r.p, fp4);
                    get_mi_vector(CC, p3->r.p, fp4);
                    for (int i=0; i < 3; i++) { fp1[i] = fp4[i] + AA[i]; fp2[i] = fp4[i] + BB[i]; fp3[i] = fp4[i] + CC[i]; }
                } else {
                    printf("Something wrong in bond/MembraneCollision.cpp: All particles in a bond are ghost particles, impossible to unfold the positions...");
                    return 0;
                }
            }
        }
    }
	
	get_n_triangle(fp2,fp3,fp4,n);
	dn=normr(n);
    if ( fabs(dn) < 0.001 )
		printf("MembraneCollision.cpp, add_bonded_four_particle_force: Length of outward vector is close to zero!\n");
	for(j=0;j<3;j++){
            p1->p.out_direction[j] = n[j]/dn;
    }
    return 0;

}

int Bond::MembraneCollision::add_bonded_four_particle_energy(Particle *p1, Particle *p2, 
							     Particle *p3,Particle *p4, 
							     double *_energy) const {

  *_energy = 0.0;
  return 0;
}

