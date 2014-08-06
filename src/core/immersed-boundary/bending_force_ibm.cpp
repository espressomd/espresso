#include "immersed-boundary/bending_force_ibm.hpp"

#ifdef BENDING_FORCE_IMMERSED_BOUNDARY
#include "communication.hpp"


int tribend_set_params(int bond_type, int ind1, int ind2, int ind3, int ind4, int boo, double max, double kb) {
    Particle p1, p2, p3, p4;
    double n1l[3], n2l[3], n1[3], n2[3];
    double dx1[3], dx2[3], dx3[3];
    double theta0, desc;
    double tmp[3];
    double sc;
    
    if(bond_type<0) {
      return ES_ERROR;
    }
    
    make_bond_type_exist(bond_type);
    
    //printf("Index: %d %d %d %d\n", ind1, ind2, ind3, ind4);
    
    get_particle_data(ind1, &p1);
    get_particle_data(ind2, &p2);
    get_particle_data(ind3, &p3);
    get_particle_data(ind4, &p4);
    
    //printf("p1:%d %lf %lf %lf\n", p1.p.identity, p1.r.p[0], p1.r.p[1], p1.r.p[2]);
    //printf("p2:%d %lf %lf %lf\n", p2.p.identity, p2.r.p[0], p2.r.p[1], p2.r.p[2]);
    //printf("p3:%d %lf %lf %lf\n", p3.p.identity, p3.r.p[0], p3.r.p[1], p3.r.p[2]);
    //printf("p4:%d %lf %lf %lf\n", p4.p.identity, p4.r.p[0], p4.r.p[1], p4.r.p[2]);
    
    //Get vectors of triangles
    get_mi_vector(dx1, p1.r.p, p3.r.p);
    get_mi_vector(dx2, p2.r.p, p3.r.p);
    get_mi_vector(dx3, p4.r.p, p3.r.p);
    
    //printf("dx1 = %lf %lf %lf\n", dx1[0], dx1[1], dx1[2]);
    //printf("dx2 = %lf %lf %lf\n", dx2[0], dx2[1], dx2[2]);
    //printf("dx3 = %lf %lf %lf\n", dx3[0], dx3[1], dx3[2]);
    
    //Get normals on triangle; pointing outwards by definition of indices sequence
    vector_product(dx1, dx2, n1l);
    vector_product(dx1, dx3, n2l);
    
    //printf("n1l: %lf %lf %lf\n", n1l[0], n1l[1], n1l[2]);
   
    
    
    if(boo == 0) {
	n2l[0]=-1*n2l[0]; n2l[1]=-1*n2l[1]; n2l[2]=-1*n2l[2];
    }
    
    //printf("n2l: %lf %lf %lf\n", n2l[0], n2l[1], n2l[2]);
    
    unit_vector(n1l,n1);
    unit_vector(n2l,n2);
    
    //printf("n1: %lf %lf %lf\n", n1[0], n1[1], n1[2]);
    //printf("n2: %lf %lf %lf\n", n2[0], n2[1], n2[2]);
    
    //calculate theta by taking the acos of the scalar n1*n2
    sc = scalar(n1,n2);
  
    if(sc>1.0) {
      sc = 1.0;
    }
    
    theta0 = acos(sc);
    vector_product(n1,n2,tmp);
    
    desc = scalar(dx1,tmp);
    
    if(desc<0) {
	theta0 = 2.0*PI-theta0;
    }
    
    //printf("%lf\n", theta0*TOANGLE);
    
    //effective springconstant = bending resistance^(1/3); Krueger2012
    bonded_ia_params[bond_type].p.tribend.kb = sqrt(3)*kb;
    bonded_ia_params[bond_type].p.tribend.theta0 = theta0;
    bonded_ia_params[bond_type].p.tribend.boo = boo;
    bonded_ia_params[bond_type].p.tribend.max = max;
    
    bonded_ia_params[bond_type].type = TRIBEND_IA;
    bonded_ia_params[bond_type].num = 3;
    
    mpi_bcast_ia_params(bond_type, -1);
    
    //free allocated paticles
    free_particle(&p1);
    free_particle(&p2);
    free_particle(&p3);
    free_particle(&p4);
    
    return ES_OK;
    
}

int tribend_reset_params(int bond_type, double bood, double theta0, double kb, double max) {
    
    int boo = bood;
  
    if(bond_type<0) {
      return ES_ERROR;
    }
    
    make_bond_type_exist(bond_type);
    
    bonded_ia_params[bond_type].p.tribend.boo = boo;
    bonded_ia_params[bond_type].p.tribend.theta0 = theta0;
    bonded_ia_params[bond_type].p.tribend.kb = kb;
    bonded_ia_params[bond_type].p.tribend.max = max;
    
    bonded_ia_params[bond_type].type = TRIBEND_IA;
    bonded_ia_params[bond_type].num = 3;
    
    mpi_bcast_ia_params(bond_type, -1);
    
    return ES_OK;
}

#endif
