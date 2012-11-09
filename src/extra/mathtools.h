#ifndef included_mathtools_h
#define included_mathtools_h

/* type conversion */
void double2double(int n, double* d, double* f);
void double2double(int n, double* f, double* d);

/* integer stuff */
int maximum(int a, int b);

/* one vectors */
void set0(int n, double* a);
void set1(int n, double* a);
void setvector3(double a, double b, double c, double* x);
double normr(int n, double* r);
void normv(int n, double* r, double* rdach);

/* vector scaling */
void btimess(int n, double s, double* b);
void bdivs(int n, double s, double* b);
void btimessx(int n, double s, double* b, double* x);
void bdivsx(int n, double s, double* b, double* x);

/* two vectors */
void copyv(int n, double* a, double* b);
void bdivsa(int n, double s, double* a, double* b); /* b_i = b_i /(s*a_i) */
double dot(int n, double* a, double* b);
void aminusb(int n, double* a, double* b, double* c);
double angleab(double* a, double* b);
void rotatecfromatob(double* m, double* h, double* c, double phi);
void cisaxb(double* a, double* b, double* c);
void cisaxaxb(double* a, double* b, double* c);
double distv(int, double* a, double* b);

/* geometry */
double volume_sphere(double r);
int meshsizesphere(double r, double dx);
int meshsphere(int n0, double r, double* center, double dx, double* x, int lx);
double area_triangle_new(double *P1, double *P2);
double area_triangle(double *P1, double *P2, double *P3);
double angle_btw_triangles_new(double *P1, double *P2, double *P3, double *P4);
double angle_btw_triangles(double *P1, double *P2, double *P3, double *P4);
void get_n_triangle(double* p1, double* p2, double* p3, double* n);  // calculate normal vector of triangle

/* magnetics */
/* physics: SI units: Magnetization Ms = Js/m0 
 *                    now use m*m0 as unit for dipole 
 * input magnetic polarization times volume */
void dipolefield(double* m, double* x, double *field);


double interpolate(int n, double* x, double *y, double x0);

#endif 
