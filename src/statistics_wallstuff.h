#ifndef WALLSTUFF_H
#define WALLSTUFF_H

/// list of the currently specified box boundaries for layerwise analysis
extern DoubleList wallstuff_boundaries;
/** particles in the bins defined by \ref wallstuff_boundaries. There
    is one less bin that boundaries. */
extern IntList *wallstuff_part_in_bin;

/// sort particles into bins along the x-coordinate
void wall_sort_particles();
/// yz-MSD in a x-bin - use \ref wall_sort_particles before
void calc_wallmsdyz(double *g, int bin);
/// x-MSD in a x-bin - use \ref wall_sort_particles before
void calc_wallmsdx(double *g, int bin);
/** bond order parameter G_6 in xy, in a x-bin -
    use \ref wall_sort_particles before */
void calc_wallbondyz(double *g, int bin, double rclocal,
		     double rmax, int rbins);
/** radial distribution function in xy, in a x-bin -
    use \ref wall_sort_particles before */
void calc_wallrdfyz(double *g, int bin, double rmin, double rmax, int rbins);
/** scaling of the bond order parameter in xy, in a x-bin - use \ref
    wall_sort_particles before */
void calc_scaling(double *g, int bin, int boxes, double rclocal);
/** scaling of the bond order parameter in xy, in a x-bin - use \ref
    wall_sort_particles before. This returns the raw values, to be
    processed later, e.g. to calculate the Binder cumulant. */
void calc_scaling2(double *g, int bin, int boxes, double rclocal);

#endif
