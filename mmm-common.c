// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
/** \file mmm-common.c
    common parts of the MMM family of methods for the electrostatic interaction, MMM1D, MMM2D and ELC.
    This file contains the code for the polygamma expansions used for the near formulas of MMM1D and MMM2D,
    and the documentation for the full family.

    The expansion of the polygamma functions is fairly easy and follows directly from Abramowitz and Stegun.
    For details, see Axel Arnold and Christian Holm, "MMM2D: A fast and accurate summation method for
    electrostatic interactions in 2D slab geometries", Comp. Phys. Comm., 148/3(2002),327-348.
*/
/** \page MMM_general The MMM family of algorithms
\section intro Introduction
In the MMM family of algorithms for the electrostatic interaction, a convergence factor approach to tackle the conditionally convergent
Coulomb sum is used (even the authors of the original MMM method have no idea what this acronym stands for). Instead of defining the
summation order, one multiplies each summand by a continuous factor \f$c(\beta,r_{ij},n_{klm})\f$ such that the sum is absolutely convergent
for \f$\beta>0\f$, but \f$c(0,.,.)=1\f$. The energy is then defined as the limit \f$\beta\rightarrow 0\f$ of the sum, i. e. \f$\beta\f$ is an
artificial convergence parameter.  For a convergence factor of \f$e^{-\beta n_{klm}^2}\f$ the limit is the same as the spherical limit,
and one can derive the classical Ewald method quite conveniently through this approach (Smith 81). To derive the formulas for MMM,
one has to use a different convergence factor, namely \f$e^{-\beta|r_{ij}+n_{klm}|}\f$, which defines the alternative energy
\f[
  \tilde{E}=\,\frac{1}{2}\lim_{\beta\rightarrow 0}\sum_{k,l,m}{\sum_{i,j=1}^N}' \frac{q_i
    q_je^{-\beta|p_{ij} + n_{klm}|}}
  {|p_{ij} + n_{klm}|}
  =:\,\frac{1}{2}\lim_{\beta\rightarrow 0}\sum_{i,j=1}^N q_iq_j\phi_\beta(x_{ij},
  y_{ij},z_{ij}).
\f]
\f$\phi_\beta\f$ is given by \f$ \phi_\beta(x,y,z)=\,\tilde\phi_\beta(x,y,z) + \frac{e^{-\beta r}}{r} \f$
for \f$(x,y,z)\neq 0\f$ and \f$\phi_\beta(0,0,0)=\,\tilde\phi_\beta(0,0,0)\f$, where
\f[
\tilde\phi_\beta(x,y,z)=\,\sum_{(k,l,m)\neq 0} \frac{e^{-\beta r_{klm}}}{r_{klm}}.
\f]

The limit \f$\tilde{E}\f$ exists, but differs for three dimensionally periodic systems by some multiple of the square of the
dipole moment from the spherical limit as obtained by the Ewald summation (Smith 81). From the physical point of view the Coulomb
interaction is replaced by a screened Coulomb interaction with screening length \f$1/\beta\f$.  \f$\tilde{E}\f$ is then the
energy in the limit of infinite screening length. But because of the conditional convergence of the electrostatic sum, this is
not necessarily the same as the energy of an unscreened system. Since the difference to the Ewald methods only depends on the dipole
moment of the system, the correction can be calculated easily in linear time and can be ignored with respect to accuracy as well as
to computation time.

For one or two dimensionally systems, however, \f$\tilde{E}=E\f$, i.e. the convergence factor approach equals the spherical summation
limit of the Ewald sum, and MMM1D and MMM2D do not require a dipole correction.

Starting from this convergence factor approach, R. Strebel and R. Sperb constructed a method of computational order \f$O(N\log N)\f$,
MMM (Strebel 99). The favourable scaling is obtained, very much like in the Ewald case, by technical tricks in the
calculation of the far formula. The far formula has a product decomposition and can be evaluated hierarchically similarly to the
fast multipole methods.

For particles sufficiently separated in the z-axis one can Fourier transform the potential along both x and y. We obtain
the far formula as
\f[
    \phi(x,y,z) =\,
      u_x u_y\sum_{p,q\neq 0}
    \frac{e^{2\pi f_{pq}z} + e^{2\pi f_{pq}(\lambda_z-z)}}{f_{pq}
      \left(e^{2\pi f_{pq}\lambda_z} - 1\right)}
    e^{2\pi i u_y q y}e^{2\pi i u_x p x} +
      2\pi u_x u_y\left(u_z z^2 - z + \frac{\lambda_z}{6}\right).
\f]
where \f$\lambda_{x,y,z}\f$ are the box dimensions,
\f$ f_{pq} =\, \sqrt{(u_x p)^2 + (u_y q)^2},\quad f_p =\, u_x p,\quad f_q =\, u_x q \f$,
\f$ \omega_p=2\pi u_x p\f$ and \f$\omega_q=2\pi u_y q\f$.
The advantage of this formula is that it allows for a product decomposition into components of the particles.  For example
\f[
  e^{2\pi f_{pq}z}=e^{2\pi f_{pq}(z_i-z_j)}=e^{2\pi f_{pq}z_i}e^{-2\pi f_{pq}z_j}
\f]
etc. Therefore one just has to calculate the sum over all these exponentials on the left side and on the right side and multiply
them together, which can be done in \f$O(N)\f$ computation time. As can be seen easily, the convergence of the series is excellent as
long as z is sufficiently large. By symmetry one can choose the coordinate with the largest distance as z to optimise the
convergence. Similar to the Lekner sum, we need a different formula if all coordinates are small, i. e. for particles close
to each other. For sufficiently small \f$u_y\rho\f$ and \f$u_xx\f$ we obtain the near formula as
\f[
  \begin{array}{rl}
    \tilde\phi(x,y,z)=\,
    &   2 u_x u_y\sum\limits_{p,q>0}
    \frac{\cosh(2\pi f_{pq}z)}{f_{pq}
      \left(e^{2\pi f_{pq}\lambda_z} - 1\right)}
    e^{2\pi i u_y q y}e^{2\pi i u_x p x} +\\
    & 4u_x\sum\limits_{l,p>0}\left(K_0(2\pi u_x p\rho_l) +
      K_N(2\pi u_x p\rho_{-l})\right)cos(2\pi u_x p x) -\\
    &  2u_x\sum\limits_{n\ge 1}\frac{b_{2n}}{2n(2n)!}\Re\bigl((2\pi u_y (z+iy))^{2n}\bigr) +\\
    &  u_x\sum\limits_{n\ge 0}\left(\begin{array}{c}-\frac{1}{2}\\ n\end{array}\right)\frac{\left(
        \psi^{(2n)}(1 + u_x x) +
        \psi^{(2n)}(1 - u_x x)\right)}{(2n)!}\rho^{2n} -\\
    &  2\log(4\pi).
  \end{array}
\f]
Note that this time we calculate \f$\tilde{\phi}\f$ instead of \f$\phi\f$, i. e. we omit the contribution of the primary simulation box.
This is very convenient as it includes the case of self energy and  makes \f$\tilde{\phi}\f$ a smooth function. To obtain \f$\phi\f$ one
has to add the \f$1/r\f$ contribution of the primary box. The self energy is given by
\f[
  \tilde\phi(0,0,0)=\,
    2 u_x u_y\sum\limits_{p,q>0} \frac{1}{f_{pq}
      \left(e^{2\pi f_{pq}\lambda_z} - 1\right)}+ 8u_x\sum\limits_{l,p>0}K_N(2\pi u_x\lambda_y p l) +
    2 u_x\psi^{(0)}(1) - 2\log(4\pi).
\f]
Both the near and far formula are derived using the same convergence factor approach, and consequently the same singularity in
\f$\beta\f$ is obtained. This is important since otherwise the charge neutrality argument does not hold.

To obtain the \f$O(N\log N)\f$ scaling, some algorithm tricks are needed, which are not used in MMM1D, MMM2D or ELC and are therefore not
discussed here. For details see Strebel,99. MMM is not implemented in espresso.

\section MMM2D MMM2D

In the case of periodicity only in the x and y directions, the far formula looks like
\f[
  \phi(x,y,z) =\,
  4 u_x u_y\sum_{p,q>0}
  \frac{e^{-2\pi f_{pq}|z|}}
  {f_{pq}}
  \cos(\omega_p x)\cos(\omega_q y)\, +
  2 u_x u_y\left(\sum_{q>0}
    \frac{e^{-2\pi f_q|z|}}{f_q}
    \cos(\omega_q y)\, + \sum_{p>0}
    \frac{e^{-2\pi f_p|z|}}{f_p}
    \cos(\omega_p x)\right)\,
  - 2\pi u_x u_y |z|,
\f]
and the near formula is
\f[
  \begin{array}{rl}
    \tilde\phi(x,y,z)=\,
    & 4u_x\sum_{l,p>0}\left(K_0(\omega_p\rho_l) +
      K_0(\omega_p\rho_{-l})\right)\cos(\omega_p x) -\\
    &  2u_x\sum_{n\ge 1}\frac{b_{2n}}{2n(2n)!}
    \Re\bigl((2\pi u_y (z+iy))^{2n}\bigr)\,+\,
    \sum_{k=1}^{N_\psi-1}\left(\frac{1}{r_{k}} +
      \frac{1}{r_{-k}}\right) -\\
    &   u_x\sum_{n\ge 0}\left(\begin{array}{c}-\frac{1}{2}\\n\end{array}\right)\frac{\left(
        \psi^{(2n)}(N_\psi + u_x x) +
        \psi^{(2n)}(N_\psi - u_x x)\right)}{(2n)!}(u_x\rho)^{2n} -\\
    &  2u_x\log\left(4\pi\frac{u_y}{u_x}\right).
  \end{array}
\f]
As said before, the energy obtained from these potentials is equal to the electrostatic energy obtained by the spherical summation
limit. The deeper reason for this is that in some sense the electrostatic sum is absolutely convergent, see Arnold,02 for a discussion.

The near formula is used for particles with a small distance along the z axis, for all other particles the far formula is used. Below
is shown, that the far formula can be evaluated much more efficiently, however, its convergence breaks down for small z distance.
To efficiently implement MMM2D, the layered cell system is required, which splits up the system in equally sized gaps along the z axis.
The interaction of all particles in a layer S with all particles in the layers S-1,S,S+1 is calculated using the near formula, for the
particles in layers 1,...,S-2, and in layers S+2,...,N, the far formula is used.

The implementation of the near formula is relatively straight forward and can be treated as any short ranged force is treated
using the link cell algorithm, here in the layered variant. The special functions in the formula are somewhat demanding, but
for the polygamma functions Taylor series can be achieved, which are implemented in \ref mmm-common.h "mmm-common.h". The Bessel
functions are calculated using a Chebychev series, see specfunc.h "specfunc.h".

The treatment of the far formula is algorithmically more complicated. For a particle i in layer \f$ S_i\f$,The formula can product
decomposed, as in
\f[
  \begin{array}{rl}
    \sum_{j\in I_S, S < S_i - 1}
    q_iq_j\frac{e^{-2\pi f_{pq}|z_i-z_j|}}{f_{pq}}
    \cos(\omega_p (x_i - x_j))\cos(\omega_q (y_i - y_j)) =
    &q_i\frac{e^{-2\pi f_{pq}z_i}}{f_{pq}}
    \cos(\omega_p x_i)\cos(\omega_q y_i)
    \sum_{j\in I_S, S < S_i - 1}q_je^{2\pi f_{pq}z_j}
    \cos(\omega_p x_j)\cos(\omega_q y_j) + \\
    &q_i\frac{e^{-2\pi f_{pq}z_i}}{f_{pq}}
    \cos(\omega_p x_i)\sin(\omega_q y_i)
    \sum_{j\in I_S, S < S_i - 1}q_je^{2\pi f_{pq}z_j}
    \cos(\omega_p x_j)\sin(\omega_q y_j)
    + \\
    &q_i\frac{e^{-2\pi f_{pq}z_i}}{f_{pq}}
    \sin(\omega_p x_i)\cos(\omega_q y_i)
    \sum_{j\in I_S, S < S_i - 1}q_je^{2\pi f_{pq}z_j}
    \sin(\omega_p x_j)\cos(\omega_q y_j)
    + \\
    &q_i\frac{e^{-2\pi f_{pq}z_i}}{f_{pq}}
    \sin(\omega_p x_i)\sin(\omega_q y_i)
    \sum_{j\in I_S, S < S_i - 1}q_je^{2\pi f_{pq}z_j}
    \sin(\omega_p x_j)\sin(\omega_q y_j).
  \end{array}
\f]
This representation has the advantage, that the contributions of the two particles are decoupled. For all particles j only the
eight terms
\f[
  \xi^{(\pm,s/c,s/c)}_j= q_je^{\pm 2\pi f_{pq}z_j} \sin/\cos(\omega_p x_j)\sin/\cos(\omega_q y_j)
\f]
are needed. The upper index describes the sign of the exponential term and whether sine or cosine is used for \f$x_j\f$ and
\f$y_j\f$ in the obvious way.  These terms can be used for all expressions on the right hand side of the product decomposition.
Moreover it is easy to see from the addition theorem for the sine function that these terms also can be used to calculate the force
information up to simple prefactors that depend only on p and q.

Every processor starts with the calculation of the terms \f$\xi^{(\pm,s/c,s/c)}_j\f$ and adds them up in each layer, so that
one obtains
\f[
  \Xi^{(\pm,s/c,s/c)}_s= \sum_{j\in S_s}\xi^{(\pm,s/c,s/c)}_j.
\f]
Now we calculate
\f[
  \Xi^{(l,s/c,s/c)}_s=\sum_{t < s - 1}\Xi^{(+,s/c,s/c)}_t
\f]
and
\f[
  \Xi^{(h,s/c,s/c)}_s=\sum_{t > s + 1}\Xi^{(-,s/c,s/c)}_t,
\f]
which are needed for the evaluation of the product decomposition. While the bottom processor can calculate \f$\Xi^{(l,s/c,s/c)}_s\f$
directly, the other processors are dependent on its results. Therefore the bottom processor starts with the calculation of its
\f$\Xi^{(l,s/c,s/c)}_s\f$ and sends up \f$\Xi^{(l,s/c,s/c)}_s\f$ and \f$\Xi^{(+,s/c,s/c)}_s\f$ of its top layer s to the next
processor dealing with the layers above. Simultaneously the top processor starts with the calculation of the
\f$\Xi^{(h,s/c,s/c)}_s\f$ and sends them down. After the communicated has been completed, every processor can use the
\f$\Xi^{(l/h,s/c,s/c)}_j\f$ and the \f$\xi^{(\pm,s/c,s/c)}_j\f$ to calculate the force rsp. energy contributions for its particles.

In pseudo code, the far formula algorithm looks like this
<ul>
  <li> for each layer s=1,...,S
  <ul>
    <li>\f$\Xi^{(\pm,s/c,s/c)}_s=0\f$</li>
    <li>for each particle j in layer s
    <ul>
      <li>calculate \f$\xi^{(\pm,s/c,s/c)}_j\f$</li>
      <li>\f$\Xi^{(\pm,s/c,s/c)}_s += \xi^{(\pm,s/c,s/c)}_j\f$</li>
    </ul></li>
  </ul></li>
  <li>\f$\Xi^{(l,s/c,s/c)}_3=\Xi^{(+,s/c,s/c)}_1\f$</li>
  <li>for each layer s=4,...,S
  <ul>
    <li>\f$\Xi^{(l,s/c,s/c)}_s=\Xi^{(l,s/c,s/c)}_{s-1} + \Xi^{(+,s/c,s/c)}_{s-2}\f$</li>
  </ul></li>
  <li>\f$\Xi^{(l,s/c,s/c)}_{S-2}=\Xi^{(-,s/c,s/c)}_S\f$</li>
  <li>for each layer s=(S-3),...,1
  <ul>
    <li>\f$\Xi^{(l,s/c,s/c)}_s=\Xi^{(l,s/c,s/c)}_{s+1} + \Xi^{(-,s/c,s/c)}_{s+2}\f$</li>
  </ul></li>
  <li>for each layer s=1,...,S
  <ul>
    <li>for each particle j in layer s
    <ul>
    <li>calculate particle interaction from \f$\xi^{(+,s/c,s/c)}_j\Xi^{(l,s/c,s/c)}_s\f$ and
        \f$\xi^{(-,s/c,s/c)}_j\Xi^{(h,s/c,s/c)}_s\f$</li>
    </ul></li>
  </ul></li>
</ul>

For further details, see the articles of Arnold and Holm, 2002.

\section MMM1D MMM1D

In one dimensionally periodic systems with z being the periodic coordinate, the far formula looks like
\f[
  \begin{array}{rl}
    \phi(\rho,z) &=\, 4 u_z\sum_{p\neq 0}
    K_0(\omega\rho)\cos(\omega z)
    - 2u_z\log(\frac{\rho}{2\lambda_z}) - 2u_z\gamma\\
    F_\rho(\rho,z) &=\, 8\pi u_z^2\sum_{p\neq 0}
    p K_1(\omega\rho)\cos(\omega z) + \frac{2 u_z}{\rho}\\
    F_z(\rho,z) &=\, 8\pi u_z^2 \sum_{p\neq 0}
    pK_0(\omega\rho)\sin(\omega z),
  \end{array}
\f]
the near formula is
\f[
  \begin{array}{rl}
    \tilde{\phi}(\rho,z) &=\, -u_z\sum_{n\ge 0} \left(\begin{array}{c}-\frac{1}{2}\\n\end{array}\right)
      \frac{\left(\psi^{(2n)}(N_\psi + u_z z) +
          \psi^{(2n)}(N_\psi - u_z z)\right)}{(2n)!}(u_z\rho)^{2n} - 2u_z\gamma + \\
    &\phantom{=\,++}
    \sum_{k=1}^{N_\psi-1}\left(\frac{1}{r_k}+\frac{1}{r_{-k}}\right)\\
    \tilde{F}_\rho(\rho,z) &=\, -u_z^3 \sum_{n\ge 0} \left(\begin{array}{c}-\frac{1}{2}\\n\end{array}\right)
    \frac{\left(\psi^{(2n)}(N_\psi + u_z z) +
        \psi^{(2n)}(N_\psi - u_z z)\right)}{(2n)!}(u_z\rho)^{2n-1} + \\
    &\phantom{=\,++}
    \sum_{k=1}^{N_\psi-1}\left(\frac{\rho}{r_k^3}+\frac{\rho}{r_{-k}^3}\right) \\
    \tilde{F}_z(\rho,z) &=\, -u_z^2 \sum_{n\ge 0} \left(\begin{array}{c}-\frac{1}{2}\\n\end{array}\right)
    \frac{\left(\psi^{(2n + 1)}(N_\psi + u_z z) +
        \psi^{(2n + 1)}(N_\psi - u_z z)\right)}{(2n)!}(u_z\rho)^{2n} + \\
    &\phantom{=\,++}
    \sum_{k=1}^{N_\psi-1}\left(\frac{z+k\lambda_z}{r_k^3}+\frac{z-k\lambda_z}{r_{-k}^3}\right),
  \end{array}
\f]
where \f$\rho\f$ denotes the xy-distance of the particles. As for the two dimensional periodic case,
the obtained energy is equal to the one dimensional Ewald sum. Algorithmically, MMM1D is uninteresting,
since neither the near nor far formula allow a product decomposition or similar tricks. MMM1D has to be
implemented as a simple NxN loop. However, the formulas can be evaluated efficiently, so that MMM1D can
still be used reasonably for up to 400 particles on a single processor.

\section ELC ELC

The ELC method differs from the other MMM algorithms in that it is not an algorithm for the calculation of the
electrostatic interaction, but rather represents a correction term which allows to use any method for threedimensionally
periodic systems with spherical summation order for twodimensional periodicity. The basic idea is to expand the two
dimensional slab system of height h in the non-periodic z-coordinate to a system with periodicity in all three dimensions,
with a period of \f$\lambda_z>h\f$, which leaves an empty gap of height \f$\delta=\lambda_z - h\f$ above the particles
in the simulation box.

Since the electrostatic potential is only finite if the total system is charge neutral, the additional image layers (those layers
above or below the original slab system) are charge neutral, too.  Now let us consider the n-th image layer which has an
offset of \f$n\lambda_z\f$ to the original layer. If \f$n\lambda_z\f$ is large enough, each particle of charge q_j at position
\f$(x_j,y_j,z_j+n\lambda_z)\f$ and its replicas in the xy-plane can be viewed as constituting a homogeneous charged sheet of charge
density \f$\sigma_j = \frac{q_j}{\lambda_x\lambda_y}\f$.  The potential of such a charged sheet at distance z is \f$2\pi \sigma_j |z|\f$.
Now we consider the contribution from a pair of image layers located at \f$\pm n\lambda_z\f$, n>0 to the energy of a charge q_i at
position \f$(x_i,y_i,z_i)\f$ in the central layer.  Since \f$|z_j - z_i| < n\lambda_z\f$, we have
\f$|z_j - z_i + n\lambda_z| = n\lambda_z + z_j - z_i\f$ and \f$|z_j - z_i - n\lambda_z|= n\lambda_z - z_j + z_i\f$, and hence the
interaction energy from those two image layers with the charge \f$q_i\f$ vanishes by charge neutrality:
\f[
  2\pi q_i \sum_{j=1}^N \sigma_j(|z_j - z_i + n\lambda_z| + |z_j - z_i - n\lambda_z|)
  = 4\pi q_i n\lambda_z \sum_{j=1}^N \sigma_j = 0.
\f]
The only errors occurring are those coming from the approximation of assuming homogeneously charged, infinite sheets instead of
discrete charges.  This assumption should become better when increasing the distance \f$n\lambda_z\f$ from the central layer.

However, in a naive implementation, even large gap sizes will result in large errors. This is due to the order of
summation for the standard Ewald sum, which is spherical, while the above approach orders the cells in layers, called
slab--wise summation. Smith has shown that by adding to the Ewald energy the term
\f[
  E_c=2\pi M_z^2 - \frac{2\pi M^2}{3},
\f]
where M is the total dipole moment, one obtains the result of a slab--wise summation instead of the spherical limit (Smith 81).
Although this is a major change in the summation order, the difference is a very simple term.  In fact, Smith shows that changes
of the summation order always result in a difference that depends only on the total dipole moment.

Using the far formula of MMM2D, one can calculate the contributions of the additional layers up to arbitrarily precision, even for
small gap sizes. This method is called electrostatic layer correction, ELC. The advantage of this approach is that for the image
layers, z is necessarily large enough, so that all interactions can be represented using the product decomposition. This allows for
an order N evaluation of the ELC term.

The electrostatic layer correction term is given by
\f[
  E_{lc}=\sum_{i,j=1}^Nq_iq_j\psi(p_i-p_j),
\f]
where
\f[
  \begin{array}{rl}
    \psi(x,y,z)=&4u_xu_y\sum_{p,q>0}\frac{\cosh(2\pi f_{pq}z)}{f_{pq}(e^{2\pi f_{pq}\lambda_z} - 1)}
    \cos(\omega_p x)\cos(\omega_q y) + \\
    &2u_xu_y\sum_{p>0}\frac{\cosh(2\pi f_p z)}{f_p(e^{2\pi f_p\lambda_z} - 1)}\cos(\omega_p x)+\\
    &2u_xu_y\sum_{q>0}\frac{\cosh(2\pi f_q z)}{f_q(e^{2\pi f_q\lambda_z} - 1)}\cos(\omega_q y).
  \end{array}
\f]
The implementation is very similar to MMM2d, except that the separation between slices closeby, and above and below is not necessary.

\section Errors

Common to all algorithms of the MMM family is that accuracy is cheap with respect to computation time. More precisely, the maximal
pairwise error, i.e. the maximal error of the \f$\psi\f$ expression, decreases exponentially with the cutoffs. In turn, the computation
time grows logarithmically with the accuracy. This is quite in contrast to the Ewald methods, for which decreasing the error bound can
lead to excessive computation time. For example, P3M cannot reach precisions above \f$10^{-5}\f$ in general. The precise form of the
error estimates is of little importance here, for details see again Arnold,02.

One important aspect is that the error estimates are also exponential in the non-periodic coordinate. Since the number of closeby and
far away particles is different for particles near the border and in the center of the system, the error distribution is highly
non--homogenous. This is unproblematic as long as the maximal error is really much smaller than the thermal energy. However, one cannot
interprete the error simply as an additional error source. Image <img src="../figs/elc_errordist.gif"> shows the error
distribution of the ELC method for a gap size of 10% of the total system height. For MMM2D and MMM1D the error distribution is less
homogenous, however, also here it is always better to have some extra precision, especially since it is computationally cheap.

\section References
<ul>
<li>E. R. Smith, "Electrostatic energy in ionic crystals",
Proc. R. Soc. Lond. A, 375 (1981), 475-505
</li><br>
<li>R. Strebel, "Pieces of software for the Coulombic m body problem",
Dissertation, ETH Zürich, 13504 (1999)
<a href="http://e-collection.ethbib.ethz.ch/show?type=diss\&nr=13504">
http://e-collection.ethbib.ethz.ch/show?type=diss\&nr=13504</a>
</li><br>
<li>A. Arnold and C. Holm, "MMM2D: A fast and accurate summation method for electrostatic interactions in 2D slab geometries",
Comp. Phys. Comm., 148/3 (2002), 327-348
</li><br>
<li>A. Arnold and C. Holm, "A novel method for calculating electrostatic interactions in 2D periodic slab geometries",
Chem. Phys. Lett.,  354 (2002), 324-330
</li><br>
<li>A. Arnold, J. de Joannis and C. Holm, "Electrostatics in Periodic Slab Geometries I",
J. Chem. Phys., 117 (2002), 2496-2502
</li><br>
<li>J. de Joannis and A. Arnold and C. Holm, "Electrostatics in Periodic Slab Geometries II",
J. Chem. Phys., 117 (2002), 2503-2512
</li>
</ul>
*/
#include "mmm-common.h"
#include "utils.h"

Polynom *modPsi = NULL;
int      n_modPsi = 0;

static void preparePolygammaEven(int n, double binom, Polynom *series)
{
  /* (-0.5 n) psi^2n/2n! (-0.5 n) and psi^(2n+1)/(2n)! series expansions
     note that BOTH carry 2n! */
  int order;
  double deriv;
  double maxx, x_order, coeff, pref;

  deriv = 2*n;
  if (n == 0) {
    // psi^0 has a slightly different series expansion
    maxx = 0.25;
    alloc_doublelist(series, 1);
    series->e[0] = 2*(1 - C_GAMMA);
    for (order = 1;; order += 1) {
      x_order = 2*order;
      coeff = -2*hzeta(x_order + 1, 2);
      if (fabs(maxx*coeff)*(4.0/3.0) < ROUND_ERROR_PREC)
	break;
      realloc_doublelist(series, order + 1);
      series->e[order] = coeff;
      maxx *= 0.25;
    }
    series->n = order;
  }
  else {
    // even, n > 0
    maxx = 1;
    pref = 2;
    init_doublelist(series);
    for (order = 0;; order++) {
      // only even exponents of x
      x_order = 2*order;
      coeff = pref*hzeta(1 + deriv + x_order, 2);
      if ((fabs(maxx*coeff)*(4.0/3.0) < ROUND_ERROR_PREC) && (x_order > deriv))
	break;
      realloc_doublelist(series, order + 1);
      series->e[order] = -binom*coeff;
      maxx *= 0.25;
      pref *= (1.0 + deriv/(x_order + 1));
      pref *= (1.0 + deriv/(x_order + 2));
    }
    series->n = order;
  }
}

static void preparePolygammaOdd(int n, double binom, Polynom *series)
{
  int order;
  double deriv;
  double maxx, x_order, coeff, pref;

  deriv  = 2*n + 1;
  maxx = 0.5;
  // to get 1/(2n)! instead of 1/(2n+1)!
  pref = 2*deriv*(1 + deriv);
  init_doublelist(series);
  for (order = 0;; order++) {
    // only odd exponents of x
    x_order = 2*order + 1;
    coeff = pref*hzeta(1 + deriv + x_order, 2);
    if ((fabs(maxx*coeff)*(4.0/3.0) < ROUND_ERROR_PREC) && (x_order > deriv))
      break;
    realloc_doublelist(series, order + 1);
    series->e[order] = -binom*coeff;
    maxx *= 0.25;
    pref *= (1.0 + deriv/(x_order + 1));
    pref *= (1.0 + deriv/(x_order + 2));
  }
  series->n = order;
}

void create_mod_psi_up_to(int new_n)
{
  int n;
  double binom;

  if (new_n > n_modPsi) {
    int old = n_modPsi;
    n_modPsi = new_n;
    modPsi = realloc(modPsi, 2*n_modPsi*sizeof(Polynom));

    binom = 1.0;
    for (n = 0; n < old; n++)
      binom *= (-0.5 - n)/(double)(n+1);

    for (; n < n_modPsi; n++) {
      preparePolygammaEven(n, binom, &modPsi[2*n]);
      preparePolygammaOdd(n, binom, &modPsi[2*n + 1]);
      binom *= (-0.5 - n)/(double)(n+1);
    }
  }
}
