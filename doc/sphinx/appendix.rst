.. _Appendix:

Appendix
========

.. _The MMM family of algorithms:

The MMM family of algorithms
----------------------------

.. _MMM_Introduction:

Introduction
~~~~~~~~~~~~

In the MMM family of algorithms for the electrostatic interaction, a
convergence factor approach to tackle the conditionally convergent
Coulomb sum is used (even the authors of the original MMM method have no
idea what this acronym stands for). Instead of defining the summation
order, one multiplies each summand by a continuous factor
:math:`c(\beta,r_{ij},n_{klm})` such that the sum is absolutely
convergent for :math:`\beta>0`, but :math:`c(0,.,.)=1`. The energy is
then defined as the limit :math:`\beta\rightarrow 0` of the sum, i. e.
:math:`\beta` is an artificial convergence parameter. For a convergence
factor of :math:`e^{-\beta n_{klm}^2}` the limit is the same as the spherical
limit, and one can derive the classical Ewald method quite conveniently through
this approach :cite:`smith81a`. To derive the formulas for MMM,
one has to use a different convergence factor, namely
:math:`e^{-\beta|r_{ij}+n_{klm}|}`, which defines the alternative energy

.. math::

   \tilde{E}=\,\frac{1}{2}\lim_{\beta\rightarrow
     0}\sum_{k,l,m}{\sum_{i,j=1}^N}' \frac{q_i q_je^{-\beta|p_{ij} +
       n_{klm}|}} {|p_{ij} + n_{klm}|}
   =:\,\frac{1}{2}\lim_{\beta\rightarrow 0}\sum_{i,j=1}^N
   q_iq_j\phi_\beta(x_{ij}, y_{ij},z_{ij}).

:math:`\phi_\beta` is given by
:math:` \phi_\beta(x,y,z)=\,\tilde\phi_\beta(x,y,z)
+ \frac{e^{-\beta r}}{r} ` for :math:`(x,y,z)\neq 0` and
:math:`\phi_\beta(0,0,0)=\,\tilde\phi_\beta(0,0,0)`, where

.. math::

   \tilde\phi_\beta(x,y,z)=\,\sum_{(k,l,m)\neq 0} \frac{e^{-\beta
       r_{klm}}}{r_{klm}}.

The limit :math:`\tilde{E}` exists, but differs for three dimensionally
periodic systems by some multiple of the square of the dipole moment
from the spherical limit as obtained by the Ewald
summation :cite:`smith81a`. From the physical point of view
the Coulomb interaction is replaced by a screened Coulomb interaction
with screening length :math:`1/\beta`. :math:`\tilde{E}` is then the
energy in the limit of infinite screening length. But because of the
conditional convergence of the electrostatic sum, this is not
necessarily the same as the energy of an unscreened system. Since the
difference to the Ewald methods only depends on the dipole moment of the
system, the correction can be calculated easily in linear time and can
be ignored with respect to accuracy as well as to computation time.

For one or two dimensionally systems, however, :math:`\tilde{E}=E`, the
convergence factor approach equals the spherical summation limit of the
Ewald sum, and MMM1D and MMM2D do not require a dipole correction.

Starting from this convergence factor approach, Strebel constructed a
method of computational order :math:`O(N\log N)`, which is called MMM
:cite:`strebel99a`. The favourable scaling is obtained,
very much like in the Ewald case, by technical tricks in the calculation
of the far formula. The far formula has a product decomposition and can
be evaluated hierarchically similarly to the fast multipole methods.

For particles sufficiently separated in the z-axis one can Fourier
transform the potential along both x and y. We obtain the far formula as

.. math::

   \phi(x,y,z) =\, u_x u_y\sum_{p,q\neq 0} \frac{e^{2\pi f_{pq}z} +
     e^{2\pi f_{pq}(\lambda_z-z)}}{f_{pq} \left(e^{2\pi f_{pq}\lambda_z}
       - 1\right)} e^{2\pi i u_y q y}e^{2\pi i u_x p x} + 2\pi u_x
   u_y\left(u_z z^2 - z + \frac{\lambda_z}{6}\right).

where :math:`\lambda_{x,y,z}` are the box dimensions, :math:` f_{pq} =\,
\sqrt{(u_x p)^2 + (u_y q)^2},\quad f_p =\, u_x p,\quad f_q =\, u_x q
`, :math:` \omega_p=2\pi u_x p` and :math:`\omega_q=2\pi u_y q`. The
advantage of this formula is that it allows for a product decomposition
into components of the particles. For example

.. math::

   e^{2\pi f_{pq}z}=e^{2\pi f_{pq}(z_i-z_j)}=e^{2\pi
     f_{pq}z_i}e^{-2\pi f_{pq}z_j}

etc. Therefore one just has to calculate the sum over all these
exponentials on the left side and on the right side and multiply them
together, which can be done in :math:`O(N)` computation time. As can be
seen easily, the convergence of the series is excellent as long as z is
sufficiently large. By symmetry one can choose the coordinate with the
largest distance as z to optimise the convergence. Similar to the Lekner
sum, we need a different formula if all coordinates are small, i. e. for
particles close to each other. For sufficiently small :math:`u_y\rho`
and :math:`u_xx` we obtain the near formula as

.. math::

   \begin{array}{rl} \tilde\phi(x,y,z)=\, & 2 u_x
     u_y\sum\limits_{p,q>0} \frac{\cosh(2\pi f_{pq}z)}{f_{pq}
       \left(e^{2\pi f_{pq}\lambda_z} - 1\right)} e^{2\pi i u_y q
       y}e^{2\pi i u_x p x} +\\ & 4u_x\sum\limits_{l,p>0}\left(K_0(2\pi
       u_x p\rho_l) + K_N(2\pi u_x p\rho_{-l})\right)cos(2\pi u_x p x)
     -\\ & 2u_x\sum\limits_{n\ge 1}\frac{b_{2n}}{2n(2n)!}\Re\bigl((2\pi
     u_y (z+iy))^{2n}\bigr) +\\ & u_x\sum\limits_{n\ge
       0}\left(\begin{array}{c}-\frac{1}{2}\\
         n\end{array}\right)\frac{\left( \psi^{(2n)}(1 + u_x x) +
         \psi^{(2n)}(1 - u_x x)\right)}{(2n)!}\rho^{2n} -\\ &
     2\log(4\pi). \end{array}

Note that this time we calculate :math:`\tilde{\phi}` instead of
:math:`\phi`, i. e. we omit the contribution of the primary simulation
box. This is very convenient as it includes the case of self energy and
makes :math:`\tilde{\phi}` a smooth function. To obtain :math:`\phi` one
has to add the :math:`1/r` contribution of the primary box. The self
energy is given by

.. math::

   \tilde\phi(0,0,0)=\, 2 u_x u_y\sum\limits_{p,q>0} \frac{1}{f_{pq}
     \left(e^{2\pi f_{pq}\lambda_z} - 1\right)}+
   8u_x\sum\limits_{l,p>0}K_N(2\pi u_x\lambda_y p l) + 2 u_x\psi^{(0)}(1)
   - 2\log(4\pi).

Both the near and far formula are derived using the same convergence
factor approach, and consequently the same singularity in :math:`\beta`
is obtained. This is important since otherwise the charge neutrality
argument does not hold.

To obtain the :math:`O(N\log N)` scaling, some algorithm tricks are
needed, which are not used in MMM1D, MMM2D or ELC and are therefore not
discussed here. For details, see :cite:`strebel99a`. 

.. _MMM2D theory:

MMM2D theory
~~~~~~~~~~~~

In the case of periodicity only in the x and y directions, the far
formula looks like

.. math::

   \begin{array}{rl} \phi(x,y,z) = \, & 4 u_x u_y\sum_{p,q>0}
     \frac{e^{-2\pi f_{pq}|z|}} {f_{pq}} \cos(\omega_p x)\cos(\omega_q y)
     +\\ & 2 u_x u_y\left(\sum_{q>0} \frac{e^{-2\pi f_q|z|}}{f_q}
       \cos(\omega_q y) + \sum_{p>0} \frac{e^{-2\pi f_p|z|}}{f_p}
       \cos(\omega_p x)\right) -\\ & 2\pi u_x u_y |z| \end{array}

,

and the near formula is

.. math::

   \begin{array}{rl} \tilde\phi(x,y,z)=\, &
     4u_x\sum_{l,p>0}\left(K_0(\omega_p\rho_l) +
       K_0(\omega_p\rho_{-l})\right)\cos(\omega_p x) -\\ & 2u_x\sum_{n\ge
       1}\frac{b_{2n}}{2n(2n)!} \Re\bigl((2\pi u_y
     (z+iy))^{2n}\bigr)\,+\, \sum_{k=1}^{N_\psi-1}\left(\frac{1}{r_{k}} +
       \frac{1}{r_{-k}}\right) -\\ & u_x\sum_{n\ge
       0}\left(\begin{array}{c}-\frac{1}{2}\\n\end{array}\right)\frac{\left(
         \psi^{(2n)}(N_\psi + u_x x) + \psi^{(2n)}(N_\psi - u_x
         x)\right)}{(2n)!}(u_x\rho)^{2n} -\\ &
     2u_x\log\left(4\pi\frac{u_y}{u_x}\right). \end{array}

As said before, the energy obtained from these potentials is equal to
the electrostatic energy obtained by the spherical summation limit. The
deeper reason for this is that in some sense the electrostatic sum is
absolutely convergent :cite:`arnold02a`.

The near formula is used for particles with a small distance along the z
axis, for all other particles the far formula is used. Below is shown,
that the far formula can be evaluated much more efficiently, however,
its convergence breaks down for small z distance. To efficiently
implement MMM2D, the layered cell system is required, which splits up
the system in equally sized gaps along the z axis. The interaction of
all particles in a layer S with all particles in the layers S-1,S,S+1 is
calculated using the near formula, for the particles in layers
:math:`1,\dots,S-2`, and in layers :math:`S+2,\dots,N`, the far formula
is used.

The implementation of the near formula is relatively straight forward
and can be treated as any short ranged force is treated using the link
cell algorithm, here in the layered variant. The special functions in
the formula are somewhat demanding, but for the polygamma functions
Taylor series can be achieved, which are implemented in mmm-common.h.
The Bessel functions are calculated using a Chebychev series.

The treatment of the far formula is algorithmically more complicated.
For a particle i in layer :math:` S_i`, the formula can product
decomposed, as in

.. math::

   \begin{array}{rl} \sum_{j\in I_S, S < S_i - 1} q_iq_j\frac{e^{-2\pi
         f_{pq}|z_i-z_j|}}{f_{pq}} \cos(\omega_p (x_i -
     x_j))\cos(\omega_q (y_i - y_j)) = \\
     q_i\frac{e^{-2\pi f_{pq}z_i}}{f_{pq}} \cos(\omega_p
     x_i)\cos(\omega_q y_i) \sum_{j\in I_S, S < S_i - 1}q_je^{2\pi
       f_{pq}z_j} \cos(\omega_p x_j)\cos(\omega_q y_j) + \\
     q_i\frac{e^{-2\pi f_{pq}z_i}}{f_{pq}} \cos(\omega_p
     x_i)\sin(\omega_q y_i) \sum_{j\in I_S, S < S_i - 1}q_je^{2\pi
       f_{pq}z_j} \cos(\omega_p x_j)\sin(\omega_q y_j) + \\
     q_i\frac{e^{-2\pi f_{pq}z_i}}{f_{pq}} \sin(\omega_p
     x_i)\cos(\omega_q y_i) \sum_{j\in I_S, S < S_i - 1}q_je^{2\pi
       f_{pq}z_j} \sin(\omega_p x_j)\cos(\omega_q y_j) + \\
     q_i\frac{e^{-2\pi f_{pq}z_i}}{f_{pq}} \sin(\omega_p
     x_i)\sin(\omega_q y_i) \sum_{j\in I_S, S < S_i - 1}q_je^{2\pi
       f_{pq}z_j} \sin(\omega_p x_j)\sin(\omega_q y_j). \end{array}

This representation has the advantage, that the contributions of the two
particles are decoupled. For all particles j only the eight terms

.. math::

   \xi^{(\pm,s/c,s/c)}_j= q_je^{\pm 2\pi f_{pq}z_j} \sin/\cos(\omega_p
   x_j)\sin/\cos(\omega_q y_j)

are needed. The upper index describes the sign of the exponential term
and whether sine or cosine is used for :math:`x_j` and :math:`y_j` in
the obvious way. These terms can be used for all expressions on the
right hand side of the product decomposition. Moreover it is easy to see
from the addition theorem for the sine function that these terms also
can be used to calculate the force information up to simple prefactors
that depend only on p and q.

Every processor starts with the calculation of the terms
:math:`\xi^{(\pm,s/c,s/c)}_j` and adds them up in each layer, so that
one obtains

.. math:: \Xi^{(\pm,s/c,s/c)}_s= \sum_{j\in S_s}\xi^{(\pm,s/c,s/c)}_j.

Now we calculate

.. math:: \Xi^{(l,s/c,s/c)}_s=\sum_{t < s - 1}\Xi^{(+,s/c,s/c)}_t

and

.. math:: \Xi^{(h,s/c,s/c)}_s=\sum_{t > s + 1}\Xi^{(-,s/c,s/c)}_t,

which are needed for the evaluation of the product decomposition. While
the bottom processor can calculate :math:`\Xi^{(l,s/c,s/c)}_s` directly,
the other processors are dependent on its results. Therefore the bottom
processor starts with the calculation of its :math:`\Xi^{(l,s/c,s/c)}_s`
and sends up :math:`\Xi^{(l,s/c,s/c)}_s` and :math:`\Xi^{(+,s/c,s/c)}_s`
of its top layer s to the next processor dealing with the layers above.
Simultaneously the top processor starts with the calculation of the
:math:`\Xi^{(h,s/c,s/c)}_s` and sends them down. After the communicated
has been completed, every processor can use the
:math:`\Xi^{(l/h,s/c,s/c)}_j` and the :math:`\xi^{(\pm,s/c,s/c)}_j` to
calculate the force rsp. energy contributions for its particles.

In pseudo code, the far formula algorithm looks like:

#. for each layer :math:`s=1,\ldots,S`

   #. :math:`\Xi^{(\pm,s/c,s/c)}_s=0`

   #. for each particle :math:`j` in layer :math:`s`

      #. calculate :math:`\xi^{(\pm,s/c,s/c)}_j`

      #. :math:`\Xi^{(\pm,s/c,s/c)}_s += \xi^{(\pm,s/c,s/c)}_j`

#. :math:`\Xi^{(l,s/c,s/c)}_3=\Xi^{(+,s/c,s/c)}_1`

#. for each layer :math:`s=4,\ldots,S`

   #. .. math::
         \Xi^{(l,s/c,s/c)}_s=\Xi^{(l,s/c,s/c)}_{s-1} +\
          \Xi^{(+,s/c,s/c)}_{s-2}

#. :math:`\Xi^{(l,s/c,s/c)}_{S-2}=\Xi^{(-,s/c,s/c)}_S`

#. for each layer :math:`s=(S-3),...,1`

   #. .. math::
         \Xi^{(l,s/c,s/c)}_s=\Xi^{(l,s/c,s/c)}_{s+1} +\
          \Xi^{(-,s/c,s/c)}_{s+2}

#. for each layer :math:`s=1,...,S`

   #. for each particle :math:`j` in layer :math:`s`

      #. calculate particle interaction from
         :math:`\xi^{(+,s/c,s/c)}_j\Xi^{(l,s/c,s/c)}_s` and
         :math:`\xi^{(-,s/c,s/c)}_j\Xi^{(h,s/c,s/c)}_s`

For further details, see :cite:`arnold02a,arnold02b,arnold02c`.

.. _Dielectric contrast:

Dielectric contrast
^^^^^^^^^^^^^^^^^^^

A dielectric contrast at the lower and/or upper simulation box boundary
can be included comparatively easy by using image charges. Apart from
the images of the lowest and topmost layer, the image charges are far
enough to be treated by the far formula, and can be included as starting
points in the calculation of the :math:`\Xi` terms. The remaining
particles from the lowest and topmost layer are treated by direct
summation of the near formula.

This means, that in addition to the algorithm above, one has to only a
few things: during the calculation of the particle and cell blocks
:math:`\xi` and :math:`\Xi`, one additionally calculates the
contributions of the image charges and puts them either in a separate
array or, for the boundary layers, into two extra :math:`\xi` cell
blocks outside the simulation box. The entries in the separate array are
then added up over all processors and stored in the :math:`\Xi`-terms of
the lowest/topmost layer. This are all modifications necessary for the
far formula part. In addition to the far formula part, there is an
additional loop over the particles at the boundary to directly calculate
their interactions with their images. For details, refer to
:cite:`icmmm2d`.

.. _MMM1D theory:

MMM1D theory
~~~~~~~~~~~~

In one dimensionally periodic systems with z being the periodic
coordinate, the far formula looks like

.. math::

   \begin{array}{rl} \phi(\rho,z) &=\, 4 u_z\sum_{p\neq 0}
     K_0(\omega\rho)\cos(\omega z) - 2u_z\log(\frac{\rho}{2\lambda_z}) -
     2u_z\gamma\\ F_\rho(\rho,z) &=\, 8\pi u_z^2\sum_{p\neq 0} p
     K_1(\omega\rho)\cos(\omega z) + \frac{2 u_z}{\rho}\\ F_z(\rho,z)
     &=\, 8\pi u_z^2 \sum_{p\neq 0} pK_0(\omega\rho)\sin(\omega z),
   \end{array}

the near formula is

.. math::

   \begin{array}{rl} \tilde{\phi}(\rho,z) &=\, -u_z\sum_{n\ge 0}
     \left(\begin{array}{c}-\frac{1}{2}\\n\end{array}\right)
     \frac{\left(\psi^{(2n)}(N_\psi + u_z z) + \psi^{(2n)}(N_\psi - u_z
         z)\right)}{(2n)!}(u_z\rho)^{2n} - 2u_z\gamma + \\
     &\phantom{=\,++}
     \sum_{k=1}^{N_\psi-1}\left(\frac{1}{r_k}+\frac{1}{r_{-k}}\right)\\
     \tilde{F}_\rho(\rho,z) &=\, -u_z^3 \sum_{n\ge 0}
     \left(\begin{array}{c}-\frac{1}{2}\\n\end{array}\right)
     \frac{\left(\psi^{(2n)}(N_\psi + u_z z) + \psi^{(2n)}(N_\psi - u_z
         z)\right)}{(2n)!}(u_z\rho)^{2n-1} + \\ &\phantom{=\,++}
     \sum_{k=1}^{N_\psi-1}\left(\frac{\rho}{r_k^3}+\frac{\rho}{r_{-k}^3}\right)
     \\ \tilde{F}_z(\rho,z) &=\, -u_z^2 \sum_{n\ge 0}
     \left(\begin{array}{c}-\frac{1}{2}\\n\end{array}\right)
     \frac{\left(\psi^{(2n + 1)}(N_\psi + u_z z) + \psi^{(2n + 1)}(N_\psi
         - u_z z)\right)}{(2n)!}(u_z\rho)^{2n} + \\ &\phantom{=\,++}
     \sum_{k=1}^{N_\psi-1}\left(\frac{z+k\lambda_z}{r_k^3}+\frac{z-k\lambda_z}{r_{-k}^3}\right),
   \end{array}

where :math:`\rho` denotes the xy-distance of the particles. As for the
two dimensional periodic case, the obtained energy is equal to the one
dimensional Ewald sum. Algorithmically, MMM1D is uninteresting, since
neither the near nor far formula allow a product decomposition or
similar tricks. MMM1D has to be implemented as a simple NxN loop.
However, the formulas can be evaluated efficiently, so that MMM1D can
still be used reasonably for up to 400 particles on a single processor
:cite:`arnold05a`.

.. _ELC theory:

ELC theory
~~~~~~~~~~

The ELC method differs from the other MMM algorithms in that it is not
an algorithm for the calculation of the electrostatic interaction, but
rather represents a correction term which allows to use any method for
threedimensionally periodic systems with spherical summation order for
twodimensional periodicity. The basic idea is to expand the two
dimensional slab system of height h in the non-periodic z-coordinate to
a system with periodicity in all three dimensions, with a period of
:math:`\lambda_z>h`, which leaves an empty gap of height
:math:`\delta=\lambda_z -
h` above the particles in the simulation box.

Since the electrostatic potential is only finite if the total system is
charge neutral, the additional image layers (those layers above or below
the original slab system) are charge neutral, too. Now let us consider
the n-th image layer which has an offset of :math:`n\lambda_z` to the
original layer. If :math:`n\lambda_z` is large enough, each particle of
charge q\_j at position :math:`(x_j,y_j,z_j+n\lambda_z)` and its
replicas in the xy-plane can be viewed as constituting a homogeneous
charged sheet of charge density
:math:`\sigma_j = \frac{q_j}{\lambda_x\lambda_y}`. The potential of such
a charged sheet at distance z is :math:`2\pi \sigma_j
|z|`. Now we consider the contribution from a pair of image layers
located at :math:`\pm n\lambda_z`, n>0 to the energy of a charge q\_i at
position :math:`(x_i,y_i,z_i)` in the central layer. Since
:math:`|z_j - z_i| <
n\lambda_z`, we have
:math:`|z_j - z_i + n\lambda_z| = n\lambda_z + z_j -
z_i` and :math:`|z_j - z_i - n\lambda_z|= n\lambda_z - z_j + z_i`, and
hence the interaction energy from those two image layers with the charge
:math:`q_i` vanishes by charge neutrality:

.. math::

   2\pi q_i \sum_{j=1}^N \sigma_j(|z_j - z_i + n\lambda_z| + |z_j -
   z_i - n\lambda_z|) = 4\pi q_i n\lambda_z \sum_{j=1}^N \sigma_j = 0.

The only errors occurring are those coming from the approximation of
assuming homogeneously charged, infinite sheets instead of discrete
charges. This assumption should become better when increasing the
distance :math:`n\lambda_z` from the central layer.

However, in a naive implementation, even large gap sizes will result in
large errors. This is due to the order of summation for the standard
Ewald sum, which is spherical, while the above approach orders the cells
in layers, called slab–wise summation. Smith has shown that by adding to
the Ewald energy the term

.. math:: E_c=2\pi M_z^2 - \frac{2\pi M^2}{3},

where M is the total dipole moment, one obtains the result of a
slab–wise summation instead of the spherical limit
:cite:`smith81a`. Although this is a major change in the
summation order, the difference is a very simple term. In fact, Smith
shows that changes of the summation order always result in a difference
that depends only on the total dipole moment.

Using the far formula of MMM2D, one can calculate the contributions of
the additional layers up to arbitrarily precision, even for small gap
sizes. This method is called electrostatic layer correction, ELC. The
advantage of this approach is that for the image layers, z is
necessarily large enough, so that all interactions can be represented
using the product decomposition. This allows for an order N evaluation
of the ELC term.

The electrostatic layer correction term is given by

.. math:: E_{lc}=\sum_{i,j=1}^Nq_iq_j\psi(p_i-p_j),

where

.. math::

   \begin{array}{rl} \psi(x,y,z)=&4u_xu_y\sum_{p,q>0}\frac{\cosh(2\pi
       f_{pq}z)}{f_{pq}(e^{2\pi f_{pq}\lambda_z} - 1)} \cos(\omega_p
     x)\cos(\omega_q y) + \\ &2u_xu_y\sum_{p>0}\frac{\cosh(2\pi f_p
       z)}{f_p(e^{2\pi f_p\lambda_z} - 1)}\cos(\omega_p x)+\\
     &2u_xu_y\sum_{q>0}\frac{\cosh(2\pi f_q z)}{f_q(e^{2\pi f_q\lambda_z}
       - 1)}\cos(\omega_q y). \end{array}

The implementation is very similar to MMM2d, except that the separation
between slices closeby, and above and below is not necessary.

.. _Errors:

Errors
~~~~~~

Common to all algorithms of the MMM family is that accuracy is cheap
with respect to computation time. More precisely, the maximal pairwise
error, i.e. the maximal error of the :math:`\psi` expression, decreases
exponentially with the cutoffs. In turn, the computation time grows
logarithmically with the accuracy. This is quite in contrast to the
Ewald methods, for which decreasing the error bound can lead to
excessive computation time. For example, P3M cannot reach precisions
above :math:`10^{-5}` in general. The precise form of the error
estimates is of little importance here, for details see
:cite:`arnold02c`.

One important aspect is that the error estimates are also exponential in
the non-periodic coordinate. Since the number of closeby and far away
particles is different for particles near the border and in the center
of the system, the error distribution is highly non–homogenous. This is
unproblematic as long as the maximal error is really much smaller than
the thermal energy. However, one cannot interpret the error simply as an
additional error source.

.. figure:: figures/elc-errordist.pdf
   :alt: Error distribution of the ELC method.

   Error distribution of the ELC method.

Figure [fig:ELC-error] shows the error distribution of the ELC method
for a gap size of :math:`10\%` of the total system height. For MMM2D and
MMM1D the error distribution is less homogenous, however, also here it
is always better to have some extra precision, especially since it is
computationally cheap.

.. _Maxwell Equations Molecular Dynamics (MEMD):

Maxwell Equations Molecular Dynamics (MEMD)
-------------------------------------------

In this chapter, we want to give a more thorough introduction to the
MEMD (or “Maggs”) algorithm for the calculation of Coulomb interactions
that is implemented in |es|. For an even more detailed description, we refer
to the publications :cite:`maggs02a,pasichnyk04a`. The method is intimately
related to the Car–Parrinello approach, while being equivalent to solving
Maxwell’s equations with freely adjustable speed of light.

.. _Equations of motion:

Equations of motion
~~~~~~~~~~~~~~~~~~~

Denoting the particle masses with :math:`m_i`, their charges with
:math:`q_i`, their coordinates and momentum with :math:`\vec r_i` and
:math:`\vec p_i` respectively, the interparticle potential (of
*non*-electromagnetic type) with :math:`U`, for the coupled system of
charges and fields we write the following equations of motion

.. math::

   \begin{aligned}
     \dot{\vec r}_i & = & \frac{1}{m_i} \vec p_i \\
     \dot{\vec p}_i & = & - \frac{\partial U}{\partial \vec r_i} + q_i \vec E (\vec r_i)- \frac{\zeta}{m_i} \vec p_i
                           + \vec f_i \\
     \dot{\vec A} & = & - \vec E \\
     \dot{\vec E} & = & 
     c^2 \vec \nabla \times \left( \vec \nabla \times \vec A \right)
     - \frac{1}{\epsilon_0} \vec j ,\end{aligned}

where :math:`\epsilon_0` is the vacuum dielectric constant, :math:`c`
the speed of light, :math:`\vec A` the vector-potential, :math:`\vec E`
the electric field, :math:`\vec j` the current density; :math:`\zeta` is
the particle friction constant, and :math:`\vec f_i` is a random force
satisfying the standard fluctuation-dissipation theorem:

.. math::

   \left< f_i^\alpha (t) f_j^\beta (t^\prime) \right> =
   2 \zeta k_B T \delta_{ij} \delta_{\alpha \beta}
   \delta (t - t^\prime),

where :math:`\alpha` and :math:`\beta` denote Cartesian indices.

If we introduce the vector :math:`\vec B=\nabla\times A` the system of
equations can be rewritten in a form similar to the usual Maxwell
equations. Currently in the version with :math:`\vec B` and
:math:`\vec E` is implemented.

.. _Discretization:

Discretization
~~~~~~~~~~~~~~

For implementation on the computer, the equations need to be discretized
with respect to both space and time.We consider a domain of physical
space as being an affine space and divide it into subdomains of
contiguous cells of cubic shape. The charges live on the vertices of our
lattice which has the spacing :math:`a`. The electric fields
:math:`E(l)` and vector potentials :math:`A(l)` live on the edges or
links and are aligned with them. We need also the operator
:math:`\nabla\times{}`. It gives the vector :math:`\vec B`, which lives
on the faces of the cube or on the plaquettes, Fig. [fig:cellstructure].

.. figure:: figures/maggs-rotation.pdf
   :alt: Spatial elements of a cell complex

   Spatial elements of a cell complex

In the implementation of the algorithm we assume that particles with
masses :math:`m_i` and charges :math:`q_i` live in the continuum
(off–lattice approach). The charges are interpolated on the lattice with
grid spacing :math:`a` using a linear interpolation scheme.

.. _Initialization of the algorithm:

Initialization of the algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The algorithm as it is implemented only calculates stepwise time updates
of the exact field solution. Therefore in order to start the simulation
for the given random distribution of charges we have to calculate the
initial electrostatic field, i. e. the exact solution of the
electrostatic problem. We find a particular solution of Gauss’ law as
the result of the following recursive procedure (see
Fig. [fig:maggs-initialization]):

#. The charge in the plane :math:`z=z_\text{plane}` is

   .. math:: q_\text{plane}=\frac{1}{N_z}\sum_iq(\vec r_i)\delta(z_i-z_\text{plane}),

   :math:`N_z` is the number of charges in plane
   :math:`z=z_\text{plane}`. Update the :math:`z`-field according to the
   formula

   .. math:: E_z^2=E_z^1+\frac{q_\text{plane}}{\epsilon_0a^2};

#. Subtract the charge :math:`q_\text{plane}` from the each charge on
   sites of :math:`z_\text{plane}`. The charge of the wire
   :math:`y=y_\text{wire}, z=z_\text{plane}` is

   .. math:: q_\text{wire}=\frac{1}{N_y}\sum_iq(\vec r_i)\delta(z_i-z_\text{plane})\delta(y_i-y_\text{wire}),

   :math:`N_y` now meaning the number of charges in the wire. Update
   :math:`y`-field

   .. math:: E_y^2=E_y^1+\frac{q_\text{wire}}{\epsilon_0a^2};

#. Subtract the charge :math:`q_\text{wire}` from the each charge on the
   sites of :math:`(y_\text{wire},z_\text{plane})`. Update :math:`x`
   field

   .. math:: E_x^2=E_x^1+\frac{q_\text{vertex}}{\epsilon_0a^2}

This scheme is repeated until the fields are completely relaxed (i. e.
the energy is minimized). During repetition, the spatial dimensions are
permutated to avoid a drift in one direction.

.. figure:: figures/maggs-initial-scheme.pdf
   :alt: Recursive solution of Gauss’ law

   Recursive solution of Gauss’ law

.. _Time integrator:

Time integrator
~~~~~~~~~~~~~~~

For the time discretization we have adopted the elegant solution which
was found by Rottler and Maggs :cite:`maggs02a` and allows
to conserve *both* time–reversibility and phase–space volume
conservation:

#. Update the particle momenta by half a time step.

#. Update the :math:`\vec B` field by half a time step.

#. Update the particle positions in :math:`x` direction by half a time
   step.

#. Update the electric field in :math:`x` direction by half a time step.

#. Update the particle positions in :math:`y` direction by half a time
   step.

#. Update the electric field in :math:`y` direction by half a time step.

#. Update the particle positions in :math:`z` direction by half a time
   step.

#. Update the electric field in :math:`z` direction by a full time step.

#. Update the particle positions in :math:`z` direction by half a time
   step.

#. Update the electric field in :math:`y` direction by half a time step.

#. Update the particle positions in :math:`y` direction by half a time
   step.

#. Update the electric field in :math:`x` direction by half a time step.

#. Update the particle positions in :math:`x` direction by half a time
   step.

#. Update the :math:`\vec B` field by half a time step.

#. Update the particle momenta by half a time step.

.. _Self–energy:

Self–energy
~~~~~~~~~~~

The interpolation of the charges onto the lattice gives rise to the
artificial force exerted on the particle by its own field. In order to
cure this remedy, the direct subtraction of the self–energy is
introduced.

For the interpolated charge cloud the self–energy can be directly
calculated. For the simple cubic lattice in three dimensions the linear
interpolation will give 8 charges which are placed at the corners of the
cube with edge length :math:`a` (see Fig. [fig:charge-assignment]).

.. figure:: figures/maggs-charge-assignment.pdf
   :alt: Linear interpolation scheme

   Linear interpolation scheme

Therefore in our case the self-energy is a symmetric bilinear form
defined by the matrix :math:`\left\{\alpha_{ij}\right\}`, the elements
of which do not depend on the position of the charge. In our algorithm
the values of the coefficients are

.. math::

   \alpha_{ij}=\frac{1}{4a\epsilon_0L^3}\sum\limits_{\vec k}
     \frac{\cos \vec k(\vec R_{\imath}-\vec R_{\jmath})}
     {\sum_{\imath=1}^3(1-\cos\vec k\vec a_{\imath})}

where :math:`L` is the number of lattice points per dimension,
:math:`\vec R_i` coordinates of the interpolated charges and
:math:`\vec k` the wave vector. Those values are calculated during the
initialization step and are used in the calculation of the self-force.
The value of the self-force which has to be subtracted from the overall
forces is given by the following ansatz

.. math::

   \vec F_{self}=-\frac{\partial \mathcal U_{self}}{\partial\vec r}
     =-\sum\limits_i\sum\limits_j\alpha_{ij}
     \left[q_i\frac{\partial q_j}{\partial\vec r}
       +q_j\frac{\partial q_i}{\partial\vec r}\right].

.. _For which systems to use the algorithm:

For which systems to use the algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Although it is not very well known by now, this algorithm is a promising
alternative to the often used Ewald-based methods. The main advantages
and disadvantages shall be named here. However, it is still best to
understand the concept of the algorithm and figure out for yourself, if
it may be an option.

-  The fields are not calculated for an arbitrary charge distribution,
   but updated from the last solution. Therefore, particles should not
   move too much between timesteps (less than a lattice cube).

-  No procedure for error tuning yet. You have to adjust the parameters
   and determine the error yourself.

-  Only 3D periodic systems are possible for now.

-  With the given interpolation scheme, the short-range part of the
   potential is highy underestimated when two particles are in the same
   lattice cube!

-  The initialization routine scales with :math:`\mathcal{O}(N^3)` and
   takes a long time for larger (and also inhomogenous) systems.

-  The algorithm is a local update scheme and spatially varying
   properties can be applied (in the future).

-  Because of the locality, the algorithm itself scales
   :math:`\mathcal{O}(N)` and has a big advantage in speed for larger
   systems.

-  Because of the locality, it is highly parallelized.

-  It is fast.

The last item is of course dependent on the system properties. But if
the charges are evenly distributed and the system is not too sparse,
this algorithm outperforms P3M easily. Especially for systems with more
than 1000 charges.

Of course, if the system is not dense enough, one will have to set the
lattice spacing in a way to avoid several particles in one cell and the
mesh will be very fine for not so many charges. Also, if you have lots
of charges but your simulation should only run for a short time, the
initialization scheme takes too long in comparison.

But, if you have dense systems with more than 1000 charges or
simulations that run for many timesteps, this method is definitely an
option.


.. _Reaction Ensemble:

Reaction Ensemble
-----------------

.. note:: The whole Reaction Ensemble module uses Monte Carlo moves which require potential energies. Therefore the Reaction Ensemble requires support for energy calculations for all interactions which are used in the simulation.

The reaction ensemble :cite:`smith94a,turner2008simulation` allows to simulate
chemical reactions which can be represented by the general equation:

.. math::

   \mathrm{\nu_1 S_1 +\ \dots\  \nu_l S_l\ \rightleftharpoons\ \nu_m S_m +\ \dots\ \nu_z S_z }
       \label{general-eq}

where :math:`\nu_i` is the stoichiometric coefficient of species
:math:`S_i`. By convention, stoichiometric coefficents of the
species on the left-hand side of the reaction (*reactants*) attain
negative values, and those on the right-hand side (*products*) attain
positive values, so that the reaction can be equivalently written as

.. math::

   \mathrm{\sum_i \nu_i S_i = 0} \,.
       \label{general-eq-sum}


The equilibrium constant of the reaction is then given as

.. math::

   K = \exp(-\Delta_{\mathrm{r}}G^{\ominus} / k_B T)
       \quad\text{with}\quad
       \Delta_{\mathrm{r}}G^{\ominus} = \sum_i \nu_i \mu_i^{\ominus}\,.
       \label{Keq}


Here :math:`k_B` is the Boltzmann constant, :math:`T` is temperature,
:math:`\Delta_{\mathrm{r}}G^{\ominus}` standard Gibbs free energy change
of the reaction, and :math:`\mu_i^{\ominus}` the standard chemical
potential (per particle) of species :math:`i`. Note that thermodynamic equilibrium is
independent of the direction in which we write the reaction. If it is
written with left and righ-hand side swapped, 
both :math:`\Delta_{\mathrm{r}}G^{\ominus}` and the stoichiometric
coefficients attain opposite signs, and the equilibrium constant attains the inverse value. 
Further, note that the equilibrium constant :math:`K` is the
dimensionless *thermodynamic, concentration-based* equilibrium constant,
defined as

.. math::

   K(c^{\ominus}) = (c^{\ominus})^{-\bar\nu} \prod_i (c_i)^{\nu_i}

wher :math:`\bar\nu=\sum_i \nu_i`, and :math:`c^{\ominus}` is the reference concentration,
at which the standard chemical potential :math:`\Delta_{\mathrm{r}}G^{\ominus}` was determined.
In practice, this constant is often used with the dimension of :math:`(c^{\ominus})^{\bar\nu}`

.. math::

   K_c(c^{\ominus}) = K(c^{\ominus})\times (c^{\ominus})^{\bar\nu}

A simulation in
the reaction ensemble consists of two types of moves: the *reaction move*
and the *configuration move*. The configuration move changes the configuration
of the system. It is not performed by the Reaction Ensemble module, and can be
performed by a suitable molecular dynamics or a Monte Carlo scheme. The
``reacton_ensemble`` command takes care only of the reaction moves.
In the *forward* reaction, the appropriate number of reactants (given by
:math:`\nu_i`) is removed from the system, and the concomitant number of
products is inserted into the system. In the *backward* reaction,
reactants and products exchange their roles. The acceptance probability
:math:`P^{\xi}` for move from state :math:`o` to :math:`n` reaction
ensemble is given by the criterion :cite:`smith94a`

.. math::

   P^{\xi} = \text{min}\biggl(1,V^{\bar\nu\xi}\Gamma^{\xi}e^{-\beta\Delta E}\prod_{i=1}\frac{N_i^0!}{(N_i^0+\nu_{i}\xi)!}
       \label{eq:Pacc}
       \biggr),

where :math:`\Delta E=E_\mathrm{new}-E_\mathrm{old}` is the change in potential energy,
:math:`V` is the simulation box volume,
and :math:`\beta=1/k_\mathrm{B}T`. 
The extent of reaction, :math:`\xi=1` for the forward, and
:math:`\xi=-1` for the backward direction. 
The parameter :math:`\Gamma` proportional to the reaction constant. It is defined as

.. math::

   \Gamma = \prod_i \Bigl(\frac{\left<N_i\right>}{V} \Bigr)^{\bar\nu} = V^{-\bar\nu} \prod_i \left<N_i\right>^{\nu_i} = K_c(c^{\ominus}=1/\sigma^3)

where :math:`\left<N_i\right>/V` is the average number density of particles of type :math:`i`.
Note that the dimension of :math:`\Gamma` is :math:`V^{\bar\nu}`, therefore its
units must be consistent with the units in which Espresso measures the box volume,
i.e. :math:`\sigma^3`.
   
It is often convenient, and in some cases even necessary, that some particles
representing reactants are not removed from or placed at randomly in the system
but their identity is changed to that of the products, or vice versa in the
backward direction.  A typical example is the ionization reaction of weak
polyelectrolytes, where the ionizable groups on the polymer have to remain on
the polymer chain after the reaction.  The replacement rule is that the identity of a given reactant type is
changed to the corresponding product type as long as the corresponding
coefficients allow for it.  Corresponding means having the same position (index) in
the python lists of reactants and products which are used to set up the
reaction.
For a description of the available methods see :mod:`espressomd.reaction_ensemble`

.. _Converting tabulated reaction constants to internal units in Espresso:

Converting tabulated reaction constants to internal units in Espresso 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The implementation in Espresso requires that the dimension of :math:`\Gamma` 
is consistent with the internal unit of volume, :math:`\sigma^3`.
The tabulated values of equilibrium constants for reactions in solution, :math:`K_c`, typically use
:math:`c^{\ominus} = 1\,\mathrm{moldm^{-3}}` as the reference concentration, 
and have the dimension of :math:`(c^{\ominus})^{\bar\nu}`.  To be used with Espresso, the
value of :math:`K_c` has to be converted as

.. math::

   \Gamma = K_c(c^{\ominus} = 1/\sigma^3) = K_c(c^{\ominus} = 1\,\mathrm{moldm^{-3}}) 
   \Bigl( N_{\mathrm{A}}\bigl(\frac{\sigma}{\mathrm{dm}}\bigr)^3\Bigr)^{\bar\nu}
   
where :math:`N_{\mathrm{A}}` is the Avogardo number.  For gas-phase reactions,
the pressure-based eaction constant, :math:`K_p` is often used, which can
be converted to :math:`K_c` as

.. math::

   K_p(p^{\ominus}=1\,\mathrm{atm}) = K_c(c^{\ominus} = 1\,\mathrm{moldm^{-3}}) \biggl(\frac{c^{\ominus}RT}{p^{\ominus}}\biggr)^{\bar\nu},

where :math:`p^{\ominus}=1\,\mathrm{atm}` is the standard pressure. 



.. The text below is commented-out because it is still an open research question how it should be used correctly.
..
.. This can be used to include water autoprotolysis in the implicit solvent simulation, 
.. by means of a reaction:
.. 
.. .. math::
.. 
..    \mathrm{2 H_2O \rightleftharpoons\ H_3O^+ + OH^- } \,,
.. 
.. 
.. add the following ex nihilo reactions to Espresso. (:math:`\emptyset`, read ex
.. nihilo). Ex nihilo means that the reaction has no reactants or products.
.. Therefore, if :math:`\emptyset` is a product, particles vanish and if
.. :math:`\emptyset` is an reactant, then particles are created ex nihilo:
.. 
.. .. math::
.. 
..    \mathrm{\emptyset \rightleftharpoons\ H_3O^+ + OH^- }  \,, 
.. 
.. with reaction constant K
.. 
.. .. math::
.. 
..    \mathrm{H_3O^+ + OH^- \rightleftharpoons\ \emptyset} \,, 
.. 
.. with reaction constant 1/K. K is given implicitly as a function of the apparent dissociation
.. constant :math:`K_w=10^{-14} \rm{mol^2/l^2}=x\cdot \rm{1/(\sigma^3)^2}` such that the dimensionless is
.. :math:`K=(x\cdot \rm{1/(\sigma^3)^2})/(\beta P^0)^{\overline{\nu}}` with
.. :math:`\overline{\nu}=2` for the dissociation reaction and where x is
.. the value of the apparent dissociation constant that is converted from
.. :math:`\rm{mol^2/l^2}` to a number density in :math:`1/(\sigma^3)^2`,
.. where :math:`\sigma` is the simulation length unit. If :math:`\beta` and
.. :math:`P^0` are provided in simulation units this will make :math:`K`
.. dimensionless. As a test for the autodissociation of water a big
.. simulation box can be set up and the autodissociation reaction can be
.. performed. Then the box should fill with the correct number of protons
.. and hydroxide ions (check for the number of protons and hydroxide ions
.. in the given simulation volume and compare this to the expected value at
.. pH 7). Further the :math:`pK_w=14` should be reproduced -also in the
.. case of an initial excess of acid or base in the simulation box. Note
.. that this only works for big enough volumes.


.. _Wang-Landau Reaction Ensemble:

Wang-Landau Reaction Ensemble
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. .. note:: Requires support for energy calculations for all used interactions since it uses Monte-Carlo moves which use energies in one way or the other.

Combination of the Reaction Ensemble with the Wang-Landau algorithm 
:cite:`wang01a`
allows for enhanced sampling of the reacting system, and
and for the determination of the density of states with respect 
to the reaction coordinate or with respect to some other collective
variable :cite:`landsgesell16a`. Here the 1/t Wang-Landau
algorithm :cite:`belardinelli07a` is implemented since it
does not suffer from systematic errors. Additionally to the above
commands for the reaction ensemble use the following commands for the
Wang-Landau reaction ensemble. For a description of the available methods see :mod:`espressomd.reaction_ensemble`:

.. _Constant pH simulation using the Reaction Ensemble:

Constant pH simulation using the Reaction Ensemble
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. .. note:: Requires support for energy calculations for all used interactions since it uses Monte-Carlo moves which use energies.

In the constant pH method due to Reed and Reed
:cite:`reed92a` it is possible to set the chemical potential
of :math:`H^{+}` ions, assuming that the simulated system is coupled to an
infinite reservoir. This value is the used to simulate dissociation
equilibrium of acids and bases. Under certain conditions, the constant
pH method can yield equivalent results as the reaction ensemble :cite:`landsgesell16b`. However, it
treats the chemical potential of :math:`H^{+}` ions and their actual
number in the simulation box as independent variables, which can lead to
serious artifacts. 
The constant pH method can be used within the reaction ensemble module by
initializing the reactions with the standard commands of the reaction ensemble. 

The dissociation constant, which is the input of the constant pH method, is the equilibrium
constant :math:`K_c` for the following reaction:

.. math::

   \mathrm{HA \rightleftharpoons\ H^+ + A^- } \,,

For an example of how to setup
a Constant pH simulation, see the file in the testsuite directory. 
For a description of the available methods see :mod:`espressomd.reaction_ensemble`:

.. _Grand canonical ensemble simulation using the Reaction Ensemble:

Grand canonical ensemble simulation using the Reaction Ensemble
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As a special case, all stoichiometric coefficients on one side of the chemical
reaction can be set to zero.  Such reaction creates particles *ex nihilo*, and 
is equivalent to exchange with a reservoir. Then the simulation in the reaction ensemble becomes equivalent with the
grandcanonical simulation. Formally, this can be expressed by the reaction

.. math::
 
    \mathrm{\emptyset \rightleftharpoons\ \nu_A A  }  \,, 

where, if :math:`\nu_A=1`, the reaction constant :math:`\Gamma` defines the chemical potential of species A.
However, if :math:`\nu_A\neq 1`, the statistics of the reaction ensemble becomes
equivalent to the grandcanonical only in the limit of large average number of species A in the box.
It the reaction contains more than one product, then the reaction constant
:math:`\Gamma` defines only the sum of their chemical potentials but not the
chemical potential of each product alone.

.. Since the Reaction Ensemble acceptance transition probability can be
.. derived from the grand canonical acceptance transition probability we
.. can use the reaction ensemble to implement grand canonical simulation
.. moves. This is done via adding reactions that only have reactants (for the
.. deletion of particles) or only have products (for the creation of
.. particles). There exists a one to one mapping of the expressions in the
.. grand canonical transition probabilities and the expressions in the
.. reaction ensemble transition probabilities.

.. _Electrokinetics:

Electrokinetics
---------------

The electrokinetics setup in allows for the description of
electro-hydrodynamic systems on the level of ion density distributions
coupled to a Lattice-Boltzmann (LB) fluid. The ion density distributions
may also interact with explicit charged particles, which are
interpolated on the LB grid. In the following paragraph we briefly
explain the electrokinetic model implemented in , before we come to the
description of the interface.

If you are interested in using the electrokinetic implementation in for
scientific purposes, please contact G. Rempfer before you start your
project.

.. _Electrokinetic Equations:

Electrokinetic Equations
~~~~~~~~~~~~~~~~~~~~~~~~

In the electrokinetics code we solve the following system of coupled
continuity, diffusion-advection, Poisson, and Navier-Stokes equations:

.. math::

   \begin{aligned}
   \label{eq:ek-model-continuity} \frac{\partial n_k}{\partial t} & = & -\, \nabla \cdot \vec{j}_k \vphantom{\left(\frac{\partial}{\partial}\right)} ; \\
   \label{eq:ek-model-fluxes} \vec{j}_{k} & = & -D_k \nabla n_k - \nu_k \, q_k n_k\, \nabla \Phi + n_k \vec{v}_{\mathrm{fl}} \vphantom{\left(\frac{\partial}{\partial}\right)} ; \\
   \label{eq:ek-model-poisson} \Delta \Phi & = & -4 \pi \, {l_\mathrm{B}}\, {k_\mathrm{B}T}\sum_k q_k n_k \vphantom{\left(\frac{\partial}{\partial}\right)}; \\
   \nonumber \left(\frac{\partial \vec{v}_{\mathrm{fl}}}{\partial t} + \vec{v}_{\mathrm{fl}} \cdot \vec{\nabla} \vec{v}_{\mathrm{fl}} \right) \rho_\mathrm{fl} & = & -{k_\mathrm{B}T}\, \nabla \rho_\mathrm{fl} - q_k n_k \nabla \Phi \\
   \label{eq:ek-model-velocity} & & +\, \eta \vec{\Delta} \vec{v}_{\mathrm{fl}} + (\eta / 3 + \eta_{\text{\,b}}) \nabla (\nabla \cdot \vec{v}_{\mathrm{fl}}) \vphantom{\left(\frac{\partial}{\partial}\right)} ; \\
   \label{eq:ek-model-continuity-fl} \frac{\partial \rho_\mathrm{fl}}{\partial t} & = & -\,\nabla\cdot\left( \rho_\mathrm{fl} \vec{v}_{\mathrm{fl}} \right) \vphantom{\left(\frac{\partial}{\partial}\right)} , \end{aligned}

which define relations between the following observables

:math:`n_k`
    the number density of the particles of species :math:`k`,

:math:`\vec{j}_k`
    the number density flux of the particles of species :math:`k`,

:math:`\Phi`
    the electrostatic potential,

:math:`\rho_{\mathrm{fl}}`
    the mass density of the fluid,

:math:`\vec{v}_{\mathrm{fl}}`
    the advective velocity of the fluid,

and input parameters

:math:`D_k`
    the diffusion constant of species :math:`k`,

:math:`\nu_k`
    the mobility of species :math:`k`,

:math:`q_k`
    the charge of a single particle of species :math:`k`,

:math:`{l_\mathrm{B}}`
    the Bjerrum length,

:math:`{k_\mathrm{B}T}`
    | the thermal energy given by the product of Boltzmann’s constant
      :math:`k_\text{B}`
    | and the temperature :math:`T`,

:math:`\eta`
    the dynamic viscosity of the fluid,

:math:`\eta_{\text{\,b}}`
    the bulk viscosity of the fluid.

The temperature :math:`T`, and diffusion constants :math:`D_k` and
mobilities :math:`\nu_k` of individual species are linked through the
Einstein-Smoluchowski relation :math:`D_k /
\nu_k = {k_\mathrm{B}T}`. The system of equations described in Eqs. -,
combining diffusion-advection, electrostatics, and hydrodynamics is
conventionally referred to as the *Electrokinetic Equations*.

The electrokinetic equations have the following properties:

-  On the coarse time and length scale of the model, the dynamics of the
   particle species can be described in terms of smooth density
   distributions and potentials as opposed to the microscale where
   highly localized densities cause singularities in the potential.

   In most situations, this restricts the application of the model to
   species of monovalent ions, since ions of higher valency typically
   show strong condensation and correlation effects – the localization
   of individual ions in local potential minima and the subsequent
   correlated motion with the charges causing this minima.

-  Only the entropy of an ideal gas and electrostatic interactions are
   accounted for. In particular, there is no excluded volume.

   This restricts the application of the model to monovalent ions and
   moderate charge densities. At higher valencies or densities,
   overcharging and layering effects can occur, which lead to
   non-monotonic charge densities and potentials, that can not be
   covered by a mean-field model such as Poisson-Boltzmann or this one.

   Even in salt free systems containing only counter ions, the
   counter-ion densities close to highly charged objects can be
   overestimated when neglecting excluded volume effects. Decades of the
   application of Poisson-Boltzmann theory to systems of electrolytic
   solutions, however, show that those conditions are fulfilled for
   monovalent salt ions (such as sodium chloride or potassium chloride)
   at experimentally realizable concentrations.

-  Electrodynamic and magnetic effects play no role. Electrolytic
   solutions fulfill those conditions as long as they don’t contain
   magnetic particles.

-  The diffusion coefficient is a scalar, which means there can not be
   any cross-diffusion. Additionally, the diffusive behavior has been
   deduced using a formalism relying on the notion of a local
   equilibrium. The resulting diffusion equation, however, is known to
   be valid also far from equilibrium.

-  The temperature is constant throughout the system.

-  The density fluxes instantaneously relax to their local equilibrium
   values. Obviously one can not extract information about processes on
   length and time scales not covered by the model, such as dielectric
   spectra at frequencies, high enough that they correspond to times
   faster than the diffusive time scales of the charged species.

.. _Setup:

Setup
~~~~~

.. _Initialization:

Initialization
^^^^^^^^^^^^^^

The command initializes the LB fluid with a given set of parameters, and
it is very similar to the Lattice-Boltzmann command in set-up. We
therefore refer the reader to Chapter [sec:lb] for details on the
implementation of LB in and describe only the major differences here.

The first major difference with the LB implementation is that the
electrokinetics set-up is a Graphics Processing Unit (GPU) only
implementation. There is no Central Processing Unit (CPU) version, and
at this time there are no plans to make a CPU version available in the
future. To use the electrokinetics features it is therefore imperative
that your computer contains a CUDA capable GPU which is sufficiently
modern.

To set up a proper LB fluid using the command one has to specify at
least the following options: , , , , , and . The other options can be
used to modify the behavior of the LB fluid. Note that the command does
not allow the user to set the time step parameter as is the case for the
command, this parameter is instead taken directly from the input of the
``t_step`` command. The LB *mass density* is set independently from the
electrokinetic *number densities*, since the LB fluid serves only as a
medium through which hydrodynamic interactions are propagated, as will
be explained further in the next paragraph. If no is specified, then our
algorithm assumes = 1.0. The two ‘new’ parameters are the temperature at
which the diffusive species are simulated and the Bjerrum length
associated with the electrostatic properties of the medium. See the
above description of the electrokinetic equations for an explanation of
the introduction of a temperature, which does not come in directly via a
thermostat that produces thermal fluctuations.

can be set to *on* or *off*. It controls whether there should be an
advective contribution to the diffusive species’ fluxes. Default is
*on*.

can be set to *friction* or *estatics*. This option determines the force
term acting on the fluid. The former specifies the force term to be the
sum of the species fluxes divided by their respective mobilities while
the latter simply uses the electrostatic force density acting on all
species. Note that this switching is only possible for the linkcentered
stencil. For all other stencils, this choice is hardcoded. The default
is *friction*.

enables the action of the electrostatic potential due to the
electrokinetics species and charged boundaries on the MD particles. The
forces on the particles are calculated by interpolation from the
electric field which is in turn calculated from the potential via finite
differences. This only includes interactions between the species and
boundaries and MD particles, not between MD particles and MD particles.
To get complete electrostatic interactions a particles Coulomb method
like Ewald or P3M has to be activate too.

.. _Diffusive Species:

Diffusive Species
^^^^^^^^^^^^^^^^^

The command followed by an integer (in the range 0 to 10) and several
options can be used to initialize the diffusive species. Here the
options specify: the number density , the diffusion coefficient , the
valency of the particles of that species , and an optional external
(electric) force which is applied to the diffusive species. As mentioned
before, the LB density is completely decoupled from the electrokinetic
densities. This has the advantage that greater freedom can be achieved
in matching the internal parameters to an experimental system. Moreover,
it is possible to choose parameters for which the LB is more stable. The
LB fluid must already be (partially) set up using the ... command,
before the diffusive species can be initialized. The variables , , and
must be set to properly initialize the diffusive species; the is
optional.

.. _Boundaries:

Boundaries
^^^^^^^^^^

The command allows one to set up (internal or external) boundaries for
the electrokinetics algorithm in much the same way as the command is
used for the LB fluid. The major difference with the LB command is given
by the option , with which a boundary can be endowed with a volume
charge density. To create a surface charge density, a combination of two
oppositely charged boundaries, one inside the other, can be used.
However, care should be taken to maintain the surface charge density
when the value of is changed. Currently, the following s are available:
wall, sphere, cylinder, rhomboid, pore, stomatocyte, hollow\_cone, and
spherocylinder. We refer to the documentation of the command
(Chapter [sec:lb]) for information on the options associated to these
shapes. In order to properly set up the boundaries, the and relevant
must be specified.

.. _Output:

Output
~~~~~~

.. _Fields:

Fields
^^^^^^

print

The print parameter of the command enables simple visualization of
simulation data. A property of the fluid field can be exported into a
file with name in one go. Currently, supported values of the parameter
are: , , , and , which give the LB fluid density, the LB fluid velocity,
the electrostatic potential, and the location and type of the
boundaries, respectively. The boundaries can only be printed when the
``EK_BOUNDARIES`` is compiled in. The additional option can be used to
directly export in the vtk format. The vtk format is readable by
visualization software such as paraview [1]_ and mayavi2 [2]_. If the
option is not specified, a gnuplot readable data file will be exported.

print

This print statement is similar to the above command. It enables the
export of diffusive species properties, namely: and , which specify the
number density and flux of species , respectively.

.. _Local Quantities:

Local Quantities
^^^^^^^^^^^^^^^^

node velocity

The option of the command allows one to output the value of a quantity
on a single LB node. The node is addressed using three integer values
which run from 0 to /, /, and /, respectively. Thus far, only the
velocity of the LB fluid can be printed in the standard electrokinetics
implementation. For other quantities the command may be used.

node density

This command can be used to output the number density of the -th
diffusive species on a single LB node.

.. _Checkpointing:

Checkpointing
~~~~~~~~~~~~~

checkpoint save checkpoint load

Variant writes the species density fields as well as all necessary LB
fields into two files called and . Variant reads these files back and
restores the fields. All algorithm parameters must be set via the
simulation script, as they are not part of the checkpointed data.

The format of the checkpoint is binary and no special care is taken with
respect to the specific binary layout of the machine.

.. _Catalytic Reaction theory:

Catalytic Reaction theory
~~~~~~~~~~~~~~~~~~~~~~~~~

.. _Concept:

Concept
^^^^^^^

The electrokinetics solver implemented in can be used to simulate a
system, for which in addition to the electrokinetic equations, there is
a (local) catalytic reaction which converts one species into another.

If you are interested in using this implementation in for scientific
purposes, please contact J. de Graaf before you start your project.

Currently, a linear reaction is implemented which converts one species
into two others, in order to model the catalytic decomposition of
hydrogen peroxide in the presence of a platinum catalyst:
:math:`2 \mathrm{H}_{2}\mathrm{O}_{2} \rightarrow 
2 \mathrm{H}_{2}\mathrm{O} + \mathrm{O}_{2}`. The decomposition of
:math:`\mathrm{H}_{2}\mathrm{O}_{2}` is in reality more complicated than
the linear reaction introduced here, since it is assumed to proceed via
several intermediate complexed-states, but our model can be thought of
as modeling the rate-limiting step. If we assume that there are three
non-ionic species with number densities :math:`n_{k}`, where
:math:`n_{0} = [ \mathrm{H}_{2}\mathrm{O}_{2} ]`,
:math:`n_{1} = [ \mathrm{H}_{2}\mathrm{O} ]`, and
:math:`n_{2} = [ \mathrm{O}_{2} ]`, then we can write the
(electro)kinetic equations for this system as

.. math::

   \begin{aligned}
   \label{eq:ek-reaction-continuity} \frac{\partial n_k}{\partial t} & = & -\, \nabla \cdot \vec{j}_k +\, f_{k} c n_{k} \vphantom{\left(\frac{\partial}{\partial}\right)} ; \\
   \label{eq:ek-reaction-fluxes} \vec{j}_{k} & = & -D_k \nabla n_k + n_k \vec{v}_{\mathrm{fl}} \vphantom{\left(\frac{\partial}{\partial}\right)} ; \\
   \nonumber \left(\frac{\partial \vec{v}_{\mathrm{fl}}}{\partial t} + \vec{v}_{\mathrm{fl}} \cdot \vec{\nabla} \vec{v}_{\mathrm{fl}} \right) \rho_\mathrm{fl} & = & -{k_\mathrm{B}T}\, \sum_{k} \nabla n_k   \\
   \label{eq:ek-reaction-velocity} & & +\, \eta \vec{\Delta} \vec{v}_{\mathrm{fl}} + (\eta / 3 + \eta_{\text{\,b}}) \nabla (\nabla \cdot \vec{v}_{\mathrm{fl}}) \vphantom{\left(\frac{\partial}{\partial}\right)} ; \\
   \label{eq:ek-reaction-continuity-fl} \frac{\partial \rho_\mathrm{fl}}{\partial t} & = & -\,\nabla\cdot\left( \rho_\mathrm{fl} \vec{v}_{\mathrm{fl}} \right) \vphantom{\left(\frac{\partial}{\partial}\right)} ,\end{aligned}

which define relations between the following observables

:math:`n_k`
    the number density of the particles of species :math:`k`,

:math:`\vec{j}_k`
    the number density flux of the particles of species :math:`k`,

:math:`\rho_{\mathrm{fl}}`
    the mass density of the fluid,

:math:`\vec{v}_{\mathrm{fl}}`
    the advective velocity of the fluid,

and input parameters

:math:`D_k`
    the diffusion constant of species :math:`k`,

:math:`{k_\mathrm{B}T}`
    | the thermal energy given by the product of Boltzmann’s constant
      :math:`k_\text{B}`
    | and the temperature :math:`T`,

:math:`\eta`
    the dynamic viscosity of the fluid,

:math:`\eta_{\text{\,b}}`
    the bulk viscosity of the fluid,

:math:`f_{k}`
    the reaction constant :math:`f_{0} \equiv -1`, :math:`f_{1} = 1` and
    :math:`f_{2} = 0.5` for the above reaction,

:math:`c`
    the reaction rate.

In this set of equations we have fully decoupled the number densities
and the fluid mass density. N.B. We have set the initial fluid mass
density is not necessarily equal to the sum of the initial species
number densities. This means that some care needs to be taken in the
interpretation of the results obtained using this feature. In
particular, the solution of the Navier-Stokes equation exclusively
models the momentum transport through the (multicomponent) fluid, while
the diffusive properties of the individual chemical species are handled
by Eqs.  and .

It is important to note that to ensure mass conservation the reaction
must satisfy:

.. math:: \label{eq:ek-mass-balance} \sum_{k} f_{k} m_{k} = 0 ,

where :math:`m_{k}` is the molecular mass of a reactive species.
Unfortunately, the current electrokinetic implementation does not
conserve mass flux locally. That is to say, the LB fluid is compressible
and the sum of the fluxes of the three species is not equal to zero in
the frame co-moving with the advective fluid velocity. It is therefore
debatable whether it is necessary to impose Eq. , since the EK algorithm
itself does not conserve mass density. However, we strived to be as
accurate as possible and in future versions of the EK algorithm the lack
of incompressiblity will be addressed.

The reaction is specified by the second term on the right-hand side of
Eq. . It is important to note that this term can be set locally, as
opposed to the other terms in the equation system Eqs. -, in our
implementation, as will become clear in the following. This has the
advantage that catalytic surfaces may be modeled.

.. _Initialization and Geometry Definition:

Initialization and Geometry Definition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The command is used to set up the catalytic reaction between three
previously defined the diffusive species, of which the i identifiers are
given by , , and , respectively. In the 1:2 reaction, these fulfill the
role of the reactant and the two products, as indicated by the naming
convention. For each species a reservoir (number) density must be set,
given by the variables , , and , respectively. These reservoir densities
correspond to the initial number densities associated with the reactive
species. The reservoir densities, in tandem with reservoir nodes, see
below, can be used to keep the reaction from depleting all the reactant
in the simulation box. The variable specifies the speed at which the
reaction proceeds. The three masses (typically given in the atomic
weight equivalent) are used to determine the total mass flux provided by
the reaction, as described above, and are also used to check whether the
reaction ratios that are given satisfy the chemical requirement of mass
conservation. Finally, the parameters and specify what fractions of the
product are generated when a given quantity of reactant is catalytically
converted. To use a chemical reaction, all options for the command must
be specified.

The option of the command allows one to set up regions in which the
reaction takes place with the help of the constraints that are available
to set up boundaries. The integer value can be used to select the
reaction: 0 no reaction takes place for this region, 1 the catalytic
reaction takes place in this region, and 2 the region functions as a
reservoir, wherein the species densities are reset to their initial (or
reservoir) concentrations. The rest of the command follows the same
format of the command. Currently, the following s are available: box,
wall, sphere, cylinder, rhomboid, pore, stomatocyte, hollow\_cone, and
spherocylinder. The box shape is a specific command, which can be used
to set the entire simulation box to a specific reaction value. To use
the command, one must first set up a reaction, as described above. To
successfully specify a region all the relevant arguments that go with
the shape constraints must be provided.

.. [1]
   http://www.paraview.org/

.. [2]
   http://code.enthought.com/projects/mayavi/
