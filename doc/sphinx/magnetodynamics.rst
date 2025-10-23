.. _Magnetodynamics:

_Magnetodynamics
==============

|es| contains methods to simulate the internal magnetisation dynamics of magnetic particles. the interactions between point dipoles. The defaout behaviour in |es| is that, once a dipole moment is assigned to a particle, it will strictly follow the rotation of the particles quaternion. This is called the "fixed point dipole" model. This model is the predominantly used one to simulate magnetic soft matter. However, there are implications to this approximation, which does not represent a realistin behavouir of a single domain magnetic nanopartilcle, in all cases expect for inifinitely high magnetic naisotropy energy. In realisy, there are several internal relaxation mechanism the dipole moment experiences and hence is not necessatily coaligned with the particle quaternion, or follows its motion. In order to incorporate this phenomenology in simualtions, |es| offers several models, als listed bellow.

.. Thermal Stoner-Wohlfarth:

Thermal Stoner-Wohlfarth
-------------------

The thermal Stoner–Wohlfarth (SW) model includes Néel relaxation in simulations of single-domain magnetic nanoparticles. This is an implementation of the algoritm as originaly presented in :cite:`mostarac2025thermal`.

The magnetic energy of a Stoner-Wohlfarth particle (a T=0) is given by:
```math
 U =  - \mu_0 \mu(\vec{e}\cdot\vec{H}) - KV (\vec{e} \cdot \vec{n})^2
```

Here we outline our **tSW algorithm**, specifically how we simulate the internal dynamics and the Néel relaxation mechanism in magnetic colloids.  
For more details, please referer to the original publication.

The algorithm requires the following input parameters: saturation magnetisation **μ**, reduced field **h**, and the anisotropy parameter **σ**.  
It can be logically separated into three main steps:

1. **Finding the extrema** of the magnetic energy for the current state of the particle.  
2. **Calculating the energy barrier** to estimate the transition probability between possible dipole moment orientation states due to thermal fluctuations.  
3. **Updating the dipole moment orientation** based on a trial move against the transition probability.

---

### Step 1 — Critical Field and Energy Minima

First, the **critical field** \( h_{cr} \) is calculated, based on the current angle between the magnetic field vector **h** and the particle’s anisotropy axis at its position, **φ**.

For a given **φ**, the algorithm finds the **θ** that minimizes the total magnetic energy, denoted as  
\( 	heta'_{	ext{min}} \), which is closest to the previous dipole moment state.  
This is ensured by initializing the state of the energy minimizer with the previous particle state.

If the field acting on the particle is less than \( h_{cr} \), the algorithm proceeds to find a **θ** that **maximizes** the total magnetic energy on both sides of \( 	heta'_{	ext{min}} \), denoted as  
\( 	heta'_{	ext{max}} \) and \( 	heta''_{	ext{max}} \).

---

### Step 2 — Energy Barrier and Néel Relaxation

Using these extrema, the algorithm calculates the **energy barriers** on both sides of \( 	heta'_{	ext{min}} \):

```math
\begin{aligned}
\Delta E' &= \frac{1}{KV} \left| U(\phi, \theta'_{\text{max}}) - U(\phi, \theta'_{\text{min}}) \right|, \\
\Delta E'' &= \frac{1}{KV} \left| U(\phi, \theta''_{\text{max}}) - U(\phi, \theta'_{\text{min}}) \right|.
\end{aligned}
```

The smaller of the two energy barriers is chosen:

```math
\Delta E = \min(\Delta E', \Delta E'').
```

This barrier is then used to estimate a **characteristic timescale** for the Néel relaxation process, using:

```math
\tau_N = \frac{\tau_D}{2\sigma} \sqrt{\frac{\pi}{\sigma}} e^{\Delta E \sigma}.
```

From this, the **transition probability** (neglecting back-switching) is obtained as:

```math
p = 1 - e^{-\delta t / \tau_N},
```

where \( \delta t \) is the integration time step.  

---

### Step 3 — Trial Move and Dipole Update

A **trial move** is made by casting a random number and comparing it with the transition probability, similar to a **Metropolis Monte Carlo step**.  

- If the trial move is **successful**, the algorithm finds a new minimum \( \theta''_{\text{min}} \) and aligns the dipole moment accordingly.  
- Otherwise, the dipole moment remains aligned with \( \theta'_{\text{min}} \).

To use this implementation of magnetodynamics, activate the feature
``THERMAL_STONER_WOHLFARTH``.

In interacting systems, the method relies on DIPOLE_FIELD_TRACKING feature. Make sure you use the method with magnetostatics actors taht support this feature.
