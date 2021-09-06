#
# Copyright (C) 2021 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""
Take in all simulation data, compute equilibrium values for every
datafile and take means for simulations with same temperature.
Plot the phase diagram and the fitted critical point.
"""

import espressomd.analyze
import numpy as np
import scipy.optimize
import pickle
import gzip
import matplotlib.pyplot as plt
import argparse
import logging


parser = argparse.ArgumentParser(
    description="""Computes mean values for every file given as input,
    computes the mean values for same temperatures and different seeds,
    tries to compute the critical point using rectilinear and scaling
    law and plots the phase diagram.""")
parser.add_argument(
    'files',
    nargs='+',
    help="(All) files of the simulated system")
args = parser.parse_args()

# Set the logging level
logging.basicConfig(
    format='%(asctime)s %(levelname)s: %(message)s',
    datefmt='%H:%M:%S',
    level=logging.INFO
)

# Warmup length to skip at beginning
WARMUP_LENGTH = 3000

# Lists for (unsorted) densities and temperature
dens_liquid = []
dens_gas = []
errors_liquid = []
errors_gas = []
kTs = []

# 3D critical exponent (See Frenkel, Smit: Understanding Molecular
# Simulation 2002, p.217)
BETA = 0.32


def autocorr(x):
    """ Compute the normalized autocorrelation function """
    _x = x - np.mean(x)
    acf = espressomd.analyze.autocorrelation(_x) * np.arange(len(_x), 0, -1)
    return acf / acf[0]


def relaxation_time(x):
    """ Compute the relaxation time """
    corr = autocorr(x)
    popt, _ = scipy.optimize.curve_fit(lambda x, a: np.exp(-a * x),
                                       range(np.size(corr)), corr, p0=(1. / 15000))
    return popt[0]


def std_error_mean(x):
    """ Compute the standard error of the mean with correlated samples """
    autocorr_time = relaxation_time(x)
    N_eff = np.size(x) / (2 * autocorr_time)
    return np.sqrt(np.var(x) / N_eff)


def scaling_law(kT, B, kT_c):
    """ Scaling law, see Frenkel, Smit, eq. (8.3.6) """
    return B * np.abs(kT - kT_c)**BETA


def rectilinear_law(kT, p_c, A, kT_c):
    """ Rectilinear law, see Frenkel, Smit, eq. (8.3.7) """
    return p_c + A * np.abs(kT - kT_c)


# Load data
for f in args.files:

    logging.info(f"Loading {f}...")
    with gzip.open(f, 'rb') as _f:
        data = pickle.load(_f)

    kT = data["temperature"]

    # Calculate densities of both phases
    key1, key2 = ("densities_box01", "densities_box02")
    density1 = np.mean(data[key1][WARMUP_LENGTH:])
    density2 = np.mean(data[key2][WARMUP_LENGTH:])

    dens_liquid.append(np.max((density1, density2)))
    dens_gas.append(np.min((density1, density2)))

    # Calculate errors
    if density1 < density2:
        key1, key2 = (key2, key1)
    errors_liquid.append(std_error_mean(data[key2]))
    errors_gas.append(std_error_mean(data[key1]))
    kTs.append(kT)

# Sort temperatures and densities respectively
order = np.argsort(kTs)
kTs = np.array(kTs)[order]
dens_liquid = np.array(dens_liquid)[order]
dens_gas = np.array(dens_gas)[order]
errors_liquid = np.array(errors_liquid)[order]
errors_gas = np.array(errors_gas)[order]

# Average over simulation results with different seeds
dens_liquid = [np.mean([dens_liquid[i] for (i, T) in enumerate(kTs) if T == kT])
               for kT in np.unique(kTs)]
dens_gas = [np.mean([dens_gas[i] for (i, T) in enumerate(kTs) if T == kT])
            for kT in np.unique(kTs)]

# Compute Errors as maximum of each simulation error
errors_liquid = [np.max([errors_liquid[i] for (i, T) in enumerate(kTs)
                         if T == kT]) for kT in np.unique(kTs)]
errors_gas = [np.max([errors_gas[i] for (i, T) in enumerate(kTs) if T == kT])
              for kT in np.unique(kTs)]
kTs = np.unique(kTs)

# Convert to arrays to do math
dens_liquid = np.array(dens_liquid)
dens_gas = np.array(dens_gas)

# Needed for scaling law and rectilinear law
y_scaling = dens_liquid - dens_gas
y_rectilinear = 0.5 * (dens_liquid + dens_gas)

# Fit using scaling law
fit_scaling, p_cov_scaling = scipy.optimize.curve_fit(
    scaling_law, kTs, y_scaling, p0=(1., 1.), maxfev=6000)

# Print fit values
logging.info(f"Fits using scaling law: B, T_c: {fit_scaling}")
logging.info(f"Fit uncertainty: {np.diag(p_cov_scaling)}")

# Critical temperature obtained via scaling law
kT_c_scaling = fit_scaling[1]

# Fit using rectilinear law
fit_rectilinear, p_cov_rectilinear = scipy.optimize.curve_fit(
    rectilinear_law, kTs, y_rectilinear, p0=(0.3, 0.2, kT_c_scaling),
    maxfev=6000)

# Print fit values
logging.info(f"Fits using rectilinear law: p_c, A, T_c: {fit_rectilinear}")
logging.info(f"Fit uncertainty: {np.diag(p_cov_rectilinear)}")

# Critical point obtained via rectilinear law
p_c = fit_rectilinear[0]
kT_c = fit_rectilinear[2]

# Save results
with gzip.open('result.dat.gz', 'wb') as f:
    pickle.dump(
        {"dens_liquid": dens_liquid,
         "dens_gas": dens_gas,
         "errors_liquid": errors_liquid,
         "errors_gas": errors_gas,
         "kTs": kTs,
         "kT_c_scaling": kT_c,
         "p_c": p_c},
        f)

# Plot results
ms = 40

plt.errorbar(dens_liquid,
             kTs,
             xerr=errors_liquid,
             fmt='o',
             c='red')

plt.errorbar(dens_gas,
             kTs,
             xerr=errors_gas,
             fmt='o',
             c='red',
             label='simulation data')

plt.scatter(p_c,
            kT_c,
            s=ms, c='black',
            marker='o',
            label='crit. point')

plt.scatter(0.5 * (dens_liquid + dens_gas),
            kTs,
            s=ms - 10, c='black',
            marker='x')

plt.plot(rectilinear_law(kTs,
                         p_c,
                         fit_rectilinear[1],
                         fit_rectilinear[2]),
         kTs)

plt.xlabel(r"Particle density $\rho*$")
plt.ylabel(r"Temperature $T*$")
plt.legend()

plt.savefig("phase_diagram.pdf")
plt.show()
