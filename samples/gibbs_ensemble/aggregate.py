import numpy as np
import sys
import pickle
import scipy.optimize

warmup_length = 55000

low_densities = {}
high_densities = {}
for f in sys.argv[1:]:
    data = pickle.load(open(f, 'rb'))
    kT = data['args'].temperature
    if not kT in low_densities: low_densities[kT] = []
    if not kT in high_densities: high_densities[kT] = []

    densities = np.array(data['densities'])
    indices = data['list_indices']
    start = np.where(np.array(indices) > warmup_length)[0][0]

    low_densities[kT].append(np.mean(np.amin(densities[:, start:], axis=0)))
    high_densities[kT].append(np.mean(np.amax(densities[:, start:], axis=0)))


for kT in low_densities:
    print(" %.3g %.4g +/- %.4g, %.4g +/- %.4g" % (kT, 
                                                  np.mean(
                                                      low_densities[kT]), np.std(
                                                      low_densities[kT]),
                                                  np.mean(high_densities[kT]), np.std(high_densities[kT])))
kTs = sorted(list(low_densities.keys()))
delta_rhos = [np.mean(high_densities[kT]) -
              np.mean(low_densities[kT]) for kT in kTs]

print(delta_rhos)

beta = 0.32


def scaling_law(kT, B, kT_c): return B * np.abs(kT - kT_c)**beta


fit, p_cov = scipy.optimize.curve_fit(
    scaling_law, kTs, delta_rhos, p0=(
        1 / 6, 2), maxfev=10000)
print(fit, np.sqrt(np.diag(p_cov)))
print(np.vstack((kTs, scaling_law(kTs, *fit), delta_rhos)).T)
