import numpy as np
import sys
import pickle
import scipy.optimize
import matplotlib.pyplot as plt
warmup_length = 100000

delta_rhos = {}
for f in sys.argv[1:]:
    data = pickle.load(open(f, 'rb'))
    kT = data['args'].temperature

    densities = np.array(data['densities'])
    indices = data['list_indices']
    start = np.where(np.array(indices) > warmup_length)[0][0]

    low_dens = np.amin(densities[:, start:], axis=0)
    high_dens = np.amax(densities[:, start:], axis=0)
    print(kT, np.mean(low_dens), np.mean(high_dens))
    delta_rho = high_dens-low_dens
    if kT not in delta_rhos:
        delta_rhos[kT] =[]
    delta_rhos[kT].append(np.mean(delta_rho))

x = []
y = []
for kT in sorted(delta_rhos):
    print(kT, np.mean(delta_rhos[kT]), np.std(delta_rhos[kT]))
    x.append(kT)
    y.append(np.mean(delta_rhos[kT]))



beta = 0.32


def scaling_law(kT, B, kT_c): return B * np.abs(kT - kT_c)**beta


fit, p_cov = scipy.optimize.curve_fit(
    scaling_law, x, y, p0=(
        1 / 6, 2), maxfev=10000)
print(fit, np.sqrt(np.diag(p_cov)))
print(np.vstack((x, scaling_law(x, *fit), y)).T)
ms=30
lw=20
plt.plot(x, y, '.k', markersize=ms)
plt.plot(x,scaling_law(x,*fit), '-k', linewidth=lw)
plt.show()
