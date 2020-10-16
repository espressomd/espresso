import numpy as np
import sys
import pickle
import scipy.optimize
import matplotlib.pyplot as plt


# Number of densities to skip at the beginning of a simulation
WARMUP_LENGTH = 100000 

delta_rhos = {} # Holds density difference between dense and dilute phase vs kT

for f in sys.argv[1:]: # Iterate over result files
    data = pickle.load(open(f, 'rb'))
    
    kT = data['args'].temperature
    densities = np.array(data['densities'])
    indices = data['list_indices']
    
    
    # First index to consider
    start = np.where(np.array(indices) > WARMUP_LENGTH)[0][0]

    # Determin which is teh dilute phase from the last density recorded
    # This way the selection is permanent and the average densities of both
    # phases will be approximately the same if T > T_c
    dilute_phase = np.argmin(densities[:,-1])
    dense_phase = np.argmax(densities[:,-1])


    low_dens = densities[dilute_phase, start:]
    high_dens = densities[dense_phase, start:]
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
