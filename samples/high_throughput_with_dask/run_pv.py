import sys
import dask.distributed
import logging
import numpy as np

import dask_espresso

logging.basicConfig(level=logging.WARN)


PYPRESSO = "/data/weeber/es/build/pypresso"
SIM_SCRIPT = "lj_pressure.py"

N_STEPS = int(2E4)
N_PARTICLES = 100
VOLUME_FRACTIONS = np.arange(0.1, 0.52, 0.01)


client = dask.distributed.Client(sys.argv[1])


futures = []


for volume_fraction in VOLUME_FRACTIONS:
    sim_params = {"volume_fraction": volume_fraction,
                  "n_particles": N_PARTICLES,
                  "n_steps": N_STEPS}
    futures.append(client.compute(dask_espresso.dask_espresso_task(
        PYPRESSO, SIM_SCRIPT, **sim_params)))

dask.distributed.progress(futures)

sim_results = client.gather(futures)

for vol_frac, sim_result in zip(VOLUME_FRACTIONS, sim_results): 
    print(vol_frac, sim_result["pressure"], sim_result["pressure_std_dev"])
