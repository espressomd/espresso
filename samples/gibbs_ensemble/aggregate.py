import numpy as np
import sys
import pickle

warmup_length = 10000


for f in sys.argv[1:]:
    data = pickle.load(open(f, 'rb'))
    densities = np.array(data['densities'])
    indices = data['list_indices']
    start = np.where(np.array(indices) > warmup_length)[0][0]

    print(data['args'].temperature, data['args'].seed, np.mean(np.amin(
        densities[:, start:], axis=0)), np.mean(np.amax(densities[:, start:], axis=0)))
