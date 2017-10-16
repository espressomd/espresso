#!/usr/bin/env python

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

vec = np.genfromtxt('ttt.dat')
vec = vec.reshape([len(vec), 19, 3])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

X = np.zeros(19)

print vec[1]

ax.quiver(X, X, X, vec[1,:,0], vec[1,:,1], vec[1,:,2])

plt.show()
