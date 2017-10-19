#!/usr/bin/env python

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import sys

vec = np.genfromtxt('ttt.dat')
print vec
print len(vec)
vec = vec.reshape([len(vec), 19, 3])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

quiver = vec[int(sys.argv[1]),:,:]

X = np.zeros(len(quiver))

print(quiver)

for i in np.arange(len(quiver)):
  ax.quiver(X[i], X[i], X[i], quiver[i,0], quiver[i,1], quiver[i,2], pivot='tail', length=np.sqrt(np.sum(quiver[i]**2)))

ax.set_xlim3d(np.min(quiver[:,0]), np.max(quiver[:,0]))
ax.set_ylim3d(np.min(quiver[:,1]), np.max(quiver[:,1]))
ax.set_zlim3d(np.min(quiver[:,2]), np.max(quiver[:,2]))

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

plt.show()
