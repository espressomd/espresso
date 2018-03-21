import numpy as np
import sys
import matplotlib.pyplot as plt
import math
import sys
from mpl_toolkits.mplot3d import Axes3D


data = np.loadtxt('snapshot.dat')

#plt.ion()
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#for i in range(260):
  #plt.clf()
  #ax.scatter(data[i,1],data[i,2],data[i,3])


  #plt.pause(0.1)
  #plt.show()
#ax.scatter(data[:,1],data[:,2],data[:,3])
ax.scatter(data[:,0],data[:,1],data[:,2])
plt.show()
