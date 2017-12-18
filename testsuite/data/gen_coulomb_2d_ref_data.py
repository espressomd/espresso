import sys
from scipy import *
import numpy.linalg as linalg

box_l=10

data = genfromtxt("coulomb_mixed_periodicity_system.data")
n=data.shape[0]

pos=data[:,1:4]
q=data[:,4]
forces=zeros((n,3))
energy=0.

q1q2=outer(q,q)
images=2000

for i in range(n):
  for x in range(-images,images+1,1):
    for y in range(-images,images+1,1):
        if x**2+y**2 > images**2:
          continue
        pos_diff=pos[i]-(pos +array((x,y,0))*box_l)
        r=sqrt(sum(pos_diff**2,1))
        r3=r**3
        qq=q1q2[i,:]
        
        tmp=qq/r
        tmp[abs(tmp)==inf] =0
        energy+=sum(tmp)
        pref=qq/r**3
        pref=pref.reshape((n,1))
        tmp=pos_diff*hstack((pref,pref,pref))
        
        forces+=nan_to_num(tmp)


ids=arange(n)



forces*=-1
savetxt("coulomb_mixed_periodicity_system.data",hstack((ids.reshape((n,1)),pos,q.reshape((n,1)),forces)))
