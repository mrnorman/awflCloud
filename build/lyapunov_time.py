
import numpy as np
import sys
import matplotlib.pyplot as plt
from netCDF4 import Dataset

nc1 = Dataset(sys.argv[1],"r")
nt = nc1.dimensions["t"].size
nx = nc1.dimensions["x"].size
ny = nc1.dimensions["y"].size
nz = nc1.dimensions["z"].size

t = nc1.variables["t"]

nc2 = Dataset(sys.argv[2],"r")

theta1 = nc1.variables["pot_temp_pert"]
theta2 = nc2.variables["pot_temp_pert"]

l1 = [0. for i in range(nt)]
for i in range(nt) :
  l1[i] = np.sum( np.abs(theta2[i,:,:,:] - theta1[i,:,:,:]) )
  # print( i , l1[i] )

slope_l1 = ( np.log(l1[2400]) - np.log(l1[40]) ) / (t[2400] - t[40])

print( 1./slope_l1 )

slope_l1 = ( np.log(l1[2]) - np.log(l1[0]) ) / (t[2] - t[0])

print( 1./slope_l1 )

# plt.plot( t , l2 )
# plt.yscale('log')
# plt.show()

