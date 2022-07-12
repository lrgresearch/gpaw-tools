'''
find3Dmin.py: Draw 3D contour plot of E vs. latticeparams and show the minimum datapoint using
              the optimize_latticeparam.py's output.
              It is not fitting yet. But it can be implemented in the future.
Usage: $ python find3Dmin
'''

import sys
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# Load from  optimize_latticeparam.py's output
chunk = np.loadtxt('Optimize-Lattice_Table-LatticeParam.txt')
data=np.array(chunk)

Xs = data[:,0]
Ys = data[:,1]
Zs = data[:,3]

# Find minimum data point and colorize it
zmin = np.min(Zs)
mask = np.array(Zs) == zmin
color = np.where(mask, 'red', 'blue')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(Xs, Ys, Zs, s=20, color=color, plotnonfinite=True)
surf = ax.plot_trisurf(Xs, Ys, Zs, cmap=cm.jet, alpha=0.7, linewidth=0)
#surf = ax.plot_surface(Xs, Ys, Zs, cmap='coolwarm',linewidth=0, antialiased=True)
#surf = ax.contour3D(Xs, Ys, Zs, 50,  cmap='coolwarm',)
fig.colorbar(surf)

ax.xaxis.set_major_locator(MaxNLocator(5))
ax.yaxis.set_major_locator(MaxNLocator(6))
ax.zaxis.set_major_locator(MaxNLocator(5))

# Axes and title
ax.set_xlabel('Lattice Param a (Ang.)')
ax.set_ylabel('Lattice Param c (Ang.)')
ax.set_zlabel('Total Energy (eV)')
ax.set_title('Minimum Point is at X:'+str(Xs[np.where(Zs == Zs.min())])+", Y:"+str(Ys[np.where(Zs == Zs.min())]))

# Draw
fig.tight_layout()
plt.show()
