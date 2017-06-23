import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import math
from Mesh import Mesh
from numpy import linalg
from ProjectedGradient import ProjectedGradient

sol = lambda x : 0

def f(x):
    result = 1  
    return result

u_d = lambda x : 0

points = np.array([[0,0],[1,0],[1,1],[0,1]])
triangles = [[0,1,2],[2,3,0]]
boundaryEdges = [[[0,1],sol],[[1,2],sol],[[2,0],None],[[2,3],sol],[[3,0],sol]]
m = Mesh(points, triangles, boundaryEdges)
error = []

m.refineMesh(4)
proj = ProjectedGradient(m,[f(x) for x in m.points], [u_d(x) for x in m.points])
proj.solve()
# m.plotTriangles()
#Initiate Finite element class, with the PDE

fig = plt.figure()
ax = Axes3D(fig)
ax.text(0,0.5,0,"Control",color="red")
ax.plot_trisurf(m.points[:,0],m.points[:,1],m.triangles.copy(),proj.control)

fig1 = plt.figure()
ax1 = Axes3D(fig1)
ax1.text(0,0.5,0,"State",color="red")
ax1.plot_trisurf(m.points[:,0],m.points[:,1],m.triangles.copy(),proj.state)
plt.show()
#plt.loglog([e[0] for e in error], [e[1] for e in error],'o')
