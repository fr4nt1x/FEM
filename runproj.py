import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import math
from Mesh import Mesh
from numpy import linalg
from ProjectedGradient import ProjectedGradient

alpha = 1

def state(x): 
    return (x[0]-x[0]**2)*( x[1]-x[1]**2)

def control(x):
    return (x[0]-x[0]**2)*( x[1]-x[1]**2)

def RHSAddendum(x):
    return 2*((x[0]-x[0]**2)+ (x[1]-x[1]**2))- (x[0]-x[0]**2)*( x[1]-x[1]**2)

def f(x):
    result = 0  
    return result

u_d = lambda x : state(x) + 2/alpha*( (x[0]-x[0]**2)+ (x[1]-x[1]**2))

def boundrayFunc(x):
    return 0;

points = np.array([[0,0],[1,0],[1,1],[0,1]])
triangles = [[0,1,2],[2,3,0]]
boundaryEdges = [[[0,1],boundrayFunc],[[1,2],boundrayFunc],[[2,0],None],[[2,3],boundrayFunc],[[3,0],boundrayFunc]]
m = Mesh(points, triangles, boundaryEdges)
error = []
error2 = []

for i in range(0,4):

    m.refineMesh(1)
    proj = ProjectedGradient(m,[f(x) for x in m.points], [u_d(x) for x in m.points], [RHSAddendum(x) for x in m.points], alpha, tol = 1e-12)
    proj.solve()
    error.append([proj.getL2ErrorControlGauss(control),proj.mesh.diam])
    error2.append([proj.getL2ErrorControl(control),proj.mesh.diam])

# m.plotTriangles()

fig = plt.figure()
ax = Axes3D(fig)
ax.text(0,0.5,0,"Control",color="red")
ax.plot_trisurf(m.points[:,0],m.points[:,1],m.triangles.copy(),proj.control)

# fig1 = plt.figure()
# ax1 = Axes3D(fig1)
# ax1.text(0,0.5,0,"State",color="red")
# ax1.plot_trisurf(m.points[:,0],m.points[:,1],m.triangles.copy(),proj.state)

# fig2 = plt.figure()
# ax4 = Axes3D(fig2)
# ax4.text(0,0.5,0,"State real",color="red")
# ax4.plot_trisurf(m.points[:,0],m.points[:,1],m.triangles.copy(),[state(p) for p in m.points])
for index, e in enumerate(error):
    if index < len(error)-1:
        print(np.log(error[index+1][0]/e[0])/np.log(error[index+1][1]/e[1]))
for index, e in enumerate(error2):
    if index < len(error2)-1:
        print(np.log(error2[index+1][0]/e[0])/np.log(error2[index+1][1]/e[1]))
# plt.show()

#plt.loglog([e[0] for e in error], [e[1] for e in error],'o')
