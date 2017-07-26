import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import math
from Mesh import Mesh
from numpy import linalg
from Boundary import PolygonalBoundary
from gaussIntegration import GaussIntegrator

"""
This runscript test the FiniteElement module via an manufactured solution.
to change the number of Points on the domain n can be changed.
"""
lam = 2/3.0
#exact solution
sol = lambda x : 0 
# sol = lambda x : np.exp(x[0]+0.2*x[1])

#right handside calculated to match exact solution
def f(x):
    return 1
points = np.array([[0,0],[1,0],[0,1]])
polBoundary= PolygonalBoundary(points, [sol]*np.shape(points)[0])
triangles = [[0,1,2]]
#the boundary edges holds all the edges of each traingle, as a  pair of an edge and 
#the index of the polygonalBoundary it belongs to.
boundaryEdges = [[[0,1],0],[[1,2],1],[[2,0],None]]
m = Mesh(points, triangles, boundaryEdges, polBoundary)
# m.refineMesh(1)
# m.plotTriangles()
GInteg = GaussIntegrator(m)
value = 0
for ele,triPoints in enumerate(m.triangles):
    g= GInteg.getFiniteElementFunctionOverTriangle([f(x) for x in m.points],ele)
    value  += GInteg.getIntegralOverTriangleGauss(g,ele,3)
print(value)

#Plots the Fem solution
# fig = plt.figure()
# ax = Axes3D(fig)
# ax.text(0.5,0.5,1,"singularPart",color="red")
# ax.plot_trisurf(m.points[:,0],m.points[:,1],m.triangles.copy(),projR.p_h_star)

# fig1 = plt.figure()
# ax1 = Axes3D(fig1)
# ax1.text(0.5,0.5,1,"singularPart",color="red")
# ax1.plot_trisurf(m.points[:,0],m.points[:,1],m.triangles.copy(),projR.r_h)
# plt.plot([x[0] for x in als],[x[1] for x in als])
#Plots the exact solution
# fig1 = plt.figure()
# ax1 = Axes3D(fig1)
# ax1.text(0.5,0.5,1,"exact - solution",color="red")
# ax1.plot_trisurf(Fem.triangulation.points[:,0],Fem.triangulation.points[:,1],Fem.triangulation.triangles.copy(),[sol(x) for x in Fem.triangulation.points])

# plt.loglog([e[0] for e in error], [e[1] for e in error],'o')
plt.show()
