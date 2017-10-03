import matplotlib.pyplot as plt
import os
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import json
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
        if np.linalg.norm(x)==0:
            return 0
        theta =  np.arctan2(x[1],x[0])
        if theta < 0:
            theta+=2*np.pi

        #negative Angles should be converted to positive Values to match the angleCOoefficient
        r = np.linalg.norm(x)
        return (r**(-lam) * np.sin(lam * theta))**2

points = np.array([[0,0],[1,0],[1,1],[0,1],[-1,1],[-1,0],[-1,-1],[0,-1]])
polBoundary= PolygonalBoundary(points, [sol]*np.shape(points)[0])
triangles = [[0,1,2],[0,2,3],[0,3,4],[0,4,5],[0,5,6],[0,6,7]] 
#the boundary edges holds all the edges of each traingle, as a  pair of an edge and 
#the index of the polygonalBoundary it belongs to.
boundaryEdges = [[[0,1],0],[[1,2],1],[[2,0],None],[[2,3],2],[[3,0],None],[[3,4],2],[[4,0],None],[[4,5],4],[[4,5],3],[[5,0],None],[[5,6],3],[[6,0],None],[[6,7],4],[[7,0],5]]
m = Mesh(points, triangles, boundaryEdges, polBoundary)
# m.refineMeshHalf(11)
# m.plotTriangles()
GInteg = GaussIntegrator(m)
value = 0

# for ele,triPoints in enumerate(m.triangles):
    # g= GInteg.getFiniteElementFunctionOverTriangle([f(x) for x in m.points],ele)
    # value  += GInteg.getIntegralOverTriangleGauss(g,ele,3)
value = GInteg.getIntegralOverDomain(f,20)
print(value)

for i in range(1,100):
    with open(os.path.join("LeggaussPoints"+str(i)),'w') as fileTriangles:
        json.dump(np.polynomial.legendre.leggauss(i)[0].tolist(),fileTriangles)
    with open(os.path.join("LeggaussWeights"+str(i)),'w') as fileTriangles:
        json.dump(np.polynomial.legendre.leggauss(i)[1].tolist(),fileTriangles)
 
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
# ax1.plot_trisurf(m.points[:,0],m.points[:,1],m.triangles.copy(),[f(x) for x in m.points])

# # plt.loglog([e[0] for e in error], [e[1] for e in error],'o')
# plt.show()
