import matplotlib.pyplot as plt
import operator
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import math
from Mesh import Mesh
from numpy import linalg
from FiniteElement import FiniteElement
from Boundary import PolygonalBoundary
from ProjectionOnR import ProjectionOnR

"""
This runscript test the FiniteElement module via an manufactured solution.
to change the number of Points on the domain n can be changed.
"""
lam = 2/3.0
#exact solution
sol = lambda x : 0 
# sol = lambda x : np.exp(x[0]+0.2*x[1])
def rhSinLam(x):
        theta =  np.arctan2(x[1],x[0])
        
        #negative Angles should be converted to positive Values to match the angleCOoefficient
        if theta <0:
            theta = 2*np.pi+theta
        r = np.linalg.norm(x)
        if r<= 1e-15:
            return 0
        else:
            return r**(-lam) * np.sin(lam*theta)

#right handside calculated to match exact solution
def u(x):
    return (x[0]**3-x[0])*(x[1]**3-x[1])

def q(x):
    return -6*(x[0]*(x[1]**3-x[1]) + x[1]*(x[0]**3-x[0]))

points = np.array([[0,0],[1,0],[1,1],[0,1],[-1,1],[-1,0],[-1,-1],[0,-1],[-0.5,0.5]])
polBoundary= PolygonalBoundary(points, [sol]*np.shape(points)[0])
triangles = [[0,1,2],[0,2,3],[0,3,8],[0,8,5],[8,3,4],[8,4,5],[0,5,6],[0,6,7]] 
#the boundary edges holds all the edges of each traingle, as a  pair of an edge and 
#the index of the polygonalBoundary it belongs to.
boundaryEdges = [[[0,1],0],[[1,2],1],[[2,0],None],[[2,3],2],[[3,0],None],[[3,8],None],[[8,0],None],[[8,5],None],[[3,4],2],[[8,4],None],[[4,5],3],[[5,0],None],[[5,6],3],[[6,0],None],[[6,7],4],[[7,0],5]]
m = Mesh(points, triangles, boundaryEdges, polBoundary)
error = []
m.refineMesh(5)
# m.plotTriangles()
projR= ProjectionOnR(angleCoefficient = lam,mesh= m,indexOfNonConvexCorner = 0,functionValuesToProject = [rhSinLam(x) for x in m.points] )

projR.functionValuesToProject = projR.p_h_tilde
projR.degreeOfGauss = 6

projR.calculateR_h()
projR.calculatePStar()
projR.calculatePTilde()
projR.functionValuesToProject = projR.p_h_tilde

maxIndex,maxValue = max(enumerate(projR.p_h_tilde.tolist()),key=operator.itemgetter(1))
print(maxIndex,maxValue)
print(projR.mesh.points[maxIndex])
projR.calculateNormPTildeSquared()
proR = projR.getProjectionOnR()
#Plots the Fem solution

fig = plt.figure()
ax = Axes3D(fig)
ax.text(0.5,0.5,1,"control",color="red")
ax.plot_trisurf(m.points[:,0],m.points[:,1],m.triangles.copy(),proR)

# fig1 = plt.figure()
# ax1 = Axes3D(fig1)
# ax1.text(0.5,0.5,1,"state",color="red")
# ax1.plot_trisurf(m.points[:,0],m.points[:,1],m.triangles.copy(),[u(x) for x in m.points])
# plt.plot([x[0] for x in als],[x[1] for x in als])
#Plots the exact solution
# fig1 = plt.figure()
# ax1 = Axes3D(fig1)
# ax1.text(0.5,0.5,1,"exact - solution",color="red")
# ax1.plot_trisurf(Fem.triangulation.points[:,0],Fem.triangulation.points[:,1],Fem.triangulation.triangles.copy(),[sol(x) for x in Fem.triangulation.points])

# plt.loglog([e[0] for e in error], [e[1] for e in error],'o')
plt.show()
