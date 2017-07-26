import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from Boundary import PolygonalBoundary
import matplotlib.pyplot as plt
import numpy as np
import math
from Mesh import Mesh
from numpy import linalg
from FiniteElement import FiniteElement

"""
This runscript test the FiniteElement module via an manufactured solution.
to change the number of Points on the domain n can be changed.
"""
lam = 2/3.0
#exact solution
def sol(x):
    return (x[1]**3-x[1]) * (x[0]**3-x[0])

# sol = lambda x : np.exp(x[0]+0.2*x[1])

def control(x):
    return -6*(x[0]*(x[1]**3-x[1]) + x[1]*(x[0]**3-x[0]))

points = np.array([[0,0],[1,0],[1,1],[0,1],[-1,1],[-1,0],[-1,-1],[0,-1],[-0.5,0.5]])
polBoundary= PolygonalBoundary(points, [sol]*np.shape(points)[0])
triangles = [[0,1,2],[0,2,3],[0,3,8],[0,8,5],[8,3,4],[8,4,5],[0,5,6],[0,6,7]] 
#the boundary edges holds all the edges of each traingle, as a  pair of an edge and 
#the index of the polygonalBoundary it belongs to.
boundaryEdges = [[[0,1],0],[[1,2],1],[[2,0],None],[[2,3],2],[[3,0],None],[[3,8],None],[[8,0],None],[[8,5],None],[[3,4],2],[[8,4],None],[[4,5],3],[[5,0],None],[[5,6],3],[[6,0],None],[[6,7],4],[[7,0],5]]
m = Mesh(points, triangles, boundaryEdges,polBoundary)
error = []

for i in range(0,4):
    m.refineMesh(1)
    # m.plotTriangles()
    #Initiate Finite element class, with the PDE
    Fem = FiniteElement(m,PDEMatrix= np.array([[1,0],[0,1]]),functionRHS = control)

    #Calculate the Stiffnessmatrix
    Fem.calculateGlobalStiffnessMatrix()
    Fem.calculateRightHandSide()
    Fem.solve()

    #Calculates the error in the infty norm 
    # error = np.linalg.norm(np.array([sol(x) for x in Fem.triangulation.points])-np.array(Fem.solution),np.infty)
    error.append([Fem.getL2Error(sol),Fem.mesh.diam])
print(error)

#Plots the Fem solution
fig = plt.figure()
ax = Axes3D(fig)
ax.text(0.5,0.5,1,"FEM - solution",color="red")
ax.plot_trisurf(Fem.triangulation.points[:,0],Fem.triangulation.points[:,1],Fem.triangulation.triangles.copy(),Fem.solution)

xs = []
ys = []
als = []
for index,x in enumerate(Fem.triangulation.points):
    if x[1] == 0:
        als.append([x[0],Fem.solution[index]])
als = sorted(als, key  = lambda x : x[0])
# plt.plot([x[0] for x in als],[x[1] for x in als])
#Plots the exact solution
# fig1 = plt.figure()
# ax1 = Axes3D(fig1)
# ax1.text(0.5,0.5,1,"exact - solution",color="red")
# ax1.plot_trisurf(Fem.triangulation.points[:,0],Fem.triangulation.points[:,1],Fem.triangulation.triangles.copy(),[sol(x) for x in Fem.triangulation.points])

# plt.loglog([e[0] for e in error], [e[1] for e in error],'o')

for index, e in enumerate(error):
    if index < len(error)-1:
        print(np.log(error[index+1][0]/e[0])/np.log(error[index+1][1]/e[1]))
plt.show()
