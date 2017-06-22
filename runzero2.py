import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
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
sol = lambda x : 0 
# sol = lambda x : np.exp(x[0]+0.2*x[1])

#right handside calculated to match exact solution
def f(x):
    if np.linalg.norm(x) >= 1e-10:
        result = np.linalg.norm(x)**(-lam) * np.sin(lam * np.arctan2(x[1],x[0]))
    else:
        result = 0
    return result

points = np.array([[0,0],[0,-1],[1,-1],[1,1],[-1,1],[-1,0]])
triangles = [[0,1,2],[0,2,3],[0,3,4],[0,4,5]]
boundaryEdges = [[[0,1],sol],[[1,2],sol],[[2,0],None],[[2,3],sol],[[3,0],None],[[3,4],sol],[[4,0],None],[[4,5],sol],[[5,0],sol]]
m = Mesh(points, triangles, boundaryEdges)
error = []
for i in range(0,5):
    m.refineMesh(1)
    # m.plotTriangles()
    #Initiate Finite element class, with the PDE
    Fem = FiniteElement(m,PDEMatrix= np.array([[1,0],[0,1]]),functionRHS = f)

    #Calculate the Stiffnessmatrix
    Fem.calculateGlobalStiffnessMatrix()
    Fem.calculateRightHandSide()
    Fem.solve()

    #Calculates the error in the infty norm 
    # error = np.linalg.norm(np.array([sol(x) for x in Fem.triangulation.points])-np.array(Fem.solution),np.infty)
    error.append([Fem.getL2Error(sol),Fem.mesh.diam])
print(error)

#Plots the Fem solution
# fig = plt.figure()
# ax = Axes3D(fig)
# ax.text(0.5,0.5,1,"FEM - solution",color="red")
# ax.plot_trisurf(Fem.triangulation.points[:,0],Fem.triangulation.points[:,1],Fem.triangulation.triangles.copy(),Fem.solution)

xs = []
ys = []
als = []
for index,x in enumerate(Fem.triangulation.points):
    if x[1] == 0:
        als.append([x[0],Fem.solution[index]])
als = sorted(als, key  = lambda x : x[0])
plt.plot([x[0] for x in als],[x[1] for x in als])
#Plots the exact solution
# fig1 = plt.figure()
# ax1 = Axes3D(fig1)
# ax1.text(0.5,0.5,1,"exact - solution",color="red")
# ax1.plot_trisurf(Fem.triangulation.points[:,0],Fem.triangulation.points[:,1],Fem.triangulation.triangles.copy(),[sol(x) for x in Fem.triangulation.points])

# plt.loglog([e[0] for e in error], [e[1] for e in error],'o')
plt.show()
