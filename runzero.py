from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import math
from Mesh import Mesh
from Boundary import PolygonalBoundary
from numpy import linalg
from FiniteElement import FiniteElement

"""
This runscript test the FiniteElement module via an manufactured solution.
to change the number of Points on the domain n can be changed.
"""
#exact solution
sol = lambda x : np.exp(x[0]+0.2*x[1])

#right handside calculated to match exact solution
def f(x):
    result = -1.04* np.exp(x[0]+0.2*x[1])
    return result

points = np.array([[0,0],[1,0],[1,1],[0,1]])
triangles = [[0,1,2],[2,3,0]]
boundaryEdges = [[[0,1],0],[[1,2],1],[[2,0],None],[[2,3],2],[[3,0],3]]
polBoundary= PolygonalBoundary(points, [sol]*np.shape(points)[0])
print(polBoundary.functionAtEdges)

m = Mesh(points, triangles, boundaryEdges,polBoundary)
error = []
for i in range(2,5):
    print("STEP: "+str(i))
    m.refineMesh(1)
    # m.plotTriangles()
    #Initiate Finite element class, with the PDE
    Fem = FiniteElement(m,PDEMatrix= np.array([[1,0],[0,1]]),functionRHS = f)

    #Calculate the Stiffnessmatrix
    print("Calculate Global Stiffness")
    Fem.calculateGlobalStiffnessMatrix()

    print("Calculate RHS:")
    Fem.calculateRightHandSide()
    print("solve")
    Fem.solve()

    #Calculates the error in the infty norm 
    # error = np.linalg.norm(np.array([sol(x) for x in Fem.triangulation.points])-np.array(Fem.solution),np.infty)
    error.append([Fem.getL2Error(sol),Fem.mesh.diam])

# print(error)

#Plots the Fem solution
fig = plt.figure()
ax = Axes3D(fig)
ax.text(0.5,0.5,1,"FEM - solution",color="red")
ax.plot_trisurf(Fem.triangulation.points[:,0],Fem.triangulation.points[:,1],Fem.triangulation.triangles.copy(),Fem.solution)

#Plots the exact solution
fig1 = plt.figure()
ax1 = Axes3D(fig1)
ax1.text(0.5,0.5,1,"exact - solution",color="red")
ax1.plot_trisurf(Fem.triangulation.points[:,0],Fem.triangulation.points[:,1],Fem.triangulation.triangles.copy(),[sol(x) for x in Fem.triangulation.points])

#plt.loglog([e[0] for e in error], [e[1] for e in error],'o')
for index, e in enumerate(error):
    if index < len(error)-1:
        print(np.log(error[index+1][0]/e[0])/np.log(error[index+1][1]/e[1]))
plt.show()
