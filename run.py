import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import math
from numpy import linalg
from FiniteElement import FiniteElement

"""
This runscript test the FiniteElement module via an manufactured solution.
to change the number of Points on the domain n can be changed.
"""
#number of points to use in x and y direction
n = 30

#exact solution
sol = lambda x : np.exp(x[0]+0.2*x[1])

#right handside calculated to match exact solution
def f(x):
    result = -1.04* np.exp(x[0]+0.2*x[1])
    return result

#boundary values  have the format list with elements [pointindex, value_at_this_point] 
boundary=[]

points = []
for x in np.linspace(0.5,1,n):
    for y in np.linspace(0,1,n):
        points.append([x,y])

#maximum value of x and y direction for the rectangular domain
maximum0 = max(np.array(points)[:,0])
maximum1 = max(np.array(points)[:,1])

#minimum value of x and y direction for the rectangular domain
minimum0 = min(np.array(points)[:,0])
minimum1 = min(np.array(points)[:,1])

#use the exact solution to add dirichlet boundary
boundaryValuesx = [[k,sol(x)] for k,x in enumerate(points) if x[0]==maximum0 ]
boundaryValuesy = [[k,sol(x)] for k,x in enumerate(points) if x[1]==maximum1 and not x[0]==maximum0 and not x[0]==minimum0 ]
boundaryValuesx += [[k,sol(x)] for k,x in enumerate(points) if x[0]==minimum0 ]
boundaryValuesy += [[k,sol(x)] for k,x in enumerate(points) if x[1]==minimum1 and not x[0]==maximum0 and not x[0]==minimum0]

boundary = np.array(boundaryValuesx + boundaryValuesy)

#Initiate Finite element class, with the PDE
Fem = FiniteElement(points,boundary,PDEMatrix= np.array([[1,0],[0,1]]),functionRHS = f)

#Calculate the Stiffnessmatrix
Fem.calculateGlobalStiffnessMatrix()
Fem.calculateRightHandSide()
Fem.solve()

#Calculates the error in the infty norm 
error = np.linalg.norm(np.array([sol(x) for x in Fem.triangulation.points])-np.array(Fem.solution),np.infty)

#Plots the Fem solution
fig = plt.figure()
ax = Axes3D(fig)
ax.text(0.5,0.5,1,"exact solution",color="red")
ax.plot_trisurf(Fem.triangulation.points[:,0],Fem.triangulation.points[:,1],Fem.triangulation.simplices.copy(),Fem.solution)

#Plots the exact solution
fig1 = plt.figure()
ax1 = Axes3D(fig1)
ax1.text(0.5,0.5,1,"FEM - solution",color="red")
ax1.plot_trisurf(Fem.triangulation.points[:,0],Fem.triangulation.points[:,1],Fem.triangulation.simplices.copy(),[sol(x) for x in Fem.triangulation.points])

print(error)
plt.show()
