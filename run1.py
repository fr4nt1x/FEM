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
n = 10

#exact solution
sol = lambda x : 1+x[0]*x[0]+2*x[1]*x[1]

#right handside calculated to match exact solution
def f(x):
    result = -6
    return result

#boundary values  have the format list with elements [pointindex, value_at_this_point] 
boundary=[]

#specifies the bounding area of the domain
minimumX = 0
maximumX = 1
minimumY = 0
maximumY = 1

points = []
i=0
for x in np.linspace(minimumX,maximumX,n):
    for y in np.linspace(minimumY,maximumY,n):
        
        points.append([x,y])

        if x ==  minimumX or x == maximumX or y == minimumY or y == maximumY:
            boundary.append([int(i),sol([x,y])])
        i+=1

boundary= np.array(boundary)

#Initiate Finite element class, with the PDE
Fem = FiniteElement(points,boundary,PDEMatrix= np.array([[1,0],[0,1]]),functionRHS = f)

#Calculate the Stiffnessmatrix
Fem.calculateGlobalStiffnessMatrix()
Fem.calculateRightHandSide()
Fem.solve()

#Calculates the error in the infty norm 
error = np.linalg.norm(np.array([sol(x) for x in points])-np.array(Fem.solution),np.infty)


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
