import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import math
from numpy import linalg
from FiniteElement import FiniteElement

nx =50
ny =50
dx = 1/(nx-1)
dy= 1/(ny-1)
points = []
for x in np.linspace(-1.5,1.5,nx):
    for y in np.linspace(-1.5,1.5,ny):
        points.append([x,y])
#points = [x  for x in points if linalg.norm(x)<=1]
boundaryValuesx = []#[([x,-1.5],0) for x in np.linspace(-1.5,1.5,nx)]
boundaryValuesy = []#[([0,y],0) for y in np.linspace(0,1,ny)]
boundaryRadial = [[k,0] for k,x in enumerate(points) if linalg.norm(x) >1]
boundaryRadial +=[]# [[k,0] for k,x in enumerate(points) if linalg.norm(x) <=0.3]
boundaryValuesx += []#[[x,1] for x in np.array(points)[:,-1]]
boundaryValuesy += []#[[x,0] for x in np.array(points)[-1,:]]

f= lambda x : x[0]

boundary = np.array(boundaryValuesx + boundaryValuesy+boundaryRadial)
Fem = FiniteElement(points,boundary,PDEMatrix= np.array([[1,0],[0,10]]),functionRHS = f)
Fem.calculateElementStiffnessMatrix(0)
Fem.calculateGlobalStiffnessMatrix()
Fem.calculateRightHandSide()
Fem.solve()
#plt.triplot(Fem.triangulation.points[:,0],Fem.triangulation.points[:,1],Fem.triangulation.simplices.copy())
#plt.plot(Fem.triangulation.points[:,0],Fem.triangulation.points[:,1],'o')
fig = plt.figure()
ax = Axes3D(fig)
ax.plot_trisurf(Fem.triangulation.points[:,0],Fem.triangulation.points[:,1],Fem.triangulation.simplices.copy(),Fem.solution)

#
#print(Fem.triangulation.points[:])
plt.show()

#rhs = lambda x,y: -0.5*x+-0.5*y+0.5 
#fig = plt.figure()
#ax = Axes3D(fig);
#xx,yy = np.meshgrid(np.arange(-1,1,0.1),np.arange(-1,1,0.1))
#ax.plot_surface(xx,yy,rhs(xx,yy))
#plt.show()

#A = np.array([[0,1,0],[0,0,1],[1,1,1]])
#B = np.array([[1,0,-1],[0,1,0],[1,1,1]])
#C=np.dot(B,np.linalg.inv(A))
#print(C)
#b= C[0:-1,0]
#A = C[0:-1,0:-1]
#print(np.dot(A,np.array([0,1]))+b)
