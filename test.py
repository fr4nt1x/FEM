import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import math
from numpy import linalg
from FiniteElement import FiniteElement

sol = lambda x : np.exp(x[0]+0.2*x[1])
nx =50
ny =50


#f= lambda x : 1 if np.linalg.norm(x)<=0.1 else 0
def f(x):
    result = 0
    """
    if np.linalg.norm(x - np.array([-0.5,0])) <= 0.1:
        result = 1
    elif np.linalg.norm(x-np.array([0.5,0]))<= 0.1:
        result= 2
    elif np.linalg.norm(x - np.array([0,0.5])) <= 0.1:
        result = 1
    elif np.linalg.norm(x-np.array([0,-0.5]))<= 0.1:
        result= 2
    else:
        result=0
    """
    result = -1.04* np.exp(x[0]+0.2*x[1])
    return result

error= []
for n in range(10,50,5):
    points = []
    for x in np.linspace(0.5,1,n):
        for y in np.linspace(0.0,1.0,n):
            points.append([x,y])

    #points = [x  for x in points if linalg.norm(x)<=1]
    maximum0 = max(np.array(points)[:,0])
    maximum1 = max(np.array(points)[:,1])
    minimum1 = min(np.array(points)[:,1])
    minimum0 = min(np.array(points)[:,0])

    boundaryValuesx = [[k,np.exp(1.0+0.2*x[1])] for k,x in enumerate(points) if x[0]==maximum0 ]
    boundaryValuesy = [[k,np.exp(x[0]+0.2)] for k,x in enumerate(points) if x[1]==maximum1 and not x[0]==maximum0 and not x[0]==minimum0 ]
    boundaryRadial = []#[[k,0] for k,x in enumerate(points) if linalg.norm(x) >1.2]
    boundaryRadial += []#[[k,1] for k,x in enumerate(points) if linalg.norm(x) <=0.7 and linalg.norm(x)>0.3]
    boundaryRadial += []#[[k,1.1] for k,x in enumerate(points) if linalg.norm(x) <=0.3]
    boundaryValuesx += [[k,np.exp(0.5+0.2*x[1])] for k,x in enumerate(points) if x[0]==minimum0 ]
    boundaryValuesy += [[k,np.exp(x[0])] for k,x in enumerate(points) if x[1]==minimum1 and not x[0]==maximum0 and not x[0]==minimum0]
    boundary = np.array(boundaryValuesx + boundaryValuesy+boundaryRadial)
    
    Fem = FiniteElement(points,boundary,PDEMatrix= np.array([[1,0],[0,1]]))#,functionRHS = f)
    Fem.calculateGlobalStiffnessMatrix()
    Fem.calculateRightHandSide()
    Fem.solve()
    error.append([Fem.maxDiam,np.linalg.norm([sol(x) for x in Fem.triangulation.points]-Fem.solution)])
error = np.array(error)
#plt.triplot(Fem.triangulation.points[:,0],Fem.triangulation.points[:,1],Fem.triangulation.simplices.copy())
#plt.plot(Fem.triangulation.points[:,0],Fem.triangulation.points[:,1],'o')
fig = plt.figure()
#ax = Axes3D(fig)
#ax.plot_trisurf(Fem.triangulation.points[:,0],Fem.triangulation.points[:,1],Fem.triangulation.simplices.copy(),abs([sol(x) for x in Fem.triangulation.points]-Fem.solution))

#fig1 = plt.figure()
#ax1 = Axes3D(fig1)
#ax1.plot_trisurf(Fem.triangulation.points[:,0],Fem.triangulation.points[:,1],Fem.triangulation.simplices.copy(),[sol(x) for x in Fem.triangulation.points])
#
#print(Fem.triangulation.points[:])
plt.loglog(error[:,0],error[:,1])
print((np.log(error[0,1])-np.log(error[-1,1]))/(np.log(error[0,0])-np.log(error[-1,0])))
plt.loglog([error[0,0],error[-1,0]],[error[0,1],error[0,1]*(error[-1,0]-error[0,0])**(-1)])
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
