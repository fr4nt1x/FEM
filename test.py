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


def f(x):
    result = -1.04* np.exp(x[0]+0.2*x[1])
    return result

error= []
pointerror= []
boundary=[]

for n in range(20,21,1):
    points = []
    for x in np.linspace(0,1,n):
        for y in np.linspace(0,1,n):
            points.append([x,y])

    maximum0 = max(np.array(points)[:,0])
    maximum1 = max(np.array(points)[:,1])
    minimum0 = min(np.array(points)[:,0])
    minimum1 = min(np.array(points)[:,1])

    boundaryValuesx = [[k,sol(x)] for k,x in enumerate(points) if x[0]==maximum0 ]
    boundaryValuesy = [[k,sol(x)] for k,x in enumerate(points) if x[1]==maximum1 and not x[0]==maximum0 and not x[0]==minimum0 ]

    boundaryValuesx += [[k,sol(x)] for k,x in enumerate(points) if x[0]==minimum0 ]
    boundaryValuesy += [[k,sol(x)] for k,x in enumerate(points) if x[1]==minimum1 and not x[0]==maximum0 and not x[0]==minimum0]

    boundary = np.array(boundaryValuesx + boundaryValuesy)

    Fem = FiniteElement(points,boundary,PDEMatrix= np.array([[1,0],[0,1]]),functionRHS = f)
    Fem.calculateGlobalStiffnessMatrix()
    Fem.calculateRightHandSide()
    Fem.solve()
    #print(Fem.solution)

    error.append([Fem.maxDiam,np.linalg.norm(np.array([sol(x) for x in Fem.triangulation.points])-np.array(Fem.solution),np.infty)])
    pointlist=[]
    for k,x in enumerate(Fem.triangulation.points):
        if abs(Fem.solution[k] - sol(x)) == error[-1][-1]:
            pointlist.append([k,x])
    pointerror.append([error[-1][-1],pointlist]) 
#for x in pointerror:
 #   print("point",x[1])

error = error
print(error[0][1])

#plt.triplot(Fem.triangulation.points[:,0],Fem.triangulation.points[:,1],Fem.triangulation.simplices.copy())
#plt.plot(Fem.triangulation.points[:,0],Fem.triangulation.points[:,1],'o')

#fig = plt.figure()
#ax = Axes3D(fig)
#ax.plot_trisurf(Fem.triangulation.points[:,0],Fem.triangulation.points[:,1],Fem.triangulation.simplices.copy(),Fem.solution)

#fig1 = plt.figure()
#ax1 = Axes3D(fig1)
#ax1.plot_trisurf(Fem.triangulation.points[:,0],Fem.triangulation.points[:,1],Fem.triangulation.simplices.copy(),[sol(x) for x in Fem.triangulation.points])

#

#plt.loglog(error[:,0],error[:,1],'o')
#print(error)
#plt.loglog([error[0,0],error[-1,0]],[error[0,1],error[0,1]*(error[-1,0]-error[0,0])**(5)])
#plt.show()
