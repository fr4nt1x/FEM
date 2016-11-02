import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import math
from numpy import linalg
from FiniteElement import FiniteElement

sol = lambda x : 1+x[0]**2+2*x[1]**2
nx =50
ny =50


def f(x):
    result = -6.
    return result

error= []
pointerror= []
boundary=[]

for n in range(5,19,1):
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
    Fem.getL2Error(sol)
    #error.append([Fem.maxDiam,np.linalg.norm(np.array([sol(x) for x in Fem.triangulation.points])-np.array(Fem.solution),np.infty)])
    error.append([Fem.maxDiam,Fem.getL2Error(sol)])
    #print(Fem.triangulation.points)
    #print(Fem.solution)
    #print(Fem.rightHandSide)
    #print([(k,f(x)) for k,x in enumerate(Fem.triangulation.points)])
    #pointlist=[]
    #for k,x in enumerate(Fem.triangulation.points):
    #    if abs(Fem.solution[k] - sol(x)) == error[-1][-1]:
     #       pointlist.append([k,x])
    #pointerror.append([error[-1][-1],pointlist]) 
#for x in pointerror:
    #print("point",x[1])
#print(Fem.solution,[sol(x) for x in Fem.triangulation.points])
#plt.triplot(Fem.triangulation.points[:,0],Fem.triangulation.points[:,1],Fem.triangulation.simplices.copy())
#plt.plot(Fem.triangulation.points[:,0],Fem.triangulation.points[:,1],'o')

#fig = plt.figure()
#ax = Axes3D(fig)
#ax.plot_trisurf(Fem.triangulation.points[:,0],Fem.triangulation.points[:,1],Fem.solution)

#fig1 = plt.figure()
#ax1 = Axes3D(fig1)
#ax1.plot_trisurf(Fem.triangulation.points[:,0],Fem.triangulation.points[:,1],[sol(x) for x in Fem.triangulation.points])

error = np.array(error)
print(error)
plt.loglog(error[:,0],error[:,1],'o')
#plt.loglog([error[0,0],error[-1,0]],[error[0,1],error[0,1]*(error[-1,0]-error[0,0])**(5)])
plt.show()
