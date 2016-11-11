import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import math
from numpy import linalg
from FiniteElement import FiniteElement
#pointcharge in the middle with radial boundary

nx =50
ny =50


#f= lambda x : 1 if np.linalg.norm(x)<=0.2 else 0
def f(x):
    result = 0
    
    if np.linalg.norm(x) <= 0.5:
        result = 1
    elif np.linalg.norm(x-np.array([0.5,0]))<= 0.1:
        result= -0
    elif np.linalg.norm(x - np.array([0,0.5])) <= 0.1:
        result = -0
    elif np.linalg.norm(x-np.array([0,-0.5]))<= 0.1:
        result= 0
    
    else:
        result=0
    #result = -1.04* np.exp(x[0]+0.2*x[1])
    return result

n= 10
error= []
pointerror= []

#for n in range(3,100,10):
points = []
for x in np.linspace(-2,2,n):
    for y in np.linspace(-2,2,n):
        points.append([x,y])

maximum0 = max(np.array(points)[:,0])
maximum1 = max(np.array(points)[:,1])
minimum0 = min(np.array(points)[:,0])
minimum1 = min(np.array(points)[:,1])

boundaryValuesx = [[k,0] for k,x in enumerate(points) if x[0]==maximum0 ]
boundaryValuesy = [[k,0] for k,x in enumerate(points) if x[1]==maximum1 and not x[0]==maximum0 and not x[0]==minimum0 ]

boundaryValuesx += [[k,0] for k,x in enumerate(points) if x[0]==minimum0 ]
boundaryValuesy += [[k,0] for k,x in enumerate(points) if x[1]==minimum1 and not x[0]==maximum0 and not x[0]==minimum0]

boundary = np.array(boundaryValuesx + boundaryValuesy)
"""
boundaryRadial = [[k,0] for k,x in enumerate(points) if linalg.norm(x) >1.2]
boundaryRadial += [[k,1] for k,x in enumerate(points) if linalg.norm(x) <=0.7 and linalg.norm(x)>0.3]
boundaryRadial += []#[[k,1.1] for k,x in enumerate(points) if linalg.norm(x) <=0.3]
boundaryValuesx = []#[[k,np.exp(0.5+0.2*x[1])] for k,x in enumerate(points) if x[0]==minimum0 ]
boundaryValuesy = []#[[k,np.exp(x[0])] for k,x in enumerate(points) if x[1]==minimum1 and not x[0]==maximum0 and not x[0]==minimum0]

boundary = np.array(boundaryValuesx + boundaryValuesy+boundaryRadial)
"""
#boundary=np.array([])

Fem = FiniteElement(points,boundary,PDEMatrix= np.array([[1,0],[0,1]]),functionRHS = f)
Fem.calculateGlobalStiffnessMatrix()
Fem.calculateRightHandSide()
Fem.solve()

#plt.triplot(Fem.triangulation.points[:,0],Fem.triangulation.points[:,1],Fem.triangulation.simplices.copy())
#plt.plot(Fem.triangulation.points[:,0],Fem.triangulation.points[:,1],'o')
print(Fem.solution)
fig = plt.figure()
ax = Axes3D(fig)
ax.plot_trisurf(Fem.triangulation.points[:,0],Fem.triangulation.points[:,1],Fem.triangulation.simplices.copy(),Fem.solution)

plt.show()


