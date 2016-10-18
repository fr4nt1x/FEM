import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from FiniteElement import FiniteElement


Fem = FiniteElement(np.array([[1.,1],[0,2],[0,1]]))
print(Fem.triangulation)
Fem.calculateElementStiffnessMatrix(0)

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
