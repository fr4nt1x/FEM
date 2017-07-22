import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import math
from Mesh import Mesh
from numpy import linalg
from ProjectedGradient import ProjectedGradient
from Boundary import PolygonalBoundary 

alpha = 1
lam = 2/3.0

def state(x):
    return (x[0]**3-x[0])*(x[1]**3-x[1])

def control(x):
    return -6*(x[0]*(x[1]**3-x[1]) + x[1]*(x[0]**3-x[0]))


def RHSAddendum(x):
    return 2*((x[0]-x[0]**2)+ (x[1]-x[1]**2))- (x[0]-x[0]**2)*( x[1]-x[1]**2)

u_d = lambda x : state(x) 

def boundaryFunc(x):
    return 0;

points = np.array([[0,0],[1,0],[1,1],[0,1],[-1,1],[-1,0],[-1,-1],[0,-1],[-0.5,0.5]])
polBoundary= PolygonalBoundary(points, [boundaryFunc]*np.shape(points)[0])
triangles = [[0,1,2],[0,2,3],[0,3,8],[0,8,5],[8,3,4],[8,4,5],[0,5,6],[0,6,7]] 
#the boundary edges holds all the edges of each traingle, as a  pair of an edge and 
#the index of the polygonalBoundary it belongs to.
boundaryEdges = [[[0,1],0],[[1,2],1],[[2,0],None],[[2,3],2],[[3,0],None],[[3,8],None],[[8,0],None],[[8,5],None],[[3,4],2],[[8,4],None],[[4,5],3],[[5,0],None],[[5,6],3],[[6,0],None],[[6,7],4],[[7,0],5]]
m = Mesh(points, triangles, boundaryEdges,polBoundary)
# m.plotTriangles()

error = []
error2 = []

for i in range(0,4):
    m.refineMesh(1)
    proj = ProjectedGradient(m,[0 for x in m.points], [u_d(x) for x in m.points],lam,  alpha, tol = 1e-10 )
    proj.solve()
    # proj.loadFromJson(".","8")
    # error.append([proj.mesh.getDiameter(),proj.getL2ErrorControl(proj.referenceControl)])

# proj.dumpToJson(".","8")
# m.plotTriangles()

# for ele,triPoints in enumerate(proj.mesh.triangles):
    # transformMatrix,translateVector = proj.calculateTransform(ele)
    # determinant = abs(np.linalg.det(transformMatrix))
    # #Last vector is the precalculated integral of the basisfunctions over a reference element
    # def control(x):
        # #onlz points at the transformed triangle should be put in
        # x = np.dot(np.linalg.inv(transformMatrix),x-translateVector)
        # return proj.control[triPoints[0]]*proj.linearBasis[0](x)+ proj.control[triPoints[1]]*proj.linearBasis[1](x)+proj.control[triPoints[2]]*proj.linearBasis[2](x)

# fig = plt.figure()
# ax = Axes3D(fig)
# ax.text(0,0.5,0,"control",color="red")
# ax.plot_trisurf(proj.mesh.points[:,0],proj.mesh.points[:,1],proj.mesh.triangles,[state(x) for x in proj.mesh.points])

# fig1 = plt.figure()
# ax1 = Axes3D(fig1)
# ax1.text(0,0.5,0,"state",color="red")
# ax1.plot_trisurf(proj.mesh.points[:,0],proj.mesh.points[:,1],proj.mesh.triangles,proj.state)
# plt.show()


# fig = plt.figure()
# ax = Axes3D(fig)
# ax.text(0,0.5,0,"Control",color="red")
# ax.plot_trisurf(m.points[:,0],m.points[:,1],m.triangles.copy(),proj.control)

# fig1 = plt.figure()
# ax1 = Axes3D(fig1)
# ax1.text(0,0.5,0,"State",color="red")
# ax1.plot_trisurf(m.points[:,0],m.points[:,1],m.triangles.copy(),proj.state)

# fig2 = plt.figure()
# ax4 = Axes3D(fig2)
# ax4.text(0,0.5,0,"State real",color="red")
# ax4.plot_trisurf(m.points[:,0],m.points[:,1],m.triangles.copy(),[state(p) for p in m.points])
print("Gauss:")
for index, e in enumerate(error):
    if index < len(error)-1:
        print(np.log(error[index+1][0]/e[0])/np.log(error[index+1][1]/e[1]))
print("WithoutGauss:")
for index, e in enumerate(error2):
    if index < len(error2)-1:
        print(np.log(error2[index+1][0]/e[0])/np.log(error2[index+1][1]/e[1]))
print("ERROR: ",error)
plt.show()

#plt.loglog([e[0] for e in error], [e[1] for e in error],'o')
