import numpy as np
import matplotlib.pyplot as plt 
from Mesh import Mesh
from Boundary import PolygonalBoundary


points = np.array([[0,0],[1,0],[1,1],[0,1]])
triangles = [[0,1,2],[2,3,0]]
sol = lambda x : np.exp(x[0]+0.2*x[1])
boundaryEdges = [[[0,1],0],[[1,2],1],[[2,0],None],[[2,3],2],[[3,0],3]]
polBoundary= PolygonalBoundary(points, [sol]*np.shape(points)[0])
m = Mesh(points, triangles, boundaryEdges,polBoundary)
def f(x) :
    return np.sin(x[0])*np.sin(x[1])

m.refineMesh(3)
m.setValuesAtPoints([f(x) for x in m.points])
m.plotValues()
m.refineMesh(3)
m.plotValues()
plt.show()
# m.dumpToJson("./meshtest","Ref9")

# m2 = Mesh(points, triangles, boundaryEdges,polBoundary)
# m2.loadFromJson("./meshtest","RefUniform01_3")
# print(m.points)
# print(m.triangles)
# print(m.boundaryValues)
# m2.plotTriangles()

