import numpy as np
import matplotlib.pyplot as plt 
from Mesh import Mesh
from Boundary import PolygonalBoundary


sol = lambda x : np.exp(x[0]+0.2*x[1])
# points = np.array([[0,0],[1,0],[1,1],[0,1]])
# triangles = [[0,1,2],[0,2,3]]
# boundaryEdges = [[[0,1],0],[[1,2],1],[[2,0],None],[[2,3],None],[[3,0],2]]
points = np.array([[0,0],[1,0],[1,1],[0,1]])
triangles = [[0,1,2],[0,2,3]]
boundaryEdges = [[[0,1],0],[[1,2],1],[[2,0],None],[[2,3],2],[[3,0],3]]
polBoundary= PolygonalBoundary(points, [sol]*np.shape(points)[0])
m = Mesh(points, triangles, boundaryEdges,polBoundary)
def f(x) :
    return np.sin(x[0])*np.sin(x[1])
print("Before")
print("edges",m.edges)
m.refineMeshHalf(9)
# print("refineHalf")
# print("edge",m.edges)
# print("points",m.points)
# print("triangles",m.triangles)
# print("trianglesedges", m.trianglesWithEdges)
m.plotTriangles()

# m.dumpToJson("./meshtest","Ref9")

# m2 = Mesh(points, triangles, boundaryEdges,polBoundary)
# m2.loadFromJson("./meshtest","RefUniform01_3")
# print(m.points)
# print(m.triangles)
# print(m.boundaryValues)
# m2.plotTriangles()

