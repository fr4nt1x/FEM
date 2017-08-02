import numpy as np
import matplotlib.pyplot as plt 
from Mesh import Mesh
from Boundary import PolygonalBoundary


sol = lambda x : np.exp(x[0]+0.2*x[1])
# points = np.array([[0,0],[1,0],[1,1],[0,1]])
# triangles = [[0,1,2],[0,2,3]]
# boundaryEdges = [[[0,1],0],[[1,2],1],[[2,0],None],[[2,3],None],[[3,0],2]]
points = np.array([[0,0],[1,0],[1,1],[0,1],[-1,1],[-1,0],[-1,-1],[0,-1]])
polBoundary= PolygonalBoundary(points, [sol]*np.shape(points)[0])
triangles = [[0,1,2],[0,2,3],[0,3,4],[0,4,5],[0,5,6],[0,6,7]] 
#the boundary edges holds all the edges of each traingle, as a  pair of an edge and 
#the index of the polygonalBoundary it belongs to.
boundaryEdges = [[[0,1],0],[[1,2],1],[[2,0],None],[[2,3],2],[[3,0],None],[[3,4],2],[[4,0],None],[[4,5],3],[[5,0],None],[[5,6],3],[[6,0],None],[[6,7],4],[[7,0],5]]
m = Mesh(points, triangles, boundaryEdges,polBoundary)
def f(x) :
    return np.sin(x[0])*np.sin(x[1])
print("Before")
print("edges",m.edges)
m.refineMeshHalf(5)
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

