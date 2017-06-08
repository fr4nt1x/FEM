import numpy as np
from Mesh import Mesh


points = np.array([[0,0],[1,0],[1,1],[0,1]])
triangles = [[0,1,2],[2,3,0]]
sol = lambda x : np.exp(x[0]+0.2*x[1])
boundaryEdges = [[[0,1],sol],[[1,2],sol],[[2,0],None],[[2,3],sol],[[3,0],sol]]
m = Mesh(points, triangles, boundaryEdges)

m.refineMesh(9)
# print(m.points)
# print(m.triangles)
# print(m.boundaryValues)
m.plotTriangles()

