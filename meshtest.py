import numpy as np
from Mesh import Mesh


points = np.array([[0,0],[1,0],[1,1],[0,1]])
triangles = [[0,1,2],[2,3,0]]

m = Mesh(points, triangles)

m.refineMesh(6)
print(len(m.triangles))
m.plotTriangles()
