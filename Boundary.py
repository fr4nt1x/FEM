import numpy as np

class PolygonalBoundary:
    """
    Class that holds the polygonal boundary of the domain
    It has points of the polygonal domain, its edges and the function to calculate
    values at a edg
    """
    def __init__(self, points,funcAtEdge, edges=None):
        self.points = points
        self.functionAtEdges = funcAtEdge
        
        if edges== None:
            self.generateEdges()
        else:
            self.edges =edges

    def generateEdges(self):
        self.edges = []
        numberOfPoints = np.shape (self.points)[0]
        for pointIndex in range(0,numberOfPoints):
            self.edges.append([pointIndex,(pointIndex+1)%numberOfPoints]) 

        print(self.edges)

