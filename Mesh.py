import numpy as np
import matplotlib.pyplot as plt

class Mesh:
    """
    Class that holds the triangulation of the domain.
    Is initiated with a list of points, and a list of triangles which hold indices.
    """

    def __init__(self,points,triangles):
        self.points = points
        self.triangles = triangles
    
    def plotTriangles (self):
        plt.triplot(self.points[:,0],self.points[:,1],self.triangles)
        plt.show()
    
    def getEdges(self):
        # listOfEdges holds a list of all the edges and which triangle they belong to  
        listOfEdges = []

        for index, triangle in enumerate(self.triangles):
            edge1 = [triangle[0],triangle[1]]
            edge2 = [triangle[1],triangle[2]]
            edge3 = [triangle[2],triangle[0]]
            
            newEdges = [edge1,edge2,edge3]
            
            for edge in newEdges:
                reversedEdge = [x for x in reversed(edge)]
                if reversedEdge not in listOfEdges:
                    listOfEdges.append(edge)
        return listOfEdges

    def createNewTriangles(self,newPointIndex, edge, oldTriangle):
        """newPointIndex holds index of the New point in self.points 
           edge pair of point indices, that specifies the edge on which the new point was created
           oldTriangleIndexR  the index of triangle right of the edge
           oldTriangleIndexL the index of the triangle left of the edge
        """
        pointNotOnEdge = list(set(oldTriangle)- set(edge))
        pointNotOnEdge = pointNotOnEdge[0]
        newTriangle1 = [edge[0],newPointIndex,pointNotOnEdge]
        newTriangle2 = [newPointIndex,edge[1],pointNotOnEdge]
        self.triangles.append(newTriangle1)
        self.triangles.append(newTriangle2)

    def findIndexOfEdge(self,oldEdges,edge):
        for index,edg in enumerate(oldEdges):
            if edg == edge or [x for x in reversed(edg)]== edge:
                return index

    def generateNewTriangles(self,oldEdges, newPoints):
        newTriangles = []
        for triangle in self.triangles:
            triangleEdges= []
            newTrianglePoints = []

            for i in range(0,3):
                triangleEdges.append([triangle[i],triangle[(i+1)%3]])

            for edge in triangleEdges:
                newTrianglePoints.append(newPoints[self.findIndexOfEdge(oldEdges,edge)])
            #bottom left riangle
            newTriangles.append([triangle[0],newTrianglePoints[0],newTrianglePoints[2]])
            #middle triangle
            newTriangles.append([newTrianglePoints[0],newTrianglePoints[1],newTrianglePoints[2]])
            #right bottom triangle
            newTriangles.append([newTrianglePoints[0],triangle[1],newTrianglePoints[1]])
            #top triangle
            newTriangles.append([newTrianglePoints[1],triangle[2],newTrianglePoints[2]])
    
        self.triangles = newTriangles

    def refineMesh(self,refinementSteps):
        for step in range(1,refinementSteps+1):
            listOfEdges = self.getEdges()
            print("Step : " +str(step) + "/" +str(refinementSteps))
            numberOfNewPoints = np.shape(listOfEdges)[0]
            numberOfOldPoints =  np.shape(self.points)[0]
            newPoints = np.zeros([numberOfNewPoints + numberOfOldPoints,2])
            newPoints[0:numberOfOldPoints,:] = self.points 
            newTriangles = [] 
            newEdges = []
            indexOfNewPoints = []

            for index,edge in enumerate(listOfEdges[:]):
                # print("edge: ",edge)
                firstPoint = self.points[edge[0]]
                secondPoint = self.points[edge[1]]
                newPoint = 0.5 * (firstPoint + secondPoint) 
                newIndex = numberOfOldPoints + index
                newPoints[newIndex,:] = newPoint 
                indexOfNewPoints.append(newIndex)
                newEdges.append([edge[0],newIndex,edge[1]])

            self.points = newPoints
            self.generateNewTriangles(listOfEdges,indexOfNewPoints)
            
            triangles = []
                
