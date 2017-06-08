import numpy as np
import matplotlib.pyplot as plt
import time
import pprint

class Mesh:
    """
    Class that holds the triangulation of the domain.
    Is initiated with a list of points, and a list of triangles which hold indices.
    """

    def __init__(self,points,triangles,boundaryEdges):
        self.points = points
        self.triangles = triangles
        self.edges = boundaryEdges
        self.boundaryValues = self.generateBoundaryValues()
        self.trianglesWithEdges = self.initiateTrianglesWithEdges() 
    
    def initiateTrianglesWithEdges(self):
        trianglesWithEdges = []
        edges = [e[0] for e in self.edges]
        for triangle in self.triangles:
            triangleEdges = []
            for i in range(0,3):
                triangleEdges.append([triangle[i],triangle[(i+1)%3]])
            trianglesWithEdges.append([self.findIndexOfEdge(edges,triangleEdges[0]), self.findIndexOfEdge(edges,triangleEdges[1]), self.findIndexOfEdge(edges,triangleEdges[2])])
        return trianglesWithEdges
            
    def generateBoundaryValues(self):
        boundaryValues = []
        for edgeSol in self.edges:
            for pointIndex in edgeSol[0]:
                if edgeSol[1] != None and pointIndex not in [x[0] for x in boundaryValues]:
                    boundaryValues.append([pointIndex,edgeSol[1](self.points[pointIndex])])
        return np.array(boundaryValues)

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

    def findIndexOfEdge(self,oldEdges,edge):
        #1 says orientation in triangles == in self.edges 
        #0 not
        for index,edg in enumerate(oldEdges):
            if edg == edge:
                return [index,1]
            elif edg == [x for x in reversed(edge)]:
                return [index,0]

    def generateNewTriangles(self,oldEdges, newPoints):
        #newPoints hold the index of the new points
        newTriangles = []
        newTrianglesWithEdges = []
        for triIndex, triangle in enumerate(self.triangles):
            #hold the indices of the edges for the oldedges
            triangleEdges= self.trianglesWithEdges[triIndex]

            newTrianglePoints = []
            # triangleEdges = []
            # for i in range(0,3):
                # triangleEdges.append([triangle[i],triangle[(i+1)%3]])

            for edge in triangleEdges:
                newTrianglePoints.append(newPoints[edge[0]])

            # triangleEdges = [edge[0] for edge in triangleEdges]
            #append the edges that are in the middle of the old triangles and remember their index in self.edges
            lengthOfEdges = len(self.edges)
            newEdgeIndices = [lengthOfEdges+i for i in range(0,3)]

            self.edges.append([[newTrianglePoints[0],newTrianglePoints[1]],None])

            self.edges.append([[newTrianglePoints[1],newTrianglePoints[2]],None])
            self.edges.append([[newTrianglePoints[2],newTrianglePoints[0]],None])

            # print("Oldedges:", oldEdges)
            # print("edges:",[e[0] for e in self.edges])
            # print("newEdge",newEdgeIndices)

            #bottom left triangle

            newTriangles.append([triangle[0],newTrianglePoints[0],newTrianglePoints[2]])
            orient0 = triangleEdges[0][1]
            orient2 = triangleEdges[2][1]
            index0 = 2*triangleEdges[0][0] +(1+orient0)%2
            index2 = 2*triangleEdges[2][0] +(orient2)%2
            newTrianglesWithEdges.append([[index0,orient0],[newEdgeIndices[2],0],[index2,orient2]])
            

            #middle triangle
            newTriangles.append([newTrianglePoints[0],newTrianglePoints[1],newTrianglePoints[2]])
            newTrianglesWithEdges.append([[newEdgeIndices[0],1],[newEdgeIndices[1],1],[newEdgeIndices[2],1]])

            #right bottom triangle
            orient0 = triangleEdges[0][1]
            orient1 = triangleEdges[1][1]
            index0 = 2*triangleEdges[0][0] +(orient0)%2
            index1 = 2*triangleEdges[1][0] +(1+orient1)%2
            newTriangles.append([newTrianglePoints[0],triangle[1],newTrianglePoints[1]])
  
            newTrianglesWithEdges.append([[index0,orient0],[index1,orient1],[newEdgeIndices[0],0]])

            #top triangle
            orient1 = triangleEdges[1][1]
            orient2 = triangleEdges[2][1]
            index1 = 2*triangleEdges[1][0] +(orient1)%2
            index2 = 2*triangleEdges[2][0] +(1+orient2)%2
            newTriangles.append([newTrianglePoints[1],triangle[2],newTrianglePoints[2]])
            newTrianglesWithEdges.append([[index1,orient1],[index2,orient2],[newEdgeIndices[1],0]])
    
        self.triangles = newTriangles
        self.trianglesWithEdges = newTrianglesWithEdges
        # print("POINTS:")
        # print(self.points)
        # print("EDGES:")
        # print([e[0] for e in self.edges])
        # print("TRIANGLES:")
        # pprint.pprint(self.triangles)
        # print("TRIANGLESWITHEDGES")
        # pprint.pprint(self.trianglesWithEdges)

    def refineMesh(self,refinementSteps):
        for step in range(1,refinementSteps+1):
            listOfEdges = [x[0] for x in self.edges]
            print("Step : " +str(step) + "/" +str(refinementSteps))
            numberOfNewPoints = np.shape(listOfEdges)[0]
            numberOfOldPoints =  np.shape(self.points)[0]
            newPoints = np.zeros([numberOfNewPoints + numberOfOldPoints,2])
            newPoints[0:numberOfOldPoints,:] = self.points 
            newTriangles = [] 
            indexOfNewPoints = []
            newEdges = [None]*(2*numberOfNewPoints)  
            newBoundaryValues = []

            for index,edge in enumerate(listOfEdges):
                # print("edge: ",edge)
                firstPoint = self.points[int(edge[0])]
                secondPoint = self.points[int(edge[1])]
                newPoint = 0.5 * (firstPoint + secondPoint) 
                newIndex = numberOfOldPoints + index
                newEdges[2*index] = [[int(edge[0]),newIndex],self.edges[index][1]] 
                newEdges[2*index+1] = [[newIndex,int(edge[1])],self.edges[index][1]] 
                newPoints[newIndex,:] = newPoint 
                indexOfNewPoints.append(newIndex)
                    
            self.points = newPoints
            self.edges = newEdges
            self.generateNewTriangles(listOfEdges,indexOfNewPoints)
            self.boundaryValues = self.generateBoundaryValues()

                
