from scipy.spatial import Delaunay
import numpy as np

class FiniteElement:
    """Initiates the triangulation and calculates the stiffness matrix and righthand side"""
    def __init__(self,points,PDEMatrix=np.eye(2)):
        self.referenceElement = np.array([[0,0],[1.,0],[0,1]])
        self.triangulation = Delaunay(points)
        self.linearBasis = ()
        self.linearBasis[0] = lambda x,y : 1-x-y
        self.linearBasis[1] = lambda x,y : x
        self.linearBasis[2] = lambda x,y : y
        self.gradBasis = ()
        self.gradBasis[0] = np.array([-1,-1]) 
        self.gradBasis[1] = np.array([1,0])
        self.gradBasis[2] = np.array([0,1])
        
        self.PDEMatrix= PDEMatrix
    def calculateTransform(self,triangleIndex):
        """ Calculates the affine linear transform from the reference Element to the given Triangle (Lx+b) 
            with matrix 2x2 L and vector b, and returns both of them in a tuple
        """
        trianglePoints = self.triangulation.points[self.triangulation.simplices.copy()[triangleIndex]]
        A = np.array([self.referenceElement[:,0],self.referenceElement[:,1],np.array([1,1,1])])
        B = np.array([trianglePoints[:,0],trianglePoints[:,1],np.array([1,1,1])])
        C = np.dot(B,np.linalg.inv(A))
        
        return (C[0:-1,0:-1],C[:-1,-1])
    
    def calculateElementStiffnessMatrix(self,triangleIndex):
        transformMatrix,translateVector = self.calculateTransform(trianlgeIndex)
        transformMatrixInv = np.linalg.inv(transformMatrix)
        print(transformMatrix)
        elementStiffnessMatrix = np.zeros((3,3))
        for row,column in range(3):
            elementStiffnessMatrix[row,column] = np.dot(np.dot(np.dot(self.PDEMatrix,transformMatrixInv.T),self.grad[row]),np.dot(transformMatrixInv.T,self.grad[column]))
        print(elementStiffnessMatrix)
