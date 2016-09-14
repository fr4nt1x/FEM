from scipy.spatial import Delaunay
import numpy as np

class FiniteElement:
    """Initiates the triangulation and calculates the stiffness matrix and righthand side"""
    def __init__(self,points):
        self.referenceElement = np.array([[0,0],[1.,0],[0,1]])
        self.triangulation = Delaunay(points)
        self.linearBasis0 = lambda x,y : 1-x-y
        self.linearBasis1 = lambda x,y : x
        self.linearBasis2 = lambda x,y : y
        self.gradBasis0 = np.array([-1,-1]) 
        self.gradBasis1 = np.array([1,0])
        self.gradBasis2 = np.array([0,1])

    def calculateTransform(self,triangleIndex):
        """ Calculates the affine linear transform from the reference Element to the given Triangle (Lx+b) 
            with matrix 2x2 L and vector b
        """
        trianglePoints = self.triangulation.points[self.triangulation.simplices.copy()[triangleIndex]]
        A = np.array([self.referenceElement[:,0],self.referenceElement[:,1],np.array([1,1,1])])
        B = np.array([trianglePoints[:,0],trianglePoints[:,1],np.array([1,1,1])])
        C = np.dot(B,np.linalg.inv(A))
        self.transformMatrix = C[0:-1,0:-1]
        self.translateVector = C[:-1,-1]
    
