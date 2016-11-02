from scipy.spatial import Delaunay
from scipy import sparse
from scipy.sparse import linalg
from itertools import product
import math

import numpy as np

class FiniteElement:
    """Initiates the triangulation and calculates the stiffness matrix and righthand side"""
    def __init__(self,points,prescribedValues,PDEMatrix=np.eye(2),functionRHS=lambda x: 0):
        #referenceElement holds the points of the reference element from which all other elements
        #are calculated
        self.functionRHS= functionRHS
        self.referenceElement = np.array([[0,0],[1.,0],[0,1.]])
        self.triangulation = Delaunay(points)
        self.numberDOF = np.size(self.triangulation.points[:,0])
        self.maxDiam = 0
        self.prescribedValues = [] 
        if self.checkPrescribedValues(prescribedValues):
            self.prescribedValues = prescribedValues
        else:
            print("Error: Prescribed Value index not an integer")
        self.linearBasis = []
        self.linearBasis.append(lambda x : 1-x[0]-x[1])
        self.linearBasis.append(lambda x : x[0])
        self.linearBasis.append(lambda x : x[1])

        #Holds integral of two basisfunctons in one reference triangle
        self.elementaryBasisMatrix = 1*1.0/12*np.array([[1.,0.5,0.5],[0.5,1.,0.5],[0.5,0.5,1.]])
        
        self.gradBasis = []
        self.gradBasis.append(np.array([-1.,-1])) 
        self.gradBasis.append(np.array([1.,0]))
        self.gradBasis.append(np.array([0,1.]))
        self.rightHandSide = np.zeros(self.numberDOF)

        self.PDEMatrix= PDEMatrix

    def checkPrescribedValues(self,prescribedValues):
        for value in prescribedValues:
            if not float(int(value[0])) == value[0]:
                print("false")
                return False
        else:
            print("true")
            return True 

    def calculateTransform(self,triangleIndex):
        """ Calculates the affine linear transform from the reference
            Element to the given Triangle (Lx+b) with matrix 2x2 L and vector b
            and returns both of them in a tuple.
        """
        
        trianglePoints = self.triangulation.points[self.triangulation.simplices.copy()[triangleIndex]]
        self.maxDiam = max(np.linalg.norm(trianglePoints[0]-trianglePoints[1]),np.linalg.norm(trianglePoints[0]-trianglePoints[2]),np.linalg.norm(trianglePoints[2]-trianglePoints[1]),self.maxDiam)
        #hold the coordinates with 1 appended as last row
        referenceCoord = np.array([self.referenceElement[:,0],self.referenceElement[:,1],np.array([1,1,1])])       
        transformedCoord = np.array([trianglePoints[:,0],trianglePoints[:,1],np.array([1,1,1])])
        C = np.dot( transformedCoord,np.linalg.inv(referenceCoord)) 
        #L is in the first n x n submatrix of (n+1) x (n+1) Matrix C, b the last column without the last entry 
        return (C[0:-1,0:-1],C[:-1,-1])

    
    def calculateElementStiffnessMatrix(self,triangleIndex):
        """Returns the Elemnt Stiffness Matrix for the given Triangle, calcluated
            from the reference Element via the Transform and the reference Gradients"""
        transformMatrix,translateVector = self.calculateTransform(triangleIndex)
        transformMatrixInv = np.linalg.inv(transformMatrix)

        #comes from the Integraltransform and is of the operator from reference to transformed
        determinant = abs(np.linalg.det(transformMatrix))
        elementStiffnessMatrix = np.zeros((3,3))
        
        for row,column in product(range(3),range(3)):
            # 0.5 is the area of the reference triangle
            elementStiffnessMatrix[row,column] =0.5*determinant *np.dot(np.dot(np.dot(self.PDEMatrix,transformMatrixInv.T),self.gradBasis[row]),np.dot(transformMatrixInv.T,self.gradBasis[column]))
        return elementStiffnessMatrix

    def calculateGlobalStiffnessMatrix(self):
        """Loops over all elements, calculates the element Stiffness and assembles 
            it in CSC to a global Matrix"""
        
        globalRowIndices = []
        globalColumnIndices=[]
        globalData = []
        for ele,triPoints in enumerate(self.triangulation.simplices):
            #print("element:"+str(ele)+"/"+ str(len(self.triangulation.simplices)))
            eleStiffness = self.calculateElementStiffnessMatrix(ele)
            for row in range(3):
                if self.prescribedValues != [] and triPoints[row] in self.prescribedValues[:,0]:
                    continue
                else:
                    for column in range(3):
                        globalRowIndices.append(triPoints[row])
                        globalColumnIndices.append(triPoints[column])
                        globalData.append(eleStiffness[row,column])
        for v in self.prescribedValues:
            globalRowIndices.append(v[0])
            globalColumnIndices.append(v[0])
            globalData.append(1)
        self.GlobalStiffness = sparse.coo_matrix((globalData,(globalRowIndices,globalColumnIndices))).tocsc() 

        #for k,v in enumerate(self.GlobalStiffness.toarray()):
          #  print(k,v)

    def calculateRightHandSide(self):
        for ele,triPoints in enumerate(self.triangulation.simplices):
            elementRHS = self.calculateElementRightHandSide(ele)
            for index,entry in enumerate(elementRHS):
                self.rightHandSide[triPoints[index]] += entry
        for indexValue in self.prescribedValues:
            self.rightHandSide[int(indexValue[0])] = indexValue[1]
        #for k,v in enumerate(self.rightHandSide):
         #   print(k,v)
    
    def calculateElementRightHandSide(self,triangleIndex):
        transformMatrix,translateVector = self.calculateTransform(triangleIndex)
        determinant = abs(np.linalg.det(transformMatrix))
        trianglePoints =self.triangulation.points[self.triangulation.simplices[triangleIndex]]
        elementRHS = determinant*np.dot(self.elementaryBasisMatrix,np.array([self.functionRHS(x) for x in trianglePoints] ))
        midpoint = np.dot(transformMatrix,np.array([0.25,0.5])) + translateVector
        elementRHS2=[]
        for i in range(0,3):
            #0.5 is because of the area of the reference triangle
            elementRHS2.append(determinant*0.5*self.functionRHS(midpoint)*self.linearBasis[i](np.array([0.25,0.5])))
        #print(self.triangulation.simplices[triangleIndex], elementRHS)
        #print(elementRHS2)
        return elementRHS

    def solve(self):
        self.solution =sparse.linalg.spsolve(self.GlobalStiffness,self.rightHandSide)
        
        

    def getL2Error(self,exactSolution):
        #
        value = 0
        error = np.array(self.solution)-np.array([exactSolution(x) for x in self.triangulation.points])
        for ele,triPoints in enumerate(self.triangulation.simplices):
            transformMatrix,translateVector = self.calculateTransform(ele)
            determinant = abs(np.linalg.det(transformMatrix))
            value+=determinant*np.dot(error[triPoints]**2,np.array([1/6.,1/3.,1/3.]))
        return(math.sqrt(value))
    def Global2Local(self,globalIndex):
        """returns for a given global index the pair of local ones"""
        ele = int(math.floor(globalIndex / 3))
        return ele,globalIndex-3*ele
    def Local2Global(self,element,localVert):
        """returns for a given pair of local indices the global one"""
        return element*3 +localVert
