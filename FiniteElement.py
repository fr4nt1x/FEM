from Mesh import Mesh
from scipy import sparse
from scipy.sparse import linalg
from itertools import product
from numpy.polynomial.legendre import leggauss
import math

import numpy as np

class FiniteElement:
    """
    Class for solving an Elliptic equation of the form div(A dot grad u) = f for given Matrix A(PDEMatrix).
    
    It's using the weak form of the equation and solves this with a Finite Element approach.
    
    The Ansatzfunctions are linear Lagrangefunctions 
    on these Triangles.

    All the integrals that are needed for Stiffnessmatrix and Righthandside get calculated via a reference triangle.
    """

    def __init__(self,mesh,PDEMatrix=np.eye(2),functionRHS = None, RHSEvaluatedAtTrianglePoints = None):
        """
        Initiates the triangulation and calculates the stiffness matrix and righthand side

        INPUT: 
                
                points             points from which the triangulation of the domain will be calculated

                prescribedValues   boundary values as a list of pairs with format [index,value], where index 
                                   specifies the  point the given value is prescribed, such that
                                   u(points[index]) = value for u the solution of the problem
                
                PDEMatrix          The Matrix which specifies the PDE, should be positive definite to guarante the PDE
                                   is elliptic
                
                functionRHS        a function handle which specifies the Righthandside of the PDE, Should be an L2 
                                   function

        """
        
        self.functionRHS= functionRHS


        #Calculate a delaunay triangulation of the input points
        self.mesh = mesh
        self.triangulation = self.mesh
        self.RHSEvaluatedAtTrianglePoints = RHSEvaluatedAtTrianglePoints
        if self.functionRHS != None: 
            self.RHSEvaluatedAtTrianglePoints = [self.functionRHS(x) for x in self.triangulation.points]
            
        #Uses to initiate the stiffness matrix and the Rhs with the correct size
        self.numberDOF = np.size(self.triangulation.points[:,0])

        self.prescribedValues = np.array(mesh.boundaryValues)

        #Holds integral of two basisfunctons over one reference triangle
        self.elementaryBasisMatrix = 1.0/12*np.array([[1.,0.5,0.5],[0.5,1.,0.5],[0.5,0.5,1.]])

        #initiate Righthandside with zeros
        self.rightHandSide = np.zeros(self.numberDOF)
        
        #strong form of PDE is: div(A dot grad(u)) = f, where A is PDEMatrix
        self.PDEMatrix= PDEMatrix
        print("FEM init finished")

    def checkPrescribedValues(self,prescribedValues):
        """
        Simple check if for the given index,value pairs, the index is an integer
        """

        for value in prescribedValues:
            if not float(int(value[0])) == value[0]:
                print("false")
                return False
        else:
            return True 

    
    def calculateElementStiffnessMatrix(self,triangleIndex):
        """
        Returns the element stiffness matrix for the given triangle, calculated
        from the reference Element with transformation formula and the exact integrals of two  reference Gradients
        """

        transformMatrix,translateVector = self.mesh.calculateTransform(triangleIndex)
        transformMatrixInv = np.linalg.inv(transformMatrix)

        #needed for the Integraltransformation
        #is Jacobian of the mapping from reference to transformed element
        determinant = abs(np.linalg.det(transformMatrix))

        elementStiffnessMatrix = np.zeros((3,3))
        
        for row,column in product(range(3),range(3)):
            # 0.5 is the area of the reference triangle, since all functions are constant for
            #linear elements this yields the integral
            elementStiffnessMatrix[row,column] = 0.5*determinant *np.dot(np.dot(np.dot(self.PDEMatrix,transformMatrixInv.T),self.mesh.gradBasis[row]),np.dot(transformMatrixInv.T,self.mesh.gradBasis[column]))
        return elementStiffnessMatrix

    def calculateGlobalStiffnessMatrix(self):
        """
        loops over all elements, calculates the element stiffness matrix and assemble them into a global 
        matrix in csc format
        """
        
        globalRowIndices = []
        globalColumnIndices=[]
        globalData = []
        

        for ele,triPoints in enumerate(self.triangulation.triangles):
            eleStiffness = self.calculateElementStiffnessMatrix(ele)

            #append the entries at the correct position in the global matrix, only if no Boundary
            #condition is prescribed at this degree of freedom 
            for row in range(3):
                if self.prescribedValues != [] and triPoints[row] in self.prescribedValues[:,0]:
                    continue
                else:
                    for column in range(3):
                        globalRowIndices.append(triPoints[row])
                        globalColumnIndices.append(triPoints[column])
                        globalData.append(eleStiffness[row,column])
        #Dirichlet boundaries at position j are enforced, by adding a row (0,0,..,1,...,0) (j-th entry),
        #and the value in the right hand side at the position j   
        for v in self.prescribedValues:

            globalRowIndices.append(v[0])
            globalColumnIndices.append(v[0])
            globalData.append(1)

        #create global matrix in coordinate format and convert it to csc for faster solving
        self.GlobalStiffness = sparse.coo_matrix((globalData,(globalRowIndices,globalColumnIndices))).tocsc() 
        # print(self.GlobalStiffness.toarray())

    def calculateRightHandSide(self):
        """
        Calculates the Righthandside, by calling the entries for each element and assembling the result
        into an global vector
        """

        #each element triangulation.simplices is a list with the point indices of the triangle
        for ele,triPoints in enumerate(self.triangulation.triangles):
            elementRHS = self.calculateElementRightHandSide(ele)

            for index,entry in enumerate(elementRHS):
                self.rightHandSide[triPoints[index]] += entry

        #overwrite the values at the prescribed boundaries condition
        for indexValue in self.prescribedValues:
            self.rightHandSide[int(indexValue[0])] = indexValue[1]
        #print(self.rightHandSide)
    
    def calculateElementRightHandSide(self,triangleIndex):
        """
        Calculates for the specified element the Righthandside via transformation to the reference triangle.
        Expands the given function in the same linear basis f(x) = Sum(f_i phi_i) and uses the precalculated
        integrals of two basisfunctions multiplied with each other. 
        """
        transformMatrix,translateVector = self.mesh.calculateTransform(triangleIndex)
        determinant = abs(np.linalg.det(transformMatrix))

        trianglePoints =self.triangulation.triangles[triangleIndex]
        elementRHS = determinant*np.dot(self.elementaryBasisMatrix,np.array([self.RHSEvaluatedAtTrianglePoints[x] for x in trianglePoints] ))
        #print("fun",np.array([[x,self.functionRHS(x)] for x in trianglePoints] ))
        #print(self.triangulation.simplices[triangleIndex])
        return elementRHS

    def calculateElementRightHandSideGaussLeg(self,triangleIndex,degree):
        """
        Calculates the Integral int{f phi_j} with transformation formula, and Gauss Legendre quadrature
        integral over a transformed element --> reference element --> standard square element
        Not used for linear elements, but can be used for higher order elements
        """

        transformMatrix,translateVector = self.mesh.calculateTransform(triangleIndex)
        determinant = abs(np.linalg.det(transformMatrix))

        
        trianglePoints =self.triangulation.points[self.triangulation.triangles[triangleIndex]]
        gPoints,gWeigths = leggauss(degree)
        elementRHS = []
        #calculate each entry per Gaussintegration
        for i in range(0,3):
            entry = 0.0
            #sum over weights evaluated at Gauss points,
            for indexX,pointX in enumerate(gPoints):
                for indexY,pointY in enumerate(gPoints):
                    transformedPoint = np.dot(transformMatrix,np.array([(1+pointX)*0.5,(1-pointX)*(1+pointY)*0.25])) +translateVector
                    entry+= determinant*0.5*(1-pointX)*0.125*self.functionRHS(transformedPoint)*gWeigths[indexX]*gWeigths[indexY]*self.mesh.linearBasis[i](np.array([(1+pointX)*0.5,(1-pointX)*(1+pointY)*0.25]))
                    #0.5 comes from transformation of reference triangle to standard square [-1,1] x [-1,1]
            elementRHS.append(entry)
        return elementRHS

    def solve(self):
        self.solution =sparse.linalg.spsolve(self.GlobalStiffness,self.rightHandSide)
        
        

    def getL2Error(self,exactSolution):
        """
        Calculates the L2-Error for given exact solution, it projects the error to an 
        elementwise linear function and integrate this exact via transformation formula
        """
        value = 0
        error = np.array(self.solution)-np.array([exactSolution(x) for x in self.triangulation.points])
        for ele,triPoints in enumerate(self.triangulation.triangles):
            transformMatrix,translateVector = self.mesh.calculateTransform(ele)
            determinant = abs(np.linalg.det(transformMatrix))
            #Last vector is the precalculated integral of the basisfunctions over a reference element
            value+=determinant*np.dot(error[triPoints]**2,np.array([1/6.,1/3.,1/3.]))
        return(math.sqrt(value))
