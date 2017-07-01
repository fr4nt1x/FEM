
from Mesh import Mesh
from FiniteElement import FiniteElement
from numpy.polynomial.legendre import leggauss
import math

import numpy as np




class ProjectedGradient:
    """Class for using a Projected Gradient Mehtod to solve the elliptic Optimiyation Problem
        min J(Sq,q).  
        UDesired, startFunction hold u_d  and q_0 at the mesh points
        
    """

    def __init__(self,mesh,startControl, uDesired, RHSAddendum=0, alpha= 1,tol= 1e-12, maxSteps = 500):
        self.mesh = mesh
        self.startControl = startControl
        self.control = self.startControl
        self.alpha = alpha
        self.uDesired = uDesired
        self.tolerance = tol
        self.state = None
        self.adjointState = None 
        self.maxSteps= maxSteps
        self.RHSAddendum = RHSAddendum
        self.maxDiam = 0
        #referenceElement holds the points of the reference element from which all other elements
        #are calculated
        self.referenceElement = np.array([[0,0],[1.,0],[0,1.]])

        #the 3 linear Basis funtctions on the reference triangle
        #each has the value 1 at one points and 0 at the other points
        #Numbering of the vertices according to self.referenceElement
        self.linearBasis = []
        self.linearBasis.append(lambda x : 1-x[0]-x[1])
        self.linearBasis.append(lambda x : x[0])
        self.linearBasis.append(lambda x : x[1])

    def solveState(self):
        Fem = FiniteElement(self.mesh,PDEMatrix= np.array([[1,0],[0,1]]),RHSEvaluatedAtTrianglePoints = np.array(self.control)+self.RHSAddendum)
        Fem.calculateGlobalStiffnessMatrix()
        Fem.calculateRightHandSide()
        Fem.solve()
        self.state = Fem.solution 

    def solveAdjoint(self):
        Fem = FiniteElement(self.mesh,PDEMatrix= np.array([[1,0],[0,1]]),RHSEvaluatedAtTrianglePoints = self.state - self.uDesired)
        Fem.calculateGlobalStiffnessMatrix()
        Fem.calculateRightHandSide()
        Fem.solve()
        self.adjointState = Fem.solution 

    def calculateOneStep (self):
        self.solveState()
        self.solveAdjoint()
        # print("state",self.state)
        # print("adjoint",self.adjointState)
        self.control = (-1/self.alpha) * self.adjointState

    def solve(self):
        for step in range(0,self.maxSteps):
            oldControl = self.control
            self.calculateOneStep()
            if np.linalg.norm(self.control - oldControl, np.infty) <= self.tolerance:
                print("Tolerance reached in "+ str(step)+ ".")
                break

    def calculateTransform(self,triangleIndex):
        """ Calculates the affine linear transform from the reference
            Element to the given Triangle (Lx+b) with matrix 2x2 L and vector b
            and returns both of them in a tuple.
        """
        
        #get the pointcoordinates via the point indices of the specified element
        trianglePoints = self.mesh.points[self.mesh.triangles.copy()[triangleIndex]]

        #calculate the Diameter and save it if a new maximum arises
        self.maxDiam = max(np.linalg.norm(trianglePoints[0]-trianglePoints[1]),np.linalg.norm(trianglePoints[0]-trianglePoints[2]),np.linalg.norm(trianglePoints[2]-trianglePoints[1]),self.maxDiam)

        #hold the coordinates both elements with ones appended as last row
        referenceCoord = np.array([self.referenceElement[:,0],self.referenceElement[:,1],np.array([1,1,1])])       
        transformedCoord = np.array([trianglePoints[:,0],trianglePoints[:,1],np.array([1,1,1])])
        C = np.dot( transformedCoord,np.linalg.inv(referenceCoord)) 
        
        #L is in the first n x n submatrix of (n+1) x (n+1) Matrix C, b the last column without the last entry 
        return (C[0:-1,0:-1],C[:-1,-1])

    def calculateIntegralOverTriangleGauss(self,f,triangleIndex,degree):
        """
        Calculates the Integral int{f} with transformation formula, and Gauss Legendre quadrature
        integral over a transformed element --> reference element --> standard square element
        Not used for linear elements, but can be used for higher order elements
        """

        transformMatrix,translateVector = self.calculateTransform(triangleIndex)
        determinant = abs(np.linalg.det(transformMatrix))

        
        trianglePoints =self.mesh.points[self.mesh.triangles[triangleIndex]]
        gPoints,gWeigths = leggauss(degree)
        #calculate each entry per Gaussintegration
        entry = 0.0

        #sum over weights evaluated at Gauss points,
        for indexX,pointX in enumerate(gPoints):
            for indexY,pointY in enumerate(gPoints):
                transformedPoint = np.dot(transformMatrix,np.array([(1+pointX)*0.5,(1-pointX)*(1+pointY)*0.25])) +translateVector
                entry+= determinant*0.5*(1-pointX)*0.125*f(transformedPoint)*gWeigths[indexX]*gWeigths[indexY]
                #0.5 comes from transformation of reference triangle to standard square [-1,1] x [-1,1]
        return entry

    def getL2ErrorControl(self,exactSolution):
        """
        Calculates the L2-Error for given exact solution, it projects the error to an 
        elementwise linear function and integrate this exact via transformation formula
        """
        value = 0
        error = np.array(self.control)-np.array([exactSolution(x) for x in self.mesh.points])
        for ele,triPoints in enumerate(self.mesh.triangles):
            transformMatrix,translateVector = self.calculateTransform(ele)
            determinant = abs(np.linalg.det(transformMatrix))
            #Last vector is the precalculated integral of the basisfunctions over a reference element
            value+=determinant*np.dot(error[triPoints]**2,np.array([1/6.,1/3.,1/3.]))
        return(math.sqrt(value))

    def getL2ErrorControlGauss(self,exactSolution):
        """
        Calculates the L2-Error for given exact solution,using Gauss Legendre quadrature
        """
        value = 0
        for ele,triPoints in enumerate(self.mesh.triangles):
            transformMatrix,translateVector = self.calculateTransform(ele)
            invTransformMatrix = np.linalg.inv(transformMatrix)
            determinant = abs(np.linalg.det(transformMatrix))
            #Last vector is the precalculated integral of the basisfunctions over a reference element
            def error(x):
                #points at the transformed triangle as call parameter, need to be transformed to the basis triangle
                xRef = np.dot(invTransformMatrix,x-translateVector)
                return (self.control[triPoints[0]]*self.linearBasis[0](xRef)+ self.control[triPoints[1]]*self.linearBasis[1](xRef)+self.control[triPoints[2]]*self.linearBasis[2](xRef)-exactSolution(x))**2
                

            value+=self.calculateIntegralOverTriangleGauss(error,ele,9)
    
        return(math.sqrt(value))
