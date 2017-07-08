
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

    def calculateIntegralOverTriangleGauss(self,f,triangleIndex,degree):
        """
        Calculates the Integral int{f} with transformation formula, and Gauss Legendre quadrature
        integral over a transformed element --> reference element --> standard square element
        Not used for linear elements, but can be used for higher order elements
        """

        transformMatrix,translateVector = self.mesh.calculateTransform(triangleIndex)
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
            transformMatrix,translateVector = self.mesh.calculateTransform(ele)
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
            transformMatrix,translateVector = self.mesh.calculateTransform(ele)
            invTransformMatrix = np.linalg.inv(transformMatrix)
            determinant = abs(np.linalg.det(transformMatrix))
            #Last vector is the precalculated integral of the basisfunctions over a reference element
            def error(x):
                #points at the transformed triangle as call parameter, need to be transformed to the basis triangle
                xRef = np.dot(invTransformMatrix,x-translateVector)
                return (self.control[triPoints[0]]*self.mesh.linearBasis[0](xRef)+ self.control[triPoints[1]]*self.mesh.linearBasis[1](xRef)+self.control[triPoints[2]]*self.mesh.linearBasis[2](xRef)-exactSolution(x))**2
                

            value+=self.calculateIntegralOverTriangleGauss(error,ele,9)
    
        return(math.sqrt(value))
