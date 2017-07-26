
from Mesh import Mesh
from FiniteElement import FiniteElement
from numpy.polynomial.legendre import leggauss
import json
import os
from gaussIntegration import GaussIntegrator
from ProjectionOnR import ProjectionOnR
import math

import numpy as np




class ProjectedGradient:
    """Class for using a Projected Gradient Mehtod to solve the elliptic Optimiyation Problem
        min J(Sq,q).  
        UDesired, startFunction hold u_d  and q_0 at the mesh points
        
    """

    def __init__(self,mesh,startControl, uDesired, angleCoefficient, RHSAddendum=0, alpha= 1,tol= 1e-12, maxSteps = 500):
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
        self.integrator = GaussIntegrator(mesh)
        self.ProjectorOnR = ProjectionOnR(angleCoefficient = angleCoefficient,mesh= self.mesh, indexOfNonConvexCorner = 0,functionValuesToProject = self.control)
        self.ProjectorOnR.calculateR_h()
        self.ProjectorOnR.calculatePStar()
        self.ProjectorOnR.calculatePTilde()
        self.ProjectorOnR.calculateNormPTildeSquared()
        print("Norm:",self.ProjectorOnR.normPTildeSquared)

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
        self.ProjectorOnR.functionValuesToProject = (-1/self.alpha)* self.adjointState
        # self.control =  self.ProjectorOnR.getProjectionOnR()
        self.control = (-1/self.alpha)* self.adjointState

    def dumpToJson(self,FolderName,SuffixName):

        with open(os.path.join(FolderName,"State_"+SuffixName),'w') as fileState:
            json.dump(self.state.tolist(),fileState)

        with open(os.path.join(FolderName,"Control_"+SuffixName),'w') as fileControl:
            json.dump(self.control.tolist(),fileControl)

    def loadFromJson(self,FolderName,SuffixName):

        state = None
        control = None
        with open(os.path.join(FolderName,"State_"+SuffixName),'r') as fileState:
            state= np.array(json.load(fileState))

        with open(os.path.join(FolderName,"Control_"+SuffixName),'r') as fileControl:
            control= np.array(json.load(fileControl))

        self.referenceControl = control
        self.referenceState =state 

    def solve(self):
        for step in range(0,self.maxSteps):
            oldControl = self.control
            self.calculateOneStep()
            if np.linalg.norm(self.control - oldControl, np.infty) <= self.tolerance:
                print("Tolerance reached in "+ str(step)+ ".")
                break

    def getL2ErrorControl(self,exactSolutionEvaluatedAtPoints):
        """
        Calculates the L2-Error for given exact solution, it projects the error to an 
        elementwise linear function and integrate this exact via transformation formula
        """
        value = 0
        newMesh = Mesh(self.mesh.points,self.mesh.triangles,self.mesh.edges,self.mesh.polygonalBoundary) 
        newMesh.setValuesAtPoints(self.control)
        print(np.shape(newMesh.points))
        while (np.shape(newMesh.points)[0] != np.shape(exactSolutionEvaluatedAtPoints)[0]):
            newMesh.refineMesh(1)
            print("refinement ",np.shape(newMesh.points))
        error = exactSolutionEvaluatedAtPoints-newMesh.valuesAtPoints
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
            #Last vector is the precalculated integral of the basisfunctions over a reference element
            def error(x):
                #points at the transformed triangle as call parameter, need to be transformed to the basis triangle
                xRef = np.dot(invTransformMatrix,x-translateVector)
                return (self.control[triPoints[0]]*self.mesh.linearBasis[0](xRef)+ self.control[triPoints[1]]*self.mesh.linearBasis[1](xRef)+self.control[triPoints[2]]*self.mesh.linearBasis[2](xRef)-exactSolution(x))**2
                

            value+=self.integrator.getIntegralOverTriangleGauss(error,ele,9)
    
        return(math.sqrt(value))
