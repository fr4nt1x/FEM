
from Mesh import Mesh
from FiniteElement import FiniteElement

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

