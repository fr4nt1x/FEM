import numpy as np
from FiniteElement import FiniteElement




class ProjectionOnR():
    """
    L^2 (Omega) can be split into two orthogonal Subspaces (for non-convex domains):
    First Subspace: 
    The image of the Laplace-Operator for H^2 functions (functions that have weak 
    weak derivatives up to order 2 and are zero (in the trace sense) at the domain 
    boundary)  
    Second Subspace:
    Dimension for Dirichlet Boundary Conditions is equal to the number of Non-Convex 
    Corners. This Case (One Corner at [0,0] , interior angle phi)
    span(ps) : ps := p_tilde + r^(-angleCoefficient)*sin(angleCoefficient*theta) , where
                r .... distance to corner
                theta ... angle in polarcoordinates
                angleCoefficient ... pi/phi
                p_tilde ... H^1 function for orthogonality of the subspaces

    This Class calculates the Projection onto R for a given Function.
    The Given function should be a FEM Function on a mesh.
    For discretizing p_s  the idea is to discretize p_tilde via: 
    p_h_tilde = p*_h - r_h , where
    r_h .... Finite Element Function for the given Mesh, where only Values at the 
            boundary are given. These are calcutalted via r^(-theta)*sin(lambda*theta)
            where this is defined and 0 for the non convex vertex (function is zero 
            at the two nonconvex edges of the boundary)
    p*_h ... is chosen such that p*_h - r_h fullfils -Laplace p*_h - r_h = 0. in a 
             weak sense with p*_h is zero at boundary.
    """

    def __init__(self,angleCoefficient,mesh,indexOfNonConvexCorner):
        """
        """
        self.mesh = mesh
        self.numberOfPoints = np.shape(self.mesh.points)[0]
        self.indexOfNonConvexCorner= indexOfNonConvexCorner
        self.angleCoefficient = angleCoefficient
        self.r_h = np.zeros(self.numberOfPoints)
        self.p_h_tilde = np.zeros(self.numberOfPoints)
        self.p_h_star = np.zeros(self.numberOfPoints)

    def calculateR_h(self):
        """
        """
        boundaryIndices = [x[0] for x in self.mesh.boundaryValues if x[0] != self.indexOfNonConvexCorner ]
        points = self.mesh.points[boundaryIndices] 

        theta =  np.arctan2(points[:,1],points[:,0])
        
        #negative Angles should be converted to positive Values to match the angleCOoefficient
        theta = np.array([x if x >=0 else 2*np.pi + x for x in theta])
        r = np.linalg.norm(points,axis=1)
        self.r_h[boundaryIndices] = r**(-self.angleCoefficient) * np.sin(self.angleCoefficient * theta)

    def calculatePStar(self):
        """
        """

        Fem = FiniteElement(self.mesh , PDEMatrix= np.array([[1,0],[0,1]]),functionRHS = lambda x : 0)
        Fem.calculateGlobalStiffnessMatrix(calculateWholeMatrix = True)
        Fem.rightHandSide = Fem.GlobalStiffness.dot(self.r_h) 
        print("RHS: ")
        print(Fem.rightHandSide)
        Fem.overwriteBoundaryValuesRightHandSide()
        print("RHS2: ")
        print(Fem.rightHandSide)

    def calculatePTilde(self):
        """
        """

