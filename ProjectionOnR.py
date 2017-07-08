imort numpy as np




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
    span(ps) : ps := p_tilde + r^(-lambda)*sin(lambda*theta) , where
                r .... distance to corner
                theta ... angle in polarcoordinates
                lambda ... pi/phi
                p_tilde ... H^1 function for orthogonality of the subspaces

    This Class calculates the Projection onto R for a given Function.
    The Given function should be a FEM Function on a mesh.
    For discretizing p_s  the idea is to discretize p_tilde via: 
    p_h_tilde = p*_h - r_h , where
    r_h .... Finite Element Function for the given Mesh, where only Values at the 
            boundary are given. These are calcutalted via r^(-lambda)*sin(lambda*theta)
            where this is defined and 0 for the non convex vertex (function is zero 
            at the two nonconvex edges of the boundary)
    p*_h ... 
    """
