
from Mesh import Mesh

import math

import numpy as np




class ProjectedGradient:
    """Class for using a Projected Gradient Mehtod to solve the elliptic Optimiyation Problem
        min J(Sq,q).  
    """

    def __init__(self,mesh,startFunction,alpha= 1, UDesired):
        self.mesh = mesh
        self.startFunction = startFunction
        self.alpha = alpha
        self.UDesired = UDesired
