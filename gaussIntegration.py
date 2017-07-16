import numpy as np
from numpy.polynomial.legendre import leggauss 

class GaussIntegrator():
    """
    Used for Calculating integrals with Gauss integration formula over the given Mesh
    """
    def __init__(self, mesh):
        self.mesh = mesh

    def getIntegralOverTriangleGauss(self,functionToIntegrate,triangleIndex,degree):
        """
        Calculates the Integral int{f} with transformation formula, and Gauss Legendre quadrature
        integral over a transformed element --> reference element --> standard square element
        """

        transformMatrix,translateVector = self.mesh.calculateTransform(triangleIndex)
        print("TriangleIndex: ",triangleIndex)
        determinant = abs(np.linalg.det(transformMatrix))
        print(transformMatrix)
        
        trianglePoints =self.mesh.points[self.mesh.triangles[triangleIndex]]
        print('trianglge', trianglePoints)
        gPoints,gWeigths = leggauss(degree)
        #calculate each entry per Gaussintegration
        entry = 0.0

        #sum over weights evaluated at Gauss points,
        for indexX,pointX in enumerate(gPoints):
            for indexY,pointY in enumerate(gPoints):
                transformedPoint = np.dot(transformMatrix,np.array([(1+pointX)*0.5,(1-pointX)*(1+pointY)*0.25])) +translateVector
                #TODO 0.5 not sure if right
                print("fuctionValue",transformedPoint,functionToIntegrate(transformedPoint),pointX,pointY)
                entry+= determinant*(1-pointX)*0.125*functionToIntegrate(transformedPoint)*gWeigths[indexX]*gWeigths[indexY]
                #0.5 comes from transformation of reference triangle to standard square [-1,1] x [-1,1]
        return entry

    def getIntegralOverDomain(self,functionToIntegrate,degree): 
        """
        """
        value = 0
        for ele, triPoints in enumerate(self.mesh.triangles):
            value += self.getIntegralOverTriangleGauss(functionToIntegrate, ele, degree)

        return value 

    def getFiniteElementFunctionOverTriangle(self,valuesAtMeshPoints,element):
            transformMatrix,translateVector = self.mesh.calculateTransform(element)
            invTransformMatrix = np.linalg.inv(transformMatrix)
            triPoints = self.mesh.triangles[element]
            print('Values: ',valuesAtMeshPoints[triPoints])

            #Last vector is the precalculated integral of the basisfunctions over a reference element
            def func(x):
                #points at the transformed triangle as call parameter, need to be transformed to the basis triangle
                xRef = np.dot(invTransformMatrix,x-translateVector)
                return valuesAtMeshPoints[triPoints[0]]*self.mesh.linearBasis[0](xRef)+ valuesAtMeshPoints[triPoints[1]]*self.mesh.linearBasis[1](xRef)+valuesAtMeshPoints[triPoints[2]]*self.mesh.linearBasis[2](xRef)
            return func
