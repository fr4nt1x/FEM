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
        # print("TriangleIndex: ",triangleIndex)
        determinant = abs(np.linalg.det(transformMatrix))
        # print(transformMatrix)
        
        trianglePoints =self.mesh.points[self.mesh.triangles[triangleIndex]]
        # print('trianglge', trianglePoints)
        gPointsTriangle,gPointsSquare,gWeigths = self.getPointsAndWeightsOverReferenceTriangle(degree)
        #calculate each entry per Gaussintegration
        entry = 0.0

        #sum over weights evaluated at Gauss points,
        for index,point in enumerate(gPointsTriangle):
            # print("PW",point,gWeigths[index])
            transformedPoint = np.dot(transformMatrix,point) +translateVector
            
            #The term 1-gPointsSquare need to be a point of the suare gauss points. this comes from the Transform of the trinangle to the square
            entry+= determinant*(1-gPointsSquare[index][0])*0.125*functionToIntegrate(transformedPoint)*gWeigths[index]
            # print("entri",entry)
        return entry

    def getPointsAndWeightsOverReferenceTriangle(self,degree):
        """
        Returns 3 Numpy Arrays, One with the Points inside the reference Triangle, 
        one with the points inside the square [-1,1]x[-1,1].
        The other holds the weights at this points
        """

        gPoints,gWeights = leggauss(degree)
        numberOfPoints = len(gPoints)**2
        resultPointsTriangle = np.zeros((numberOfPoints,2))
        resultPointsSquare = np.zeros((numberOfPoints,2))
        resultWeights = np.zeros((numberOfPoints,1))
        index2D = 0
        for indexX,pointX in enumerate(gPoints):
            for indexY,pointY in enumerate(gPoints):
                resultWeights[index2D] = gWeights[indexX]*gWeights[indexY] 
                resultPointsTriangle[index2D,:] = np.array([(1+pointX)*0.5,(1-pointX)*(1+pointY)*0.25])
                resultPointsSquare[index2D,:] = np.array([pointX,pointY])
                index2D += 1
        return resultPointsTriangle,resultPointsSquare,resultWeights

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
            # print('Values: ',valuesAtMeshPoints[triPoints])

            #Last vector is the precalculated integral of the basisfunctions over a reference element
            def func(x):
                #points at the transformed triangle as call parameter, need to be transformed to the basis triangle
                xRef = np.dot(invTransformMatrix,x-translateVector)
                return valuesAtMeshPoints[triPoints[0]]*self.mesh.linearBasis[0](xRef)+ valuesAtMeshPoints[triPoints[1]]*self.mesh.linearBasis[1](xRef)+valuesAtMeshPoints[triPoints[2]]*self.mesh.linearBasis[2](xRef)
            return func
