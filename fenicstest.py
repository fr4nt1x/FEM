from fenics import *
error= []
for i in range(5,19,1):
    mesh = UnitSquareMesh(i-1,i-1)
    V=FunctionSpace(mesh,'P',1)

    u_D = Expression('x[0]*x[0]+x[1]*x[1]',degree=2)

    def boundary(x,on_boundary):
        return on_boundary

    bc = DirichletBC(V,u_D,boundary)


    u=TrialFunction(V)
    v=TestFunction(V)
    f= Constant(-4)
    a= dot(grad(u),grad(v))*dx
    L= f*v*dx

    u=Function(V)
    solve(a==L,u,bc)

    error_L2=errornorm(u_D,u,'L2')
    error.append(error_L2)
print(error)

