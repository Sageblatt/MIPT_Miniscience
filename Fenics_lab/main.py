from fenics import *
from ufl import nabla_grad
from ufl import nabla_div
import mshr
from tr_prism import *

# Scaled variables
FLAG_PRISM = True # or cylinder
Len = 60
W = 0.3
mu = 1
rho = 1
delta = W / Len
gamma = 0.4 * delta ** 2
beta = 1.25
lambda_ = beta
g = gamma

# Create mesh and define function space
if FLAG_PRISM:
    create_off_file(Len, 2*W)
    triangular_prism = mshr.Surface3D("base2.off")
    mesh = mshr.generate_mesh(triangular_prism, 128)
else:
    cyl = mshr.Cylinder(Point(0, 0, 0), Point(Len, 0, 0), W, W, 256)
    mesh = mshr.generate_mesh(cyl, 256)

accel = 20*g

V = VectorFunctionSpace(mesh, 'P', 1)

# Define boundary condition
tol = 1E-14


def clamped_boundary_left_right(x, on_boundary):
    return on_boundary and (x[0] < tol or x[0] > Len - tol)


def clamped_boundary_center(x, on_boundary):
    return on_boundary and (Len/2 - 0.5 < x[0] < Len/2 + 0.5)


bc1 = DirichletBC(V, Constant((0, 0, 0)), clamped_boundary_left_right)
bc2 = DirichletBC(V, Constant((0, 0, 0.5*W)), clamped_boundary_center)
bc = [bc1, bc2]

# Define strain and stress


def epsilon(func):
    return 0.5 * (nabla_grad(func) + nabla_grad(func).T)
    # return sym(nabla_grad(u))


def sigma(func):
    return lambda_ * nabla_div(func) * Identity(d) + 2 * mu * epsilon(func)


# Define variational problem
u = TrialFunction(V)
d = u.geometric_dimension()  # space dimension
v = TestFunction(V)
f = Constant((0, 0, -rho * accel))
T = Constant((0, 0, 0))
a = inner(sigma(u), epsilon(v)) * dx
L = dot(f, v) * dx + dot(T, v) * ds

# Compute solution
u = Function(V)
solve(a == L, u, bcs=bc)

# Plot stress
s = sigma(u) - (1. / 3) * tr(sigma(u)) * Identity(d)  # deviatoric stress
von_Mises = sqrt(3. / 2 * inner(s, s))
V = FunctionSpace(mesh, 'P', 1)
von_Mises = project(von_Mises, V)

# Compute magnitude of displacement
u_magnitude = sqrt(dot(u, u))
u_magnitude = project(u_magnitude, V)

# Save solution to file in VTK format
if FLAG_PRISM:
    File('Bending/displacement_tr.pvd') << u
    File('Bending/von_mises_tr.pvd') << von_Mises
    File('Bending/magnitude_tr.pvd') << u_magnitude
else:
    File('Bending/displacement_cyl.pvd') << u
    File('Bending/von_mises_cyl.pvd') << von_Mises
    File('Bending/magnitude_cyl.pvd') << u_magnitude
