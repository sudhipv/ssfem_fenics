## Modified from FEniCS demo:
#
# Copyright (C) 2017 Ajit Desai, PhD Candidate, ajit.ndesai@gmail.com
#
# Code to test s simple FEniCS code

import numpy as np
import os as os
from dolfin import *

print("**************** Running FEniCS... *****************")

# Create mesh and define function space
mesh = UnitSquareMesh(2, 2)
V = FunctionSpace(mesh, "Lagrange", 1)

# Define Dirichlet boundary (x = 0 or x = 1)
def boundary(x):
    return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS

# Define boundary condition
u0 = Constant(0.0)
bc = DirichletBC(V, u0, boundary)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(1.0)
g = Constant(0.0)
a = inner(grad(u), grad(v))*dx
L = f*v*dx + g*v*ds

# Compute solution
u = Function(V)
solve(a == L, u, bc)
nodal_values = u.vector().get_local()
print(nodal_values)

# Save solution in VTK format
file = File("demo_poisson.pvd")
file << u

# Plot solution
##plot(u, interactive=True)

print("========================================================")
print("======================== Success =======================")
print("========================================================")
