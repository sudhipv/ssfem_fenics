
# Copyright (C) 2017 Ajit Desai, Ph.D. Candidate ajit.ndesai@gmail.com

################################################################
## AJIT DESAI: Modified Oct-2017
## from terminal go to the folder where you have this code
## from "Docker" start "fenicsproject start myFenics"
## go to the working directory: Eg $cd det-ddm/pwd/
## after that run this code using "python ddm_**.py"
################################################################

#from petsc4py import PETSc
import numpy as np
import scipy.io as sio
from dolfin import *

print("========================================================")
print("Running FEniCS...")

parameters['reorder_dofs_serial'] = False
    
nParts = np.genfromtxt("../../../data/meshData//num_partition.dat",delimiter=',')
nParts = nParts.astype(int)

for i in range(nParts):
    #mesh = Mesh("../inputs/gmsh3D.xml")         ## GMSH converted
    pwd1 = "../data/foo"+str(i+1)+".xml"
    mesh = Mesh(pwd1)
    print("Input Mesh:", pwd1)

    # initialize the connectivity between facets and cells
    tdim = mesh.topology().dim()
    mesh.init(tdim-1, tdim)

    V = FunctionSpace(mesh, "Lagrange", 1)
    uh= Function(V)

    # Mark boundary subdomians (Rectangular beam x(0,1), y(0,1), z(0,4)
    left =  CompiledSubDomain("near(x[2], side) && on_boundary", side = 0.0) #z=0
    right = CompiledSubDomain("near(x[2], side) && on_boundary", side = 4.0) #z=4

    # Define boundary condition
    u0 = Constant(0.0)
    bcl = DirichletBC(V, u0, left)
    bcr = DirichletBC(V, u0, right)
    #bc = [bcl, bcr]  ## fixing both ends
    #bc = bcl         ## Fixing left-side only

    # Test/Trial space & constants
    u, v = TrialFunction(V), TestFunction(V)
    f = Constant(-6.0)
    g = Constant(0.0)

    ## Variation formulation
    # Define forms to be assembled (crucial to turn into Forms in parallel)
    a = Form(inner(grad(u), grad(v))*dx)
    L = Form(f*v*dx + g*v*ds)

    ## FEniCS assembly
    # Dummy problem to define LAYOUT of global problem (the easy way)
    A_g = assemble( Constant(0.)*inner(grad(u), grad(v))*dx )
    f_g = assemble( Constant(0.)*v*dx )

    # Get dofmap to construct cell-to-dof connectivity
    dofmap = V.dofmap()

    # Perform assembling
    for cell in cells(mesh):
        dof_idx = dofmap.cell_dofs(cell.index())
        
        # Assemble local rhs and lhs
        a_local  = assemble_local(a, cell)
        L_local  = assemble_local(L, cell)
        
        # Assemble into global system
        A_g.add_local(a_local,dof_idx, dof_idx)
        f_g.add_local(L_local,dof_idx)

    # Finalize assembling
    A_g.apply("add"), f_g.apply("add")

    # APPLY BC
    bcl.apply(A_g,f_g)
    bcr.apply(A_g,f_g)

    ## To sove assembly FEniCS assembled matrix and vector
    pwd3 = "../data/outputs/Ab"+str(i+1)+".mat"
    sio.savemat(pwd3, {'An':A_g.array(), 'bn':f_g.array()})
    print'Local assembly After-BC  :', pwd3

#    # Compute solution
#    solve(A_g, uh.vector(), f_g)
#
#    # Save solution in VTK format
#    pwd4 ="../outputs/poisson"+str(i+1)+".pvd"
#    File(pwd4) << uh
#    print'Local assembly solution  :', pwd4
print("========================================================")
