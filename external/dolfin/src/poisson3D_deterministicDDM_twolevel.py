
################################################################
#################### SUDHI P V ################################
################################################################

#from petsc4py import PETSc
import numpy as np
import scipy.io as sio
import os as os
from dolfin import *

print("===================================================")
print("Running FEniCS for 3D stochastic Poisson/twoLevel..")

## For avoiding FEniCS internal node reordering optimized for speed
parameters['reorder_dofs_serial'] = False

##---------------------------------------------------------------------------
## Global deterministic assembly for each KLE mode
def detAssembly(a, L):
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

        # Define Dirichlet boundary (x = 0 or x = 1)
    left =  CompiledSubDomain("near(x[2], side) && on_boundary", side = 0.0)
    right = CompiledSubDomain("near(x[2], side) && on_boundary", side = 1.0)

        # Define boundary condition
    u0 = Constant(0.0)
    bcl = DirichletBC(V, u0, left)
    bcr = DirichletBC(V, u0, right)
    bcl.apply(A_g,f_g)
    bcr.apply(A_g,f_g)

    #A = np.zeros([V.dim(), V.dim()])
    A = A_g.array()
    b = f_g.get_local()
    #A = A_g.mat()
    #Adense = as_backend_type(A_g).mat().convert('dense')
    #A = A_dense.getDenseArray()
    return A, b

##---------------------------------------------------------------------------
## OneLevel Deterministic assembly matrix decomposer

def OneLevelDDMat(A, b, nNodes, nBnodes):
    ADii = A[0:(nNodes-nBnodes),0:(nNodes-nBnodes)]
    ADig = A[0:(nNodes-nBnodes),(nNodes-nBnodes):nNodes]
    #A[(nNodes-nBnodes):nNodes,0:(nNodes-nBnodes)]   ## A[0:(nNodes-nBnodes),(nNodes-nBnodes):nNodes]
    ADgg = A[(nNodes-nBnodes):nNodes,(nNodes-nBnodes):nNodes]
    bi   = b[0:(nNodes-nBnodes)]
    bg   = b[(nNodes-nBnodes):nNodes]
    return ADii, ADig, ADgg, bi, bg

##---------------------------------------------------------------------------
## OneLevel Deterministic assembly matrix decomposer
def TwoLevelDDMat(A, b, nNodes, nBnodes, nCnodes):
    #    ADii = A[0:(nNodes-nBnodes),0:(nNodes-nBnodes)]
    #    ADgi = A[0:(nNodes-nBnodes),(nNodes-nBnodes):nNodes]
    #    ADgg = A[(nNodes-nBnodes):nNodes,(nNodes-nBnodes):nNodes]

    ADic = A[0:(nNodes-nBnodes),(nNodes-nCnodes):nNodes]
    ADir = A[0:(nNodes-nBnodes),(nNodes-nBnodes):(nNodes-nCnodes)]
    ADcc = A[(nNodes-nCnodes):nNodes,(nNodes-nCnodes):nNodes]
    ADrr = A[(nNodes-nBnodes):(nNodes-nCnodes),(nNodes-nBnodes):(nNodes-nCnodes)]
    ADrc = A[(nNodes-nBnodes):(nNodes-nCnodes),(nNodes-nCnodes):nNodes]

    #    bi   = b[0:(nNodes-nBnodes)]
    #    bg   = b[(nNodes-nBnodes):nNodes]
    #    bc   = b[(nNodes-nCnodes):nNodes]
    #    br   = b[(nNodes-nBnodes):(nNodes-nCnodes)]
    return ADic, ADir, ADcc, ADrr, ADrc   # ADii, ADgi, ADgg, bi, bg  #, br, bc

##---------------------------------------------------------------------------
#################### MAIN CODE #####################
## Path to the global meshdims.dat
meshdimpath = '../../../data/meshData/meshdim.dat'
nParts = np.genfromtxt(meshdimpath)
nParts = nParts[5].astype(int)   ## For 3D
print("number of partitions:", nParts)

PDE_Dim = [1]
pdeDimPath = '../../../data/meshData/PDE_Dim.txt'
np.savetxt(pdeDimPath,PDE_Dim,fmt='%.d')

##---------------------------------------------------------------------------
## Fore each subdomain invokes assembly procedure for deterministic DD blocks
for ip in range(nParts):

    ## Get DDM mesh dimension data  ## NOTE: Zfill=4 for both Matlab & Fortran DDM data
    mdpath = "../../../data/meshData/meshdim"+str(ip+1).zfill(4)+".dat"
    meshdim = np.genfromtxt(mdpath)
    nNodes = meshdim[0].astype(int)     ## For 3D
    nBnodes = meshdim[4].astype(int)    ## For 3D

    mcrdpath = "../../../data/meshData/dimcrn"+str(ip+1).zfill(4)+".dat"
    meshdim = np.genfromtxt(mcrdpath)
    nCnodes = meshdim[0].astype(int)   ## For 2D
    nRnodes = meshdim[1].astype(int)   ## For 2D

    ## Create mesh and define function space
    mpath = "../data/foo"+str(ip+1)+".xml"   ## path to mesh files
    print("Input XML-Mesh:", mpath)
    mesh = Mesh(mpath)

    # Initialize the connectivity between facets and cells
    tdim = mesh.topology().dim()
    mesh.init(tdim-1, tdim)

    # Functional/Test/Trial space & constants
    V = FunctionSpace(mesh, "Lagrange", 1)
    u, v = TrialFunction(V), TestFunction(V)
    f = Constant(1.0)
    g = Constant(0.0)


    a = Form(inner(grad(u), grad(v))*dx)
    L = Form(f*v*dx + g*v*ds)

        ## Save deterministic assembly matrices in ".mat" or ".txt" format
    subpath = "../data/Amats/subdom000".rstrip("0")+str(ip+1)

    if not os.path.exists(subpath):
        os.makedirs(subpath)


        ## Invoke Deterministic Assembly procedure for each KLE/PCE mode
    A, b = detAssembly(a, L) ## New assembly procedure

        ##if k==0:b = b0  ## Decompose assembly matrices intor DD-blocks
    ADii, ADig, ADgg, bi, bg = OneLevelDDMat(A, b, nNodes, nBnodes)
    ADic, ADir, ADcc, ADrr, ADrc = TwoLevelDDMat(A, b, nNodes, nBnodes, nCnodes)

        #dpath = subpath+"/AD"+str(k+1)+".dat",  np.savetxt(dpath, A)
        #print("Out Deterministic Assembly:", dpath)
        ## To save output matrices in ".dat" for Fortran
    k = 0
    ADiipath = subpath+"/ADii000".rstrip("0")+str(k+1)+".dat"
    np.savetxt(ADiipath, ADii)
    ADigpath = subpath+"/ADig000".rstrip("0")+str(k+1)+".dat"
    np.savetxt(ADigpath, ADig)
    ADggpath = subpath+"/ADgg000".rstrip("0")+str(k+1)+".dat"
    np.savetxt(ADggpath, ADgg)
    fipath = subpath+"/bi000".rstrip("0")+str(k+1)+".dat"
    np.savetxt(fipath, bi)
    fgpath = subpath+"/bg000".rstrip("0")+str(k+1)+".dat"
    np.savetxt(fgpath, bg)

    ADicpath = subpath+"/ADic000".rstrip("0")+str(k+1)+".dat"
    np.savetxt(ADicpath, ADic)
    ADirpath = subpath+"/ADir000".rstrip("0")+str(k+1)+".dat"
    np.savetxt(ADirpath, ADir)
    ADccpath = subpath+"/ADcc000".rstrip("0")+str(k+1)+".dat"
    np.savetxt(ADccpath, ADcc)
    ADrrpath = subpath+"/ADrr000".rstrip("0")+str(k+1)+".dat"
    np.savetxt(ADrrpath, ADrr)
    ADrcpath = subpath+"/ADrc000".rstrip("0")+str(k+1)+".dat"
    np.savetxt(ADrcpath, ADrc)

    if (ip == 0):
        print("output Mat-Vec:",subpath)


print("==========================Success============================")
###---------------------------------------------------------------------------
##### TO ASSEMBLE & SAVE FULL STOCHATIC MATRIX FOR ENTIRE DOMAIN
#indexi = 0
#N = V.dim()              ## deterministic matrix size
#Ns = (PCE_u*N)           ## stochastic matrix size
#As = np.zeros([Ns,Ns])   ## number of elements
#
###---------------------------------------------------------------------------
### Fore each KLE mode invokes assembly procedure to build stochastic matrices
#for k in range(PCE_A):
#    ## Call stochastic vartional form for each input PCE mode
#    params = ([k, dim_L])
#    a, L = stoVartional(params, u, v, g, f)  ## new
#    A, b = detAssembly(a, L, k) ## new
#    #if k==0:b = b0
#
#    for i in range(PCE_u):
#        for j in range(PCE_u):
#            if k==ijk[indexi,0] and i==ijk[indexi,1] and j==ijk[indexi,2]:
#                As[i*N:(i+1)*N, j*N:(j+1)*N] += cijk[indexi]*A
#                indexi = indexi+1
#
### Save assembly matrices in ".mat" or ".txt" format
#pwd1 = '../data/outputs/As.mat'
#sio.savemat(pwd1, {'As':As, 'b':b})
#print('Stochastic Mat-Vec: New', pwd1)
###---------------------------------------------------------------------------
##################################
# Compute solution
# solve(A_g, uh.vector(), f_g)

## Save solution in VTK format
#pwd4 ='../outputs/poisson.pvd'
#File(pwd4) << uh
#print'Local assembly solution  :', pwd4
