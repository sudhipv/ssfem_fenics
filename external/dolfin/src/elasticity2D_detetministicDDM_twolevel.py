
# Copyright (C) 2017 Ajit Desai, Ph.D. Candidate ajit.ndesai@gmail.com

################################################################
## AJIT DESAI: Modified Mar-2018
## from terminal go to the folder where you have this code
## from "Docker" start "fenicsproject start dolfin"
## after that run this code using "python ddm_**.py"
################################################################

#from petsc4py import PETSc
import numpy as np
import scipy.io as sio
import os as os
from fenics import *
import ufl

print("==========================================================")
print("Running FEniCS for 2D deterministic elasticity/oneLevel...")

## For avoiding FEniCS internal node reordering optimized for speed
parameters['reorder_dofs_serial'] = False

##---------------------------------------------------------------------------
## One-level deterministic assembly matrix decomposer
def OneLevelDDMat(A, b, nP, nB, nComp):
    nI = nP - nB
    bi = np.zeros([nI*nComp])
    bg = np.zeros([nB*nComp])
    ADii = np.zeros([nI*nComp,nI*nComp])
    ADgg = np.zeros([nB*nComp,nB*nComp])
    ADgi = np.zeros([nI*nComp,nB*nComp])
    
    for j in range(nComp):
        bi[j*nI:(j+1)*nI] = b[nP*j:(nP*j)+nI]
        bg[j*nB:(j+1)*nB] = b[(nP*j)+nI:(nP*j)+nP]
        for k in range(nComp):
            ADii[j*nI:(j+1)*nI, k*nI:(k+1)*nI] = A[nP*j:(nP*j)+nI, nP*k:(nP*k)+nI]
            ADgg[j*nB:(j+1)*nB, k*nB:(k+1)*nB] = A[(nP*j)+nI:(nP*j)+nP, (nP*k)+nI:(nP*k)+nP];
            ADgi[j*nI:(j+1)*nI, k*nB:(k+1)*nB] = A[nP*j:(nP*j)+nI, (nP*k)+nI:(nP*k)+nP];

    return  ADii, ADgi, ADgg, bi, bg

##---------------------------------------------------------------------------
## Two-level deterministic assembly matrix decomposer
def TwoLevelDDMat(A, b, nP, nB, nC, nComp):
    nI = nP - nB
    nR = nB - nC
    nIR = nP - nC
    
    ADci = np.zeros([nI*nComp,nC*nComp])
    ADri = np.zeros([nI*nComp,nR*nComp])
    ADcc = np.zeros([nC*nComp,nC*nComp])
    ADrr = np.zeros([nR*nComp,nR*nComp])
    ADcr = np.zeros([nR*nComp,nC*nComp])
    
    for j in range(nComp):
        for k in range(nComp):
            ADri[j*nI:(j+1)*nI, k*nR:(k+1)*nR] = A[nP*j:(nP*j)+nI, (nP*k)+nI:(nP*k)+nIR]
            ADci[j*nI:(j+1)*nI, k*nC:(k+1)*nC] = A[nP*j:(nP*j)+nI, (nP*k)+nIR:(nP*k)+nP]
            ADrr[j*nR:(j+1)*nR, k*nR:(k+1)*nR] = A[(nP*j)+nI:(nP*j)+nIR, (nP*k)+nI:(nP*k)+nIR]
            ADcc[j*nC:(j+1)*nC, k*nC:(k+1)*nC] = A[(nP*j)+nIR:(nP*j)+nP, (nP*k)+nIR:(nP*k)+nP]
            ADcr[j*nR:(j+1)*nR, k*nC:(k+1)*nC] = A[(nP*j)+nI:(nP*j)+nIR, (nP*k)+nIR:(nP*k)+nP]

    return ADci, ADri, ADcc, ADrr, ADcr

##---------------------------------------------------------------------------
## Global deterministic assembly for each KLE mode (here it's only mean term)
def detAssembly(a, L):
    ## Dummy problem to define LAYOUT of global problem (the easy way)
    A_g = assemble( Constant(0.)*inner(sigma(u), epsilon(v))*dx)
    f_g = assemble( Constant(0.)*(dot(f, v)*dx + dot(T, v)*ds))
    
    ## Get dofmap to construct cell-to-dof connectivity
    dofmap = V.dofmap()
    
    ## Perform assembling
    for cell in cells(mesh):
        dof_idx = dofmap.cell_dofs(cell.index())
        
        ## Assemble local rhs and lhs
        a_local  = assemble_local(a, cell)
        L_local  = assemble_local(L, cell)
        
        ## Assemble into global system
        A_g.add_local(a_local,dof_idx, dof_idx)
        f_g.add_local(L_local,dof_idx)

    ## Finalize assembling
    A_g.apply("add"), f_g.apply("add")
    
    ## Define boundary conditions
    tol = 1E-14   ## Use this kind of tol for assembly as well
    def clamped_boundary(x, on_boundary):
        return on_boundary and x[1] < tol
    
    bc = DirichletBC(V, Constant((0, 0)), clamped_boundary)
    
    ## Apply boundary conditions
    bc.apply(A_g,f_g)

    ## A = np.zeros([V.dim(), V.dim()])
    A = A_g.array()
    b = f_g.get_local()
    ##A = A_g.mat()
    ##Adense = as_backend_type(A_g).mat().convert('dense')
    ##A = A_dense.getDenseArray()
    
    return A, b


##---------------------------------------------------------------------------
################################ MAIN CODE ##################################
## Scaled variables: System parameters (similar to FEniCS demo_elasticity.py)
# Lenth = 1; Width = 0.2
# mu = 1              ## Scaled to 1
# beta = 1.25         ## Scaled to 1.25
# lambda_ = beta      ## Scaled to 1.25

rho = 1
delta = 0.2         ## Width(1)/Length(5)
g = 0.04*delta**2   ## Use 0.04 for beam(1,1,5) and 0.4 for beam(0.2,0.2,1)
nu = 0.2778         ## Calculated
E = 2.5556
mu = 1/(2.0*(1.0+nu))
lambda_ = (1*nu)/((1+nu)*(1-(2.0*nu)))

## Actual variables
##L = 5; W = 1
#E = 200e9    #PA
#rho = 7860   #(Kg/m3) ## Aplified to increase deflaction Steel_rho=7860 Kg/m3
#g  = 9.8     #(m/s2)
#nu = 0.3
#mu = 1/(2.0*(1.0+nu))
#lambda_ = (1*nu)/((1+nu)*(1-(2.0*nu)))

nComp = 2 # 3 for 3D, 2 for 2D

## Path to the global meshdims.dat
meshdimpath = '../../../data/meshData/meshdim.dat'
nParts = np.genfromtxt(meshdimpath)
nParts = nParts[4].astype(int)     ## For 2D
print("number of partitions:", nParts)

##---------------------------------------------------------------------------
## Fore each subdomain invokes assembly procedure for deterministic DD blocks
for ip in range(nParts):

    ## Get DDM mesh dimension data  ## NOTE: Zfill=4 for both Matlab & Fortran DDM data
    mdpath = "../../../data/meshData/meshdim"+str(ip+1).zfill(4)+".dat"
    meshdim = np.genfromtxt(mdpath)
    nNodes = meshdim[0].astype(int)    ## For 2D
    nBnodes = meshdim[3].astype(int)   ## For 2D
    
    mcrdpath = "../../../data/meshData/dimcrn"+str(ip+1).zfill(4)+".dat"
    meshdim = np.genfromtxt(mcrdpath)
    nCnodes = meshdim[0].astype(int)
    nRnodes = meshdim[1].astype(int)
    
    ## Create mesh and define function space
    mpath = "../data/foo"+str(ip+1)+".xml"   ## path to mesh files
    if ip == 0:
        print("Input XML-Mesh:", mpath)
    mesh = Mesh(mpath)
    V = VectorFunctionSpace(mesh, 'P', 1)

    # Initialize the connectivity between facets and cells
    tdim = mesh.topology().dim()
    mesh.init(tdim-1, tdim)
    
    # Define strain and stress
    def epsilon(u):
        return sym(nabla_grad(u))
        #return 0.5*(nabla_grad(u) + nabla_grad(u).T)
    
    #def sigma(u):
    #    return 1 * (lambda_*nabla_div(u)*Identity(d) + 2*mu*epsilon(u))

    def sigma(u):
        return E * (lambda_*nabla_div(u)*Identity(d) + 2*mu*epsilon(u))

    # Define variational problem
    u = TrialFunction(V)
    d = ufl.domain.find_geometric_dimension(u)# space dimension #
    v = TestFunction(V)
    f = Constant((0, -rho*g))
    T = Constant((0, 0))

    # Variation Form
    a = Form(inner(sigma(u), epsilon(v))*dx)
    L = Form(dot(f, v)*dx + dot(T, v)*ds)

    ## Save deterministic assembly matrices in ".mat" or ".txt" format
    subpath = "../data/Amats/subdom000".rstrip("0")+str(ip+1)
    
    if not os.path.exists(subpath):
        os.makedirs(subpath)

    ## Invoke Deterministic Assembly procedure for each KLE/PCE mode
    A, b = detAssembly(a, L) ## New assembly procedure

    ## Decompose assembly matrices intor DD-blocks
    ADii, ADig, ADgg, bi, bg = OneLevelDDMat(A, b, nNodes, nBnodes, nComp)
    ADci, ADri, ADcc, ADrr, ADcr = TwoLevelDDMat(A, b, nNodes, nBnodes, nCnodes, nComp)


    ## To save output matrices in ".dat" for Fortran
    k = 0
    ADiipath = subpath+"/ADii000".rstrip("0")+str(k+1)+".dat"
    np.savetxt(ADiipath, ADii)
    ADgipath = subpath+"/ADgi000".rstrip("0")+str(k+1)+".dat"
    np.savetxt(ADgipath, ADig)
    ADggpath = subpath+"/ADgg000".rstrip("0")+str(k+1)+".dat"
    np.savetxt(ADggpath, ADgg)
    fipath = subpath+"/bi000".rstrip("0")+str(k+1)+".dat"
    np.savetxt(fipath, bi)
    fgpath = subpath+"/bg000".rstrip("0")+str(k+1)+".dat"
    np.savetxt(fgpath, bg)

    ADcipath = subpath+"/ADci000".rstrip("0")+str(k+1)+".dat"
    np.savetxt(ADcipath, ADci)
    ADripath = subpath+"/ADri000".rstrip("0")+str(k+1)+".dat"
    np.savetxt(ADripath, ADri)
    ADccpath = subpath+"/ADcc000".rstrip("0")+str(k+1)+".dat"
    np.savetxt(ADccpath, ADcc)
    ADrrpath = subpath+"/ADrr000".rstrip("0")+str(k+1)+".dat"
    np.savetxt(ADrrpath, ADrr)
    ADcrpath = subpath+"/ADcr000".rstrip("0")+str(k+1)+".dat"
    np.savetxt(ADcrpath, ADcr)

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
#print("========================================================")

