#
# Copyright (C) 2017 Ajit Desai, Ph.D. Candidate ajit.ndesai@gmail.com
#
####################################################################
## AJIT DESAI: Original code written in Nov-2017
## from terminal go to the folder where you have this code
## from "Docker" start "fenicsproject start fenics" at "dolfin" dir
## run this code using "python demo_stoDDM_poissonModified.py"
####################################################################
#
# This is demo for stochastic ddm using FEniCS assembly
# This code only extract subdomain level sto-dd matricies
#
# Begin demo

# from petsc4py import PETSc
import numpy as np
import scipy.io as sio
import os as os
from dolfin import *

print("========================================================")
print("Running FEniCS...")

## For avoiding FEniCS internal node reordering optimized for speed
parameters['reorder_dofs_serial'] = False

## Spatialy varying KLE/PCE modes
class MyExpression(Expression):
    def __init__(self, params, **kwargs):
        self.index = params[0]
        self.ndim = params[1]
    
    def eval(self, values, x):
        multipliers = [0.92184, 0.49248, 0.29374]
        omegas = [1.30654, 3.67319, 6.58462]
    
        # Mean and Standard Deviation of Underlaying Gaussian
        meang = 0.1      # similar to meanc in Fortran package
        sigma = 0.3      # similar to sigma in Fortran package
        
        if self.ndim == 2:
            g = []
            g.append(sigma*multipliers[0]*multipliers[0])
            g.append(sigma*multipliers[0]*multipliers[1])
            
            if self.index == 0:
                values[0] = meang*1.0
            elif self.index == 1:
                values[0] = meang*(g[0]*cos(omegas[0]*(x[0]-0.5))*cos(omegas[0]*(x[1]-0.5)))
            elif self.index == 2:
                values[0] = meang*(g[1]*cos(omegas[0]*(x[0]-0.5))*sin(omegas[1]*(x[1]-0.5)))
            elif self.index == 3:
                values[0] = meang*(g[0]*cos(omegas[0]*(x[0]-0.5))*cos(omegas[0]*(x[1]-0.5)))**2/np.sqrt(2)
            elif self.index == 4:
                values[0] = meang*(g[0]*cos(omegas[0]*(x[0]-0.5))*cos(omegas[0]*(x[1]-0.5))) * \
                            (g[1]*cos(omegas[0]*(x[0]-0.5))*sin(omegas[1]*(x[1]-0.5)))
            else:
                values[0] = meang*(g[1]*cos(omegas[0]*(x[0]-0.5))*sin(omegas[1]*(x[1]-0.5)))**2/np.sqrt(2)

        if self.ndim ==3:
            g = []
            g.append(sigma*multipliers[0]*multipliers[0])
            g.append(sigma*multipliers[0]*multipliers[1])
            g.append(sigma*multipliers[1]*multipliers[0])
            if self.index == 0:
                values[0] = meang*1.0
            elif self.index == 1:
                values[0] = meang*(g[0]*cos(omegas[0]*(x[0]-0.5))*cos(omegas[0]*(x[1]-0.5)))
            elif self.index == 2:
                values[0] = meang*(g[1]*cos(omegas[0]*(x[0]-0.5))*sin(omegas[1]*(x[1]-0.5)))
            elif self.index == 3:
                values[0] = meang*(g[2]*sin(omegas[1]*(x[0]-0.5))*cos(omegas[0]*(x[1]-0.5)))
            elif self.index == 4:
                values[0] = meang*(g[0]*cos(omegas[0]*(x[0]-0.5))*cos(omegas[0]*(x[1]-0.5)))**2/np.sqrt(2)
            elif self.index == 5:
                values[0] = meang*(g[0]*cos(omegas[0]*(x[0]-0.5))*cos(omegas[0]*(x[1]-0.5))) * \
                (g[1]*cos(omegas[0]*(x[0]-0.5))*sin(omegas[1]*(x[1]-0.5)))
            elif self.index == 6:
                values[0] = meang*(g[0]*cos(omegas[0]*(x[0]-0.5))*cos(omegas[0]*(x[1]-0.5))) * \
                (g[2]*sin(omegas[1]*(x[0]-0.5))*cos(omegas[0]*(x[1]-0.5)))
            elif self.index == 7:
                values[0] = meang*(g[1]*cos(omegas[0]*(x[0]-0.5))*sin(omegas[1]*(x[1]-0.5)))**2/np.sqrt(2)
            elif self.index == 8:
                values[0] = meang*(g[1]*cos(omegas[0]*(x[0]-0.5))*sin(omegas[1]*(x[1]-0.5))) * \
                (g[2]*sin(omegas[1]*(x[0]-0.5))*cos(omegas[0]*(x[1]-0.5)))
            else:
                values[0] = meang*(g[2]*sin(omegas[1]*(x[0]-0.5))*cos(omegas[0]*(x[1]-0.5)))**2/np.sqrt(2)


## Stochastic variational form
def stoVartional(params, u, v, gp, gd, gn, f):
    ## Method-3: Similar to puffin/Matlab package (more non-zeros):
    ## Boundary conditions embedded in the variational form
    kleID = params[0]
    if kleID == 0:
        cd = MyExpression(params, degree=0)
        a = cd*inner(grad(u), grad(v))*dx + gp*u*v*ds
        L = f*v*dx + ((gp*gd) - gn)*v*ds
    else:
        cd = MyExpression(params, degree=0)
        a = cd*inner(grad(u), grad(v))*dx
        L = f*v*dx + ((gp*gd) - gn)*v*ds
    return a, L

## Stochastic variational form
def newStoVartional(params, u, v, g, f):
    ## Method-3: Similar to puffin/Matlab package (more non-zeros):
    ## Boundary conditions applied externally
    kleID = params[0]
    cd = MyExpression(params, degree=0)
    a = Form(cd*inner(grad(u), grad(v))*dx)
    L = Form(f*v*dx + g*v*ds)
    
    return a, L

## Global deterministic assembly for each KLE mode
def detAssembly(a, L):
    
    b = np.zeros([V.dim()])            ## number of elements
    b_g = assemble(L)                  ## OLD Assembly
    b = b_g.array()
    
    A = np.zeros([V.dim(), V.dim()])   ## number of elements
    for cell in cells(V.mesh()):
        A_local = assemble_local(a, cell)
        local_to_global_map = V.dofmap().cell_dofs(cell.index())
        # Copy local matrix into global matrix
        for (idx_loc_1, idx_glob_1) in enumerate(local_to_global_map):
            for (idx_loc_2, idx_glob_2) in enumerate(local_to_global_map):
                A[idx_glob_1, idx_glob_2] += A_local[idx_loc_1, idx_loc_2]
    return A, b

def newDetAssembly(a, L):
    # Define Dirichlet boundary (x = 0 or x = 1)
    def boundary(x):
        return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS or x[1] < DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS
    # Define boundary condition
    u0 = Constant(0.0)
    bc = DirichletBC(V, u0, boundary)
    
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
    if k==0:
        bc.apply(A_g,f_g)

    #A = np.zeros([V.dim(), V.dim()])
    A = A_g.array()
    b = f_g.array()
    return A, b


## OneLevel Deterministic assembly matrix decomposer
def OneLevelDDMat(A, b, nNodes, nBnodes):
    ADii = A[0:(nNodes-nBnodes),0:(nNodes-nBnodes)]
    ADgi = A[0:(nNodes-nBnodes),(nNodes-nBnodes):nNodes]
    ADgg = A[(nNodes-nBnodes):nNodes,(nNodes-nBnodes):nNodes]
    bi   = b[0:(nNodes-nBnodes)]
    bg   = b[(nNodes-nBnodes):nNodes]
    return ADii, ADgi, ADgg, bi, bg


#################### MAIN CODE #####################
## Path to the global meshdims.dat       ## Fortran Decomposed
meshdimpath = '../../../data/meshData/meshdim.dat'
nParts = np.genfromtxt(meshdimpath)
nParts = nParts[4].astype(int)
print("number of partitions:", nParts)

## Input/output PCE terms
dim_L   = 3      ## number of KLE modes or number of RVs
order_A = 2      ## order of input PCE
order_u = 3      ## order of output PCE
# Function to get factorial to get number of PCE terms
def nfactorial(nf):
    if nf == 0:
        return 1
    else:
        return nf * nfactorial(nf-1)
# Input/output number of PCE terms
PCE_A = nfactorial((order_A+dim_L))/(nfactorial(order_A)*nfactorial(dim_L))
PCE_u = nfactorial((order_u+dim_L))/(nfactorial(order_u)*nfactorial(dim_L))

## Load ijk and Cijk files
cijkpath = "../data/inputs/cijk"+str(dim_L)+str(order_u)+".mat"   ## path to mesh files
cijkMat = sio.loadmat(cijkpath)
ijk = cijkMat['ijk']
cijk = cijkMat['cijk']
print('input Cijk:', cijkpath)


for ip in range(nParts):
    
    ## Get DDM meshdims data
    meshdim = np.genfromtxt("../../../data/meshData/meshdim000"+str(ip+1)+".dat")
    nNodes = meshdim[0].astype(int)
    nBnodes = meshdim[3].astype(int)

    ## Create mesh and define function space
    #mesh = Mesh("inputs/foo.xml")  # info(mesh, True) ## print mesh details
    mpath = "../data/foo"+str(ip+1)+".xml"   ## path to mesh files
    print("Input Mesh:", mpath)
    mesh = Mesh(mpath)

    V = FunctionSpace(mesh, "Lagrange", 1)
    u = TrialFunction(V)
    v = TestFunction(V)

    ## Assemble without applying boundary condition
    gp = Constant(1e50)         ## Conductivity: penalty factor
    gd = Constant(0.0)          ## Dirichlet boundary condition
    gn = Constant(0.0)          ## Neumann boundary condition
    f = Constant(1.0)           ## Right-hand source term
    g = Constant(0.0)           ##

    ## Call stochastic vartional form for each input PCE mode
    #a, L = stoVartional(0,dim_L u, v, gp, gd, gn, f)
    ## Call deteministic assembly for each stochastic vartional form
    #b = assemble(L)             ## Need to be replaced by detAssembly(L)
    #A = detAssembly(a)

    #indexi = 0
    #N = V.dim()              ## deterministic matrix size
    #Ns = (PCE_u*N)           ## stochastic matrix size
    #As = np.zeros([Ns,Ns])   ## number of elements

    for k in range(int(float(PCE_A))):
        ## Call stochastic vartional form for each input PCE mode
        params = ([k, dim_L])
        a_old, L_old = stoVartional(params, u, v, gp, gd, gn, f)      ## OLD
        a, L = newStoVartional(params, u, v, g, f)                    ## NEW
        
        ## Call deteministic assembly for each stochastic vartional form
        A_old, b_old = detAssembly(a_old, L_old)
        
        A, b = newDetAssembly(a, L)     ## New Assembly
        #if k==0:b = b0
        
        ADii_old, ADig_old, ADgg_old, bi_old, bg_old = OneLevelDDMat(A_old, b_old, nNodes, nBnodes)
        ADii, ADig, ADgg, bi, bg = OneLevelDDMat(A, b, nNodes, nBnodes)
        
        ## Save deterministic assembly matrices in ".mat" or ".txt" format
        if (ip < 9):
            subpath = "../data/Amats/subdom000"+str(ip+1)
        else:
            subpath = "../data/Amats/subdom00"+str(ip+1)

        if not os.path.exists(subpath):
            os.makedirs(subpath)
        
        dpath = subpath+"/AD"+str(k+1)+".dat"
        np.savetxt(dpath, A)
        #print("Out Deterministic Assembly:", dpath)
        ## To save output matrices in ".dat" for Fortran
        if (k < 9):
            ADiipath = subpath+"/ADii000"+str(k+1)+".dat"
            np.savetxt(ADiipath, ADii)
            ADgipath = subpath+"/ADgi000"+str(k+1)+".dat"
            np.savetxt(ADgipath, ADig)
            ADggpath = subpath+"/ADgg000"+str(k+1)+".dat"
            np.savetxt(ADggpath, ADgg)
            fipath = subpath+"/bi000"+str(k+1)+".dat"
            np.savetxt(fipath, bi)
            fgpath = subpath+"/bg000"+str(k+1)+".dat"
            np.savetxt(fgpath, bg)
        else:
            ADiipath = subpath+"/ADii00"+str(k+1)+".dat"
            np.savetxt(ADiipath, ADii)
            ADgipath = subpath+"/ADgi00"+str(k+1)+".dat"
            np.savetxt(ADgipath, ADig)
            ADggpath = subpath+"/ADgg00"+str(k+1)+".dat"
            np.savetxt(ADggpath, ADgg)
            fipath = subpath+"/bi00"+str(k+1)+".dat"
            np.savetxt(fipath, bi)
            fgpath = subpath+"/bg00"+str(k+1)+".dat"
            np.savetxt(fgpath, bg)

        ## OLD Assemlby
        dpath = subpath+"/AD"+str(k+1)+"_old.dat"
        np.savetxt(dpath, A)
        ## To save output matrices in ".dat" for Fortran
        if (k < 9):
            ADiipath = subpath+"/ADii000"+str(k+1)+"_old.dat"
            np.savetxt(ADiipath, ADii_old)
            ADgipath = subpath+"/ADgi000"+str(k+1)+"_old.dat"
            np.savetxt(ADgipath, ADig_old)
            ADggpath = subpath+"/ADgg000"+str(k+1)+"_old.dat"
            np.savetxt(ADggpath, ADgg_old)
            fipath = subpath+"/bi000"+str(k+1)+"_old.dat"
            np.savetxt(fipath, bi_old)
            fgpath = subpath+"/bg000"+str(k+1)+"_old.dat"
            np.savetxt(fgpath, bg_old)
        else:
            ADiipath = subpath+"/ADii00"+str(k+1)+"_old.dat"
            np.savetxt(ADiipath, ADii_old)
            ADgipath = subpath+"/ADgi00"+str(k+1)+"_old.dat"
            np.savetxt(ADgipath, ADig_old)
            ADggpath = subpath+"/ADgg00"+str(k+1)+"_old.dat"
            np.savetxt(ADggpath, ADgg_old)
            fipath = subpath+"/bi00"+str(k+1)+"_old.dat"
            np.savetxt(fipath, bi_old)
            fgpath = subpath+"/bg00"+str(k+1)+"_old.dat"
            np.savetxt(fgpath, bg_old)
            #np.savetxt(ADiipath, ADii) #np.savetxt("b.txt", b.array())

            print("output paths:",subpath)

### This section is to build-save the subdomain level stochastic matrices AS^s
#        ## Build stochastic assembly matrix from deterministic calls
#        for i in range(PCE_u):
#            for j in range(PCE_u):
#                if k==ijk[indexi,0] and i==ijk[indexi,1] and j==ijk[indexi,2]:
#                    #print(indexi)
#                    #print(cijk[indexi])
#                    #fromi = i*N toi = (i+1)*N
#                    #fromj = j*N toj = (j+1)*N
#                    As[i*N:(i+1)*N, j*N:(j+1)*N] += cijk[indexi]*A
#                    indexi = indexi+1
#
#    ## Save assembly matrices in ".mat" or ".txt" format
#    spath = "outputs/As"+str(ip+1)+".mat"
#    print("Out Stochastic Assembly:", spath)
#    sio.savemat(spath, {'As':As, 'b':b.array()})
#################################### END #########################################
