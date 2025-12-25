
# Copyright (C) 2017 Ajit Desai, Ph.D. Candidate ajit.ndesai@gmail.com

################################################################
## AJIT DESAI: Modified Oct-2017
## from terminal go to the folder where you have this code
## from "Docker" start "fenicsproject start myFenics"
## after that run this code using "python ddm_**.py"
################################################################

## from petsc4py import PETSc
import numpy as np
import scipy.io as sio
import os as os
from dolfin import *
#import mpi4py                          ## pip3 install mpi4py
#mpi4py.rc(initialize=False, finalize=False)
from mpi4py import MPI

## Initialize MPI for python
#mpi_comm = MPI.comm_world()
#ip = MPI.rank(mpi_comm)
#comm_mpi4py = mpi_comm_world().tompi4py()
#ip = comm_mpi4py.Get_rank()
#MPI.Init()
comm = MPI.COMM_WORLD
ip = comm.Get_rank()
## print('My rank is ',ip)

if ip==0:
    print("===================================================")
    print("Running FEniCS for 3D stochastic Poisson/twoLevel..")

## For avoiding FEniCS internal node reordering optimized for speed
parameters['reorder_dofs_serial'] = False

## Spatialy varying KLE/PCE modes
class MyExpression(Expression):
    def __init__(self, params, **kwargs):
        self.index = params[0]
        self.ndim = params[1]
        self.sIndex = params[2].astype(int)
        self.mIndex = params[3].astype(int)
    
    def eval(self, values, x):
        ## For bz=1, a=1 good for upto 20 RVs cases (unit square)
        multipliersZ = [0.92184, 0.49248, 0.29374, 0.20437, 0.15576, 0.1256, 0.1051, 0.0903, 0.0791, 0.0704]
        omegasZ = [1.30654, 3.67319, 6.58462, 9.63168, 12.72324, 15.8341, 18.9550, 22.0817, 25.2120, 28.345]
        ## For bx=by=1, a=0.2 good for upto 20 RVs cases
        multipliers = [0.9835, 0.2685, 0.1402, 0.0942, 0.0709, 0.0569, 0.0474]
        omegas = [3.1105, 16.3199, 31.7310, 47.3351, 62.9906, 78.5398, 94.3538]

        ### For bx=by=1 good for upto 15 RVs cases (unit square)
        #multipliers = [0.92184, 0.49248, 0.29374, 0.20437, 0.15576]
        #omegas = [1.30654, 3.67319, 6.58462, 9.63168, 12.72324]
        ### For bz=1 good for upto 15 RVs cases (L=5)
        #multipliersZ = [0.60662, 0.18584, 0.10193, 0.090852, 0.077738]
        #omegasZ = [0.62832, 2.0653, 3.87103, 4.39823, 5.10394]
        
        # Mean and Standard Deviation of Underlaying Gaussian
        meang = 1.046     # exp(0+1/2*(sigma^2)), meanc in Fortran package
        sigma = 0.3      # mul=1.03 for mug=0,sigma=0.3 | mul=1.27 for mug=0,sigma=0.7
        
### Trunctaed PCE of lognormal process
## Automation to n-RVs : ----------------------------------------------------------------
        ## KLE: obtain KLE terms
        g = []
        for i in range(self.ndim):
            Xindex = self.sIndex[i,0]
            if (Xindex % 2) == 0:
            ##print("even")
                gg1 = multipliers[Xindex-1] * (sin(omegas[Xindex-1]*(x[0]-0.1)))
            else:
                ##print("odd")
                gg1 = multipliers[Xindex-1] * (cos(omegas[Xindex-1]*(x[0]-0.1)))

            Yindex = self.sIndex[i,1]
            if (Yindex % 2) == 0:
                ##print("even")
                gg2 = multipliers[Yindex-1] * (sin(omegas[Yindex-1]*(x[1]-0.1)))
            else:
                ##print("odd")
                gg2 = multipliers[Yindex-1] * (cos(omegas[Yindex-1]*(x[1]-0.1)))

            Zindex = self.sIndex[i,2]
            if (Zindex % 2) == 0:
                ##print("even")
                gg3 = multipliersZ[Zindex-1] * (sin(omegasZ[Zindex-1]*(x[2]-0.5)))
            else:
                ##print("odd")
                gg3 = multipliersZ[Zindex-1] * (cos(omegasZ[Zindex-1]*(x[2]-0.5)))

            g.append(sigma*gg1*gg2*gg3)

        ## PCE: obtain PCE terms
        newY = 1.0
        i = self.index
        for j in range(self.ndim):
            idx = mIndex[i,j]
            if idx == 0:
                yy = 1.0
            elif idx == 1:
                yy = g[j]
            else:
                varFacto = nfactorial(idx)
                yy = (g[j]**(idx))/np.sqrt(varFacto)

            ## Keep all indentations properly (there was error here due to indentation)
            newY =  newY*yy

        values[0] = meang * newY

##---------------------------------------------------------------------------
## Stochastic variational formic variational form
def stoVartional(params, u, v, g, f):
    ## Method-3: Similar to puffin/Matlab package (more non-zeros):
    ## Boundary conditions embedded in the variational form
    
    cd = MyExpression(params, degree=0)
    a = Form(cd*inner(grad(u), grad(v))*dx)
    L = Form(f*v*dx + g*v*ds)
    
    return a, L

##---------------------------------------------------------------------------
## Global deterministic assembly for each KLE mode
def detAssembly(a, L, k):
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
        # Define Dirichlet boundary (x = 0 or x = 1)
        # Mark boundary subdomians
        left =  CompiledSubDomain("near(x[2], side) && on_boundary", side = 0.0)
        right = CompiledSubDomain("near(x[2], side) && on_boundary", side = 1.0)
        #left2 =  CompiledSubDomain("near(x[2], side) && on_boundary", side = 0.0)
        #right2 = CompiledSubDomain("near(x[2], side) && on_boundary", side = 2.0)
        
        # Define boundary condition
        u0 = Constant(0.0)
        bcl = DirichletBC(V, u0, left)
        bcr = DirichletBC(V, u0, right)
        bcl.apply(A_g,f_g)
        bcr.apply(A_g,f_g)
        #print("Appying BCs")

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
    ADgi = A[0:(nNodes-nBnodes),(nNodes-nBnodes):nNodes]
    #A[(nNodes-nBnodes):nNodes,0:(nNodes-nBnodes)]   ## A[0:(nNodes-nBnodes),(nNodes-nBnodes):nNodes]
    ADgg = A[(nNodes-nBnodes):nNodes,(nNodes-nBnodes):nNodes]
    bi   = b[0:(nNodes-nBnodes)]
    bg   = b[(nNodes-nBnodes):nNodes]
    return ADii, ADgi, ADgg, bi, bg

##---------------------------------------------------------------------------
## OneLevel Deterministic assembly matrix decomposer
def TwoLevelDDMat(A, b, nNodes, nBnodes, nCnodes):
    #    ADii = A[0:(nNodes-nBnodes),0:(nNodes-nBnodes)]
    #    ADgi = A[0:(nNodes-nBnodes),(nNodes-nBnodes):nNodes]
    #    ADgg = A[(nNodes-nBnodes):nNodes,(nNodes-nBnodes):nNodes]
    
    ADci = A[0:(nNodes-nBnodes),(nNodes-nCnodes):nNodes]
    ADri = A[0:(nNodes-nBnodes),(nNodes-nBnodes):(nNodes-nCnodes)]
    ADcc = A[(nNodes-nCnodes):nNodes,(nNodes-nCnodes):nNodes]
    ADrr = A[(nNodes-nBnodes):(nNodes-nCnodes),(nNodes-nBnodes):(nNodes-nCnodes)]
    ADcr = A[(nNodes-nBnodes):(nNodes-nCnodes),(nNodes-nCnodes):nNodes]
    
    #    bi   = b[0:(nNodes-nBnodes)]
    #    bg   = b[(nNodes-nBnodes):nNodes]
    #    bc   = b[(nNodes-nCnodes):nNodes]
    #    br   = b[(nNodes-nBnodes):(nNodes-nCnodes)]
    return ADci, ADri, ADcc, ADrr, ADcr   # ADii, ADgi, ADgg, bi, bg  #, br, bc

##-------------------------------------------------------------------------
############################### MAIN CODE #################################

## Path to the global meshdims.dat
meshdimpath = '../../../data/meshData/meshdim.dat'
nParts = np.genfromtxt(meshdimpath)
nParts = nParts[5].astype(int)   ## For 3D

pcedatapath = "../../../data/klePceData/pcedata.dat"
pcedata = np.genfromtxt(pcedatapath)
order_u = pcedata[0].astype(int)     ## order of output PCE (input==2)
dim_L   = pcedata[1].astype(int)     ## number of KLE modes (RVs)
PCE_u   = pcedata[2].astype(int)     ## number of PCE outputs
PCE_A   = pcedata[3].astype(int)     ## number of PCE inputs

### Load ijk and Cijk files
#cijkpath = "../data/inputs/cijk"+str(dim_L)+str(order_u)+".mat"   ## path to mesh files
#cijkMat = sio.loadmat(cijkpath)
#ijk = cijkMat['ijk']
#cijk = cijkMat['cijk']
#print("input Cijk:", cijkpath)

## Load KLE/PCE Indices for Automation
#sIndexPath = "../../../data/klePceData/sortIndex3D.dat"
sIndexPath = "../../../data/klePceData/sortIndex3D_a02.dat"
sIndex = np.genfromtxt(sIndexPath)
mIndexPath = "../../../data/klePceData/multiIndex.dat"
mIndex = np.genfromtxt(mIndexPath)

if ip == 0:
    print("number of dimensions:", dim_L)
    print("number of pceInputss:", PCE_A)
    print("number of pceOutputs:", PCE_u)
    print("number of partitions:", nParts)

# Function to get factorial to get number of PCE terms
def nfactorial(nf):
    if nf == 0:
        return 1
    else:
        return nf * nfactorial(nf-1)

##---------------------------------------------------------------------------
## Fore each subdomain invokes assembly procedure for deterministic DD blocks
## for ip in range(nParts):

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
if ip == 0:
    print("Input XML-Mesh:", mpath)
mesh = Mesh(mpi_comm_self(), mpath)

# Initialize the connectivity between facets and cells
tdim = mesh.topology().dim()
mesh.init(tdim-1, tdim)

# Functional/Test/Trial space & constants
V = FunctionSpace(mesh, "Lagrange", 1)
u, v = TrialFunction(V), TestFunction(V)
f = Constant(1.0)
g = Constant(0.0)

for k in range(int(float(PCE_A))):
    ## Call stochastic vartional form for each input PCE mode
    
    ## Save deterministic assembly matrices in ".mat" or ".txt" format
    subpath = "../data/Amats/subdom000".rstrip("0")+str(ip+1)
    
    if not os.path.exists(subpath):
        os.makedirs(subpath)
    
    params = ([k, dim_L, sIndex, mIndex])

    ##--------------------------------------------------------------------
    ## FEniCS (genearalize) assembly procedure ## NEW
    ## Invoke variation fomrulaiton for each KLE/PCE mode
    a, L = stoVartional(params, u, v, g, f)

    ## Invoke Deterministic Assembly procedure for each KLE/PCE mode
    A, b = detAssembly(a, L, k) ## New assembly procedure

    ##if k==0:b = b0  ## Decompose assembly matrices intor DD-blocks
    ADii, ADig, ADgg, bi, bg = OneLevelDDMat(A, b, nNodes, nBnodes)
    ADci, ADri, ADcc, ADrr, ADcr = TwoLevelDDMat(A, b, nNodes, nBnodes, nCnodes)

    #dpath = subpath+"/AD"+str(k+1)+".dat",  np.savetxt(dpath, A)
    #print("Out Deterministic Assembly:", dpath)
    ## To save output matrices in ".dat" for Fortran
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

    if ((k == 1) and (ip ==0)):
        print("output Mat-Vec:",subpath)

if ip == 0:
    print("==========================Success============================")

#################################### END ######################################
####---------------------------------------------------------------------------------------
#### Manual PCE temrs definition:: Keep it for varification
#### 2RVs case:---------------------------------------------------------------------------
#        g = []
#        g.append(sigma*multipliers[0]*multipliers[0]*multipliersZ[0])
#        g.append(sigma*multipliers[0]*multipliers[1]*multipliersZ[0])
#        g.append(sigma*multipliers[1]*multipliers[0]*multipliersZ[0])
#        if self.ndim == 2:
#            if self.index == 0:
#                values[0] = meang*1.0
#            elif self.index == 1:
#                values[0] = meang*(g[0]*cos(omegas[0]*(x[0]-0.5))*cos(omegas[0]*(x[1]-0.5))*cos(omegasZ[0]*(x[2]-2.5)))
#            elif self.index == 2:
#                values[0] = meang*(g[1]*cos(omegas[0]*(x[0]-0.5))*sin(omegas[1]*(x[1]-0.5))*cos(omegasZ[0]*(x[2]-2.5)))
#            elif self.index == 3:
#                values[0] = meang*(g[0]*cos(omegas[0]*(x[0]-0.5))*cos(omegas[0]*(x[1]-0.5))*cos(omegasZ[0]*(x[2]-2.5)))**2/np.sqrt(2)
#            elif self.index == 4:
#                values[0] = meang*((g[0]*cos(omegas[0]*(x[0]-0.5))*cos(omegas[0]*(x[1]-0.5))*cos(omegasZ[0]*(x[2]-2.5))) * \
#                    (g[1]*cos(omegas[0]*(x[0]-0.5))*sin(omegas[1]*(x[1]-0.5))*cos(omegasZ[0]*(x[2]-2.5))))
#            else:
#                values[0] = meang*(g[1]*cos(omegas[0]*(x[0]-0.5))*sin(omegas[1]*(x[1]-0.5))*cos(omegasZ[0]*(x[2]-2.5)))**2/np.sqrt(2)
#
### 3RVs case:---------------------------------------------------------------------------
#        if self.ndim ==3:
#            if self.index == 0:
#                values[0] = meang*1.0
#            elif self.index == 1:
#                values[0] = meang*(g[0]*cos(omegas[0]*(x[0]-0.5))*cos(omegas[0]*(x[1]-0.5)))
#            elif self.index == 2:
#                values[0] = meang*(g[1]*cos(omegas[0]*(x[0]-0.5))*sin(omegas[1]*(x[1]-0.5)))
#            elif self.index == 3:
#                values[0] = meang*(g[2]*sin(omegas[1]*(x[0]-0.5))*cos(omegas[0]*(x[1]-0.5)))
#            elif self.index == 4:
#                values[0] = meang*(g[0]*cos(omegas[0]*(x[0]-0.5))*cos(omegas[0]*(x[1]-0.5)))**2/np.sqrt(2)
#            elif self.index == 5:
#                values[0] = meang*(g[0]*cos(omegas[0]*(x[0]-0.5))*cos(omegas[0]*(x[1]-0.5))) * \
#                (g[1]*cos(omegas[0]*(x[0]-0.5))*sin(omegas[1]*(x[1]-0.5)))
#            elif self.index == 6:
#                values[0] = meang*(g[0]*cos(omegas[0]*(x[0]-0.5))*cos(omegas[0]*(x[1]-0.5))) * \
#                (g[2]*sin(omegas[1]*(x[0]-0.5))*cos(omegas[0]*(x[1]-0.5)))
#            elif self.index == 7:
#                values[0] = meang*(g[1]*cos(omegas[0]*(x[0]-0.5))*sin(omegas[1]*(x[1]-0.5)))**2/np.sqrt(2)
#            elif self.index == 8:
#                values[0] = meang*(g[1]*cos(omegas[0]*(x[0]-0.5))*sin(omegas[1]*(x[1]-0.5))) * \
#                (g[2]*sin(omegas[1]*(x[0]-0.5))*cos(omegas[0]*(x[1]-0.5)))
#            else:
#                values[0] = meang*(g[2]*sin(omegas[1]*(x[0]-0.5))*cos(omegas[0]*(x[1]-0.5)))**2/np.sqrt(2)
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
