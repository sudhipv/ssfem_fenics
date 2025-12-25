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

print("=============================================================")
print("Running FEniCS for 2D stochastic Poisson equation/twolevel...")

## For avoiding FEniCS internal node reordering optimized for speed
parameters['reorder_dofs_serial'] = False

## Spatialy varying KLE/PCE modes
class MyExpression(Expression):
    def __init__(self, params, **kwargs):
        self.index  = params[0]
        self.ndim   = params[1]
        self.sIndex = params[2].astype(int)
        self.mIndex = params[3].astype(int)
    
    def eval(self, values, x):
        ## For bx=by=1 good for upto 15 RVs cases (unit square)
        multipliers = [0.92184, 0.49248, 0.29374, 0.20437, 0.15576]
        omegas = [1.30654, 3.67319, 6.58462, 9.63168, 12.72324]
    
        # Mean and Standard Deviation of Underlaying Gaussian
        meang = 0.1      # similar to meanc in Fortran package
        sigma = 0.3      # similar to sigma in Fortran package
        
### Trunctaed PCE of lognormal process
## Automation to n-RVs : ----------------------------------------------------------------
## KLE: obtain KLE terms
        g = []
        for i in range(self.ndim):
            Xindex = self.sIndex[i,0]
            if (Xindex % 2) == 0:
                ##print("even")
                gg1 = multipliers[Xindex-1] * (sin(omegas[Xindex-1]*(x[0]-0.5)))
            else:
                ##print("odd")
                gg1 = multipliers[Xindex-1] * (cos(omegas[Xindex-1]*(x[0]-0.5)))

            Yindex = self.sIndex[i,1]
            if (Yindex % 2) == 0:
                ##print("even")
                gg2 = multipliers[Yindex-1] * (sin(omegas[Yindex-1]*(x[1]-0.5)))
            else:
                ##print("odd")
                gg2 = multipliers[Yindex-1] * (cos(omegas[Yindex-1]*(x[1]-0.5)))

            g.append(sigma*gg1*gg2)

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

###---------------------------------------------------------------------------------------
### Manual PCE temrs definition:: Keep it for varification
### 2RVs case: ---------------------------------------------------------------------------
#        if self.ndim == 2:
#            g = []
#            g.append(sigma*multipliers[0]*multipliers[0])
#            g.append(sigma*multipliers[0]*multipliers[1])
#
#            if self.index == 0:
#                values[0] = meang*1.0
#            elif self.index == 1:
#                values[0] = meang*(g[0]*cos(omegas[0]*(x[0]-0.5))*cos(omegas[0]*(x[1]-0.5)))
#            elif self.index == 2:
#                values[0] = meang*(g[1]*cos(omegas[0]*(x[0]-0.5))*sin(omegas[1]*(x[1]-0.5)))
#            elif self.index == 3:
#                values[0] = meang*(g[0]*cos(omegas[0]*(x[0]-0.5))*cos(omegas[0]*(x[1]-0.5)))**2/np.sqrt(2)
#            elif self.index == 4:
#                values[0] = meang*(g[0]*cos(omegas[0]*(x[0]-0.5))*cos(omegas[0]*(x[1]-0.5))) * \
#                            (g[1]*cos(omegas[0]*(x[0]-0.5))*sin(omegas[1]*(x[1]-0.5)))
#            else:
#                values[0] = meang*(g[1]*cos(omegas[0]*(x[0]-0.5))*sin(omegas[1]*(x[1]-0.5)))**2/np.sqrt(2)
#
### 3RVs case: ---------------------------------------------------------------------------
#        if self.ndim ==3:
#            g = []
#            g.append(sigma*multipliers[0]*multipliers[0])
#            g.append(sigma*multipliers[0]*multipliers[1])
#            g.append(sigma*multipliers[1]*multipliers[0])
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
#
##---------------------------------------------------------------------------
##---------------------------------------------------------------------------
## Stochastic variational form
def stoVartional(params, u, v, g, f):
    ## Method-3: Similar to puffin/Matlab package (more non-zeros):
    ## Boundary conditions applied externally
    kleID = params[0]       ##Check if we need kleID here??
    cd = MyExpression(params, degree=0)
    a = Form(cd*inner(grad(u), grad(v))*dx)
    L = Form(f*v*dx + g*v*ds)
    
    return a, L

##---------------------------------------------------------------------------
def detAssembly(a, L):    ##Check if we need to pass 'k' here??
    
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
    
    # APPLY BC (on physical boundary at Oth PCE only)
    if k==0:
        bc.apply(A_g,f_g)

    #A = np.zeros([V.dim(), V.dim()]) #b = np.zeros([V.dim()])
    A = A_g.array()
    b = f_g.get_local()
    #b = f_g.array() #DeprecationWarning:
    return A, b

##---------------------------------------------------------------------------
## OneLevel Deterministic assembly matrix decomposer
def OneLevelDDMat(A, b, nNodes, nBnodes):
    ADii = A[0:(nNodes-nBnodes),0:(nNodes-nBnodes)]
    ADgi = A[0:(nNodes-nBnodes),(nNodes-nBnodes):nNodes]
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

##---------------------------------------------------------------------------
#################### MAIN CODE #####################
## Path to the global meshdims.dat       ## Fortran Decomposed
meshdimpath = "../../../data/meshData/meshdim.dat"
nParts = np.genfromtxt(meshdimpath)
nParts = nParts[4].astype(int)     ## For 2D
print("number of partitions:", nParts)

pcedatapath = "../../../data/klePceData/pcedata.dat"
pcedata = np.genfromtxt(pcedatapath)
order_u = pcedata[0].astype(int)     ## order of output PCE (input==2)
dim_L   = pcedata[1].astype(int)     ## number of KLE modes (RVs)
PCE_u   = pcedata[2].astype(int)     ## number of PCE outputs
PCE_A   = pcedata[3].astype(int)    ## number of PCE inputs

print("number of dimensions:", dim_L)
print("number of pceInputss:", PCE_A)
print("number of pceOutputs:", PCE_u)

### Input/output PCE terms
def nfactorial(nf):
    if nf == 0:
        return 1
    else:
        return nf * nfactorial(nf-1)
#dim_L = input("Enter PCE Dimension (RVs):  ")  ## number of KLE modes or number of RVs
#order_A = 2      ## order of input PCE
#order_u = 3      ## order of output PCE
## Function to get factorial to get number of PCE terms
## Input/output number of PCE terms
#PCE_A = nfactorial((order_A+dim_L))/(nfactorial(order_A)*nfactorial(dim_L))
#PCE_u = nfactorial((order_u+dim_L))/(nfactorial(order_u)*nfactorial(dim_L))

## Load ijk and Cijk files
cijkpath = "../data/inputs/cijk"+str(dim_L)+str(order_u)+".mat"   ## path to mesh files
cijkMat = sio.loadmat(cijkpath)
ijk = cijkMat['ijk']
cijk = cijkMat['cijk']
print('input Cijk:', cijkpath)

## Load KLE/PCE Indices for Automation
sIndexPath = "../../../data/klePceData/sortIndex.dat"
sIndex = np.genfromtxt(sIndexPath)
#print(sIndex)
#print(sIndex[0,0])
#print(sIndex[1,1])

mIndexPath = "../../../data/klePceData/multiIndex.dat"
mIndex = np.genfromtxt(mIndexPath)
#print(mIndex)
#print(mIndex[0,0])
#print(mIndex[9,1])

##---------------------------------------------------------------------------
## Fore each subdomain invokes assembly procedure for deterministic DD blocks
for ip in range(nParts):

    ## Get DDM meshdims data  ## NOTE: Zfill=3 >Matlab & 4 >Fortran DDM data
    mdpath = "../../../data/meshData/meshdim"+str(ip+1).zfill(4)+".dat"
    meshdim = np.genfromtxt(mdpath)
    nNodes = meshdim[0].astype(int)    ## For 2D
    nBnodes = meshdim[3].astype(int)   ## For 2D

    mcrdpath = "../../../data/meshData/dimcrn"+str(ip+1).zfill(4)+".dat"
    meshdim = np.genfromtxt(mcrdpath)
    nCnodes = meshdim[0].astype(int)   ## For 2D
    nRnodes = meshdim[1].astype(int)   ## For 2D

    ## Create mesh and define function space
    #mesh = Mesh("inputs/foo.xml")  # info(mesh, True) ## print mesh details
    mpath = "../data/foo"+str(ip+1)+".xml"   ## path to mesh files
    print("Input XML-Mesh:", mpath)
    mesh = Mesh(mpath)

    # Initialize the connectivity between facets and cells
    tdim = mesh.topology().dim()
    mesh.init(tdim-1, tdim)

    V = FunctionSpace(mesh, "Lagrange", 1)
    u = TrialFunction(V)
    v = TestFunction(V)

    ## Assemble without applying boundary condition
    gp = Constant(1e50)         ## Conductivity: penalty factor
    gd = Constant(0.0)          ## Dirichlet boundary condition
    gn = Constant(0.0)          ## Neumann boundary condition
    f = Constant(1.0)           ## Right-hand source term
    g = Constant(0.0)           ##

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
        A, b = detAssembly(a, L)     ## New Assembly

        ##if k==0:b = b0  ## Decompose assembly matrices intor DD-blocks
        ADii, ADgi, ADgg, bi, bg = OneLevelDDMat(A, b, nNodes, nBnodes)
        ADci, ADri, ADcc, ADrr, ADcr = TwoLevelDDMat(A, b, nNodes, nBnodes, nCnodes)

        #dpath = subpath+"/AD"+str(k+1)+".dat",  np.savetxt(dpath, A)
        #print("Out Deterministic Assembly:", dpath)
        ## To save output matrices in ".dat" for Fortran
        ADiipath = subpath+"/ADii000".rstrip("0")+str(k+1)+".dat"
        np.savetxt(ADiipath, ADii)
        ADgipath = subpath+"/ADgi000".rstrip("0")+str(k+1)+".dat"
        np.savetxt(ADgipath, ADgi)
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

        if (k == 1):
            print("output Mat-Vec:",subpath)

print("==========================Success============================")
##--------------------------------------------------------------------
#meshdata = "../../../data/meshData/"
#listing = os.listdir(meshdata)
#for infile in listing:
#    if infile.startswith("meshdim0") :
#        print "current file is: " + infile


## Call stochastic vartional form for each input PCE mode
#a, L = stoVartional(0,dim_L u, v, gp, gd, gn, f)
## Call deteministic assembly for each stochastic vartional form
#b = assemble(L)             ## Need to be replaced by detAssembly(L)
#A = detAssembly(a)

#indexi = 0
#N = V.dim()              ## deterministic matrix size
#Ns = (PCE_u*N)           ## stochastic matrix size
#As = np.zeros([Ns,Ns])   ## number of elements
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
