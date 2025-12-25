#
# Copyright (C) 2017 Ajit Desai, PhD Candidate, ajit.ndesai@gmail.com
#
# parallel code to convert mesh data in points and cells format to dolfin-xml format
# Using meshio and lxml module https://github.com/nschloe/meshio
#
# First added:  2007-10-08
# meshio can be directly use to convert any mesh: meshio-convert input.msh output.vtu

##from petsc4py import PETSc
import mpi4py                          ## pip3 install mpi4py
mpi4py.rc(initialize=False, finalize=False)
from mpi4py import MPI
import numpy as np
import meshio as ms                    ## pip3 install meshio
import lxml                            ## pip3 install lxml

### Initialize MPI for python
MPI.Init()
comm = MPI.COMM_WORLD
i = comm.Get_rank()
## print('My rank is ',i)

if i == 0:
    print("============================================================")
    print("Parallel: Convert GMSH msh-mesh to FEniCS xml-mesh")

## Path to the global points*.dat         ## Fortran Decomposed
ppath = "../../../data/meshData/points"+str(i+1).zfill(4)+".dat"
pts = np.genfromtxt(ppath)
nPoints = len(pts)                     ## Number of points (used for 2D)
points = np.c_[pts, np.zeros(nPoints)] ## adding zeros in 3rd column for 2D mesh
if (i==0):
    print("Input GMSH: points:    ", ppath)

## Path to the global triangles*.dat      ## Fortran Decomposed
rpath = "../../../data/meshData/triangles"+str(i+1).zfill(4)+".dat" ## Global
elements = np.genfromtxt(rpath)
cells = elements.astype(int)           ## converted all elements to integer
cells = cells - 1                      ## index start with zero in python
cells = {'triangle': cells}            ## meshio required input cells format
if (i==0):
    print("Input GMSH: triangles: ", rpath)

## Path to write foo*.xml files
mpath = "../data/foo"+str(i+1)+".xml"  ## output file path
mesh = ms.Mesh(points,cells)
ms.write(mpath,mesh)           ## write mesh.xml for dolfin
if (i==0):
    print("Output FEniCS-XML mesh:" , mpath)
    print("==========================Success============================")

MPI.Finalize()

################========================================################
#data = None
#if i == 0:
#    print("========================================================")
#    print("Parallel: Convert GMSH msh-mesh to FEniCS xml-mesh")
#
#    ## Path to the local meshdim.dat         ## Matlab Decomposed
#    #nParts = np.genfromtxt("../data/meshData/num_partition.dat",delimiter=',')
#    #nParts = nParts.astype(int)
#
#    ## Path to the global meshdims.dat       ## Fortran Decomposed
#    nParts = np.genfromtxt("../../../data/meshData/meshdim.dat")
#    nParts = nParts[4].astype(int)
#    print("number of partitions:", nParts)
#else:
#    nParts = None
#
#data = comm.bcast(nParts, root=0)
#
#for i in range(nParts):
## Provide path to the input points/cells files in ".dat/.txt" formats
## Path to the local points*.dat     ## Matlab Decomposed
#ppath = "../data/meshData/rpoints00"+str(i+1)+".dat"
#pts = np.genfromtxt(ppath,delimiter=',')
#
## Path to the local traingles*.dat  ## Matlab Decomposed
#rpath = "../data/meshData/rtriangles00"+str(i+1)+".dat"
#elements = np.genfromtxt(rpath, delimiter=',')
#cell = np.delete(elements, 3, 1)    ## deleted last column of all ones (tags)
#cells = cell.astype(int)
