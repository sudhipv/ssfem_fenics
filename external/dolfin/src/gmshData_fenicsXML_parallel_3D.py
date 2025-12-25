#
# Copyright (C) 2017 Ajit Desai, Ph.D. Candidate ajit.ndesai@gmail.com
#
####################################################################
## AJIT DESAI: Original code written in Dec-2017
## This code convert points.dat, elements.dat to FEniCS-XML mesh
## run this code using "python3 gmshPointElements_fenicsXML"
## IMP: Need meshio and lxml packages installed to use this code
####################################################################
#
# This is demo for stochastic ddm using FEniCS assembly
# This code only extract subdomain level sto-dd matricies
# To "pip install mpi4py numpy meshio: first "source fenics_activate.sh"

##from petsc4py import PETSc
import mpi4py                          ## pip3 install mpi4py
#mpi4py.rc(initialize=False, finalize=False)
from mpi4py import MPI
import numpy as np
import meshio as ms                    ## pip3 install meshio
import lxml                            ## pip3 install lxml

### Initialize MPI for python
#MPI.Init()
comm = MPI.COMM_WORLD
i = comm.Get_rank()
## print('My rank is ',i)

if i == 0:
    print("======================================================")
    print("Convert Local 3D GMSH msh-mesh data to FEniCS xml-mesh")

## NOTE: zfill=4 for both matlab & fortran
## NOTE: rpoints0*/points0* are same in matlab for only one-level decomposition
## rtetrahedrans=>a=matlab, tetrahedrans=fortran
ppath = "../../../data/meshData/points"+str(i+1).zfill(4)+".dat"
tpath = "../../../data/meshData/tetrahedrons"+str(i+1).zfill(4)+".dat"

if (i==0):
    print("Input GMSH: points: ", ppath)
    print("Input GMSH: cells : ", tpath)

## Read co-ordinates (points) and convert to FEniCS foramt
#pts = np.genfromtxt(ppath)
pts = np.genfromtxt(ppath) ## Works with both Matlab/Fortran
nPoints = len(pts)             ## Number of points (used for 2D)
points = np.c_[pts]            ## adding zeros in 3rd column for 2D mesh

## Read co-ordinates connectivity (elements/cells) and convert to FEniCS foramt
#elements = np.genfromtxt(tpath)
elements = np.genfromtxt(tpath)
cells = elements.astype(int)   ## converted all elements to integer
cells = cells - 1              ## index start with zero in python
cells = {'tetra': cells}       ## meshio required input cells format

## Write output to the FEniCS readable XML format
mpath = "../data/foo"+str(i+1)+".xml"   ## output file path
mesh = ms.Mesh(points,cells)
ms.write(mpath,mesh)   ## write mesh.xml for dolfin
if (i==0):
    print("Output FEniCS-mesh: " , mpath)
    print("==========================Success============================")

#MPI.Finalize()
################========================================################
### Need for validation and postProcessing (to write VTK in Matlab)
#print("========================================================")
#print("Convert Global GMSH data to FEniCS mesh")
#ppath = '../../../data/meshData/points.dat'
#tpath = '../../../data/meshData/tetrahedrons4c.dat' ## Only 4 columns
#
### Read co-ordinates (points) and convert to FEniCS foramt
##pts = np.genfromtxt(ppath,delimiter=',')
#pts = np.genfromtxt(ppath)
#nPoints = len(pts)                      ## Number of points (used for 2D)
#points = np.c_[pts]
#print('P:Inputs  :',ppath)
#
### Read co-ordinates connectivity (elements/cells) and convert to FEniCS foramt
#elements = np.genfromtxt(tpath)##elements = elements[:,0:3]
#cells = elements.astype(int)   ## converted all elements to integer
#cells = cells - 1              ## index start with zero in python
#cells = {'tetra': cells}       ## meshio required input cells format
#print('T:Inputs  :', tpath)
#
### Write output to the FEniCS readable XML format
#mpath = '../data/foo.xml'    ## output file path
#ms.write(mpath,points,cells)   ## write mesh.xml for dolfin
#print('XML:Output:', mpath)


