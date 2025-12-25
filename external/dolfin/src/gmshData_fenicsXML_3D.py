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
#

# Begin demo

#from petsc4py import PETSc
import numpy as np
import meshio as ms
import lxml

print("======================================================")
print("Convert Local 3D GMSH msh-mesh data to FEniCS xml-mesh")
## Path to the global meshdims.dat       ## Fortran Decomposed
nParts = np.genfromtxt("../../../data/meshData/num_partition.dat",delimiter=',')
nParts = nParts.astype(int)
print("number of partitions:", nParts)

for i in range(nParts):
    ## Specify path to the input meshdata files
    ## NOTE: zfill=4 for both matlab & fortran
    ## NOTE: rpoints0*/points0* are same in matlab for only one-level decomposition
    ppath = "../../../data/meshData/points"+str(i+1).zfill(4)+".dat"
    ## rtetrahedrans=>a=matlab, tetrahedrans=fortran
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

print("==========================Success============================")
