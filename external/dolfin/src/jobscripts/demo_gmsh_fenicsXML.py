#
# Copyright (C) 2017 Ajit Desai, PhD Candidate, ajit.ndesai@gmail.com
#
# Code to convert mesh data in points and cells format to dolfin-xml format
# Using meshio and lxml module https://github.com/nschloe/meshio

# First added:  2007-10-08
# meshio can be directly use to convert any mesh: meshio-convert input.msh output.vtu

import numpy as np
import meshio as ms                   ## pip3 install meshio
import lxml                           ## pip3 install lxml

print("========================================================")
print("****** Demo: Convert GMSH msh to FEniCS xml *******")

## Path to the global points*.dat    ## Fortran Decomposed
ppath = "points.dat"
pts = np.genfromtxt(ppath)
nPoints = len(pts)
points = np.c_[pts, np.zeros(nPoints)]

## Path to the global triangles*.dat ## Fortran Decomposed
rpath = "triangles.dat"
elements = np.genfromtxt(rpath)
cells = elements.astype(int)

cells = cells - 1
cells = {'triangle': cells}

mpath = "foo.xml" ## output file path
ms.write(mpath,points,cells)   ## write mesh.xml for dolfin
print("========================================================")
print("======================== Success =======================")
print("========================================================")
