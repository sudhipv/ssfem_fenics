## Copyright (C) Oct-2017 Ajit Desai
# Code to convert mesh data in points and cells format to dolfin-xml format
# Using meshio and lxml module https://github.com/nschloe/meshio

# First added:  2007-10-08
# meshio can be directly use to convert any mesh: meshio-convert input.msh output.vtu

import numpy as np
import meshio as ms  ## pip3 install meshio
import lxml          ## pip3 install lxml

nParts = np.genfromtxt("../../../data/meshData/num_partition.dat",delimiter=',')
nParts = nParts.astype(int)

for i in range(nParts):
    # Provide path to the input points/cells files in ".dat",".txt" formats
    ppath = "../../../data/meshData/rpoints00"+str(i+1)+".dat"   ## points file path
    print(ppath)
    pts = np.genfromtxt(ppath,delimiter=',')
    nPoints = len(pts)   ## Number of points (used for 2D)
    points = np.c_[pts, np.zeros(nPoints)]  ## adding zeros in 3rd column for 2D mesh
    # print(points)

    rpath = "../../../data/meshData/rtriangles00"+str(i+1)+".dat" ## elements file path
    print(rpath)
    elements = np.genfromtxt(rpath, delimiter=',')
    cell = np.delete(elements, 3, 1)   ## deleted last column of all ones (tags)
    cells = cell.astype(int)           ## converted all elements to integer
    cells = cells - 1                  ## index start with zero in python
    cells = {'triangle': cells}        ## meshio required input cells format
    # print(cells)

    mpath = "../data/foo"+str(i+1)+".xml" ## output file path
    ms.write(mpath,points,cells)   ## write mesh.xml for dolfin

