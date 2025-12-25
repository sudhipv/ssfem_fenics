#!/bin/bash
## Adjust the mesh density parameter "lc"
## lc=old;/lc=new for foo*.geo
sed -i -e 's/lc=0.25;/lc=0.2;/g' foo3D.geo

## select the number of partitions
NP=8

## create mesh file::
gmsh -3 foo3D.geo -part $NP -o foo3D.msh

## create mesh data::
gfortran preprocmesh3D1_AD.F90 -O2 -o ./a.out;./a.out
gfortran preprocmesh3D2_AD.F90 -O2 -o ./a.out;./a.out

## create measurement data:: we don't have this for 3D
#gfortran preprocmesh_meas1.F90 -O2 -o ./a.out;./a.out
#gfortran preprocmesh_meas2.F90 -O2 -o ./a.out;./a.out

## mesh visualization::
## gmsh gmsh.msh

## remove unnecessory files created by "sed"
rm *.geo-e
