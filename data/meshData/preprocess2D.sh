#!/bin/bash
##cd ../../jobscripts/
##sh preprocessKLE.sh
##cd -

cp geofiles/foo.geo foo.geo

## Adjust the mesh density parameter "lc"
## lc=old;/lc=new for foo*.geo
sed -i -e 's/lc=0.09;/lc=0.05;/g' foo.geo

## select the number of partitions
NP=4

## create mesh file::
gmsh -2 foo.geo -part $NP -o foo.msh

## create mesh data::
gfortran preprocmesh1_AD.F90 -O2 -o ./a.out;./a.out
gfortran preprocmesh2_AD.F90 -O2 -o ./a.out;./a.out

## create measurement data::
gfortran preprocmesh_meas1.F90 -O2 -o ./a.out;./a.out
gfortran preprocmesh_meas2.F90 -O2 -o ./a.out;./a.out

## mesh visualization::
#gmsh gmsh.msh 

## remove unnecessory files created by "sed"
rm *.geo-e
