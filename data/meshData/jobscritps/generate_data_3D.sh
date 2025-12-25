#!/bin/bash
##SBATCH --time=0-00:30
##SBATCH --mem-per-cpu=3700M
##SBATCH --nodes=1
##SBATCH --tasks-per-node=1

## select the number of partitions
NP=32

## create mesh file::
/home/sudhipv/project/sudhipv/packages/gmsh/gmsh-3.0.4-Linux/bin/./gmsh -3 foo.geo -part $NP -o foo3D.msh

## create mesh data::
gfortran preprocmesh3D1_AD.F90 -O2 -o ./a.out;./a.out
gfortran preprocmesh3D2_AD.F90 -O2 -o ./a.out;./a.out

## mesh visualization::
# gmsh gmsh.msh
