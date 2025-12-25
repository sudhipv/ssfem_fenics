#!/usr/bin/env bash

module load gcc/4.9.1

## select the number of partitions
NP=176

## create mesh file::
/home/ajitd/softwares/gmsh-2.6.1-source/build/./gmsh -2 square.geo -part $NP -o gmsh.msh

## create mesh data::
gfortran preprocmesh1_AD.F90 -O2 -o ./a.out;./a.out
gfortran preprocmesh2_AD.F90 -O2 -o ./a.out;./a.out
gfortran preprocmesh_meas1.F90 -O2 -o ./a.out;./a.out
gfortran preprocmesh_meas2.F90 -O2 -o ./a.out;./a.out

## mesh visualization::
# gmsh gmsh.msh


