#!/bin/bash
#SBATCH --time=0-00:30
#SBATCH --mem-per-cpu=3700M
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1

### first lc=current lc and second lc=desiredlc
### for cluster use lc=0.8 for 10K, lc=0.5 for 30K
#sed -i -e 's/lc=0.5;/lc=0.08;/g' fooServer3D.geo

## select the number of partitions
NP=32

module load arch/avx2
## create mesh file::
/home/sudhipv/packages/gmsh/gmsh-3.0.4-Linux/bin/./gmsh -2 foo.geo -part $NP -o foo.msh

## create mesh data::
gfortran preprocmesh1_AD.F90 -O2 -o ./a.out;./a.out
gfortran preprocmesh2_AD.F90 -O2 -o ./a.out;./a.out
gfortran preprocmesh_meas1.F90 -O2 -o ./a.out;./a.out
gfortran preprocmesh_meas2.F90 -O2 -o ./a.out;./a.out
