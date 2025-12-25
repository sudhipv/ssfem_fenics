#!/bin/bash
#SBATCH --account=def-asarkar
#SBATCH --nodes=1
##SBATCH --ntasks=40
#SBATCH --time=0-00:15
#SBATCH --mem-per-cpu=7700M
#SBATCH --tasks-per-node=32
#SBATCH --job-name=nRV_nNodes_nParts_3DE
#SBATCH --output=%j-%x.txt

## Import the required libraries to run fenics
source fenics_activate.sh

## Compile and execute the preprocessor code
## This will convert meshDataFile (*.dat) to fenicsDataFiles (*.xml)
mpiexec -n 32 python gmshData_fenicsXML_parallel_3D.py

## Compile and execute the FEniCS/dolfin assembly routines
## This will create assemlby matrices-vectors in ../data/Amats/
mpiexec -n 32 python elasticity3D_stochasticDDM_parallel_twolevel.py

deactivate

cd ../../../src/3Delasticity_twolevel/

## Import the required libraries to run fenics
source petsc_activate.sh

make all
mpiexec -n 32 ./a.out -log_view

exit
