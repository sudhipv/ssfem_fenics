#!/bin/bash
#SBATCH --time=0-00:10
#SBATCH --mem-per-cpu=7700M
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1

## Source the fenics activation file
source ~/testDDM/sgfem_ddm_hpc/external/dolfin/src/fenics_activate.sh

## Compile and execute the preprocessor code
## This will convert meshDataFile (*.dat) to fenicsDataFiles (*.xml) in ../data/
python pointElements_xml.py

## Import the required libraries to run fenics
module load hdf5-mpi/1.8.18 boost eigen python/3 python35-scipy-stack/2017a petsc/3.7.5 fftw-mpi/3.3.6
source ~/software/dolfin/share/dolfin/dolfin.conf
source ~/fenics/bin/activate

## Compile and execute the FEniCS/dolfin assembly routines
## This will create assemlby matrices-vectors in ../data/Amats/
python poisson3D_stoDDMAssembly_twoLevel

exit
