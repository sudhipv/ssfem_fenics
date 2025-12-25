#!/bin/bash
#SBATCH --time=0-00:15
#SBATCH --mem-per-cpu=7700M
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1

### Import the required libraries to run fenics: graham
#module load hdf5-mpi/1.8.18 boost eigen python/3 python35-scipy-stack/2017a petsc/3.7.5 fftw-mpi/3.3.6
#source ~/software/dolfin/share/dolfin/dolfin.conf
#source ~/fenics/bin/activate
### Import the required libraries to run fenics: cedar
#module load hdf5-mpi/1.8.18 boost eigen python/3.5 scipy-stack/2017b petsc/3.7.5 fftw-mpi/3.3.6
#source ~/software/dolfin/share/dolfin/dolfin.conf
#source ~/fenics/bin/activate
source fenics_activate.sh

## Compile and execute the preprocessor code
## This will convert meshDataFile (*.dat) to fenicsDataFiles (*.xml) in ../data/
python gmshData_fenicsXML_3D.py

## Compile and execute the FEniCS/dolfin assembly routines
## This will create assemlby matrices-vectors in ../data/Amats/
python poisson3D_stochasticDDM_onelevel.py

exit
