#!/bin/bash
#SBATCH --account=def-asarkar
#SBATCH --nodes=1
##SBATCH --ntasks=32
#SBATCH --time=0-00:15
#SBATCH --mem-per-cpu=7700M
#SBATCH --tasks-per-node=32
#SBATCH --output=%j-3DPfenics.out
#SBATCH --mail-user=sudhipv@cmail.carleton.ca
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


## Import the required libraries to run fenics
source fenics_activate.sh
#module load hdf5-mpi/1.8.18 boost eigen python/3.5 scipy-stack/2017b petsc/3.7.5 fftw-mpi/3.3.6
#source ~/software/dolfin_17_1/share/dolfin/dolfin.conf
#source ~/fenics_17_1/bin/activate

#PYTHONPATH=/home/sudhipv/fenics_17_r/lib/python3.5/site-packages:$PYTHONPATH


#echo $PATH

#echo $PYTHONPATH	


#PATH=$HOME/fenics_17_1/lib/python3.5/site-packages:$PATH
## Compile and execute the preprocessor code
## This will convert meshDataFile (*.dat) to fenicsDataFiles (*.xml) in ../data/
mpiexec -n 32 python -m mpi4py  gmshData_fenicsXML_parallel_3D.py

## Compile and execute the FEniCS/dolfin assembly routines
## This will create assemlby matrices-vectors in ../data/Amats/
mpiexec -n 32 python -m mpi4py elasticity3D_stochasticDDM_parallel_twolevel.py

exit
