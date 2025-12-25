#!/bin/bash

##### Assuming NiaEnv/2018a loaded by default 

module load intel/2018.2
module load intelmpi/2018.2
#module load intelpython3/2018.2
module load python/3.6.5
module load anaconda3/5.1.0
module load cmake/3.11.1
module load boost/1.66.0
module load eigen/3.3.4
module load hdf5-mpi/1.8.20
module load swig/3.0.12
module load netcdf-mpi/4.6.1
module load trilinos/12.12.1
module load gmp/6.1.2
module load mpfr/4.0.1
module load petsc/3.8.4

source ~/software/dolfin_intel/share/dolfin/dolfin.conf
source ~/packages/fenics_intel/bin/activate



#module --force purge CCEnv
#module list
#module load NiaEnv
#module load intel/2018.3 intelmpi intelpython3





