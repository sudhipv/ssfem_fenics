#!/bin/bash


module purge --force
module load arch/avx2
module load nixpkgs/16.09 imkl/11.3.4.258 intel/2016.4 openmpi/2.1.1 StdEnv/2016.4


module load hdf5-mpi/1.8.18 boost eigen python/3.5 scipy-stack/2017b petsc/3.7.5 fftw-mpi/3.3.6
### Loading cmake to detect Boost
module load cmake
### Loading ipp - intel parallel processosing library
module load ipp/9.0.4

source $HOME/software/dolfin_17_cedarcopy/share/dolfin/dolfin.conf
source $HOME/fenics_17_cedarcopy/bin/activate
