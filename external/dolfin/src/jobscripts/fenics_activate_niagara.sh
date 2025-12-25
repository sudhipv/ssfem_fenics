#!/bin/bash
module load CCEnv
module load StdEnv
module load hdf5-mpi/1.8.18 boost eigen python/3.5 scipy-stack/2017b petsc/3.7.5 fftw-mpi/3.3.6

source ~/software/dolfin_17_1/share/dolfin/dolfin.conf
source ~/packages/fenics_17_1/bin/activate
