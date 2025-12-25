#!/bin/bash



# Assuming all the below modules are loaded by default ##

##Currently Loaded Modules:
#  1) nixpkgs/16.09   (S)      3) gcccore/.7.3.0  (H)   5) ifort/.2018.3.222 (H)   7) openmpi/3.1.2 (m)
#    2) imkl/2018.3.222 (math)   4) icc/.2018.3.222 (H)   6) intel/2018.3      (t)   8) StdEnv/2018.3 (S)

#      Where:
#         S:     Module is Sticky, requires --force to unload or purge
#	    m:     MPI implementations / Implémentations MPI
#	       math:  Mathematical libraries / Bibliothèques mathématiques
#	          t:     Tools for development / Outils de développement
#		     H:                Hidden Module

###############

module load intel/2018.3
module load fftw-mpi/3.3.8
module load boost
module load hdf5-mpi/1.10.3
module load petsc/3.7.5
module load eigen python/3.5 scipy-stack
module load cmake

source $HOME/software/dolfin_17_r/share/dolfin/dolfin.conf
source $HOME/fenics_17_r/bin/activate
