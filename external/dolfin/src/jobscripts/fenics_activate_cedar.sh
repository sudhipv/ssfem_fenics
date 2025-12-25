#!/bin/bash

##### Assuming all the below modules are loaded by default ############

#################
#Currently Loaded Modules:
#  1) nixpkgs/16.09   (S)      3) gcccore/.5.4.0  (H)   5) ifort/.2016.4.258 (H)   7) openmpi/2.1.1 (m)
#    2) imkl/11.3.4.258 (math)   4) icc/.2016.4.258 (H)   6) intel/2016.4      (t)   8) StdEnv/2016.4 (S)

#      Where:
#         S:     Module is Sticky, requires --force to unload or purge
#	    m:     MPI implementations / Implémentations MPI 
#              math:  Mathematical libraries / Bibliothèques mathématiques
#	          t:     Tools for development / Outils de développement
#		     H:                Hidden Module


#################

module load StdEnv/2016.4


module load hdf5-mpi/1.8.18 boost eigen python/3.5 scipy-stack/2017b petsc/3.7.5 fftw-mpi/3.3.6
source ~/software/dolfin_17_1/share/dolfin/dolfin.conf
source ~/fenics_17_1/bin/activate
