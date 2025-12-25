#!/bin/bash
echo "====================================================================="
echo "FOLLOW THESE STEPS: '$ cd external/dolfin/src/' and then run "
echo "Eg1: For 2D-onelevel: '$ python poisson2D_stochasticDDM_twolevel.py' "
echo "Eg2: For 3D-twolevel: '$ python poisson3D_stochasticDDM_twolevel.py' "
echo "Do the same for other python executables available in same folder "
echo "============================== PARALLE =============================="
echo "Eg1: mpiexec -n 4 python poisson2D_stochasticDDM_parallel_twolevel.py"
echo "Eg2: mpiexec -n 8 python poisson3D_stochasticDDM_parallel_twolevel.py"
echo "====================================================================="

###*** Poisson2D/StochasticDDM/Twolevel
#python poisson3D_stochasticDDM_twolevel.py

## KLE/PCE data generator and importer
fenicsproject start dolfin
