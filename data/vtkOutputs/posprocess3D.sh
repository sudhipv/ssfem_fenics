#!/bin/bash

#Get data from Server (Graham/Cedar)

##datapath="/scratch/sudhipv/poisson3D"

echo "Enter 1 for Graham / 2 for Cedar/ 3 for Niagara/ 4 for Beluga"
read input1
if [ $input1 == 1 ]
then
    datapath="/scratch/sudhipv/ssfem_fenics"
    echo "Copying data from" ${datapath}
    scp -r sudhipv@graham.computecanada.ca:${datapath}/\{data/vtkOutputs/*.dat,data/vtkOutputs/*.m,data/meshData/points.dat,data/meshData/tetrahedrons4c.dat,data/klePceData/pcedata.dat,data/slurm/*.out\} .
elif [ $input1 == 2 ]
then
    datapath="/scratch/sudhipv/sudhipv/ssfem_fenics"
    echo "Copying data from" ${datapath}
    scp -r sudhipv@cedar.computecanada.ca:${datapath}/\{data/vtkOutputs/*.dat,data/vtkOutputs/*.m,data/meshData/points.dat,data/meshData/tetrahedrons4c.dat,data/klePceData/pcedata.dat,data/slurm/*.out\} .
elif [ $input1 == 3 ]
then
    datapath="/scratch/a/asarkar/sudhipv/ssfem_fenics"
    echo "Copying data from" ${datapath}
    scp -r sudhipv@niagara.computecanada.ca:${datapath}/\{data/vtkOutputs/*.dat,data/vtkOutputs/*.m,data/meshData/points.dat,data/meshData/tetrahedrons4c.dat,data/klePceData/pcedata.dat,data/slurm/*.out\} .
elif [ $input1 == 4 ]
then
    datapath="/scratch/sudhipv/ssfem_fenics"
    echo "Copying data from" ${datapath}
    scp -r sudhipv@beluga.computecanada.ca:${datapath}/\{data/vtkOutputs/*.dat,data/vtkOutputs/*.m,data/meshData/points.dat,data/meshData/tetrahedrons4c.dat,data/klePceData/pcedata.dat,data/slurm/*.out\} .
else
    echo "STOP: Wrong Input"
fi

##matlabTerminal -r writeVtk   # not working, need to sort this issue
				    ##matlabTerminal -r writeVecVtk
echo "open 'matlabTerminal' and run 'writeVtk' for 3DPoisson or 'writeVecVtk' for 3D elasticity"

