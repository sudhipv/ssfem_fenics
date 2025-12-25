#!/bin/bash
echo "Enter 2 for 2D / 3 for 3D"
read input1

if [ $input1 == 2 ]
then
    echo "Entered in 2D"
    echo "Select one of the following"
    echo "1-Poisson onelevel"
    echo "2-Poisson twolevel"
    echo "3-Elasticity onelevel-NOT READY"
    echo "4-Elasticity twolevel-NOT READY"
    read input2
    if [ $input2 == 1 ]
    then
        echo "Selected: 2D Poisson onelevel"
        cd ../src/2Dpoisson_onelevel/clusterMakefiles/
        cp makefile_mac ../makefile
    elif [ $input2 == 2 ]
    then
        echo "Selected: 2D Poisson twolevel"
        cd ../src/2Dpoisson_twolevel/
        cp clusterMakefiles/makefile_mac ../makefile
    elif [ $input2 == 3 ]
    then
        echo "Selected: 2D Elasticity onelevel"
    elif [ $input2 == 4 ]
    then
        echo "Selected: 2D Elasticity twolevel"
    else
    echo "STOP: Wrong Input"
    fi

elif [ $input1 ==  3 ]
then
    echo "Entered in 3D"
    echo "Select one of the following"
    echo "1-Poisson onelevel"
    echo "2-Poisson twolevel"
    echo "3-Elasticity onelevel"
    echo "4-Elasticity twolevel"
    read input2
    if [ $input2 == 1 ]
    then
echo "Selected: 3D Poisson onelevel"
        cd ../src/3Dpoisson_onelevel/
        cp clusterMakefiles/makefile_mac ../makefile
    elif [ $input2 == 2 ]
    then
        echo "Selected: 3D Poisson twolevel"
        cd ../src/3Dpoisson_twolevel/
        cp clusterMakefiles/makefile_mac ../makefile
    elif [ $input2 == 3 ]
    then
        echo "Selected: 3D Elasticity onelevel"
        cd ../src/3Delasticity_onelevel/
        cp clusterMakefiles/makefile_mac ../makefile
    elif [ $input2 == 4 ]
    then
        echo "Selected: 3D Elasticity twolevel"
        cd ../src/3Delasticity_twolevel/
        cp clusterMakefiles/makefile_mac ../makefile
    else
    echo "STOP: Wrong Input"
    fi
else
    echo "STOP: Wrong Input"
fi

echo "Enter 0 for Nothing / 1 for VTK/Dat-outputs "
read input4

if [ $input4 == 1 ]
then
    echo "Selected VTK/Dat-outputs"
    sed -i -e 's/outputFlag=0/outputFlag=1/g' main.F90
elif [ $input4 == 0 ]
then
    echo "No output selected"
    sed -i -e 's/outputFlag=1/outputFlag=0/g' main.F90
else
    echo "STOP: Wrong input"
    exit
fi

rm *-e

####*** Submit-Job method: same for all case
NP=$(<../../data/meshData/num_partition.dat)
echo "Using $NP processors"
make all
/Users/sudhipv/documents/PETSc/petsc-3.7.5/arch-darwin-c-debug/bin/mpiexec -np $NP ./a.out

echo "====================================================================="
echo "FOR 3D Poisson or 3D elasticity"
echo "$ sh postprocessMac3D.sh"
echo "FOR 2D check *.vtk files in vtkOutputs"
echo "====================================================================="

