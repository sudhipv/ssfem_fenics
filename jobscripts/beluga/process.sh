#!/bin/bash
cd ../../external/dolfin/src/

echo "FEniCS: Select one of the following options: "
echo " 1 for 2D-Poisson two-level parallel "
echo " 2 for 3D-Poisson two-level parallel "
echo " 3 for 3D-Elasticity two-level parallel "
read input2
if [ $input2 == 1 ]
then
    echo "** FEniCS: 2D-Poisson Two-level Parallel Selected **"
    cp jobscripts/process_poisson2D_stochasticDDM_parallel_twolevel.sh run_fenics.sh

elif [ $input2 == 2 ]
then
    echo "** FEniCS: 3D-Poisson Two-level Parallel Selected **"
    cp jobscripts/process_poisson3D_stochasticDDM_parallel_twolevel.sh run_fenics.sh
elif [ $input2 == 3 ]
then
    echo "** FEniCS: 3D-Elasticity Two-level Parallel Selected **"
    cp jobscripts/process_elasticity3D_stochasticDDM_parallel_twolevel.sh run_fenics.sh
else
    echo "STOP: Wrong input"
exit
fi

echo "FEniCS: Enter / 1 for Graham / 2 for Cedar / 3 for Niagara "
read input1

if [ $input1 == 1 ]
then
    echo "FEniCS: Selected Graham"
    cp jobscripts/fenics_activate_graham.sh fenics_activate.sh
elif [ $input1 == 2 ]
then
    echo "FEniCS: Selected Cedar"
    cp jobscripts/fenics_activate_cedar.sh fenics_activate.sh
elif [ $input1 == 3 ]
then
    echo "FEniCS: Selected Niagara"
    cp jobscripts/fenics_activate_niagara.sh fenics_activate.sh
    #cp jobscripts/fenics_activate_niagaraModule.sh fenics_activate.sh
    sed -i -e 's/##SBATCH --ntasks/#SBATCH --ntasks/g' run_fenics.sh
    sed -i -e 's/#SBATCH --mem-per-cpu/##SBATCH --mem-per-cpu/g' run_fenics.sh
    sed -i -e 's/#SBATCH --tasks-per-node/##SBATCH --tasks-per-node/g' run_fenics.sh
else
    echo "STOP: Wrong input"
exit
fi

### Use next line to increase run-time for FEniCS code
sed -i -e 's/-n 32/-n 120/g' run_fenics.sh
sed -i -e 's/nodes=1/nodes=3/g' run_fenics.sh
sed -i -e 's/ntasks=40/ntasks=120/g' run_fenics.sh
sed -i -e 's/time=0-00:15/time=0-03:20/g' run_fenics.sh
sed -i -e 's/mem-per-cpu=7700M/mem-per-cpu=3700M/g' run_fenics.sh
sed -i -e 's/tasks-per-node=32/tasks-per-node=40/g' run_fenics.sh
sed -i -e 's/job-name=nRV_nNodes_nParts/job-name=3RV_30KNodes_120Parts/g' run_fenics.sh

### To avoid recreating FEniCS-xml files uncomment next line gmshData_fenicsXML_3D.py
#sed -i -e 's/python gmshData_fenicsXML_3D.py/#python gmshData_fenicsXML_3D.py/g' run_fenics3D.sh

cd -
cd ../src

if [ $input2 == 1 ]
then
    echo "2D-Poisson Two-level Selected"
    cd 2Dpoisson_twolevel/clusterMakefiles/
elif [ $input2 == 2 ]
then
    echo "3D-Poisson Two-level Selected"
    cd 3Dpoisson_twolevel/clusterMakefiles/
elif [ $input2 == 3 ]
then
    echo "3D-Elasticity Two-level Selected"
    cd 3Delasticity_twolevel/clusterMakefiles/
else
    echo "STOP: Wrong input"
    exit
fi


if [ $input1 == 1 ]
then
    echo "Selected Graham"
    cp makefile_graham ../makefile
    cp petsc_activate_graham.sh ../petsc_activate.sh
elif [ $input1 == 2 ]
then
    echo "Selected Cedar"
    cp makefile_cedar ../makefile
    cp petsc_activate_graham.sh ../petsc_activate.sh
elif [ $input1 == 3 ]
then
    echo "Selected Niagara"
    cp makefile_niagara ../makefile
    cp petsc_activate_niagara.sh ../petsc_activate.sh
else
    echo "STOP: Wrong input"
    exit
fi

cd ..
echo "MPI-Job: Enter 0 for No-outputs / 1 for VTK/Dat-outputs "
read input3

if [ $input3 == 1 ]
then
    echo "** MPI-Job: VTK/Dat-outputs selected **"
    sed -i -e 's/outputFlag=0/outputFlag=1/g' main.F90
elif [ $input3 == 0 ]
then
    echo "** MPI-Job: No VTK/Dat outputs selected **"
    sed -i -e 's/outputFlag=1/outputFlag=0/g' main.F90
else
    echo "STOP: Wrong input"
    exit
exit
fi

cd ../../external/dolfin/src/
mv run_fenics.sh run_fenics_mpi_job.sh
echo "verify and submit using * sbatch run_fenics_mpi_job.sh * "

####*** Submit-Job method: same for all cases
#sbatch run_fenics_mpi_job.sh
