#!/bin/bash
cd ../../external/dolfin/src/jobscripts/

echo "FEniCS: Enter 1 for Niagara"
read input1

if [ $input1 == 1 ]
then
    echo "FEniCS: Selected Niagara"
    cp fenics_activate_niagara_intel2017.sh ../fenics_activate.sh
else
    echo "STOP: Wrong input"
    exit
fi

echo "FEniCS: Select one of the following options: "
echo " 1 for 2D-Poisson two-level parallel "
echo " 2 for 3D-Poisson two-level parallel "
echo " 3 for 3D-Elasticity two-level parallel "
read input2
if [ $input2 == 1 ]
then
    echo "** FEniCS: 2D-Poisson Two-level Parallel Selected **"
    cp run_fenics_poisson2D_stochasticDDM_parallel_twolevel_NIAGARA.sh ../run_fenics.sh
elif [ $input2 == 2 ]
then
    echo "** FEniCS: 3D-Poisson Two-level Parallel Selected **"
    cp run_fenics_poisson3D_stochasticDDM_parallel_twolevel_NIAGARA.sh ../run_fenics.sh
elif [ $input2 == 3 ]
then
    echo "** FEniCS: 3D-Elasticity Two-level Parallel Selected **"
    cp run_fenics_elasticity3D_stochasticDDM_parallel_twolevel_NIAGARA.sh ../run_fenics.sh
else
    echo "STOP: Wrong input"
    exit
fi


cd ..
### Use next line to increase run-time for FEniCS code
sed -i -e 's/time=0-00:15/time=0-00:30/g' run_fenics.sh
sed -i -e 's/#SBATCH --mem-per-cpu=7700M/##SBATCH --mem-per-cpu=3700M/g' run_fenics.sh
sed -i -e 's/nodes=1/nodes=2/g' run_fenics.sh
sed -i -e 's/--tasks-per-node=32/--ntasks-per-node=40/g' run_fenics.sh
sed -i -e 's/-n 32/-n 80/g' run_fenics.sh

### To avoid recreating FEniCS-xml files uncomment next line gmshData_fenicsXML_3D.py
#sed -i -e 's/python gmshData_fenicsXML_3D.py/#python gmshData_fenicsXML_3D.py/g' run_fenics3D.sh

###*** Submit-Job method: same for all cases
sbatch run_fenics.sh
