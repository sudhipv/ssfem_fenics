#!/bin/bash
cd ../external/dolfin/src/jobscripts/

echo "Enter 1 for Graham / 2 for Cedar /3 for Niagara "
read input1

if [ $input1 == 1 ]
then
    echo "Selected 'Graham'"
    cp fenics_activate_graham.sh ../fenics_activate.sh
elif [ $input1 == 2 ]
then
    echo "Selected 'Cedar'"
    cp fenics_activate_cedar.sh ../fenics_activate.sh
elif [ $input1 == 3 ]
then
    echo "Selected 'Niagara'"
    cp fenics_activate_niagara.sh ../fenics_activate.sh
else
 
    echo "STOP: Wrong input"
fi

echo "Select one of the following options: "
echo " 1 for 2D Poisson two-level"
echo " 2 for 3D Poisson two-level"
echo " 3 for 3D Elasticity two-level"
read input2
if [ $input2 == 1 ]
then
    echo "Selected 2D Poisson Two-level"
    cp run_fenics_poisson2D_stochasticDDM_twolevel.sh ../run_fenics3D.sh
elif [ $input2 == 2 ]
then
    echo "Selected 3D Poisson Two-level"
    cp run_fenics_poisson3D_stochasticDDM_twolevel.sh ../run_fenics3D.sh
elif [ $input2 == 3 ]
then
    echo "Selected 3D Elasticity Two-level"
    cp run_fenics_elasticity3D_stochasticDDM_twolevel.sh ../run_fenics3D.sh
else
    echo "STOP: Wrong input"
fi


cd ..
### Use next line to increase run-time for FEniCS code



if [ $input1 == 3 ]
then
sed -i -e 's/#SBATCH --mem-per-cpu=7700M/##SBATCH --mem-per-cpu=7700M/g' run_fenics3D.sh
sed -i -e 's/--tasks-per-node=1/--ntasks-per-node=40/g' run_fenics3D.sh
sed -i -e 's#output=%j-3Dfenics.out#output=/scratch/a/asarkar/sudhipv/ssfem_fenics/data/slurm/3Dfenics-%j.out#g' run_fenics3D.sh
fi

sed -i -e 's/time=0-00:15/time=0-02:00/g' run_fenics3D.sh
sed -i -e 's/mem-per-cpu=7700M/mem-per-cpu=3700M/g' run_fenics3D.sh
sed -i -e 's/--tasks-per-node=1/--tasks-per-node=32/g' run_fenics3D.sh
sed -i -e 's#output=%j-3Dfenics.out#output=$SCRATCH/sfem_fenics/data/slurm/3Dfenics-%j.out#g' run_fenics3D.sh

### To avoid recreating FEniCS-xml files uncomment next line gmshData_fenicsXML_3D.py
#sed -i -e 's/python gmshData_fenicsXML_3D.py/#python gmshData_fenicsXML_3D.py/g' run_fenics3D.sh

###*** Submit-Job method: same for all cases
sbatch run_fenics3D.sh


###*** Poisson3D/StochasticDDM/Twolevel
#cp run_fenics_poisson3D_stochasticDDM_twolevel.sh ../run_fenics3D.sh

###*** Poisson3D/StochasticDDM/Onelevel
#cp run_fenics_poisson3D_stochasticDDM_onelevel.sh ../run_fenics3D.sh

###*** Poisson2D/StochasticDDM/Twolevel
#cp run_fenics_poisson2D_stochasticDDM_twolevel.sh ../run_fenics2D.sh

###*** Poisson2D/StochasticDDM/Onelevel
#cp run_fenics_poisson2D_stochasticDDM_onelevel.sh ../run_fenics2D.sh
