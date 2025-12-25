cd ../../src

echo " MPI-Job: Select one of the following options: "
echo " 1 for 2D-Poisson two-level"
echo " 2 for 3D-Poisson two-level"
echo " 3 for 3D-Elasticity two-level"
read input2
if [ $input2 == 1 ]
then
    echo "MPI-Job: 2D-Poisson Two-level Selected"
    cd 2Dpoisson_twolevel/clusterMakefiles/
elif [ $input2 == 2 ]
then
    echo "MPI-Job: 3D-Poisson Two-level Selected"
    cd 3Dpoisson_twolevel/clusterMakefiles/
elif [ $input2 == 3 ]
then
    echo "MPI-Job: 3D-Elasticity Two-level Selected"
    cd 3Delasticity_twolevel/clusterMakefiles/
else
    echo "STOP: Wrong input"
    exit
fi


echo "MPI-Job: Enter 1 for Graham"
read input1

if [ $input1 == 1 ]
then
    echo "Selected Graham"
    cp makefile_graham ../makefile
    cp run_ddm_graham.sh ../run_ddm.sh
else
    echo "STOP: Wrong input"
    exit
fi


###*** Personalize each job based on requirements
cd ..
sed -i -e 's/time=0-00:10/time=0-9:00/g' run_ddm.sh
sed -i -e 's/mem-per-cpu=7700M/mem-per-cpu=7700M/g' run_ddm.sh
sed -i -e 's/tasks-per-node=32/tasks-per-node=32/g' run_ddm.sh
sed -i -e 's/nodes=1/nodes=10/g' run_ddm.sh
sed -i -e 's/-np 32/-np 320/g' run_ddm.sh
## Below d->dims, lc-> mesh density parameter, p->partitions
sed -i -e 's/job-name=DP_nRV_nNodes_nParts/job-name=3DE_5RV_30KNodes_320Parts/g' run_ddm.sh
sed -i -e 's#output=%x-%j.out#output=/scratch/sudhipv/ssfem_fenics/data/slurm/%j-%x.out#g' run_ddm.sh
#sed -i -e 's/a.out/a.out -ksp_converged_reason/g' run_ddm.sh

echo "MPI-Job: Enter 0 for Nothing / 1 for VTK/Dat-outputs "
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

###*** Submit-Job method: same for all cases
sbatch run_ddm.sh
