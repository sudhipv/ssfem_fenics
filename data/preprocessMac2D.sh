## KLE/PCE data generator and importer
## This take one command line input after number of dimensions (RVs)
cd klePceData/
gfortran KLE_PCE_Data_commandLineArg.F90
./a.out
cd -

## GMSH data generator
## Change to respective folder

cd meshData/
cp geofiles/foo.geo foo.geo


## Adjust the mesh density parameter "lc"
## lc=old;/lc=new for foo*.geo

sed -i -e 's/lc=0.09;/lc=0.05;/g' foo.geo

## select the number of partitions

NP=4

## create mesh file::
gmsh -2 foo.geo -part $NP -o foo.msh

## create mesh data::
gfortran preprocmesh1_AD.F90 -O2 -o ./a.out;./a.out
gfortran preprocmesh2_AD.F90 -O2 -o ./a.out;./a.out

## create measurement data::
#gfortran preprocmesh_meas1.F90 -O2 -o ./a.out;./a.out
#gfortran preprocmesh_meas2.F90 -O2 -o ./a.out;./a.out

## remove unnecessory files created by "sed"
rm *.geo-e
cd -

## FEniCS data generator
cd ../external/dolfin/src/

echo "Enter 1 for Parallel / 2 for Serial"
read input1
if [ $input1 == 1 ]
then
    export EVENT_NOKQUEUE=1
    ##var=`cat ./../../data/meshData/num_partition.dat`
    NP=$(<../../../data/meshData/num_partition.dat)
    echo "Using $NP processors"
    mpiexec -n $NP --oversubscribe python3 gmshData_fenicsXML_parallel_2D.py
elif [ $input1 == 2 ]
then
    python3 gmshData_fenicsXML_2D.py
else
    echo "STOP: Wrong Input"
fi
echo "NEXT: sh processFEniCS_Mac.sh"
echo "=========================================================="
