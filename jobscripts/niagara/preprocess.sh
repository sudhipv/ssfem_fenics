## KLE/PCE data
#module load gcc

source ./loadmodules.sh

echo "KLE/PCE Data: "
cd ../../data/klePceData/
gfortran KLE_PCE_Data_commandLineArg.F90
./a.out
cd -

## GMSH/DDM data
cd ../../data/meshData/

echo "GMSH: Select one of the following options: "
echo " 1 for 2D-mesh "
echo " 2 for 3D-mesh "

read input1

if [ $input1 == 1 ]
then
    echo "GMSH: 2D-mesh Selected"
    cp geofiles/fooServer2D.geo foo.geo
    cp jobscritps/generate_data_niagara2D.sh generate_data.sh
elif [ $input1 == 2 ]
then
    echo "GMSH: 3D-mesh Selected"
    cp geofiles/slenderBeam.geo foo.geo
    cp jobscritps/generate_data_niagara3D.sh generate_data.sh
else
    echo "STOP: Wrong input"
exit
fi

###*** 3D MESH :: submit a job to create mesh file and mesh decomposition data
## change the lc parameter to control mesh desnity
## first lc=current lc and second lc=desiredlc
## for cluster use lc=0.01 with (0.2,0.2,1.0) and lc=0.05 with (1,1,5) for 30K
sed -i -e 's/lc=0.1/lc=0.01/g' foo.geo
sed -i -e 's/NP=32/NP=320/g' generate_data.sh
## select wb=1 for wirebasket and wb=0 for vertex grid
sed -i -e 's/wb=1/wb=1/g' preprocmesh3D2_AD.F90

## Job submission method
bash generate_data.sh
#module unload gcc
#sbatch generate_data_graham.sh
cd -
