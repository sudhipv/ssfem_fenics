## Change to respective folder
cd ../data/meshData/

echo "GMSH: Select one of the following options: "
echo " 1 for 2D-Poisson "
echo " 2 for 3D-Poisson "
echo " 3 for 3D-Elasticity "
read input1

if [ $input1 == 1 ]
then
    echo "GMSH: 2D-Poisson Selected"
    cp geofiles/fooServer2D.geo foo.geo
    cp jobscritps/generate_data_beluga2D.sh generate_data_graham.sh
elif [ $input1 == 2 ]
then
    echo "GMSH: 3D-Poisson Selected"
    cp geofiles/slenderBeam.geo foo.geo
    cp jobscritps/generate_data_beluga3D.sh generate_data_graham.sh
elif [ $input1 == 3 ]
then
    echo "GMSH: 3D-Elasticity Selected"
    cp geofiles/slenderBeam.geo foo.geo
    cp jobscritps/generate_data_beluga3D.sh generate_data_graham.sh
else
echo "STOP: Wrong input"
exit
fi

###*** 3D MESH :: submit a job to create mesh file and mesh decomposition data
## change the lc parameter to control mesh desnity
## first lc=current lc and second lc=desiredlc
## for cluster use lc=0.01 with (0.2,0.2,1.0) and lc=0.05 with (1,1,5) for 30K
sed -i -e 's/lc=0.1/lc=0.01/g' foo.geo
sed -i -e 's/time=0-00:30/time=0-00:30/g' generate_data_graham.sh
sed -i -e 's/NP=32/NP=960/g' generate_data_graham.sh
## select wb=1 for wirebasket and wb=0 for vertex grid
sed -i -e 's/wb=1/wb=1/g' preprocmesh3D2_AD.F90

## This is currently required as GMSH is not working in Cedar under job-script
echo "GMSH: Enter 1 for Graham/ 2 for Cedar"
read input2
if [ $input2 == 2 ]
then
    cp jobscritps/generate_mesh_cedar.sh generate_gmsh.sh
    sed -i -e 's/$NP/32/g' generate_gmsh.sh
    sh generate_gmsh.sh
    echo "** GMSH: Cedar is selected **"
elif [ $input2 == 1 ]
then
    echo "** GMSH: Graham is selected **"
else
    echo "STOP: Wrong input"
    exit
fi

## Job submission method
sbatch generate_data_graham.sh
