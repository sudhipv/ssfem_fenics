## Call mesh generator & decomposer
## Call KLE/PCE data generator and importer
cd ../../data/klePceData/
gfortran KLE_PCE_Data.F90
./a.out

cd ../../data/meshData/
sh preprocess.sh

## Use meshio to convet GMSH mesh data to FEniCS-XML data
cd ../../external/dolfin/src/
python3 pointElements_xml.py

### select the number of partitions
#NP=4
#
### create mesh file::
##gmsh -2 square.geo -part $NP -o foo.msh
##gmsh -2 foo.geo -part $NP -o foo.msh
#
### create mesh data::
#gfortran preprocmesh1_AD.F90 -O2 -o ./a.out;./a.out
#gfortran preprocmesh2_AD.F90 -O2 -o ./a.out;./a.out
#
### create measurement data::
#gfortran preprocmesh_meas1.F90 -O2 -o ./a.out;./a.out
#gfortran preprocmesh_meas2.F90 -O2 -o ./a.out;./a.out

## mesh visualization::
#gmsh gmsh.msh
