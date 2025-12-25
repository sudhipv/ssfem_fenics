gmsh -2 square.geo -part 8 -o gmsh.msh

gfortran preprocmesh1_AD.F90 -O2 -o ./a.out;./a.out
gfortran preprocmesh2_AD.F90 -O2 -o ./a.out;./a.out
gfortran preprocmesh_meas1.F90 -O2 -o ./a.out;./a.out
gfortran preprocmesh_meas2.F90 -O2 -o ./a.out;./a.out

#cp klePceData/multiIndex0003.dat ../VerBeta2S/

#gmsh gmsh.msh
