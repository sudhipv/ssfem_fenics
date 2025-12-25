## SRC:

### poisson3D_stoDDMAssembly.py (For EACH PCE TERM):
Written using new assemlby procedure: Here we apply BC on FEniCS assembly mat-vec as oppose to old where we used penalty factor (Fortran style)
input Cijk: ../data/inputs/cijk33.mat
Input Mesh: ../data/foo*.xml
output Assembly: ../data/Amats/subdom*

The primary difference is we apply BC using FEniCS procedure instead in new comared with old where we use panulty factor (like in Fortran package)


### poisson3D_pointsCells_fenicsXML.py (Need python3)
This code builds FEniCS(XML) mesh using meshdata from Fortran preprocessor:
input: ../../../data/meshData/rpoints00*.dat
input: ../../../data/meshData/rtetrahedrons00*.dat
output: '../data/foo*.xml

#### stoDDM_poisson.py & pointElements_xml.py
For stochastic DDM assembly (each PCE mode) for 2D Poisson problem
pointElements_xml.py for 2D GMSH mesh to FEniCS mesh converter

#### poisson3D_detDDMAssembly.py
For deterministic DDM matrices using FEniCS assembly procedure

## data:
1. FEniCS readable meshfile (xml) built using Fortran preprocessory data
2. Amats: Deterministic DD matrices for each subdomain for each input PCE.
3. Bvecs: Should contains DD vectors for each subdomain for each input PCE.
4. meshData: code can modified to read meshdata from this folder as well.
5. inputs: other input data such as cijkOD.mat is stored and read from here.

NOTE: cijkOD.mat: these files are generated using UQTk based KLE/PCE package in Matlab.
They store non-zero tripple product and their ijk-indices for given order(O) and dimension(D).


## Instuctions for Graham/Cedar:
NOTE: for first time > source fenics_activate.sh (Then install meshio, xml)
(fenics) [ajitd@cedar5 sgfem_ddm_hpc]$ pip3 install meshio
(fenics) [ajitd@cedar5 src]$ pip3 install lxml

1.  from HPCSstoDDM/data/meshData: $ generate_data_graham.sh
This read input "foo.geo" from the same folder and generate foo.msh
and all DDM meshdata in the same folder

2. for dolfin/src: $ run_fenics.sh
This code first read the meshdata from sgfem_ddm_hpc/data/meshData and create "xml" for each subdomain, then invokes main code to build stochastic DDM matrices using FEniCS
