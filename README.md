# Scalable Solver for Spectral Stochastic FEM using Domain Decomposition with One/Two-Level Preconditioner

This solver have capabilities to solve stochastic steady-state diffusion equation and static linear elasticity (2D & 3D) using Domain Decomposition Methods. 

This SSFEM code has different branches master, cedar, graham, niagara, beluga etc. 
Master branch is just the master copy of the whole code. 


### Components
1. FEniCS: Employed for deterministic FE matrix-vector assembly
2. PETSc: Employed for (local) sparse storage of stochastic FE matrix-vector and associated linear algebraic calculations
3. MPI: Employed for parallel processing of PETSc local routines

### Required packages: (Version Tested) ##
* FORTRAN    :: GNU GCC 
* MPI        :: Open MPI 
* GMSH       :: 3.0.4 (Use only this version)
  * Binary : http://gmsh.info/bin/
  * Source : https://gitlab.onelab.info/gmsh/gmsh/blob/master/README.txt
* PETSc      :: 3.7.5
    * Local Install  : https://www.mcs.anl.gov/petsc/documentation/installation.html
    * Remote Install : Refer to 'PETSC_Install.pdf' under docs folder
* FEniCS     :: 2017.1.0 (Use only this version)
    * Local installation : https://fenics.readthedocs.io/en/latest/installation.html#
      Remote Install:                    
    *  Follow instructions from "fenics_install.pdf" available in docs/fenics folder. One could also check corresponding script files for installation and modify accordingly.
    *  For FEniCS installation issues and updates refer to [computeCanada wiki page for FEniCS](https://docs.computecanada.ca/wiki/FEniCS)
    *  You also need to install following small packages: "meshio" (pip3 install meshio), "lxml" (pip3 install lxml), "mpi4py"(pip3 install mpi4py)

* PARAVIEW   :: 5.4.0
* PYTHON     :: 3
* MATLAB     :: 2017


## Running in Local machine (Tested in MacBookPro)

### Step-1:: preprocess :: KLE/PCE & MeshData Generation
* from top level project directoy open terminal (1)
* $ cd data/klePceData
* $ gfortran KLE_PCE_Data_commandLineArg.F90
* $ ./a.out                :: enter the number of nDim (RVs) intend to use
* $ cd ../data/meshData                     :: to prepare mesh data

### 3D Domain
* preprocess3D.sh                      :: here you can adjust mesh density (LC) and partitions (NP)
* $ sh preprocess3D.sh              :: should create 3D mesh data in same folder (.dat)
* $ cd ../../external/dolfin/src/
* $ python3 gmshData_fenicsXML_3D.py   :: should create FEniCS-xml data using GMSH-msh data

### 2D Domain
* preprocess2D.sh                        :: here you can adjust mesh density (LC) and partitions (NP)
* $ sh preprocess2D.sh                   :: should create 2D GMSH-mesh data in same folder (.dat)
* $ cd ../../external/dolfin/src/
* $ python3 gmshData_fenicsXML_2D.py   :: should create FEniCS-xml data using GMSH-msh data

### Step-2:: process:: FEM using FEniCS
* from top level project directoy open terminal (2)
* $ fenicsproject start dolfin
* $ cd external/dolfin/src/

### 3D Domain: one level & two level DDM
* $ python poisson3D_stochasticDDM_onelevel.py     :: create onelevel 3D DD-matrices in /data/Amats/
* $ python poisson3D_stochasticDDM_twolevel.py     :: create twolevel 3D DD-matrices in /data/Amats/

### 2D Domain: onelevel & twolevel DDM
* $ python poisson2D_stochasticDDM_onelevel.py     :: create onelevel 2D DD-matrices in /data/Amats/
* $ python poisson2D_stochasticDDM_twolevel.py     :: create twolevel 2D DD-matrices in /data/Amats/

### Step-3:: process:: SSFEM using Fortran/MPI/PETSc
* from top level project directoy open terminal (3)
* $ cd src/2Dpoisson_*  or  cd src/3Dpoisson_*  or  cd src/3Delasticity_*     :: go to intended solvers-scripts
* $ cp clusterMakefiles/makefile_mac makefile
* $ make all
* $ petscexec -np $NP ./a.out      :: should print number of iterations & write outputs in data/vtkOutputs/

### Step-4:: postprocess:: Visualize output using Matlab/ParaView
* from top level project directory terminal (4)
* $ cd data/vtkOutputs    ::
#### For poisson3D
* $ matlabTerminal -r writeVtk
#### For elasticity3D
* $ matlabTerminal -r writeVecVtk
* $ paraview out_*.vtk          :: open *.vtk files with paraview to check outputs (Note: For 2D vtk files are created directly)


### Special case: 2D Domain (Comparison): All in one steps
* $ vi preprocess2D.sh    :: to adjust number of partitions & mesh density
* $ sh preprocess2D.sh                   :: should create DD-mesh data in same folder (.dat)
* $ cd ../../external/dolfin/src/
* $ python3 pointElements_xml.py :: should create DD-mesh data in /data (.xml)
* $ fenicsproject start dolfin
* $ cd external/dolfin/src/
* $ python stoDDM_poisson.py      :: should create DD-matrices in /data/Amats/subdom000* (.dat)
* $ cd src/onelevel or  cd src/twolevel
* $ make all
* $ petscexec -np $NP ./a.out        :: should print number of iterations & write outputs in data/vtkOutputs/


## Running Code in REMOTE CLUSTERS (CEDAR, GRAHAM, NIAGARA, BELUGA)

### Step-1 :: Get the repository from bitbucket (Note: private repo, needs permission to access)
* $ git clone https://sudhipv@bitbucket.org/sudhipv/ssfem_fenics.git
* $ git fetch && git checkout "branch" (by default it is master branch).
* $ cd ssfem_fenics/


### Step-2 :: preprocess & process:: GMSH/FEniCS/MPI/PETSc 
### NOTE : Please ensure to change the paarmeters inside each script file according to the instructions numbered below it. There are places one needs to adapt the installation paths or output file paths.
* $ cd jobscripts/        :: from top level project directory
* $ cd machine-name/ (cedar, niagara, beluga,graham) 
* $ sh preprocess.sh      :: This script file takes input from the user for number of random variable to use for stochastic expansion as well as the physical dimension for the problem.     
                             One could edit out the parameters inside this script file as in steps 2, 3,4.
      
    1. It promt for one command line input for nDim (number of RVs): enter the desired integer (Eg: 2 for 2RVs)
     
     generates KLE/PCE data for selected nDim (RVs) case. Data generated inside ssfem_fenics/data/klePceData/
       
    2. Control mesh density using lc parameter: Eg: lc=0.1/lc=0.09  , here lc=0.1 changes to lc=0.09 (NOTE: keep first "lc" always fixed to "0.1")
    3. Control type of coarse grid using "wb" parameter: Eg: wb=1/wb=0 , changes coarse grid from wire-basket(wb=1) to vertex grid(wb=0)
    4. Control number of partitions using "NP" parameter: Eg: NP=32/NP=64 , changes numper of partitions from 32 to 64 (default=32)

    generates GMSH/DDM data for selected LC and NP inside ssfem_fenics/data/meshData/
        
* $ sh processFEniCS.sh/ processFEniCS_MPI.sh    ::  Generates FEniCS mesh data and DDM assembly Mat-Vec data
    1. Control request time for job using "time" parameter: Eg: time=0-00:10/time=0-00:15 , changes job time from 10mins to 15mins (default=10mins)
    2. Be sure to change the output path to your desired directory
    
## Note : By default mail comes only when the process ends or fails. Changes to mail id and default value have to be changed inside other script files inside : ssfem_fenics/external/dolfin/src/jobscripts.
## Depending upon the problem you have selected this script file changes. The name of this file can be found from commands inside processFEniCS_MPI.sh. ex : run_fenics_poisson2D_stochasticDDM_parallel_twolevel.sh, 
## run_fenics_elasticity3D_stochasticDDM_parallel_twolevel.sh

## In many cases the process fails because of an issue with the way MPI Job is initialized. Please ensure your job is successfull by looking at the output folder you selected while running the job
## folder -default : /scratch/sudhipv/sudhipv/ssfem_fenics/data/slurm/.
## In general a successful job prints out details of random variable, applying BCs, and finally success message.
## Running FEniCS for 3D stochastic Poisson/twoLevel..
## number of partitions: 320
## number of dimensions: 3
## number of pceInputss: 10
## number of pceOutputs: 20
## ==========================Success============================



* $ sh processMPI.sh          :: compiles FORTRAN executables and run MPI/PETSc simulation
    1. Control request time for job using "time" parameter: Eg: time=0-00:10/time=0-00:15 , changes job time from 10mins to 15mins (default=10mins)
    2. Control request tasks per node using "tasks-per-node" parameter: Eg: tasks-per-node=32/tasks-per-node=16 , changes from 32 to 16 (default=32)
    3. Control request number of nodes using "nodes" parameter: Eg: nodes=1/nodes=2 , changes requested node from 1 to 2 (default=1)
    4. Contorl number of MPI processor to compile and run using "np" parameter: Eg: -np 32/-np 16 , changes from 32 to 16 (Note: np=nodes*tasks-per-node)
    5. Option to control job-name and output-file name eg: job-name=test/job-name=d2_lc05_p32 ,  (here d->dims, lc-> mesh density parameter, p->partitions)
    
### Output will be produced inside ./ssfem_fenics/data/vtkOutputs/
    
### Step-3 :: postprocess3D :: Copy data to HomePC and use Matlab/ParaView for postprocessing
* from top level project directory @HomePC
* $ cd data/vtkOutputs
* $ copy "posprocess3D.sh" to your local pc and edit the file locations    
* $ Then run 'sh posprocess3D.sh'  :: should do all of the steps outlined below (or follow step by step guide)

### Step-3 :: step-by-step
* $ scp -r sudhipv@graham.computecanada.ca:/~/scratch/sudhipv/ssfem_fenics/data/vtkOutputs/*.dat .
### For 3D mesh :: Mesh-connectivity for ParaView
* $ scp -r sudhipv@graham.computecanada.ca:/~/scratch/sudhipv/ssfem_fenics/data/meshData/points.dat .
* $ scp -r sudhipv@graham.computecanada.ca:/~/scratch/sudhipv/ssfem_fenics/data/meshData/tetrahedrons4c.dat .
* $ scp -r sudhipv@graham.computecanada.ca:/~/scratch/sudhipv/ssfem_fenics/data/klePceData/pcedata.dat .

#### For poisson3D
* $ matlabTerminal -r writeVtk
#### For elasticity3D
* $ matlabTerminal -r writeVecVtk
#### Use ParaView for visualization
* $ paraview out_*.vtk      :: open *.vtk files with paraview to check outputs


## NOTE:: Additional Useful Commands
* $ sacct  :: to check status for all of your recent jobs submitted using slrum
* $ more slrum-JobID.out  :: for outputs info of perticular JobID
* $ git fetch && git checkout graham
* $ git checkout -f   :: this will discard any local changes (useful to pull latest changes and discard old one)


### Folder structure : ###
* src    :: executable (all stochastic cases)
    - 3Delasticity_twolevel 
    - 3Delasticity_onelevel
    - 3Dpoisson_twolevel
    - 3Dpoisson_onelevel
    - 2Dpoisson_twolevel
    - 2Dpoisson_onelevel
    - onelevel (2D original-old)
    - twolevel (2D original-old)
* jobscripts   ::
    - preprocessKLE.sh     :: generates KLE/PCE data
    - preprocessMESH.sh :: generates GMSH/DDM data
    - preprocess.sh           ::  Generates KLE/PCE data and GMSH/DDM data
    
    - processFEniCS.sh    :: generates FEniCS/DDM assembly
    - processFEniCS_MPI.sh :: generates FEniCS/DDM assembly simultaneously in each core for each subdomain by parallel code
     
    - processMPI.sh          :: MPI/PETSc compile and execute
    
* external ::
    - puffin     :: external matlab FEM package
    - dolfin     :: external python/C++ FEM package
* data   ::
    - klePceData :: KLE and PCE data
    - mallocData :: malloc allocation data
    - meshData   :: preprocessed mesh data
    - vtkOutputs  :: outputs of simulation
* docs ::
    - contains project related documents

### Submit the job ###
* $ qsub run_ddm.sh     !! submit the job
* $ qstat -u userID     !! check status
* $ qdel jobID          !! kill the job

### Once status is complete then check ###
* $ more outputfile*
* $ more errorfile*

### Output visualization : ###
* $ paraview out_deterministic.vtk               

### Clean compilation/output data ###
* $ sh clean.sh             :: cleans everything

### Who do I talk to? ###
* Sudhi Sharma P V : sudhi.pv@cmail.carleton.ca


### Reference ##
* CMAME-Ajit: [Scalable Domain Decomposition Solvers for Stochastic PDEs in High Performance Computing](http://www.sciencedirect.com/science/article/pii/S0045782516313056)
