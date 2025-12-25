#!/bin/bash
#SBATCH --time=0-00:10
#SBATCH --mem-per-cpu=7700M
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1

source ~/testDDM/sgfem_ddm_hpc/external/dolfin/src/fenics_activate.sh
python pointElements_xml.py

## cd /home/ajitd/petscTest/PETSc_helloworld_Guillimin
module load hdf5-mpi/1.8.18 boost eigen python/3 python35-scipy-stack/2017a petsc/3.7.5 fftw-mpi/3.3.6
source ~/software/dolfin/share/dolfin/dolfin.conf
source ~/fenics/bin/activate

## execute using PETSc-MPIEXEC :: named here as 'petscexec'
python stoDDM_poisson.py

##  -log_summary                    ## PETSc Log summary
##  -ksp_monitor                    ## PETSc KSP iteration
##  -mat_view ::ascii_info          ## PETSc Mat mallocs
##  -mat_view draw -draw_pause 10   ## Mat sparsity pattern
##  -ksp_converged_reason	        ## print reason for converged or diverged
##  -ksp_monitor_solution	        ## plot solution at each iteration
##  -ksp_max_it                     ## maximum number of linear iterations
##  -ksp_rtol rtol	                ## default relative tolerance used for convergence

exit
