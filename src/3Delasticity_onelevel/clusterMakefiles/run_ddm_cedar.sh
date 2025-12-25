#!/bin/bash
#SBATCH --account=def-asarkar
#SBATCH --time=0-00:10
#SBATCH --mem-per-cpu=7700M
#SBATCH --nodes=1
#SBATCH --tasks-per-node=32
#SBATCH --job-name=DP_nRV_nNodes_nParts
#SBATCH --output=%x-%j.out

module load nixpkgs/16.09  intel/2017.1  openmpi/2.0.2  petsc/3.7.5
## cd /home/ajitd/petscTest/PETSc_helloworld_Guillimin

## execute using PETSc-MPIEXEC :: named here as 'petscexec'
make all
mpiexec -np 64 ./a.out

##  -log_summary                    ## PETSc Log summary
##  -ksp_monitor                    ## PETSc KSP iteration
##  -mat_view ::ascii_info          ## PETSc Mat mallocs
##  -mat_view draw -draw_pause 10   ## Mat sparsity pattern
##  -ksp_converged_reason	        ## print reason for converged or diverged
##  -ksp_monitor_solution	        ## plot solution at each iteration
##  -ksp_max_it                     ## maximum number of linear iterations
##  -ksp_rtol rtol	                ## default relative tolerance used for convergence

exit
