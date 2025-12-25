#!/bin/bash
#SBATCH --account=def-asarkar
#SBATCH --time=0-00:30
#SBATCH --mem-per-cpu=3700M
#SBATCH --nodes=2
#SBATCH --tasks-per-node=40
#SBATCH --job-name=3DP_3RV_30KNodes_80Parts_VB
#SBATCH --output=/scratch/sudhipv/ssfem_fenics/data/slurm/%j-%x.out

module load nixpkgs/16.09  intel/2018.3  openmpi/3.1.2  petsc/3.7.5

## cd /home/ajitd/petscTest/PETSc_helloworld_Guillimin
## cd $PBS_O_WORKDIR

## execute using PETSc-MPIEXEC :: named here as 'petscexec'
make all
mpiexec -np 80 ./a.out --mca pml ob1 -log_view

##  -log_summary                    ## PETSc Log summary
##  -ksp_monitor                    ## PETSc KSP iteration
##  -mat_view ::ascii_info          ## PETSc Mat mallocs
##  -mat_view draw -draw_pause 10   ## Mat sparsity pattern
##  -ksp_converged_reason	          ## print reason for converged or diverged
##  -ksp_monitor_solution	          ## plot solution at each iteration
##  -ksp_max_it                     ## maximum number of linear iterations
##  -ksp_rtol rtol	              ## default relative tolerance used for convergence

exit
