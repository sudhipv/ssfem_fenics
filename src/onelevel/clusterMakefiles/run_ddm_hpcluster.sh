#!/usr/bin/env bash
source  /scratch/software/use_openmpi.sh

## select the number of partitions
NP=176

## execute using PETSc-MPIEXEC :: named here as 'petscexec'
mpirun -np $NP -machinefile  hostfile-ib2  ./a.out

##  -log_summary                    ## PETSc Log summary
##  -ksp_monitor                    ## PETSc KSP iteration
##  -mat_view ::ascii_info          ## PETSc Mat mallocs
##  -mat_view draw -draw_pause 10   ## Mat sparsity pattern
##  -ksp_converged_reason	        ## print reason for converged or diverged
##  -ksp_monitor_solution	        ## plot solution at each iteration
##  -ksp_max_it                     ## maximum number of linear iterations
##  -ksp_rtol rtol	                ## default relative tolerance used for convergence
