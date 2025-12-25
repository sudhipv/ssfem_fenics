#!/bin/bash
#PBS -l pmem=6400m
#PBS -l walltime=00:30:00
#PBS -l nodes=1:ppn=8
#PBS -o outputfile
#PBS -e errorfile
## PBS -l procs=16
## PBS -N jobname

module load gcc/4.9.1 openmpi/1.8.3-gcc cmake/2.8.5

## cd /home/ajitd/petscTest/PETSc_helloworld_Guillimin

cd $PBS_O_WORKDIR

## execute using PETSc-MPIEXEC :: named here as 'petscexec'
make all
mpiexec -np 8 ./a.out

##  -log_summary                    ## PETSc Log summary
##  -ksp_monitor                    ## PETSc KSP iteration
##  -mat_view ::ascii_info          ## PETSc Mat mallocs
##  -mat_view draw -draw_pause 10   ## Mat sparsity pattern
##  -ksp_converged_reason	        ## print reason for converged or diverged
##  -ksp_monitor_solution	        ## plot solution at each iteration
##  -ksp_max_it                     ## maximum number of linear iterations
##  -ksp_rtol rtol	                ## default relative tolerance used for convergence

exit