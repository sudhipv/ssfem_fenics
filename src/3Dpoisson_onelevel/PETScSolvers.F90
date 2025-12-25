!! This MODULE is a part of HPUQ package : Ajit Desai : Dec 2015 to Current
!! This is the solver MODULE contains Stochastic-DDM solvers using PETSc
!! This is mainly written to make use of PETSc for solving FEM or DDM problems
!! If you want to include any other PETSc Solvers, then this is the place
!!
!input:
!    DDM-Data  : Descretized blocks of matrices & vectors
!    Mesh-Data : Global & local mesh data
!    MPI-Data  : For MPI routines
!    PCE-Data  : For stochastic formulation
!
!output:
!    Ui, Ub    : Local solution vectors (interior, interface)
!    Ub_g      : Global solution vectors (interface)
!
!The following include statements are required for KSP Fortran programs:
!    petscsys.h: base PETSc routines
!    petscvec.h: vectors
!    petscmat.h: matrices
!    petscksp.h: Krylov subspace methods
!    petscpc.h : preconditioners
!
!! Contains: Stochastic Solvers
!! 1. PETSc_stolpcgm  : One-Level Lumped-PCGM Solver using PETSc
!!
!!-----------------------------------------------------------------------------------

MODULE PETScSolvers

use PETScommon
use PETScAssembly
use common

IMPLICIT NONE


CONTAINS


!!!---------------------------------------------------------------------------------------------
!!!Lumped-Preconditioned PCGM Solver : Algorithm:2 from the report
!!!---------------------------------------------------------------------------------------------

SUBROUTINE PETSc_stolpcgm(pid,np,nb,nbg,npcein,npceout,ncijk,ijk,cijk,Ui,Ub,Ub_g, &
                          mallocsCals,maxiter,tol)

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscpc.h>

!!!---------------------------------------------------------------------------------------------
!! PETSc-Data Declaration
    Mat              Asii,Asgg,Asgi                  !! PETSc matrix objects
    Mat              RsMat

    Vec              PetVecb,SolVecb                 !! PETSc vectors objects
    Vec              PetVeci,SolVeci
    Vec              PetVecbg,Fsi,Fsg
    Vec              PetVecbx

    KSP              kspAi, kspAg                    !! PETSc Solver Context
    PC               pcAi, pcAg                      !! PETSc Preconditioner

!!!---------------------------------------------------------------------------------------------
!! Fortran-Data Declaration
    integer :: i, pid, npceout, ncijk, maxiter, ierr, npcein
    integer :: np, nb, nbg, nip, nbp, nbgp, mallocsCals

    integer, dimension(nb*npceout)     :: tempb
    integer, dimension((np-nb)*npceout):: tempi
    integer, dimension(nbg*npceout)    :: tempbg
    integer, dimension(ncijk,3)        :: ijk

    double precision :: rho_next, rho_curr, alpha, beta, err, tol, NegOne

    double precision, dimension(nbg*npceout)     :: Ub_g, rb_g, Qb_g
    double precision, dimension(nbg*npceout)     :: Pb_g, Zb_g, RZb

    double precision, dimension(ncijk)           :: cijk
    double precision, dimension(nb*npceout)      :: Ub
    double precision, dimension((np-nb)*npceout) :: Ui

!!!-------------------------------------------------------------------------------------------
!! PETSc-Initialize
    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

    nip = (np-nb)*npceout      !! sto-interior nodes
    nbp = nb*npceout           !! sto-local-boundary nodes
    nbgp= nbg*npceout          !! sto-global-boundary nodes

!! Initiate Array's
    NegOne   = -1.0d0
    rho_next = 0.0d0
    rho_curr = 0.0d0

!!!-------------------------------------------------------------------------------------------
!! Assemble only Subdomain-Level Deterministic matrices to calculate Mallocs

    IF (mallocsCals .eq. 1) then

        if (pid .eq. 0) print*, '--------------------------------------------------'
        if (pid .eq. 0) print*, 'Initializing One-Level-PCGM Mallocs Calculation...'
        if (pid .eq. 0) print*, '--------------------------------------------------'

        call GetMallocs(pid,np,nb,npceout,nip,nbp,ierr)

        !Ui(:)  = 0.0d0
        !Ub(:)  = 0.0d0
        !Ub_g(:)= 0.0d0
    END IF

!! ELSE

!!!-------------------------------------------------------------------------------------------
    if (pid .eq. 0) print*, '------------------------------------------------'
    if (pid .eq. 0) print*, 'Initializing Stochastic-Lumped One-Level-PCGM...'
    if (pid .eq. 0) print*, '------------------------------------------------'


!!!-------------------------------------------------------------------------------------------
!! Assemble Subdomain-Level Matrices directly as PETSc Matrices
!! These routine has both (in-house Fortran or FEniCS-Dolfin) capabilities inherited
    !call StoMatSeqOneLevel(pid,Asii,Asgg,Asgi,p,e,t,np,ne,nt,nb,ndim,npcein,npceout,nomga,nip,nbp, &
    !                       casep,ncijk,ijk,cijk,dbounds,const_diff,omegas,multipliers,sigma, &
    !                       mIndex,sIndex,ierr)
    call StoMatSeqOneLevel(pid,Asii,Asgg,Asgi,np,nb,npcein,npceout,nip,nbp,ncijk,ijk,cijk,ierr)


    !call StoVecSeqOneLevel(pid,Fsi,Fsg,p,e,t,np,ne,nt,nb,ndim,npceout,nomga,nip,nbp,amp,dbounds, &
    !                        mIndex,sIndex,ierr)
    call StoVecSeqOneLevel(pid, Fsi, Fsg, np, nb, npceout, nip, nbp, ierr)

    call GetVecSeqTemp2(nip,PetVeci,Solveci,tempi,ierr)
    call GetVecSeqTemp2(nbp,PetVecb,Solvecb,tempb,ierr)
    call GetVecSeqDummy(nbgp,PetVecbg,tempbg,ierr)
    call GetVecSeqTemp1(nbp,PetVecbx,ierr)
    !!call GetVecSeqTemp1(nbp,PQb,ierr)
    !!call GetVecSeqTemp1(nbp,PPb,ierr)
    !!call GetVecSeqTemp2(nbp,PetVecbx,Solvecbx,tempbx,ierr)


!!!------------------------------------------------------------------------------------------
    if (pid .eq. 0) print*, '----------------------------------------------------------------'
    if (pid .eq. 0) print*, 'Using Sparse-Iterative-PETSc-KSP Solver For Interior Problems...'
    if (pid .eq. 0) print*, '----------------------------------------------------------------'

!!---------------------------------------------------------------------------------------------
!!  Uncomment this section for Execution Time Calculations
    !integer :: c1,c2,cr
    !double precision :: time1
    !double precision :: time2
    !call system_clock(count_rate=cr)

!!---------------------------------------------------------------------------------------------
! Pre-Iteration : 0th Iteration : Step 1 to 9
!!!---------------------------------------------------------------------------------------------
    !! Step 3 : Solve
    call PETScKSP(kspAi,pcAi,Asii,Fsi,SolVeci,ierr)       !! PETSc-KSP-Solver
    !!call PETScMUMPS(kspAi,pcAi,Asii,Fsi,SolVeci,ierr)       !! PETSc-Mumps-Solver

    !! Step 4 : Compute
    call MatMult(Asgi,SolVeci,SolVecb,ierr)               !! A(1) x(2) = b(3)
    call VecAYPX(SolVecb,NegOne,Fsg,ierr)                 !! (y,a,x):y = x + a*y

    !! Step 5 : Gather
    !call getubg(pid,nb,nbg,npceout,RZb,Zb)
    call GetRs(pid,nb,nbg,npceout,RsMat)                  !! Sparse Rs Operator
    call MatMultTranspose(RsMat,SolVecb,PetVecbg,ierr)     !! A(1)'x(2) = b(3)
    call VecGetValues(PetVecbg, nbgp, tempbg, RZb, ierr)
    call MPI_ALLREDUCE(RZb,rb_g,nbgp,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

!!!---------------------------------------------------------------------------------------------
    !! Parallel-Preconditioning effect on initial residual : Step 2 to 9
    !! Step 6 : Scatter
    !! call getub(pid,nb,nbg,npceout,rb_g,Zb)
    call SetVecSeq(nbgp,rb_g,PetVecbg,ierr)                !! Fortran to PETSc Vec
    call MatMult(RsMat,PetVecbg,PetVecb,ierr)              !! Sparse PETSc Mat*Vec

    !! Step 7 : Solve
    !! call SetVecSeq(nbp,Zb,PetVecb,ierr)
    call PETScKSP(kspAg,pcAg,Asgg,PetVecb,SolVecb,ierr)
    !call PETScMUMPS(kspAg,pcAg,Asgg,PetVecb,SolVecb,ierr)

    !! Step 8 : Gather
    !! call getubg(pid,nb,nbg,npceout,RZb,Zb)
    !! call VecGetValues(SolVecb, nbp, tempb, Zb, ierr)
    call MatMultTranspose(RsMat,SolVecb,PetVecbg,ierr)     !! A(1)'x(2) = b(3)
    call VecGetValues(PetVecbg, nbgp, tempbg, RZb, ierr)
    CALL MPI_REDUCE(RZb,Zb_g,nbgp,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

    !! Step 10 & 11 : Compute
    if (pid .eq. 0) Pb_g = Zb_g
    if (pid .eq. 0) rho_next = dot_product(rb_g,Zb_g)

!!!---------------------------------------------------------------------------------------------
!! PCGM Iteration : Main For Loop for each iteration: Step 12 to 32
!!!---------------------------------------------------------------------------------------------
    Ub_g(:)  = 0.0d0
    DO i = 1,maxiter

        !! Parallel Matrix-Vector product: Step 13 to 18
        !! Step 14 : Scatter
        CALL MPI_BCAST(Pb_g,nbgp,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        !call getub(pid,nb,nbg,npceout,Pb_g,Zb)
        !call SetVecSeq(nbp,Zb,PetVecb,ierr)
        call SetVecSeq(nbgp,Pb_g,PetVecbg,ierr)
        call MatMult(RsMat,PetVecbg,PetVecb,ierr)

        !! Step 15 : Solve
        call MatMultTranspose(Asgi,PetVecb,PetVeci,ierr)
        call KSPSolve(kspAi,PetVeci,SolVeci,ierr)

        !! Step 16 : Compute
        call MatMult(Asgi,SolVeci,SolVecb,ierr)
        call MatMult(Asgg,PetVecb,PetVecbx,ierr)
        call VecAXPY(PetVecbx,NegOne,SolVecb,ierr)

        !! Step 17 : Gather
        !call VecGetValues(PQb, nbp, tempb, Zb, ierr)
        !call getubg(pid,nb,nbg,npceout,RZb,Zb)
        call MatMultTranspose(RsMat,PetVecbx,PetVecbg,ierr)
        call VecGetValues(PetVecbg, nbgp, tempbg, RZb, ierr)

        CALL MPI_REDUCE(RZb,Qb_g,nbgp,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

        !!!--------------------------------------------------------------------------------------
        !! Step 19,20, 21 & 22 : Compute & Update
        if (pid .eq. 0) then
            rho_curr = rho_next
            alpha = rho_curr/dot_product(Qb_g,Pb_g)
            err = alpha*alpha*dot_product(Pb_g,Pb_g)/dot_product(Ub_g,Ub_g)
            Ub_g = Ub_g + alpha*Pb_g
            rb_g = rb_g - alpha*Qb_g
            print*, '----------------'
            print*, 'Main Iteration #',i, ', relative error of ',err
            print*, '----------------'
            !!print*, 'and sum of Ub of ', sum(Ub_g)
        end if

        !! Step 23 : Exit
        CALL MPI_BCAST(err,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        if (err .lt. tol) exit

        !!!----------------------------------------------------------------------------------------
        !! Parallel-Preconditioning effect for each iteration: Step 24 to 28
        !! Step 25 : Scatter
        CALL MPI_BCAST(rb_g,nbgp,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call SetVecSeq(nbgp,rb_g,PetVecbg,ierr)
        call MatMult(RsMat,PetVecbg,PetVecb,ierr)

        !! Step 26 : Solve
        !call getub(pid,nb,nbg,npceout,rb_g,Zb)
        !call SetVecSeq(nbp,Zb,PetVecb,ierr)
        call KSPSolve(kspAg,PetVecb,SolVecb,ierr)

        !! Step 27 : Gather
        !call VecGetValues(SolVecb, nbp, tempb, Zb, ierr)
        !call getubg(pid,nb,nbg,npceout,RZb,Zb)
        call MatMultTranspose(RsMat,SolVecb,PetVecbg,ierr)
        call VecGetValues(PetVecbg, nbgp, tempbg, RZb, ierr)
        CALL MPI_REDUCE(RZb,Zb_g,nbgp,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

        !! Step 29, 30 & 31 : Compute & Update
        if (pid .eq. 0) then
        rho_next = dot_product(rb_g,Zb_g)
        beta = rho_next/rho_curr
        Pb_g = Zb_g + beta*Pb_g
        end if

    END DO
    !! Step 33 : Output : Global interface solution vector

!!!---------------------------------------------------------------------------------------------
!! Post Iteration : To calculate local interior/interface solutions : Step 34 to 38
!!!---------------------------------------------------------------------------------------------
    !! Step 35 : Scatter
    CALL MPI_BCAST(Ub_g,nbgp,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    !call getub(pid,nb,nbg,npceout,Ub_g,Ub)
    call SetVecSeq(nbgp,Ub_g,PetVecbg,ierr)
    call MatMult(RsMat,PetVecbg,PetVecb,ierr)
    call VecGetValues(PetVecb, nbp, tempb, Ub, ierr)

    !! Step 36 : Solve
    !call SetVecSeq(nbp,Ub,PetVecb,ierr)
    call MatMultTranspose(Asgi,PetVecb,SolVeci,ierr)
    call VecAXPY(Fsi,NegOne,SolVeci,ierr)
    call KSPSolve(kspAi,Fsi,SolVeci,ierr)
    call VecGetValues(SolVeci, nip, tempi, Ui, ierr)
    !! Step 38 : Output : Local interface & interior solution vectors

!!!---------------------------------------------------------------------------------------------
!!! PETSc-Destroy & Finalize
    call VecDestroy(PetVeci,ierr)
    call VecDestroy(SolVeci,ierr)
    call VecDestroy(PetVecb,ierr)
    call VecDestroy(SolVecb,ierr)
    call VecDestroy(PetVecbg,ierr)
    call VecDestroy(Fsi,ierr)
    call VecDestroy(Fsg,ierr)
    call VecDestroy(PetVecbx,ierr)
    Call MatDestroy(RsMat,ierr)
    Call MatDestroy(Asii,ierr)
    Call MatDestroy(Asgg,ierr)
    Call MatDestroy(Asgi,ierr)
    Call KspDestroy(kspAi,ierr)
    Call KspDestroy(kspAg,ierr)


!! END IF
!!-----------------------------------------------------------------------------------------------
!!! PETSc-Finalize
    call PetscFinalize(ierr)


END SUBROUTINE PETSc_stolpcgm


END MODULE PETScSolvers
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


