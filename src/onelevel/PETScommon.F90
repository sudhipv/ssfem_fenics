!! This MODULE is a part of HPSFEM package : Ajit Desai : Dec 2015
!! This is the common MODULE contains various subrountines for PETSc
!! This is mainly written to make use of PETSc for solving FEM or DDM problems
!! If you want to include any other PETSc subroutines, then this is the place
!! Contains: 
!! 1. PETScKSP     : Construct PETSc KSP Set-up and solve using iterative solver
!! 2. PETScMUMPS   : Construct PETSc MUMPS set-up and slover using direct solver 
!! 3. GetPetVec    : Directly assemble Sparse-PETSc Vector : No dense vector
!! 4. GetPetMat    : Directly assemble Sparse-PETSc Matrix : No dense Matrix
!! 5. GetPetVecSelf:
!! 6. GetPetMatSelf:
!! 5. GetVec       : Construct PETSc Vector : General Way : Serial/Parallel
!! 6. GetMat       : Construct PETSc Matrix : General Way : Serial/Parallel
!! 7. GetVecSeq    : Construct PETSc Vector : Sequential  : Serial
!! 8. GetMatSeq    : Construct PETSc Matrix : Sequential  : Serial
!! 9. GetVecMPI    : Construct PETSc Vector : MPI format  : Parallel
!!10. PETScSolver  : Call PETSc assembled Mat, Vec and solver Ax=b using ksp
!!-----------------------------------------------------------------------------------

MODULE PETScommon

 use PETScAssembly
 use assembly
 use common

IMPLICIT NONE


CONTAINS


!!!*********************************************************
!!! Subroutine: To Call PETSc KSP solver in most general way
SUBROUTINE PETScKSP(ksp,pc,PetMat,PetVec,Solvec,ierr)

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscpc.h>

    Vec              PetVec,SolVec
    Mat              PetMat
    KSP              ksp
    PC               pc

    PetscErrorCode   ierr

    call KSPCreate(PETSC_COMM_SELF,ksp,ierr)
    call KSPSetOperators(ksp,PetMat,PetMat,ierr)
    call KSPSetFromOptions(ksp,ierr)   !! position changed : moved down

    call KSPSetType(ksp,KSPCG,ierr)
    call KSPGetPC(ksp,pc,ierr)
    call PCSetType(pc,PCBJACOBI,ierr)
    !call PCSetType(pc,PCKSP,ierr)

    call KSPSolve(ksp,PetVec,SolVec,ierr)

END SUBROUTINE PETScKSP


!!!!*********************************************************
!!!! Subroutine
!SUBROUTINE PETScMUMPS(ksp,pc,PetMat,PetVec,Solvec,ierr)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!#include <petsc/finclude/petscmat.h>
!#include <petsc/finclude/petscksp.h>
!#include <petsc/finclude/petscpc.h>
!
!    Vec              PetVec,SolVec
!    Mat              PetMat
!    KSP              ksp
!    PC               pc
!
!    Mat              F
!    PetscInt         ival,icntl
!
!    PetscErrorCode   ierr
!
!    call KSPCreate(PETSC_COMM_SELF,ksp,ierr)
!    call KSPSetOperators(ksp,PetMat,PetMat,ierr)
!
!    call KSPSetType(ksp,KSPPREONLY,ierr)
!    call KSPGetPC(ksp,pc,ierr)
!    !call PCSetType(pc,PCLU,ierr)           !! LU Factorization
!    call PCSetType(pc,PCCHOLESKY,ierr)    !! Cholesky Factorization
!    call PCFactorSetMatSolverPackage(pc,MATSOLVERMUMPS,ierr)
!    call PCFactorSetUpMatSolverPackage(pc,ierr)
!    call PCFactorGetMatrix(pc,F,ierr)
!
!    !! sequential ordering
!    icntl = 7
!    ival  = 2
!    call MatMumpsSetIcntl(F,icntl,ival,ierr)
!
!    call KSPSetFromOptions(ksp,ierr)
!    call KSPGetPC(ksp,pc,ierr)
!
!    call KSPSolve(ksp,PetVec,SolVec,ierr)
!
!
!END SUBROUTINE PETScMUMPS

!!!*********************************************************
!!! Subroutine: To construc PETSc Vector in Sequential Form
SUBROUTINE GetVecSeq(n,fvec,PetVec,Solvec,temp1,ierr)

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    Vec              PetVec,SolVec
    PetscInt         i,ii
    PetscInt         n,one
    PetscScalar      ivec
    PetscErrorCode   ierr

    double precision, dimension(n) :: fvec
    integer, dimension(n)          :: temp1

    call VecCreateSeq(PETSC_COMM_SELF, n, PetVec, ierr)

    call VecDuplicate(PetVec,SolVec,ierr)

    one = 1
    do i = 1,n
        ii = i-1
        temp1(i) = i-1
        ivec = fvec(i)
        call VecSetValues(PetVec,one,ii,ivec,INSERT_VALUES,ierr)
    end do

    call VecAssemblyBegin(PetVec,ierr)
    call VecAssemblyEnd(PetVec,ierr)


END SUBROUTINE GetVecSeq


!!!*********************************************************
!!! Subroutine: To construc PETSc Vector in Sequential Form
SUBROUTINE SetVecSeq(n,fvec,PetVec,ierr)

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    Vec              PetVec
    PetscInt         i,ii
    PetscInt         n,one
    PetscScalar      ivec
    PetscErrorCode   ierr

    double precision, dimension(n) :: fvec
    !integer, dimension(n)          :: temp1
    !call VecCreateSeq(PETSC_COMM_SELF, n, PetVec, ierr)
    !call VecDuplicate(PetVec,SolVec,ierr)

    one = 1
    do i = 1,n
        ii = i-1
        !temp1(i) = i-1
        ivec = fvec(i)
        call VecSetValues(PetVec,one,ii,ivec,INSERT_VALUES,ierr)
    end do

    call VecAssemblyBegin(PetVec,ierr)
    call VecAssemblyEnd(PetVec,ierr)

END SUBROUTINE SetVecSeq


!!!*********************************************************
!!! Subroutine: To construc PETSc Matrix in Sequential Form
SUBROUTINE GetMatSeq(n,m,Amat,PetMat,ierr)

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>


    Mat              PetMat
    PetscInt         i,j,ii,jj
    PetscInt         n,m,nz,one
    PetscScalar      imat
    PetscErrorCode   ierr

    double precision, dimension(n,m) :: Amat

    nz = 12  !! number of non-zeros per row for matrix pre-allocation
    call MatCreateSeqAIJ(PETSC_COMM_SELF,n,m,nz,PETSC_NULL_INTEGER,PetMat,ierr)

    one = 1
    do i = 1,n
    do j = 1,m
        if (Amat(i,j) .ne. 0) then
            ii = i-1
            jj = j-1
            imat = Amat(i,j)
            call MatSetValues(PetMat,one,ii,one,jj,imat,INSERT_VALUES,ierr)
        end if
    end do
    end do

    call MatAssemblyBegin(PetMat,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(PetMat,MAT_FINAL_ASSEMBLY,ierr)


END SUBROUTINE GetMatSeq

!!!*********************************************************
!!! Subroutine: To construc PETSc Vector in Sequential Form
SUBROUTINE GetVecSeqDummy(n,PetVec,temp1,ierr)

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    Vec              PetVec
    PetscInt         i, n
    PetscErrorCode   ierr

    integer, dimension(n) :: temp1

    call VecCreateSeq(PETSC_COMM_SELF, n, PetVec, ierr)

    do i = 1,n
       temp1(i) = i-1
    end do


END SUBROUTINE GetVecSeqDummy


!!!*********************************************************
!!! Subroutine: To construc PETSc Vector in Sequential Form
SUBROUTINE GetVecSeqTemp1(n,PetVec,ierr)

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    Vec              PetVec
    PetscInt         n
    PetscErrorCode   ierr

    call VecCreateSeq(PETSC_COMM_SELF, n, PetVec, ierr)


END SUBROUTINE GetVecSeqTemp1


!!!*********************************************************
!!! Subroutine: To construc PETSc Vector in Sequential Form
SUBROUTINE GetVecSeqTemp2(n,PetVec,Solvec,temp1,ierr)

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    Vec              PetVec,SolVec
    PetscInt         i, n
    PetscErrorCode   ierr

    integer, dimension(n) :: temp1

    call VecCreateSeq(PETSC_COMM_SELF, n, PetVec, ierr)
    call VecDuplicate(PetVec,SolVec,ierr)

    do i = 1,n
        temp1(i) = i-1
    end do

END SUBROUTINE GetVecSeqTemp2


!!!*********************************************************
!! PETSc-Vec assembly using Fortrao assembly (PETScAssembly.F90) # OLD
SUBROUTINE StoVecSeqOneLevel_Fortron(Fsi, Fsg, p, e, t, np,ne,nt,nb, ndim, npceout,nomga,&
                             nip,nbp,amp,dbounds,mIndex,sIndex,ierr)

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    Vec              Fsi, Fsg
    PetscInt         nip, nbp
    PetscErrorCode   ierr

    integer :: np, ne, nt, nb, ndim, npceout, nomga
    double precision                  :: amp
    integer, dimension(3,ne)          :: e
    integer, dimension(3,nt)          :: t
    double precision, dimension(2,np) :: p
    double precision, dimension(2,2)  :: dbounds

    integer, dimension(ndim,npceout)  :: mIndex
    integer, dimension(2,ndim)        :: sIndex

    !!integer, dimension(nip) :: tempi
    !!integer, dimension(nbp) :: tempb
    !!nip = ((np-nb)*npceout)
    !!nbp = nb*npceout

    call VecCreateSeq(PETSC_COMM_SELF, nip, Fsi, ierr)
    call VecCreateSeq(PETSC_COMM_SELF, nbp, Fsg, ierr)
    !!call VecDuplicate(Fsi,SolVeci2,ierr)
    !!call VecDuplicate(Fsg,SolVecb2,ierr)

    !! Fortran Vec-Assembly procedur
    call PETscSubVecAssembly(Fsi, Fsg, p, e, t, np, ne, nt, nb, ndim, npceout,nomga,amp,dbounds,mIndex,sIndex)

    call VecAssemblyBegin(Fsi,ierr)
    call VecAssemblyBegin(Fsg,ierr)
    call VecAssemblyEnd(Fsi,ierr)
    call VecAssemblyEnd(Fsg,ierr)


END SUBROUTINE StoVecSeqOneLevel_Fortron


!!!*********************************************************
!! PETSc-Vec assembly using FEniCS assembled ddm Vecs # NEW
SUBROUTINE StoVecSeqOneLevel(pid, Fsi, Fsg, p, e, t, np,ne,nt,nb, ndim, npceout,nomga,&
                             nip,nbp,amp,dbounds,mIndex,sIndex,ierr)

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    Vec              Fsi, Fsg
    PetscInt         nip, nbp
    PetscScalar      ivec
    PetscErrorCode   ierr

    integer :: np, ne, nt, nb, ndim, npceout, nomga
    double precision                  :: amp
    integer, dimension(3,ne)          :: e
    integer, dimension(3,nt)          :: t
    double precision, dimension(2,np) :: p
    double precision, dimension(2,2)  :: dbounds

    integer, dimension(ndim,npceout)  :: mIndex
    integer, dimension(2,ndim)        :: sIndex

    double precision, allocatable, dimension(:) :: bi,bg
    character(len=255)      :: str2, str3
    integer :: pid, nid, id, ii, one, selectPackage

!!-----------------------------------------------------------------------------------------------
    !! Allocate memory to read pre-assembled vectors
    allocate(bi(np-nb), bg(nb))

    call VecCreateSeq(PETSC_COMM_SELF, nip, Fsi, ierr)
    call VecCreateSeq(PETSC_COMM_SELF, nbp, Fsg, ierr)

!!-----------------------------------------------------------------------------------------------
    selectPackage = 1  !! 1-Dolfin, 2-Puffin, 3-Fortran Assembly
!!-----------------------------------------------------------------------------------------------

    if (selectPackage == 3) then
        !! in-house Fortran Vec-Assembly procedure
        call PETscSubVecAssembly(Fsi, Fsg, p, e, t, np, ne, nt, nb, ndim,                npceout,nomga,amp,dbounds,mIndex,sIndex)

    else
        !! FEniCS based Vec-Assembly procedure (using preassembled vecs)
        nid= np-nb
        one= 1

        bi = 0.0d0
        bg = 0.0d0

        call int2str(str2,pid+1,1)
        str3 = '../../external/dolfin/data/Amats/subdom' // trim(str2) // '/bi1.dat'
        open(unit=1,file=str3,status='old')
        read(unit=1,fmt=*) bi
        close(1)

        str3 = '../../external/dolfin/data/Amats/subdom' // trim(str2) // '/bg1.dat'
        open(unit=2,file=str3,status='old')
        read(unit=2,fmt=*) bg
        close(2)

        !---------------------------------------->fsi
        do id = 1,nid
            if (bi(id) .ne. 0) then
                ii = (id-1)
                ivec = bi(id)
                call VecSetValues(Fsi,one,ii,ivec,ADD_VALUES,ierr)
            end if
        end do

        !!---------------------------------------->fsg
        do id = 1,nb
            if (bg(id) .ne. 0) then
                ii = (id-1)
                ivec = bg(id)
                call VecSetValues(Fsg,one,ii,ivec,ADD_VALUES,ierr)
            end if
        end do

    end if

    call VecAssemblyBegin(Fsi,ierr)
    call VecAssemblyBegin(Fsg,ierr)
    call VecAssemblyEnd(Fsi,ierr)
    call VecAssemblyEnd(Fsg,ierr)

    !!-----------------------------------------------------------------------------------------------
    if (pid==0) then
        if (selectPackage == 3) then
            print*, 'PETSc-Vec assembly using: In-house Fortran assembly'
        else
            print*, 'PETSc-Vec assembly using: FEniCS-dolfin assemlby'
        end if
    end if

    DEALLOCATE(bi,bg)

END SUBROUTINE StoVecSeqOneLevel


!!!*********************************************************
!!! Subroutine: To construc PETSc Matrix in most general way
SUBROUTINE StoMatSeqOneLevel(pid,Asii,Asgg,Asgi,p,e,t,np,ne,nt,nb,ndim,npcein,npceout,nomga,nip,nbp,&
                 casep,ncijk,ijk,cijk,dbounds,const_diff,omegas,multipliers,sigma,mIndex,sIndex,ierr)

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>


    Mat              Asii,Asgg,Asgi
    PetscInt         nip, nbp, one !! iNZ, nzi, nzb,
    PetscInt         ii, jj, id, jd, nni, nnj
    PetscScalar      imat
    PetscErrorCode   ierr

    PetscInt,ALLOCATABLE    :: nnzi(:), nnzb(:), nnzbi(:)
    character(len=255)      :: str1

    integer :: i, j, k, indexi, nid, ncijk, casep, pid
    integer :: np, ne, nt, nb, ndim, nomga, npceout, npcein
    double precision :: const_diff, sigma

    integer, dimension(ncijk,3)       :: ijk
    integer, dimension(3,ne)          :: e
    integer, dimension(3,nt)          :: t
    double precision, dimension(2,np) :: p
    double precision, dimension(ncijk):: cijk
    double precision, dimension(2,2)  :: dbounds
    double precision, dimension(nomga):: omegas, multipliers

    integer, dimension(ndim,npceout)  :: mIndex
    integer, dimension(2,ndim)        :: sIndex

    double precision, allocatable, dimension(:,:) :: Adii,Adgg,Adgi

    character(len=255)      :: str2, str3

    integer :: selectPackage

!!-----------------------------------------------------------------------------------------------
    !! Allocate memory to read pre-assembled matrices
    allocate(Adii((np-nb),(np-nb)), Adgg(nb,nb), Adgi(nb,(np-nb)))
    allocate(nnzi(nip),nnzb(nbp),nnzbi(nbp))

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzi' // trim(str1) // '.dat'
    open(unit=1,file=str1,status='old')
    read(unit=1,fmt=*) nnzi
    close(1)

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzb' // trim(str1) // '.dat'
    open(unit=2,file=str1,status='old')
    read(unit=2,fmt=*) nnzb
    close(2)

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzbi' // trim(str1) // '.dat'
    open(unit=3,file=str1,status='old')
    read(unit=3,fmt=*) nnzbi
    close(3)

    call MatCreateSeqAIJ(PETSC_COMM_SELF, nip, nip, 0, nnzi, Asii, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, nbp, nip, 0, nnzbi, Asgi, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, nbp, nbp, 0, nnzb, Asgg, ierr)

!!-----------------------------------------------------------------------------------------------
!! Manual (approximate) memory PETSc memory allocation procedure
!! number of non-zeros per row for Aii&Agi/Agg  matrix mallocs
!! For nOrd = 1; nDim = 50; nMesh=15K; nzi = 600,  nzb = 300
!! For nOrd = 2; nDim = 15; nMesh=15K; nzi = 1500, nzb = 900
!! For nOrd = 2; nDim = 20; nMesh=15K; nzi = 2200, nzb = 1200
!! For nOrd = 2; nDim = 25; nMesh=15K; nzi = 3200, nzb = 2000
!
!    nzi = 300    !! ndim=50: 6000
!    nzb = 200    !! ndim=50; 3500
!
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, nip, nip, nzi, PETSC_NULL_INTEGER, Asii, ierr)
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, nbp, nip, nzb, PETSC_NULL_INTEGER, Asgi, ierr)
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, nbp, nbp, nzb, PETSC_NULL_INTEGER, Asgg, ierr)
!    !!call MatCreateSeqAIJ(PETSC_COMM_SELF, nip, nbp, nzb, PETSC_NULL_INTEGER, Asig, ierr)
!
!    call MatSetOption(Asii,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE,ierr)
!    call MatSetOption(Asgg,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE,ierr)
!    call MatSetOption(Asgi,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE,ierr)
!    call MatSetOption(Asii,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,ierr)
!    call MatSetOption(Asgg,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,ierr)
!    call MatSetOption(Asgi,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,ierr)
!
!!-----------------------------------------------------------------------------------------------
!! Method : 3
nid    = np-nb
one    = 1
indexi = 1

!!-----------------------------------------------------------------------------------------------
selectPackage = 1  !! 1-Dolfin, 2-Puffin, 3-Fortran  !! 1 and 3 validated
!!-----------------------------------------------------------------------------------------------

DO k = 1,npcein !! npceout

    !! Deterministic Matices
    Adii = 0.0d0
    Adgg = 0.0d0
    Adgi = 0.0d0

!!-----------------------------------------------------------------------------------------------
!! For this option we use Fortran assemble routines (in-house)
    if (selectPackage == 3) then
        if (k .eq. 1) then
            !! Det-Advection Matrix
            call SubAssembleMatrix(Adii,Adgg,Adgi,p,e,t,np,ne,nt,nb,ndim,npceout,nomga,1,0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],&
                    [0.0d0,0.0d0,0.0d0,0.0d0],0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],dbounds,mIndex,sIndex,0,0)
            !! Det-Diffusion Matrix
            call SubAssembleMatrix(Adii,Adgg,Adgi,p,e,t,np,ne,nt,nb,ndim,npceout,nomga,2,const_diff,omegas,multipliers,&
                    sigma,[0.0d0,0.0d0,0.0d0,0.0d0],dbounds,mIndex,sIndex,casep,1)
        else
            !! Det-Advection Matrix  !! Included on Sun, Aug, 27, 2017 !! on stochastic dof we don't apply Penulty factor
            ! call SubAssembleMatrix(Adii,Adgg,Adgi,p,e,t,np,ne,nt,nb,ndim,npceout,nomga,1,0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],&
            !         [0.0d0,0.0d0,0.0d0,0.0d0],0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],dbounds,mIndex,sIndex,0,0)
            !! Sto-Diffusion Matrices
            call SubAssembleMatrix(Adii,Adgg,Adgi,p,e,t,np,ne,nt,nb,ndim,npceout,nomga,2,const_diff,omegas,multipliers,&
                    sigma,[0.0d0,0.0d0,0.0d0,0.0d0],dbounds,mIndex,sIndex,casep,k)
        end if
    end if

!!!-----------------------------------------------------------------------------------------------
!! Deterministic(Aii):: ADii
    call int2str(str1,pid+1,1)
    call int2str(str2,k,1)

!!! Python-dolfin
!!! The pressembled FEniCS DD Matrices for each subdomain for each PCE mode
    if (selectPackage == 1) then
    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/ADii' // trim(str2) // '.dat'
    open(unit=1,file=str3,status='old')
    read(unit=1,fmt=*) Adii
    end if

!!! Matlab-puffin
!!! The pressembled Matlab DD Matrices for each subdomain for each PCE mode
    if (selectPackage == 2) then
    str3 = '../../external/puffin/data/Amats/subdom' // trim(str1) // '/ADii' // trim(str2) // '.dat'
    open(unit=1,file=str3,status='old')
    read(unit=1,fmt=*) Adii
    end if

!!! Fortran-assembly:: Fortran assemble routines (in-house)
!!! This part only stores assembled matrices to compare with FEniCS assembled matrices
    if (selectPackage == 3) then
    str3 = '../../external/puffin/data/Amats/subdom' // trim(str1) // '/OrigADii' // trim(str2) // '.dat'
    open(unit=1,file=str3,status='replace')
    write(unit=1,fmt=*) Adii
    end if

    close(1)

!!!-----------------------------------------------------------------------------------------------
!! Deterministic(Agg):: ADgg
    call int2str(str1,pid+1,1)
    call int2str(str2,k,1)

!! Python-dolfin
    if (selectPackage == 1) then
    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/ADgg' // trim(str2) // '.dat'
    open(unit=2,file=str3,status='old')
    read(unit=2,fmt=*) Adgg
    end if

!! Matlab-puffin
    if (selectPackage == 2) then
    str3 = '../../external/puffin/data/Amats/subdom' // trim(str1) // '/ADgg' // trim(str2) // '.dat'
    open(unit=2,file=str3,status='old')
    read(unit=2,fmt=*) Adgg
    end if

!!! Fortran-assembly
    if (selectPackage == 3) then
    str3 = '../../external/puffin/data/Amats/subdom' // trim(str1) // '/OrigADgg' // trim(str2) // '.dat'
    open(unit=2,file=str3,status='replace')
    write(unit=2,fmt=*) Adgg
    end if

    close(2)

!!!-----------------------------------------------------------------------------------------------
!! Deterministic(Agi):: ADgi        !! Note we don't need ADig (due to symmetry)
    call int2str(str1,pid+1,1)
    call int2str(str2,k,1)

!! Python-dolfin
    if (selectPackage == 1) then
    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/ADgi' // trim(str2) // '.dat'
    open(unit=3,file=str3,status='old')
    read(unit=3,fmt=*) Adgi
    end if

!! Matlab-puffin
    if (selectPackage == 2) then
    str3 = '../../external/puffin/data/Amats/subdom' // trim(str1) // '/ADgi' // trim(str2) // '.dat'
    open(unit=3,file=str3,status='old')
    read(unit=3,fmt=*) Adgi
    end if

!! Fortran-assembly
    if (selectPackage == 3) then
    str3 = '../../external/puffin/data/Amats/subdom' // trim(str1) // '/OrigADgi' // trim(str2) // '.dat'
    open(unit=3,file=str3,status='replace')
    write(unit=3,fmt=*) Adgi
    end if

    close(3)

!!-----------------------------------------------------------------------------------------------
    if ((pid==0) .and. (k==1)) then
        if (selectPackage == 3) then
            print*, 'PETSc-Mat assembly using: In-house Fortran assembly'
        else
            print*,'PETSc-Mat assembly using:', str3
        end if
    end if

!!-----------------------------------------------------------------------------------------------
!! PETSc-Mat setting: Stochastic DD block matrices using deteministic DD blocks
    do i = 1,npceout
    do j = 1,npceout

        if ((k .eq. ijk(indexi,1)) .and. (i .eq. ijk(indexi,2)) .and. (j .eq. ijk(indexi,3))) then

    !!---------------------------------------->Asii
        nni = ((i-1)*nid)
        nnj = ((j-1)*nid)

        do id = 1,nid
        do jd = 1,nid
            if (Adii(id,jd) .ne. 0) then
                ii = (nni+(id-1))
                jj = (nnj+(jd-1))
                imat = cijk(indexi)*Adii(id,jd)
                call MatSetValues(Asii,one,ii,one,jj,imat,ADD_VALUES,ierr)
            end if
        end do
        end do

    !!---------------------------------------->Asgg
        nni = ((i-1)*nb)
        nnj = ((j-1)*nb)

        do id = 1,nb
        do jd = 1,nb
            if (Adgg(id,jd) .ne. 0) then
                ii = (nni+(id-1))
                jj = (nnj+(jd-1))
                imat = cijk(indexi)*Adgg(id,jd)
                call MatSetValues(Asgg,one,ii,one,jj,imat,ADD_VALUES,ierr)
            end if
        end do
        end do

    !!---------------------------------------->Asgi
        nni = ((i-1)*nb)
        nnj = ((j-1)*nid)

        do id = 1,nb
        do jd = 1,nid
            if (Adgi(id,jd) .ne. 0) then
                ii = (nni+(id-1))
                jj = (nnj+(jd-1))
                imat = cijk(indexi)*Adgi(id,jd)
                call MatSetValues(Asgi,one,ii,one,jj,imat,ADD_VALUES,ierr)
            end if
        end do
        end do

    !!---------------------------------------->Asig !! don't need due to symmetry(Asgi==Asig)
        !nni = ((i-1)*nid)
        !nnj = ((j-1)*nb)
        !
        !do id = 1,nid
        !do jd = 1,nb
        !    if (Adig(id,jd) .ne. 0) then
        !        ii = (nni+(id-1))
        !        jj = (nnj+(jd-1))
        !        imat = cijk(indexi)*Adig(id,jd)
        !        call MatSetValues(Asig,one,ii,one,jj,imat,ADD_VALUES,ierr)
        !    end if
        !end do
        !end do

        indexi = indexi+1
        end if

    end do
    end do

END DO

    call MatAssemblyBegin(Asii,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyBegin(Asgg,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyBegin(Asgi,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(Asii,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(Asgg,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(Asgi,MAT_FINAL_ASSEMBLY,ierr)
    !!call MatAssemblyBegin(Asig,MAT_FINAL_ASSEMBLY,ierr)
    !!call MatAssemblyEnd(Asig,MAT_FINAL_ASSEMBLY,ierr)

    DEALLOCATE(Adii,Adgg,Adgi)
    DEALLOCATE(nnzb,nnzi,nnzbi)


END SUBROUTINE StoMatSeqOneLevel


!!!*********************************************************
!!! Subroutine: To extract non-zero elements from PETSc Matrix Assembly
SUBROUTINE GetMallocs(pid, p, e, t,np,ne,nt,nb,ndim,npcein,npceout,nomga,nip,nbp, &
                      casep,ncijk,ijk,cijk,dbounds,const_diff,omegas,multipliers, &
                      sigma,mIndex,sIndex,ierr)

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>


    PetscInt         nip, nbp, one, iNZ  !!nzi,nzb,
    PetscInt         id, jd !!ii, jj, nni, nnj
    PetscErrorCode   ierr

    !!PetscInt                :: m, n!
    PetscInt,ALLOCATABLE    :: idxm(:), idxn(:)!! nnzb(:) !!nnzbi(:), nnzi(:)
    PetscInt,ALLOCATABLE    :: idxmi(:), idxni(:), idxmbi(:), idxnbi(:)

    integer :: k, kk, indexi, nid, ncijk, casep, pid, npcein
    integer :: np, ne, nt, nb, ndim, npceout, nomga
    double precision :: const_diff, sigma

    integer, dimension(ncijk,3)       :: ijk
    integer, dimension(3,ne)          :: e
    integer, dimension(3,nt)          :: t
    double precision, dimension(2,np) :: p
    double precision, dimension(ncijk):: cijk
    double precision, dimension(2,2)  :: dbounds
    double precision, dimension(nomga):: omegas, multipliers

    integer, dimension(ndim,npceout)  :: mIndex
    integer, dimension(2,ndim)        :: sIndex
    integer, dimension(npceout)       :: nZcijk
    character(len=255)                :: str1, str2, str3

    double precision, allocatable, dimension(:,:) :: Adii,Adgg,Adgi

    integer :: selectPackage

!!-----------------------------------------------------------------------------------------------
    nid    = np-nb
    one    = 1
    indexi = 1

!!-----------------------------------------------------------------------------------------------
    selectPackage = 1 !! 1-Dolfin, 2-Puffin, 3-Fortran  !! 1 and 3 validated
!!-----------------------------------------------------------------------------------------------

    ALLOCATE (idxm(nb),idxn(nbp))
    ALLOCATE (idxmi(nid),idxni(nip))
    ALLOCATE (idxmbi(nb),idxnbi(nbp))
    allocate (Adii(nid,nid),Adgg(nb,nb),Adgi(nb,nid))

!!-----------------------------------------------------------------------------------------------
!! Precalculated nZijk for each PCE order and Dimension
    open(unit=2,file='../../data/klePceData/nZijk.dat',status='old')
    read(unit=2,fmt='(I8)') nZcijk
    close(2)

!!-----------------------------------------------------------------------------------------------
!! non-zeros(ADii)
    kk = 1   !! kk=1, because the non-zero structure is assumed to be same for each PCE mode
    call int2str(str1,pid+1,1)
    call int2str(str2,kk,1)

    !! For Dolfin we use ADii***2.dat: because 1-matrix has different non-zero structure
    if (selectPackage == 1) then
    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/ADii2.dat'
    open(unit=1,file=str3,status='old')
    read(unit=1,fmt=*) Adii
    end if

    if (selectPackage == 2) then
    str3 = '../../external/puffin/data/Amats/subdom' // trim(str1) // '/ADii0001.dat'
    open(unit=1,file=str3,status='old')
    read(unit=1,fmt=*) Adii
    end if

    close(1)

!!-----------------------------------------------------------------------------------------------
!! non-zeros(ADgg)
    call int2str(str1,pid+1,1)
    call int2str(str2,kk,1)

    if (selectPackage == 1) then
    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/ADgg2.dat'
    open(unit=2,file=str3,status='old')
    read(unit=2,fmt=*) Adgg
    end if

    if (selectPackage == 2) then
    str3 = '../../external/puffin/data/Amats/subdom' // trim(str1) // '/ADgg0001.dat'
    open(unit=2,file=str3,status='old')
    read(unit=2,fmt=*) Adgg
    end if

    close(2)

!!-----------------------------------------------------------------------------------------------
!! non-zeros(ADgi)
    call int2str(str1,pid+1,1)
    call int2str(str2,kk,1)

    if (selectPackage == 1) then
    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/ADgi2.dat'
    open(unit=3,file=str3,status='old')
    read(unit=3,fmt=*) Adgi
    end if

    if (selectPackage == 2) then
    str3 = '../../external/puffin/data/Amats/subdom' // trim(str1) // '/ADgi0001.dat'
    open(unit=3,file=str3,status='old')
    read(unit=3,fmt=*) Adgi
    end if

    close(3)

!!-----------------------------------------------------------------------------------------------
    if (pid==0) then
        if (selectPackage == 3) then
            print*, 'PETSc-mallocs using: In-house Fortran assembly'
        else
            print*,'PETSc-mallocs using:', str3
        end if
    end if

!!-----------------------------------------------------------------------------------------------
!! For this option we use Fortran assemble routines (in-house)
    if (selectPackage == 3) then
    Adii = 0.0d0
    Adgg = 0.0d0
    Adgi = 0.0d0
    kk = 1
        if (kk .eq. 1) then
            !! Det-Advection Matrix
            call SubAssembleMatrix(Adii,Adgg,Adgi,p,e,t,np,ne,nt,nb,ndim,npceout,nomga,1,0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],&
            [0.0d0,0.0d0,0.0d0,0.0d0],0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],dbounds,mIndex,sIndex,0,0)
            !! Det-Diffusion Matrix
            call SubAssembleMatrix(Adii,Adgg,Adgi,p,e,t,np,ne,nt,nb,ndim,npceout,nomga,2,const_diff,omegas,multipliers,&
            sigma,[0.0d0,0.0d0,0.0d0,0.0d0],dbounds,mIndex,sIndex,casep,1)
            else
            !! Sto-Diffusion Matrices
            call SubAssembleMatrix(Adii,Adgg,Adgi,p,e,t,np,ne,nt,nb,ndim,npceout,nomga,2,const_diff,omegas,multipliers,&
            sigma,[0.0d0,0.0d0,0.0d0,0.0d0],dbounds,mIndex,sIndex,casep,k)
        end if
    end if

!!-----------------------------------------------------------------------------------------------
!! Here, using first level non-Zero & cijk structure we extrac second level non-Zero structure
    DO k = 1,npceout

        do id = 1,nb
            iNZ = 0
            do jd = 1,nb
            if (Adgg(id,jd) .ne. 0) then
            !if (abs(Adgg(id,jd)) .gt. 0.000000001) then
            iNZ = iNZ+1
            end if
            end do
            idxm(id) = (iNZ)
        end do
        !print*,'-----------------'
        !print*,nb
        !print*,'--------$$-------'
        !print*,idxm
        idxn(((k-1)*nb+1):(k*nb)) = nZcijk(k)*idxm
        !print*,'-------**--------'
        !print*,idxn


        do id = 1,nid
            iNZ = 0
            do jd = 1,nid
            if (Adii(id,jd) .ne. 0) then
            !if (abs(Adii(id,jd)) .gt. 0.000000001) then
            iNZ = iNZ+1
            end if
            end do
            idxmi(id) = (iNZ)
        end do
        idxni(((k-1)*nid+1):(k*nid)) = nZcijk(k)*idxmi


        do id = 1,nb
            iNZ = 0
            do jd = 1,nid
            if (Adgi(id,jd) .ne. 0) then
            iNZ = iNZ+1
            end if
            end do
            idxmbi(id) = (iNZ)
        end do
        idxnbi(((k-1)*nb+1):(k*nb)) = nZcijk(k)*idxmbi

    END DO

!!-----------------------------------------------------------------------------------------------
!! Write mallocs files for respective matrices  !! Symmetry structure exploited
    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzb' // trim(str1) // '.dat'    !!**
    open(unit=1,file=str1,status='replace')
    write(1,*) idxn
    close(1)

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzi' // trim(str1) // '.dat'    !!**
    open(unit=1,file=str1,status='replace')
    write(1,*) idxni
    close(1)

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzbi' // trim(str1) // '.dat'    !!**
    open(unit=1,file=str1,status='replace')
    write(1,*) idxnbi
    close(1)

    DEALLOCATE(Adii,Adgg,Adgi)
    DEALLOCATE(idxn,idxm,idxni,idxmi,idxnbi,idxmbi)

END SUBROUTINE GetMallocs


!!!*********************************************************
!!! Subroutine: To construc PETSc Matrix in most general way
!SUBROUTINE GetMatSeqOneLevel(n,m,PetMat, pg, eg, tg, npg, neg, ntg, amp, dbounds, meanc, ierr)
SUBROUTINE StoMatSeqOneLevelFullMallocs(pid,Asii,Asgg,Asgi, p, e, t, np, ne, nt,nb,ndim,npceout,nomga,nip,nbp, &
             casep,ncijk,ijk,cijk,dbounds,const_diff,omegas,multipliers,sigma,mIndex,sIndex,ierr)

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>


    Mat              Asii,Asgg,Asgi
    PetscInt         nip, nbp, one !!nzi, nzb
    PetscInt         ii, jj, id, jd, nni, nnj
    PetscScalar      imat
    PetscErrorCode   ierr

    PetscInt,ALLOCATABLE    :: nnzi(:), nnzb(:), nnzbi(:)
    character(len=255)      :: str1

    integer :: i, j, k, indexi, nid, ncijk, casep, pid
    integer :: np, ne, nt, nb, ndim, npceout, nomga
    double precision :: const_diff, sigma

    integer, dimension(ncijk,3)       :: ijk
    integer, dimension(3,ne)          :: e
    integer, dimension(3,nt)          :: t
    double precision, dimension(2,np) :: p
    double precision, dimension(ncijk):: cijk
    double precision, dimension(2,2)  :: dbounds
    double precision, dimension(nomga):: omegas, multipliers

    integer, dimension(ndim,npceout)  :: mIndex
    integer, dimension(2,ndim)        :: sIndex

    double precision, allocatable, dimension(:,:)  :: Adii,Adgg,Adgi  

!!-----------------------------------------------------------------------------------------------
    allocate(Adii((np-nb),(np-nb)), Adgg(nb,nb), Adgi(nb,(np-nb)))

!! number of non-zeros per row for Aii&Agi/Agg  matrix mallocs
!! For nOrd = 1; nDim = 50; nMesh=15K; nzi = 600,  nzb = 300
!! For nOrd = 2; nDim = 15; nMesh=15K; nzi = 1500, nzb = 900
!! For nOrd = 2; nDim = 20; nMesh=15K; nzi = 2200, nzb = 1200
!! For nOrd = 2; nDim = 25; nMesh=15K; nzi = 3200, nzb = 2000

!    nzi = 1   !! ndim=50: 6000
!    nzb = 1   !! ndim=50; 3500
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, nip, nip, nzi, PETSC_NULL_INTEGER, Asii, ierr)
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, nbp, nbp, nzb, PETSC_NULL_INTEGER, Asgg, ierr)
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, nbp, nip, nzb, PETSC_NULL_INTEGER, Asgi, ierr)
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, nip, nbp, nzb, PETSC_NULL_INTEGER, Asig, ierr)


allocate(nnzi(nip),nnzb(nbp),nnzbi(nbp))

call int2str(str1,pid+1,4)
str1 = '../../data/mallocData/nnzi' // trim(str1) // '.dat'    !!**
open(unit=1,file=str1,status='old')
read(unit=1,fmt=*) nnzi
close(1)

call int2str(str1,pid+1,4)
str1 = '../../data/mallocData/nnzb' // trim(str1) // '.dat'    !!**
open(unit=2,file=str1,status='old')
read(unit=2,fmt=*) nnzb
close(2)

call int2str(str1,pid+1,4)
str1 = '../../data/mallocData/nnzbi' // trim(str1) // '.dat'   !!**
open(unit=3,file=str1,status='old')
read(unit=3,fmt=*) nnzbi
close(3)

call MatCreateSeqAIJ(PETSC_COMM_SELF, nip, nip, 0, nnzi, Asii, ierr)
call MatCreateSeqAIJ(PETSC_COMM_SELF, nbp, nbp, 0, nnzb, Asgg, ierr)
call MatCreateSeqAIJ(PETSC_COMM_SELF, nbp, nip, 0, nnzbi, Asgi, ierr)


call MatSetOption(Asii,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE,ierr)
call MatSetOption(Asgg,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE,ierr)
call MatSetOption(Asgi,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE,ierr)
call MatSetOption(Asii,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,ierr)
call MatSetOption(Asgg,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,ierr)
call MatSetOption(Asgi,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,ierr)
!!-----------------------------------------------------------------------------------------------
!! Method : 3
    nid    = np-nb
    one    = 1
    indexi = 1

    do k = 1,npceout

        !! Deterministic Matices
        Adii = 0.0d0
        Adgg = 0.0d0
        Adgi = 0.0d0

        if (k .eq. 1) then
            !! Det-Advection Matrix
            call SubAssembleMatrix(Adii,Adgg,Adgi,p,e,t,np,ne,nt,nb,ndim,npceout,nomga,1,0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],&
                    [0.0d0,0.0d0,0.0d0,0.0d0],0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],dbounds,mIndex,sIndex,0,0)
            !! Det-Diffusion Matrix
            call SubAssembleMatrix(Adii,Adgg,Adgi,p,e,t,np,ne,nt,nb,ndim,npceout,nomga,2,const_diff,omegas,multipliers,&
                    sigma,[0.0d0,0.0d0,0.0d0,0.0d0],dbounds,mIndex,sIndex,casep,1)
        else
            !! Sto-Diffusion Matrices
            call SubAssembleMatrix(Adii,Adgg,Adgi,p,e,t,np,ne,nt,nb,ndim,npceout,nomga,2,const_diff,omegas,multipliers,&
                    sigma,[0.0d0,0.0d0,0.0d0,0.0d0],dbounds,mIndex,sIndex,casep,k)
        end if

        do i = 1,npceout
        do j = 1,npceout

            if ((k .eq. ijk(indexi,1)) .and. (i .eq. ijk(indexi,2)) .and. (j .eq. ijk(indexi,3))) then

            nni = ((i-1)*nid)
            nnj = ((j-1)*nid)

            do id = 1,nid
            do jd = 1,nid
                if (Adii(id,jd) .ne. 0) then
                    ii = (nni+(id-1))
                    jj = (nnj+(jd-1))
                    imat = cijk(indexi)*Adii(id,jd)
                    call MatSetValues(Asii,one,ii,one,jj,imat,ADD_VALUES,ierr)
                end if
            end do
            end do

            nni = ((i-1)*nb)
            nnj = ((j-1)*nb)

            do id = 1,nb
            do jd = 1,nb
                if (Adgg(id,jd) .ne. 0) then
                    ii = (nni+(id-1))
                    jj = (nnj+(jd-1))
                    imat = cijk(indexi)*Adgg(id,jd)
                    call MatSetValues(Asgg,one,ii,one,jj,imat,ADD_VALUES,ierr)
                end if
            end do
            end do

!            nni = ((i-1)*nid)
!            nnj = ((j-1)*nb)
!
!            do id = 1,nid
!            do jd = 1,nb
!                if (Adig(id,jd) .ne. 0) then
!                    ii = (nni+(id-1))
!                    jj = (nnj+(jd-1))
!                    imat = cijk(indexi)*Adig(id,jd)
!                    call MatSetValues(Asig,one,ii,one,jj,imat,ADD_VALUES,ierr)
!                end if
!            end do
!            end do

            nni = ((i-1)*nb)
            nnj = ((j-1)*nid)

            do id = 1,nb
            do jd = 1,nid
                if (Adgi(id,jd) .ne. 0) then
                    ii = (nni+(id-1))
                    jj = (nnj+(jd-1))
                    imat = cijk(indexi)*Adgi(id,jd)
                    call MatSetValues(Asgi,one,ii,one,jj,imat,ADD_VALUES,ierr)
                end if
            end do
            end do

            indexi = indexi+1
            end if

        end do
        end do
    end do

    call MatAssemblyBegin(Asii,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyBegin(Asgg,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyBegin(Asgi,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(Asii,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(Asgg,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(Asgi,MAT_FINAL_ASSEMBLY,ierr)
    !!call MatAssemblyBegin(Asig,MAT_FINAL_ASSEMBLY,ierr)
    !!call MatAssemblyEnd(Asig,MAT_FINAL_ASSEMBLY,ierr)

    deallocate(Adii,Adgg,Adgi)

deallocate(nnzi,nnzb,nnzbi)

END SUBROUTINE StoMatSeqOneLevelFullMallocs

!!!!*********************************************************
subroutine GetRs(pid,nb,nbg,npceout,RsMat)

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>


    Mat              RsMat
    PetscScalar      imat
    PetscInt         one, nbp, nbgp, ii, jj
    PetscErrorCode   ierr

    character(len=255) :: str1
    integer :: nb, nbg, i, pid, temp1, npceout, j
    !double precision, dimension(:,:) :: Rmat

    one  = 1
    nbp  = nb*npceout
    nbgp = nbg*npceout

    call MatCreateSeqAIJ(PETSC_COMM_SELF, nbp, nbgp, one, PETSC_NULL_INTEGER, RsMat,ierr)

    call int2str(str1,pid+1,4)
    str1 = 'bnodes' // trim(str1) // '.dat'
    open(unit=1,file=str1,status='old')
    !Rmat(:,:) = 0.0d0
    do i = 1,nb
    read(unit=1,fmt=*) temp1

        do j = 1,npceout
            ii = (i+(j-1)*nb)-1
            jj = (temp1+(j-1)*nbg)-1
            imat = 1.0
            call MatSetValues(RsMat,one,ii,one,jj,imat,INSERT_VALUES,ierr)
        end do

    end do
    close(1)

    call MatAssemblyBegin(RsMat,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(RsMat,MAT_FINAL_ASSEMBLY,ierr)

end subroutine GetRs


!!!!*********************************************************
!!!! Subroutine: To call PETSc Assembly Vector in most general way
!SUBROUTINE GetPetVec(n,PetVec,SolVec,pg,eg,tg,npg,neg,ntg,amp,dbounds,ierr)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!
!    Vec              PetVec,SolVec
!    PetscInt         i,ii
!    PetscInt         n,one
!    PetscScalar      ivec
!    PetscErrorCode   ierr
!
!    integer :: npg, neg, ntg
!    double precision :: amp
!    integer, dimension(3,neg) :: eg
!    integer, dimension(3,ntg) :: tg
!    double precision, dimension(2,npg) :: pg
!    double precision, dimension(2,2) :: dbounds
!
!    call VecCreate(PETSC_COMM_WORLD,PetVec,ierr)
!    call VecSetFromOptions(PetVec,ierr)
!    call VecSetSizes(PetVec,PETSC_DECIDE,n,ierr)
!    call VecDuplicate(PetVec,SolVec,ierr)
!
!    call PETScAssemblevector(PetVec, SolVec, pg, eg, tg, npg, neg, ntg, amp, dbounds)
!
!    call VecAssemblyBegin(PetVec,ierr)
!    call VecAssemblyEnd(PetVec,ierr)
!
!END SUBROUTINE GetPetVec
!
!
!!!!*********************************************************
!!!! Subroutine: To construc PETSc Vector in Sequential Form
!SUBROUTINE GetPetVecSeq(n,fvec,PetVec,ierr)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!
!    Vec              PetVec
!    PetscInt         i,ii
!    PetscInt         n,one
!    PetscScalar      ivec
!    PetscErrorCode   ierr
!
!    double precision, dimension(n) :: fvec
!    !integer, dimension(n)          :: temp1
!
!    call VecCreateSeq(PETSC_COMM_SELF, n, PetVec, ierr)
!    !call VecDuplicate(PetVec,SolVec,ierr)
!
!    one = 1
!    do i = 1,n
!        ii = i-1
!!        temp1(i) = i-1
!        ivec = fvec(i)
!        call VecSetValues(PetVec,one,ii,ivec,INSERT_VALUES,ierr)
!    end do
!
!    call VecAssemblyBegin(PetVec,ierr)
!    call VecAssemblyEnd(PetVec,ierr)
!
!END SUBROUTINE GetPetVecSeq
!
!
!!!!*********************************************************
!!!! Subroutine: To construc PETSc Matrix in most general way
!SUBROUTINE GetPetMat(n,PetMat, pg, eg, tg, npg, neg, ntg, amp, dbounds, meanc, ierr)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!#include <petsc/finclude/petscmat.h>
!
!
!    Mat              PetMat
!    PetscInt         i,j,ii,jj
!    PetscInt         n,one
!    PetscScalar      imat
!    PetscErrorCode   ierr
!
!    integer :: npg, neg, ntg
!    double precision :: amp, meanc
!    integer, dimension(3,neg) :: eg
!    integer, dimension(3,ntg) :: tg
!    double precision, dimension(2,npg) :: pg
!    double precision, dimension(2,2) :: dbounds
!
!    call MatCreate(PETSC_COMM_WORLD,PetMat,ierr)
!    call MatSetFromOptions(PetMat,ierr)
!    call MatSetSizes(PetMat,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr)
!    call MatSetUp(PetMat,ierr)
!
!    call PETScAssembleADVmatrix(PetMat,pg,eg,tg,npg,neg,ntg,1,0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],&
!                               [0.0d0,0.0d0,0.0d0,0.0d0],0.0d0,[0.0d0,0.0d0,0.0d0],dbounds,0,0)
!    call PETScAssembleDIFmatrix(PetMat,pg,eg,tg,npg,neg,ntg,2,meanc,[0.0d0,0.0d0,0.0d0,0.0d0],&
!                               [0.0d0,0.0d0,0.0d0,0.0d0], 0.0d0,[0.0d0,0.0d0,0.0d0],dbounds,1,1)
!
!    call MatAssemblyBegin(PetMat,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyEnd(PetMat,MAT_FINAL_ASSEMBLY,ierr)
!
!END SUBROUTINE GetPetMat
!
!
!!!!*********************************************************
!!!! Subroutine: To call PETSc Assembly Vector in most general way
!SUBROUTINE GetPetVecSelf(n,PetVec,SolVec,pg,eg,tg,npg,neg,ntg,amp,dbounds,ierr)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!
!    Vec              PetVec,SolVec
!    PetscInt         i,ii
!    PetscInt         n,one
!    PetscScalar      ivec
!    PetscErrorCode   ierr
!
!    integer :: npg, neg, ntg
!    double precision :: amp
!    integer, dimension(3,neg) :: eg
!    integer, dimension(3,ntg) :: tg
!    double precision, dimension(2,npg) :: pg
!    double precision, dimension(2,2) :: dbounds
!
!    call VecCreateSeq(PETSC_COMM_SELF, n, PetVec, ierr)
!    call VecDuplicate(PetVec,SolVec,ierr)
!
!    call PETScAssemblevector(PetVec, SolVec, pg, eg, tg, npg, neg, ntg, amp, dbounds)
!
!    call VecAssemblyBegin(PetVec,ierr)
!    call VecAssemblyEnd(PetVec,ierr)
!
!END SUBROUTINE GetPetVecSelf
!
!!!!*********************************************************
!!!! Subroutine: To construc PETSc Matrix in most general way
!SUBROUTINE GetPetMatSelf(n,PetMat, pg, eg, tg, npg, neg, ntg, amp, dbounds, meanc, ierr)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!#include <petsc/finclude/petscmat.h>
!
!    Mat              PetMat
!    PetscInt         i,j,ii,jj
!    PetscInt         n,one
!    PetscScalar      imat
!    PetscErrorCode   ierr
!
!    integer :: npg, neg, ntg
!    double precision :: amp, meanc
!    integer, dimension(3,neg) :: eg
!    integer, dimension(3,ntg) :: tg
!    double precision, dimension(2,npg) :: pg
!    double precision, dimension(2,2) :: dbounds
!
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, n, n, n,PETSC_NULL_INTEGER, PetMat, ierr)
!
!    call PETScAssembleADVmatrix(PetMat,pg,eg,tg,npg,neg,ntg,1,0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],&
!    [0.0d0,0.0d0,0.0d0,0.0d0],0.0d0,[0.0d0,0.0d0,0.0d0],dbounds,0,0)
!    call PETScAssembleDIFmatrix(PetMat,pg,eg,tg,npg,neg,ntg,2,meanc,[0.0d0,0.0d0,0.0d0,0.0d0],&
!    [0.0d0,0.0d0,0.0d0,0.0d0], 0.0d0,[0.0d0,0.0d0,0.0d0],dbounds,1,1)
!
!    call MatAssemblyBegin(PetMat,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyEnd(PetMat,MAT_FINAL_ASSEMBLY,ierr)
!
!END SUBROUTINE GetPetMatSelf
!
!
!!!!*********************************************************
!!!! Subroutine: To construc PETSc Vector in most general way
!SUBROUTINE GetVec(n,fvec,PetVec,Solvec,ierr)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!
!    Vec              PetVec,SolVec
!    PetscInt         i,ii
!    PetscInt         n,one
!    PetscScalar      ivec
!    PetscErrorCode   ierr
!
!    double precision, dimension(n) :: fvec
!
!    call VecCreate(PETSC_COMM_WORLD,PetVec,ierr)
!    call VecSetFromOptions(PetVec,ierr)
!    call VecSetSizes(PetVec,PETSC_DECIDE,n,ierr)
!    call VecDuplicate(PetVec,SolVec,ierr)
!
!    one = 1
!    do i = 1,n
!    ii = i-1
!    ivec = fvec(i)
!    call VecSetValues(PetVec,one,ii,ivec,INSERT_VALUES,ierr)
!    end do
!
!    call VecAssemblyBegin(PetVec,ierr)
!    call VecAssemblyEnd(PetVec,ierr)
!
!END SUBROUTINE GetVec
!
!
!!!!*********************************************************
!!!! Subroutine: To construc PETSc Matrix in most general way
!SUBROUTINE GetMat(n,AMAT,PetMat,ierr)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!#include <petsc/finclude/petscmat.h>
!
!
!    Mat              PetMat
!    PetscInt         i,j,ii,jj
!    PetscInt         n,one
!    PetscScalar      imat
!    PetscErrorCode   ierr
!
!    double precision, dimension(n,n) :: Amat
!
!    call MatCreate(PETSC_COMM_WORLD,PetMat,ierr)
!    call MatSetFromOptions(PetMat,ierr)
!    call MatSetSizes(PetMat,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr)
!    call MatSetUp(PetMat,ierr)
!
!    one = 1
!    do i = 1,n
!    do j = 1,n
!    if (Amat(i,j) .ne. 0) then
!    ii = i-1
!    jj = j-1
!    imat = Amat(i,j)
!    call MatSetValues(PetMat,one,ii,one,jj,imat,INSERT_VALUES,ierr)
!    end if
!    end do
!    end do
!
!    call MatAssemblyBegin(PetMat,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyEnd(PetMat,MAT_FINAL_ASSEMBLY,ierr)
!
!END SUBROUTINE GetMat



!!!!*********************************************************
!SUBROUTINE StoVecSeqOneLevel(Fsi,Fsg,SolVeci,SolVecb,tempi,tempb,p,e,t,np,ne,nt,nb,npceout,amp,dbounds,ierr)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!
!    Vec              Fsi,Fsg,SolVeci,SolVecb
!    PetscInt         i,ii
!    PetscInt         one
!    PetscScalar      ivec
!    PetscErrorCode   ierr
!
!    integer :: np, ne, nt, nb, npceout, nip, nbp, id
!    double precision                  :: amp
!    integer, dimension(3,ne)          :: e
!    integer, dimension(3,nt)          :: t
!    double precision, dimension(2,np) :: p
!    double precision, dimension(2,2)  :: dbounds
!
!    integer, dimension(np-nb*npceout) :: tempi
!    integer, dimension(nb*npceout)    :: tempb
!
!    double precision, dimension((np-nb)*npceout)  :: Fdi
!    double precision, dimension(nb*npceout)       :: Fdg
!
!    nip = ((np-nb)*npceout)
!    nbp = nb*npceout
!
!    call VecCreateSeq(PETSC_COMM_SELF, nip, Fsi, ierr)
!    call VecCreateSeq(PETSC_COMM_SELF, nbp, Fsg, ierr)
!
!    call VecDuplicate(Fsi,SolVeci,ierr)
!    call VecDuplicate(Fsg,SolVecb,ierr)
!
!    !call PETScOneLevelAssemVec(Fdi, Fdg, p, e, t, np, ne, nt, nb, amp, dbounds,tempi,tempb)
!    call SubAssembleVector(Fdi, Fdg, p, e, t, np, ne, nt, nb, npceout, amp, dbounds)
!
!    do id = 1,nip
!        tempi(id) = (id-1)
!        if (Fdi(id) .ne. 0) then
!            ii = (id-1)
!            ivec = Fdi(id)
!            call VecSetValues(Fsi,one,ii,ivec,ADD_VALUES,ierr)
!        end if
!    end do
!
!    do id = 1,nbp
!        tempb(id) = id-1
!        if (Fdi(id) .ne. 0) then
!            ii = (id-1)
!            ivec = Fdg(id)
!            call VecSetValues(Fsg,one,ii,ivec,ADD_VALUES,ierr)
!        end if
!    end do
!
!    call VecAssemblyBegin(Fsi,ierr)
!    call VecAssemblyBegin(Fsg,ierr)
!    call VecAssemblyEnd(Fsi,ierr)
!    call VecAssemblyEnd(Fsg,ierr)
!
!
!END SUBROUTINE StoVecSeqOneLevel

!!!!*********************************************************
!SUBROUTINE GetVecSeqOneLevel(PFi,PFg,SolVeci,SolVecb,tempi,tempb,p,e,t,np,ne,nt,nb,amp,dbounds,ierr)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!
!    Vec              PFi,PFg,SolVeci,SolVecb
!    PetscInt         i,ii
!    PetscInt         ni,one
!    PetscScalar      ivec
!    PetscErrorCode   ierr
!
!    integer :: np, ne, nt, nb
!    double precision :: amp
!    integer, dimension(3,ne) :: e
!    integer, dimension(3,nt) :: t
!    double precision, dimension(2,np) :: p
!    double precision, dimension(2,2) :: dbounds
!
!    integer, dimension(np-nb) :: tempi
!    integer, dimension(nb) :: tempb
!
!    ni = np-nb
!    call VecCreateSeq(PETSC_COMM_SELF, ni, PFi, ierr)
!    call VecCreateSeq(PETSC_COMM_SELF, nb, PFg, ierr)
!
!    call VecDuplicate(PFi,SolVeci,ierr)
!    call VecDuplicate(PFg,SolVecb,ierr)
!
!    call PETScOneLevelAssemVec(PFi, PFg, p, e, t, np, ne, nt, nb, amp, dbounds,tempi,tempb)
!
!    call VecAssemblyBegin(PFi,ierr)
!    call VecAssemblyBegin(PFg,ierr)
!    call VecAssemblyEnd(PFi,ierr)
!    call VecAssemblyEnd(PFg,ierr)
!
!
!END SUBROUTINE GetVecSeqOneLevel
!
!
!!!!*********************************************************
!!!! Subroutine: To construc PETSc Matrix in most general way
!!SUBROUTINE GetMatSeqOneLevel(n,m,PetMat, pg, eg, tg, npg, neg, ntg, amp, dbounds, meanc, ierr)
!SUBROUTINE GetMatSeqOneLevel(PAii,PAgg,PAig,PAgi, p, e, t, np, ne, nt, nb, dbounds, const_diff, ierr)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!#include <petsc/finclude/petscmat.h>
!
!    Mat              PAii,PAgg,PAig,PAgi
!    PetscInt         i,j,ii,jj
!    PetscInt         ni,m,nz,one
!    PetscScalar      imat
!    PetscErrorCode   ierr
!
!    integer :: np, ne, nt, nb
!    double precision :: const_diff
!    integer, dimension(3,ne) :: e
!    integer, dimension(3,nt) :: t
!    double precision, dimension(2,np) :: p
!    double precision, dimension(2,2) :: dbounds
!
!    ni = np-nb
!    nz = 12  !! number of non-zeros per row for matrix pre-allocation
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, ni, ni, nz, PETSC_NULL_INTEGER, PAii, ierr)
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, nb, nb, nz, PETSC_NULL_INTEGER, PAgg, ierr)
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, ni, nb, nz, PETSC_NULL_INTEGER, PAig, ierr)
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, nb, ni, nz, PETSC_NULL_INTEGER, PAgi, ierr)
!
!    call PETScOneLevelAssemADVmatrix(PAii,PAgg,PAig,PAgi,p,e,t,np,ne,nt,nb,1,0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],&
!    [0.0d0,0.0d0,0.0d0,0.0d0],0.0d0,[0.0d0,0.0d0,0.0d0],dbounds,0,0)
!    call PETScOneLevelAssemDIFmatrix(PAii,PAgg,PAig,PAgi,p,e,t,np,ne,nt,nb,2,const_diff,[0.0d0,0.0d0,0.0d0,0.0d0],&
!    [0.0d0,0.0d0,0.0d0,0.0d0], 0.0d0,[0.0d0,0.0d0,0.0d0],dbounds,1,1)
!
!    call MatAssemblyBegin(PAii,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyBegin(PAgg,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyBegin(PAig,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyBegin(PAgi,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyEnd(PAii,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyEnd(PAgg,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyEnd(PAig,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyEnd(PAgi,MAT_FINAL_ASSEMBLY,ierr)
!
!
!END SUBROUTINE GetMatSeqOneLevel
!
!
!!!!*********************************************************
!SUBROUTINE GetVecSeqTwoLevel(PFi,PFg,PFc,SolVeci,SolVecb,SolVecc,tempi,tempb,tempc,p,e,t,np,ne,nt,nb,nci,amp,dbounds,ierr)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!
!    Vec              PFi,PFg,PFc
!    Vec              SolVeci,SolVecb,SolVecc
!    PetscInt         i,ii
!    PetscInt         ni,one
!    PetscScalar      ivec
!    PetscErrorCode   ierr
!
!    integer :: np, ne, nt, nb, nci
!    double precision :: amp
!    integer, dimension(3,ne) :: e
!    integer, dimension(3,nt) :: t
!    double precision, dimension(2,np) :: p
!    double precision, dimension(2,2) :: dbounds
!
!    integer, dimension(np-nb) :: tempi
!    integer, dimension(nb)    :: tempb
!    integer, dimension(nci)   :: tempc
!
!    ni = np-nb
!
!    call VecCreateSeq(PETSC_COMM_SELF, ni, PFi, ierr)
!    call VecCreateSeq(PETSC_COMM_SELF, nb, PFg, ierr)
!    call VecCreateSeq(PETSC_COMM_SELF, nci, PFc, ierr)
!
!    call VecDuplicate(PFi,SolVeci,ierr)
!    call VecDuplicate(PFg,SolVecb,ierr)
!    call VecDuplicate(PFc,SolVecc,ierr)
!
!    call PETScTwoLevelAssemVec(PFi, PFg, PFc, p, e, t, np, ne, nt, nb, nci,amp,dbounds,tempi,tempb,tempc)
!
!    call VecAssemblyBegin(PFi,ierr)
!    call VecAssemblyBegin(PFg,ierr)
!    call VecAssemblyBegin(PFc,ierr)
!    call VecAssemblyEnd(PFi,ierr)
!    call VecAssemblyEnd(PFg,ierr)
!    call VecAssemblyEnd(PFc,ierr)
!
!
!END SUBROUTINE GetVecSeqTwoLevel
!
!!!!*********************************************************
!!!! Subroutine: To construc PETSc Matrix in most general way
!!SUBROUTINE GetMatSeqOneLevel(n,m,PetMat, pg, eg, tg, npg, neg, ntg, amp, dbounds, meanc, ierr)
!SUBROUTINE GetMatSeqTwoLevel(PAmat,PAcc,PAic,PAci,PAir,PAri,PArc,PAcr,p, e, t, &
!                             np, ne, nt, nb, nci, dbounds, const_diff, ierr)
!           
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!#include <petsc/finclude/petscmat.h>
!
!    Mat              PAmat,PAcc,PAic,PAci
!    Mat              PAir,PAri,PArc,PAcr
!    PetscInt         i,j,ii,jj
!    PetscInt         ni,m,nz,one
!    PetscErrorCode   ierr
!
!    integer :: np, ne, nt, nb, nci, nc, nr
!    double precision :: const_diff
!    integer, dimension(3,ne) :: e
!    integer, dimension(3,nt) :: t
!    double precision, dimension(2,np) :: p
!    double precision, dimension(2,2) :: dbounds
!
!    ni = np-nb
!    nc = np-nci
!    nr = nb-nci
!
!    nz = 12  !! number of non-zeros per row for matrix pre-allocation
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, nc, nc, nz, PETSC_NULL_INTEGER, PAmat, ierr)
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, nci, nci, nz, PETSC_NULL_INTEGER, PAcc, ierr)
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, ni, nci, nz, PETSC_NULL_INTEGER, PAic, ierr)
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, nci, ni, nz, PETSC_NULL_INTEGER, PAci, ierr)
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, ni, nr, nz, PETSC_NULL_INTEGER, PAir, ierr)
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, nr, ni, nz, PETSC_NULL_INTEGER, PAri, ierr)
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, nr, nci, nz, PETSC_NULL_INTEGER, PArc, ierr)
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, nci, nr, nz, PETSC_NULL_INTEGER, PAcr, ierr)
!
!    call PETScTwoLevelAssemADVmatrix(PAmat,PAcc,PAic,PAci,PAir,PAri,PArc,PAcr,p,e,t,np,ne,nt,nb,nci,1,0.0d0, &
!    [0.0d0,0.0d0,0.0d0,0.0d0],[0.0d0,0.0d0,0.0d0,0.0d0],0.0d0,[0.0d0,0.0d0,0.0d0],dbounds,0,0)
!    call PETScTwoLevelAssemDIFmatrix(PAmat,PAcc,PAic,PAci,PAir,PAri,PArc,PAcr,p,e,t,np,ne,nt,nb,nci,2,const_diff, &
!    [0.0d0,0.0d0,0.0d0,0.0d0],[0.0d0,0.0d0,0.0d0,0.0d0],0.0d0,[0.0d0,0.0d0,0.0d0],dbounds,1,1)
!
!    call MatAssemblyBegin(PAmat,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyBegin(PAcc,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyBegin(PAic,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyBegin(PAci,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyBegin(PAir,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyBegin(PAri,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyBegin(PArc,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyBegin(PAcr,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyEnd(PAmat,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyEnd(PAcc,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyEnd(PAic,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyEnd(PAci,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyEnd(PAir,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyEnd(PAri,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyEnd(PArc,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyEnd(PAcr,MAT_FINAL_ASSEMBLY,ierr)
!
!END SUBROUTINE GetMatSeqTwoLevel
!
!
!!!!*********************************************************
!!!! Subroutine: To construc PETSc Vector in Parallel MPI Form
!SUBROUTINE GetVecMPI(n,fvec,PetVec,Solvec,ierr)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!
!    Vec              PetVec,SolVec
!    PetscInt         i,ii
!    PetscInt         n,one
!    PetscScalar      ivec
!    PetscErrorCode   ierr
!
!    double precision, dimension(n) :: fvec
!
!    call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,n,PetVec,ierr);
!
!    call VecDuplicate(PetVec,SolVec,ierr)
!
!    one = 1
!    do i = 1,n
!    ii = i-1
!    ivec = fvec(i)
!    call VecSetValues(PetVec,one,ii,ivec,INSERT_VALUES,ierr)
!    end do
!
!    call VecAssemblyBegin(PetVec,ierr)
!    call VecAssemblyEnd(PetVec,ierr)
!
!END SUBROUTINE GetVecMPI
!
!
!!!*********************************************************
!!! Subroutine: To construc PETSc Matrix in Sequential Form
!!! Not working : Dec 10,2015
!SUBROUTINE GetMatMPI(n,Amat,PetMat,ierr)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!#include <petsc/finclude/petscmat.h>
!
!
!    Mat              PetMat
!    PetscInt         i,j,ii,jj
!    PetscInt         n,one
!    PetscInt         Istart,Iend
!    PetscScalar      imat
!    PetscErrorCode   ierr
!
!    double precision, dimension(n,n) :: Amat
!
!    call MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,&
!     & n,n,0,PETSC_NULL_INTEGER,0,PETSC_NULL_INTEGER,PetMat,ierr);
!
!    call MatGetOwnershipRange(PetMat,Istart,Iend,ierr)
!
!    one = 1
!    do i = Istart,Iend-1
!    do j = Istart,Iend-1
!        if (Amat(i,j) .ne. 0) then
!            ii = i
!            jj = j
!            imat = Amat(i,j)
!            call MatSetValues(PetMat,one,ii,one,jj,imat,INSERT_VALUES,ierr)
!        end if
!    end do
!    end do
!
!    call MatAssemblyBegin(PetMat,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyEnd(PetMat,MAT_FINAL_ASSEMBLY,ierr)
!
!END SUBROUTINE GetMatMPI
!
!!*********************************************************
!! Subroutine : Call PETSc Matrix, Vector and Ksp
!! and solve the linear problem of the form Ax = b
!SUBROUTINE PETScSolver(n, pg, eg, tg, npg, neg, ntg, amp, dbounds, meanc, temp1, temp2)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!#include <petsc/finclude/petscmat.h>
!#include <petsc/finclude/petscksp.h>
!#include <petsc/finclude/petscpc.h>
!
!    Vec              PetVec,SolVec
!    Mat              PetMat
!    KSP              ksp
!    PC               pc
!
!    PetscInt         rank,size
!    PetscBool        flg
!    PetscReal        tol
!    PetscInt         n
!    PetscErrorCode   ierr
!
!    double precision, dimension(n)   :: temp2
!    integer, dimension(n)            :: temp1
!
!    integer :: npg, neg, ntg
!    double precision :: amp, meanc
!    integer, dimension(3,neg) :: eg
!    integer, dimension(3,ntg) :: tg
!    double precision, dimension(2,npg) :: pg
!    double precision, dimension(2,2) :: dbounds
!
!!!----------------------------------------------------------------
!!! PETSc-Initialize
!    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
!    call MPI_Comm_size(PETSC_COMM_WORLD,size,ierr)
!    call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
!
!!! PETSc-Hello World
!    WRITE(*,*) 'Hello from PETSc World, rank:',rank, ' of total',size,' processes.'
!
!!!----------------------------------------------------------------
!!!PETSc-Assembled-Vector
!    call GetPetVec(n,PetVec,SolVec,pg,eg,tg,npg,neg,ntg,amp,dbounds,ierr)
!
!!!PETSc-Assembled-Matrix
!    call GetPetMat(n,PetMat, pg, eg, tg, npg, neg, ntg, amp, dbounds, meanc, ierr)
!
!!!PETSc-KSP
!    call PETScKSP(ksp,pc,PetMat,PetVec,Solvec,ierr)
!
!!!----------------------------------------------------------------
!!!PETSc-Output to Fortran Output
!    call VecGetValues(SolVec, n, temp1, temp2, ierr)
!
!!!----------------------------------------------------------------
!!! PETSc-Finalize
!    call VecDestroy(PetVec,ierr)
!    call VecDestroy(SolVec,ierr)
!    Call MatDestroy(PetMat,ierr)
!    Call KspDestroy(ksp,ierr)
!
!    call PetscFinalize(ierr)
!
!END SUBROUTINE PETScSolver
!
!!!*********************************************************
!!! Subroutine : Call PETSc Matrix, Vector and Ksp
!!! and solve the linear problem of the form Ax = b
!SUBROUTINE PETScSolverPreAssembled(n,fvec,Amat,temp1,temp2)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!#include <petsc/finclude/petscmat.h>
!#include <petsc/finclude/petscksp.h>
!#include <petsc/finclude/petscpc.h>
!
!Vec              PetVec,SolVec
!Mat              PetMat
!KSP              ksp
!PC               pc
!
!PetscInt         rank,size
!PetscBool        flg
!PetscReal        tol
!PetscInt         n
!PetscErrorCode   ierr
!
!double precision, dimension(n,n) :: Amat
!double precision, dimension(n)   :: fvec
!double precision, dimension(n)   :: U
!double precision, dimension(n)   :: temp2
!integer, dimension(n)            :: temp1
!
!    !!----------------------------------------------------------------
!    !! PETSc-Initialize
!    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
!    call MPI_Comm_size(PETSC_COMM_WORLD,size,ierr)
!    call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
!
!    !! PETSc-Hello World
!    WRITE(*,*) 'Hello from PETSc World, rank:',rank, ' of total',size,' processes.'
!
!    !!----------------------------------------------------------------
!    !!PETSc-Vector
!    call GetVec(n,fvec,PetVec,Solvec,ierr)      !! Generalize Form
!    !call GetVecSeq(n,fvec,PetVec,Solvec,ierr)   !! Sequential Form
!    !call GetVecMPI(n,fvec,PetVec,Solvec,ierr)   !! ParallelMPI Form
!
!    !!PETSc-Matrix
!    call GetMat(n,Amat,PetMat,ierr)             !! Generalize Form
!    !call GetMatSeq(n,Amat,PetMat,ierr)          !! Sequential Form
!    !call GetMatMPI(n,Amat,PetMat,ierr)          !! ParallelMPI Form : ERROR
!
!    !!PETSc-KSP
!    call PETScKSP(ksp,pc,PetMat,PetVec,Solvec,ierr)
!
!    !!----------------------------------------------------------------
!    !!PETSc-Output to Fortran Output
!    call VecGetValues(SolVec, n, temp1, temp2, ierr)
!
!    !!----------------------------------------------------------------
!    !! PETSc-Finalize
!    call VecDestroy(PetVec,ierr)
!    call VecDestroy(SolVec,ierr)
!    Call MatDestroy(PetMat,ierr)
!    Call KspDestroy(ksp,ierr)
!
!call PetscFinalize(ierr)
!
!END SUBROUTINE PETScSolverPreAssembled




END MODULE PETScommon
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%** END **%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!!*********************************************************
!!! Subroutine :
!SUBROUTINE PETScMatrix(PetMat,pg,eg,tg,npg,neg,ntg,dbounds,const_diff,ierr)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!#include <petsc/finclude/petscmat.h>
!
!Mat              PetMat
!PetscInt         n
!PetscErrorCode   ierr
!
!double precision :: const_diff
!double precision, dimension(2,2) :: dbounds
!double precision, dimension(2,npg) :: pg
!integer, dimension(3,ntg) :: tg
!integer, dimension(3,neg) :: eg
!integer :: npg, neg, ntg
!
!n=npg
!call MatCreate(PETSC_COMM_WORLD,PetMat,ierr)
!call MatSetFromOptions(PetMat,ierr)
!call MatSetSizes(PetMat,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr)
!call MatSetUp(PetMat,ierr)
!
!call PETScAssembleADVmatrix(PetMat,pg,eg,tg,npg,neg,ntg,1,0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],&
![0.0d0,0.0d0,0.0d0,0.0d0],0.0d0,[0.0d0,0.0d0,0.0d0],dbounds,0,0)
!call PETScAssembleDIFmatrix(PetMat,pg,eg,tg,npg,neg,ntg,2,const_diff,[0.0d0,0.0d0,0.0d0,0.0d0],&
![0.0d0,0.0d0,0.0d0,0.0d0], 0.0d0,[0.0d0,0.0d0,0.0d0],dbounds,1,1)
!
!call MatAssemblyBegin(PetMat,MAT_FINAL_ASSEMBLY,ierr)
!call MatAssemblyEnd(PetMat,MAT_FINAL_ASSEMBLY,ierr)
!
!END SUBROUTINE PETScMatrix


