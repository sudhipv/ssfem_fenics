!! myCommon Module: By Ajit Desai, 8/2015
!! purpose: it contains various subrountines in addition to common.f90 used in this package 
!! Mainly written to solve 2-level problems with Coarse Grid for NNC & FETIDP
!! if you want to include more subroutine, this is the place to add
!!-----------------------------------------------------------------------------------

MODULE myCommon

USE common

IMPLICIT NONE

CONTAINS

!!**************************************************************************
!! Subroutine to constrct Block-diagonal-weighting matrix 'D_s' 
!! The diagonal entries of the Det-Weighting-Mat are equal to the reciprocal
!! of the number of subdomains which contain the same interface node 
SUBROUTINE getBDM(pid,nb,npceout,Dmat)
    
    CHARACTER(len=255) :: extension, str1
    
    INTEGER :: npg,neg,ntg,nbg,ndom,np,ne,nt,nb,npceout
    INTEGER :: i, j, k, m, Ri, Rj, Ri2, Rj2, pid
    INTEGER, DIMENSION(:), ALLOCATABLE :: GlobalBN, bnodes12, bnodes13
    
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: CountR
    DOUBLE PRECISION, DIMENSION(:,:) :: Dmat
    
    
    open(unit=1,file='meshdim.dat',status='old')
    read(unit=1,fmt='(5I8)') npg, neg, ntg, nbg, ndom
    close(1)
    
    ALLOCATE(GlobalBN(nbg),CountR(nbg))
    open(unit=1,file='boundary_nodes.dat',status='old')
    READ(1,*) GlobalBN
    close(1)
    
    CountR(:) = 0.0d0
    DO k = 1,ndom
    CALL int2str(extension,k,4)
    extension = trim(extension)
    CALL readmeshdim(np,ne,nt,nb,extension)
    ALLOCATE(bnodes12(nb))
    
    CALL int2str(str1,k,4)
    str1 = 'nbnodes' // trim(str1) // '.dat'  ! changing from nbnodes to nodes wont work
    OPEN(unit=1,file=str1,status='old')       ! Need to confirm this change
    READ(1,*) bnodes12                        ! Date Feb 1, 2014
    CLOSE(1)
    
    DO i = 1,nbg
       Ri = GlobalBN(i)
       DO j = 1,nb
          Rj = bnodes12(j)
          IF (Ri == Rj) THEN
             CountR(i) = CountR(i) + 1
          END IF
       END DO
    END DO
    
    DEALLOCATE(bnodes12)
    END DO
    
    CALL int2str(extension,pid+1,4)
    extension = trim(extension)
    CALL readmeshdim(np,ne,nt,nb,extension)
    ALLOCATE(bnodes13(nb))
    
    CALL int2str(str1,pid+1,4)
    str1 = 'nbnodes' // trim(str1) // '.dat'
    OPEN(unit=1,file=str1,status='old')
    READ(1,*) bnodes13
    CLOSE(1)
    
    Dmat(:,:) = 0.0d0
    
    Do i = 1,nb
       Ri2 = bnodes13(i)
       DO j = 1,nbg
          Rj2 = GlobalBN(j)
          IF (Ri2 == Rj2) THEN
             Dmat(i,i) = 1/CountR(j)
          END IF
       END DO
    END DO
    
    DO m = 2,npceout
        Dmat((m-1)*nb+1:m*nb,(m-1)*nb+1:m*nb) = Dmat(1:nb,1:nb)
    END DO
    
    DEALLOCATE(bnodes13, GlobalBN, CountR)
    
END SUBROUTINE getBDM
    
!!*****************************************************************
!! Subroutine to constrct Boolean operators 'Rc_S' that extract the 
!! local corner 'Uc_s' from local interface unknowns 'Ub_s'
SUBROUTINE getRcS(pid, nb, nci, npceout, RcSmat)
    
    CHARACTER(len=255) :: str1 !,extension, filename, 
    
    INTEGER :: nb,nci,npceout
    INTEGER :: i, j, m, Ri, Rj, pid
    INTEGER, DIMENSION(:), ALLOCATABLE :: cnode12, bnodes12
    DOUBLE PRECISION, DIMENSION(:,:) :: RcSmat
    
    ALLOCATE(cnode12(nci), bnodes12(nb))
    
    call int2str(str1,pid+1,4)
    str1 = 'cnodes' // trim(str1) // '.dat'
    open(unit=1,file=str1,status='old')
    READ(1,*) cnode12
    close(1)
    
    CALL int2str(str1,pid+1,4)
    str1 = 'nbnodes' // trim(str1) // '.dat'
    open(unit=2,file=str1,status='old')
    READ(2,*) bnodes12
    close(2)
    
    RcSmat(:,:) = 0.0d0
    
    Do i = 1,nci
       Ri = cnode12(i)
       DO j = 1,nb
          Rj = bnodes12(j)
          IF (Ri == Rj) THEN
             RcSmat(i,j) = 1
          END IF
       END DO
    END DO
    
    DO m = 2,npceout
        RcSmat((m-1)*nci+1:m*nci,(m-1)*nb+1:m*nb) = RcSmat(1:nci,1:nb)
    END DO
    
    DEALLOCATE(cnode12, bnodes12)
    
END SUBROUTINE getRcS
    
!!*****************************************************************
!! Subroutine to constrct Boolean operators 'Rr_s' that extract the 
!! local interface unknowns 'Ub_s' into corner 'Ur_s' 
SUBROUTINE getRrS(pid, nb, nri, npceout, RrSmat)
    
    CHARACTER(len=255) :: str1 !,extension, filename, 
    
    INTEGER :: nb,nri,npceout
    INTEGER :: i, j, m, Ri, Rj, pid
    INTEGER, DIMENSION(:), ALLOCATABLE :: rnode12, bnodes12
    DOUBLE PRECISION, DIMENSION(:,:) :: RrSmat
    
    ALLOCATE(rnode12(nri), bnodes12(nb))
    
    call int2str(str1,pid+1,4)
    str1 = 'rnodes' // trim(str1) // '.dat'
    open(unit=1,file=str1,status='old')
    READ(1,*) rnode12
    close(1)
    
    CALL int2str(str1,pid+1,4)
    str1 = 'nbnodes' // trim(str1) // '.dat'
    open(unit=2,file=str1,status='old')
    READ(2,*) bnodes12
    close(2)
    
    RrSmat(:,:) = 0.0d0
    
    Do i = 1,nri
       Ri = rnode12(i)
       DO j = 1,nb
          Rj = bnodes12(j)
          IF (Ri == Rj) THEN
             RrSmat(i,j) = 1
          END IF
       END DO
    END DO
    
    DO m = 2,npceout
        RrSmat((m-1)*nri+1:m*nri,(m-1)*nb+1:m*nb) = RrSmat(1:nri,1:nb)
    END DO
    
    DEALLOCATE(rnode12, bnodes12)
    
END SUBROUTINE getRrS
    
!!****************************************************************
!! Subroutine to constrct Block-diagonal-restriction matrix 'Bc_s' 
!! restriction operator that maps the global corner unknown 'Uc' 
!! into its local 'Uc_s' components
SUBROUTINE getBcS(pid, nci, nbgc, npceout, BcSmat)
    
    CHARACTER(len=255) :: str1 !,extension, filename
    
    INTEGER :: nci,nbgc,npceout
    INTEGER :: i, j, m, Ri, Rj, pid
    INTEGER, DIMENSION(:), ALLOCATABLE :: cnode12, GlobalCN
    
    DOUBLE PRECISION, DIMENSION(:,:) :: BcSmat
    
    ALLOCATE(GlobalCN(nbgc))
    open(unit=2,file='corner_nodes.dat',status='old')
    read(2,*) GlobalCN
    close(2)
    
    ALLOCATE(cnode12(nci))
    call int2str(str1,pid+1,4)
    str1 = 'cnodes' // trim(str1) // '.dat'
    open(unit=2,file=str1,status='old')
    READ(2,*) cnode12
    close(2)
    
    BcSmat(:,:) = 0.0d0
    
    DO i = 1,nci
       Ri = cnode12(i)
       DO j = 1,nbgc
          Rj = GlobalCN(j)
          IF (Ri == Rj) THEN
              BcSmat(i,j) = 1
          END IF
       END DO
    END DO
    
    DO m = 2,npceout
        BcSmat((m-1)*nci+1:m*nci,(m-1)*nbgc+1:m*nbgc) = BcSmat(1:nci,1:nbgc)
    END DO
    
    DEALLOCATE(cnode12, GlobalCN)
    
END SUBROUTINE getBcS
    
!!****************************************************************
!! subroutine to construct block-diagonal-jump-operator 'Br_s'
!! such that : sum_1^ns(Br_s ur_s) = 0 

SUBROUTINE getBrS(pid, nri, nbgr, npceout, BrSmat)

    CHARACTER(len=255) :: str1             !!extension, filename

    INTEGER :: nri, nbgr, npceout
    INTEGER :: i, j, m, Ri, Rj, Rc, pid
    INTEGER, DIMENSION(:), ALLOCATABLE :: rnode12, GlobalRN, BrSmatCount2
    DOUBLE PRECISION, DIMENSION(:,:) :: BrSmat

    ALLOCATE(BrSmatCount2(nbgr),GlobalRN(nbgr))
    open(unit=2,file='remaining_nodes.dat',status='old')
    read(2,*) GlobalRN
    close(2)

    open(unit=3,file='final_count.dat',status='old')
    read(3,*) BrSmatCount2
    close(3)

    ALLOCATE(rnode12(nri))
    call int2str(str1,pid+1,4)
    str1 = 'rnodes' // trim(str1) // '.dat'
    open(unit=2,file=str1,status='old')
    READ(2,*) rnode12
    close(2)

    BrSmat(:,:) = 0.0d0

    Do i = 1,nri
        Ri = rnode12(i)
        DO j = 1,nbgr
            Rj = GlobalRN(j)
            Rc = BrSmatCount2(j)
            IF ((Ri == Rj) .AND. (Rc == (pid+1))) THEN
                BrSmat(i,j) = 1
                ELSEIF ((Ri == Rj) .AND. (Rc /= (pid+1))) THEN
                BrSmat(i,j) = -1
            ELSE
                BrSmat(i,j) = 0.0d0
            END IF
        END DO
    END DO

    DO m = 2,npceout
    BrSmat((m-1)*nri+1:m*nri,(m-1)*nbgr+1:m*nbgr) = BrSmat(1:nri,1:nbgr)
    END DO
    DEALLOCATE(rnode12, GlobalRN, BrSmatCount2)

END SUBROUTINE getBrS

!!*****************************************************************
!! subroutine to construct a block-diagonal-averaging-matrix 'Dr_s'
SUBROUTINE getDrS(pid, nbgr, nri, npceout, DrSmat)

    CHARACTER(len=255) :: extension, str1

    INTEGER :: npg,neg,ntg,nbg,ndom,nbgr,nri,nci2,nri2,npceout  !!,np,ne,nt,nb
    INTEGER :: i, j, k, m, Ri, Rj, Ri2, Rj2, pid
    INTEGER, DIMENSION(:), ALLOCATABLE :: GlobalRN, rnode12, rnode13

    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: CountR2
    DOUBLE PRECISION, DIMENSION(:,:) :: DrSmat


    open(unit=1,file='meshdim.dat',status='old')
    read(unit=1,fmt='(5I8)') npg, neg, ntg, nbg, ndom
    close(1)

    ALLOCATE(GlobalRN(nbgr),CountR2(nbgr))
    open(unit=1,file='remaining_nodes.dat',status='old')
    READ(1,*) GlobalRN
    close(1)

    CountR2(:) = 0.0d0
    DO k = 1,ndom
        CALL int2str(extension,k,4)
        extension = trim(extension)
        call readcrdim(nci2,nri2,extension)
        ALLOCATE(rnode12(nri2))

        call int2str(str1,k,4)
        str1 = 'rnodes' // trim(str1) // '.dat'
        open(unit=2,file=str1,status='old')
        READ(2,*) rnode12
        close(2)

        DO i = 1,nbgr
            Ri = GlobalRN(i)
            DO j = 1,nri2
                Rj = rnode12(j)
                IF (Ri == Rj) THEN
                    CountR2(i) = CountR2(i) + 1
                END IF
            END DO 
        END DO

        DEALLOCATE(rnode12)
    END DO

    ALLOCATE(rnode13(nri))
    call int2str(str1,pid+1,4)
    str1 = 'rnodes' // trim(str1) // '.dat'
    open(unit=2,file=str1,status='old')
    READ(2,*) rnode13
    close(2)

    DrSmat(:,:) = 0.0d0

    Do i = 1,nri
        Ri2 = rnode13(i)
        DO j = 1,nbgr
        Rj2 = GlobalRN(j)
            IF (Ri2 == Rj2) THEN
                DrSmat(i,i) = 1/CountR2(j);
            END IF
        END DO
    END DO

    DO m = 2,npceout
       DrSmat((m-1)*nri+1:m*nri,(m-1)*nri+1:m*nri) = DrSmat(1:nri,1:nri)
    END DO

    DEALLOCATE(rnode13, GlobalRN, CountR2)

END SUBROUTINE getDrS
!!*********************************************************
!! Subroutine to read 2nd level local dimensions 
!! nri : number of remaining nodes for each partition
!! nci : number of corner nodes for each partition
SUBROUTINE readcrdim(nci,nri,extension)
    
    integer :: nci, nri
    character(len=255) :: filename,extension
    
    filename = 'dimcrn' // trim(extension) // '.dat'
    open(unit=1,file=filename,status='old')
    read(1,*) nci, nri
    close(1)
    
END SUBROUTINE readcrdim
    
!!*********************************************************
!! Subroutine to read 2nd level global dimensions 
!! nbgr : number of remaining nodes for whole mesh
!! nbgc : number of corner nodes for whole mesh
SUBROUTINE readgcrdim(nbgc,nbgr)
    integer :: nbgc, nbgr

    open(unit=1,file='nbgc.dat',status='old')
    READ(1,*) nbgc
    close(1)

    open(unit=2,file='nbgr.dat',status='old')
    READ(2,*) nbgr
    close(2)
    
END SUBROUTINE readgcrdim


!!*********************************************************
!! Subroutine containes 4-dimensional, 3rd-order PCE basis 
!! Try to atomate it by using Matlab UQTK like approach
SUBROUTINE evaluate_psi_sp4(ndim,term,xi,psi)

    double precision :: psi
    integer :: term,ndim
    double precision, dimension(ndim) :: xi

    SELECT CASE (term)
        CASE (1)                         !! 0th-Order (p=0)
        psi = 1.0d0
        CASE (2)                         !! 1st-Order (p=1)
        psi = xi(1)
        CASE (3)
        psi = xi(2)
        CASE (4)
        psi = xi(3)
        CASE (5)
        psi = xi(4)
        CASE (6)                           !! 2nd-Order (p=2)
        psi = (xi(1)**2 - 1.0d0)/sqrt(2.0d0)
        CASE (7)
        psi = xi(1)*xi(2)   
        CASE (8)
        psi = xi(1)*xi(3)
        CASE (9)
        psi = xi(1)*xi(4)   
        CASE (10)   
        psi = (xi(2)**2 - 1.0d0)/sqrt(2.0d0)
        CASE (11)
        psi = xi(2)*xi(3)
        CASE (12) 
        psi = xi(2)*xi(4)
        CASE (13)
        psi = (xi(3)**2 - 1.0d0)/sqrt(2.0d0)   
        CASE (14)
        psi = xi(3)*xi(4)
        CASE (15)
        psi = (xi(4)**2 - 1.0d0)/sqrt(2.0d0)
        CASE (16)                          !!3rd-Order (p=3)
        psi = (xi(1)**3 - 3.0d0*xi(1))/sqrt(6.0d0)
        CASE (17)
        psi = (xi(1)**2*xi(2) - xi(2))/sqrt(2.0d0)
        CASE (18)
        psi = (xi(1)**2*xi(3) - xi(3))/sqrt(2.0d0)
        CASE (19)
        psi = (xi(1)**2*xi(4) - xi(4))/sqrt(2.0d0)
        CASE (20)
        psi = (xi(2)**2*xi(1) - xi(1))/sqrt(2.0d0)
        CASE (21)
        psi = xi(1)*xi(2)*xi(3)
        CASE (22)
        psi = xi(1)*xi(2)*xi(4)
        CASE (23)
        psi = (xi(3)**2*xi(1) - xi(1))/sqrt(2.0d0)
        CASE (24)
        psi = xi(1)*xi(3)*xi(4)
        CASE (25)
        psi = (xi(1)*xi(4)**2 - xi(4))/sqrt(2.0d0)
        CASE (26)
        psi = (xi(2)**3 - 3.0d0*xi(2))/sqrt(6.0d0)
        CASE (27)
        psi = (xi(2)**2*xi(3) - xi(3))/sqrt(2.0d0)
        CASE (28)
        psi = (xi(2)**2*xi(4) - xi(4))/sqrt(2.0d0)
        CASE (29)
        psi = (xi(3)**2*xi(2) - xi(2))/sqrt(2.0d0)
        CASE (30)
        psi = xi(2)*xi(3)*xi(4)
        CASE (31)
        psi = (xi(4)**2*xi(2) - xi(2))/sqrt(2.0d0)
        CASE (32)
        psi = (xi(3)**3 - 3.0d0*xi(3))/sqrt(6.0d0)
        CASE (33)
        psi = (xi(3)**2*xi(4) - xi(4))/sqrt(2.0d0)
        CASE (34)
        psi = (xi(4)**2*xi(3) - xi(3))/sqrt(2.0d0)
        CASE (35)
        psi = (xi(4)**3 - 3.0d0*xi(4))/sqrt(6.0d0)
    END SELECT

END SUBROUTINE evaluate_psi_sp4



SUBROUTINE getlmindex(lmindex,ndim)

    integer, dimension(2,ndim) :: lmindex
    integer :: ndim

    open(unit=1,file='lmindex',status='old')
    READ(1,*) lmindex
    close(1)

END SUBROUTINE


!!*********************************************************
!! Subroutine for "PCKF": used to construc actual solution
SUBROUTINE create_d_pce_sp4(nmeas,np,ndim,extension,npceout,measnodes,Af,d,gamma,xi)

    double precision, dimension(np,npceout+nmeas) :: Af
    double precision, dimension(nmeas) :: d
    double precision, dimension(nmeas) :: error

    integer, dimension(nmeas) :: measnodes
    integer :: nmeas,np,i,npceout,ndim
    double precision :: gamma,psi
    double precision, dimension(ndim) :: xi

    character(len=255) :: str1,extension

    d(:) = 0.0d0
    do i = 1,npceout
    call evaluate_psi_sp4(ndim,i,xi,psi)
    d = d + Af(measnodes,i)*psi
    end do

    call random_number(error)
    d = d + sqrt(gamma)*error

    str1 = 'meas' // trim(extension) // '.dat'
    open(unit=1,file=str1,status='replace')
    do i = 1,nmeas
    write(unit=1,fmt='(1ES17.8E3)') d(i)
    end do
    close(1)

END SUBROUTINE create_d_pce_sp4

!!**********************************************************************
!! Subroutine for "PCKF": used to construc actual interior solution
SUBROUTINE create_realization_sp_i4(np,nb,ndim,npceout,nmeas,Af,uveci,xi)

    double precision, dimension(np,npceout+nmeas) :: Af
    double precision, dimension(np-nb) :: uveci
    double precision, dimension(ndim) :: xi
    double precision :: psi

    integer :: np,npceout,nmeas,nb,i,ndim

    uveci(:) = 0.0d0
    do i = 1,npceout
    call evaluate_psi_sp(i,xi,psi)
    uveci = uveci + Af(1:(np-nb),i)*psi
    end do

END SUBROUTINE create_realization_sp_i4

!!*********************************************************************
!! Subroutine for "PCKF": used to construc actual boundary solution
SUBROUTINE create_realization_sp_b4(nbg,ndim,npceout,nmeas,Afb,uvecb,xi)

    double precision, dimension(nbg,npceout+nmeas) :: Afb
    double precision, dimension(nbg) :: uvecb
    double precision, dimension(ndim) :: xi
    double precision :: psi

    integer :: nbg,npceout,nmeas,i,ndim

    uvecb(:) = 0.0d0
    do i = 1,npceout
    call evaluate_psi_sp(i,xi,psi)
    uvecb = uvecb + Afb(:,i)*psi
    end do

END SUBROUTINE create_realization_sp_b4

!!!------------------------------------------------------------------------
!!! Subroutines for FETI-DP Post-Processing 
!!!------------------------------------------------------------------------
!!*********************************************************
!! Subroutine To gather the Boundary node solutions, Ub = Ur + Uc
subroutine create_Ub_FETIddm(pid,nb,npceout,nri,nci,UrS,UcS,Ub)

    integer :: i, pid
    integer :: nb, nri, nci, npceout
    double precision, dimension(nb*npceout) :: Ub
    double precision, dimension(nri*npceout) :: UrS
    double precision, dimension(nci*npceout) :: UcS

    Ub(:) = 0.0d0

    do i = 1,npceout
    Ub((i-1)*nb+1:(((i-1)*nb)+nri)) = UrS((i-1)*nri+1:i*nri)
    Ub((((i-1)*nb)+nri)+1:i*nb)     = UcS((i-1)*nci+1:i*nci) 
    End Do

end subroutine create_Ub_FETIddm

!!***********************************************************
!! Subroutine to construct matrix of solution vectors Ui & Ub
subroutine create_Af_FETIddm(nmeasg,npceout,np,nb,Ui,Ub,Af)

    double precision, dimension(np,npceout+nmeasg) :: Af
    double precision, dimension((np-nb)*npceout) :: Ui
    double precision, dimension(nb*npceout) :: Ub

    integer :: npceout,i,np,nmeasg,nb

    do i = 1,npceout
        Af(1:(np-nb),i) = Ui((i-1)*(np-nb)+1:i*(np-nb))
        Af(((np-nb)+1):np,i) = (0.5d0)*Ub((i-1)*nb+1:i*nb)
    end do
    Af(:,(npceout+1):(npceout+nmeasg)) = 0.0d0

end subroutine create_Af_FETIddm

!!***********************************************************
!! Subroutine to construct matrix of solution vector Uc
subroutine create_Afb_FETIddm(nmeasg,npceout,nbgc,Uc,Afb)

    double precision, dimension(nbgc,npceout+nmeasg) :: Afb
    double precision, dimension(nbgc*npceout) :: Uc

    integer :: npceout,i,nmeasg,nbgc

    do i = 1,npceout
    Afb(:,i) = Uc((i-1)*nbgc+1:i*nbgc)
    end do
    Afb(:,(npceout+1):(npceout+nmeasg)) = 0.0d0

end subroutine create_Afb_FETIddm

!!***********************************************************
!! Subroutine to 
subroutine create_realization_FETIsp_i(np,nci,npceout,nmeas,Af,FETIuveci,xi)

    double precision, dimension(np,npceout+nmeas) :: Af
    double precision, dimension(np-nci) :: FETIuveci
    double precision, dimension(3) :: xi
    double precision :: psi

    integer :: np,npceout,nmeas,nci,i

    FETIuveci(:) = 0.0d0
    do i = 1,npceout
    call evaluate_psi_sp(i,xi,psi)
    FETIuveci = FETIuveci + Af(1:(np-nci),i)*psi
    end do

end subroutine create_realization_FETIsp_i

!!***********************************************************
!! Subroutine to 
subroutine construct_FETIcoln(pid,np,nb,npg,nri,nci,FETIuveci,U_gi)

character(len=255) :: str1
integer :: i, j, pid, np, nb, npg, nri, nci, npi
double precision, dimension(npg) :: U_gi
double precision, dimension(np-nci) :: FETIuveci
INTEGER, DIMENSION(:), ALLOCATABLE :: temp2

    npi = np - nb ! size(UiS)
    ALLOCATE(temp2(np))
    call int2str(str1,pid+1,4)
    str1 = 'nodes' // trim(str1) // '.dat'
    open(unit=1,file=str1,status='old')
    read(1,*) temp2
    close(1)

    U_gi(:) = 0.0d0

    do i = 1,npi
        U_gi(temp2(i)) = FETIuveci(i)
    End Do
    do j = 1,nri
        U_gi(temp2(npi+j)) = FETIuveci(npi+j)
    End Do

end subroutine construct_FETIcoln

!!***********************************************************
!! Subroutine to 
subroutine create_realization_FETIsp_b(nbgc,npceout,nmeas,Afb,FETIuvecb,xi)

    double precision, dimension(nbgc,npceout+nmeas) :: Afb
    double precision, dimension(nbgc) :: FETIuvecb
    double precision, dimension(3) :: xi
    double precision :: psi

    integer :: nbgc,npceout,nmeas,i

    FETIuvecb(:) = 0.0d0
    do i = 1,npceout
    call evaluate_psi_sp(i,xi,psi)
    FETIuvecb = FETIuvecb + Afb(:,i)*psi
    end do

end subroutine create_realization_FETIsp_b

!!***********************************************************
!! Subroutine to 
subroutine construct_FETIcoln_boundary(FETIuvecb,U_g)

character(len=255) :: str1
integer :: nbgc, i, temp1
double precision, dimension(:) :: FETIuvecb, U_g

nbgc = size(FETIuvecb)

    str1 = 'corner_nodes.dat'
    open(unit=1,file=str1,status='old')
    do i = 1,nbgc
    read(unit=1,fmt=*) temp1
    U_g(temp1) = FETIuvecb(i)
    end do
    close(1)

end subroutine construct_FETIcoln_boundary


!!*********************************************************
!! Subroutine
subroutine create_Af_pckfddm3D(nmeasg,npceout,np,nb,Ui,Ub,nVec,Af)

integer :: npceout,i,np,nmeasg,nb,nVec,nid,nbd,npd
double precision, dimension(np*nVec,npceout+nmeasg) :: Af
double precision, dimension((np-nb)*nVec*npceout) :: Ui
double precision, dimension(nb*nVec*npceout) :: Ub

nid = (np-nb)*nVec
nbd = nb*nVec
npd = np*nVec

do i = 1,npceout
!Af(1:(np-nb),i) = Ui((i-1)*(np-nb)+1:i*(np-nb))
!Af((np-nb+1):np,i) = Ub((i-1)*nb+1:i*nb)
Af(1:nid,i) = Ui((i-1)*npd+1:i*nid)
Af((nid+1):npd,i) = Ub((i-1)*npd+1:i*nbd)
end do
Af(:,(npceout+1):(npceout+nmeasg)) = 0.0d0

end subroutine create_Af_pckfddm3D

!!*********************************************************
!! Subroutine
subroutine create_Afb_pckfddm3D(nmeasg,npceout,nbg,Ub_g,nVec,Afb)

double precision, dimension(nbg*nVec,npceout+nmeasg) :: Afb
double precision, dimension(nbg*nVec*npceout) :: Ub_g

integer :: npceout,i,nmeasg,nbg,nbgd,nVec

nbgd = nbg*nVec

do i = 1,npceout
Afb(:,i) = Ub_g((i-1)*nbgd+1:i*nbgd)
end do
Afb(:,(npceout+1):(npceout+nmeasg)) = 0.0d0

end subroutine create_Afb_pckfddm3D


!!*********************************************************
!! Subroutine
subroutine construct_coln3D(pid,np,nb,npg,Ui,U_gi)

character(len=255) :: str1
integer :: npi, i, pid, temp1, np, nb, ni, npg
double precision, dimension(:) :: Ui, U_gi

npi = size(Ui)
ni = np-nb

call int2str(str1,pid+1,4)
str1 = '../../data/meshData/nodes' // trim(str1) // '.dat'
open(unit=1,file=str1,status='old')
U_gi(:) = 0.0d0
do i = 1,ni
read(unit=1,fmt=*) temp1
U_gi(temp1) = Ui(i)
U_gi(temp1+npg) = Ui(i+ni)
U_gi(temp1+(2*npg)) = Ui(i+(2*ni))
end do
close(1)

end subroutine construct_coln3D


!!*********************************************************
!! Subroutine
subroutine construct_coln_boundary3D(nbg,npg,Ub_g,U_g)

character(len=255) :: str1
integer :: nbg, nbgi, i, temp1, npg
double precision, dimension(:) :: Ub_g, U_g

nbgi = size(Ub_g)

str1 = '../../data/meshData/boundary_nodes.dat'
open(unit=1,file=str1,status='old')
do i = 1,nbg
read(unit=1,fmt=*) temp1
U_g(temp1) = Ub_g(i)
U_g(temp1+npg) = Ub_g(i+nbg)
U_g(temp1+(2*npg)) = Ub_g(i+(2*nbg))
end do
close(1)

end subroutine construct_coln_boundary3D




!!%%%%%%%%%%%%%%%%%%%  5-Dim, 3rd-Order %%%%%%%%%%%%%%%%%%%%
!!!*********************************************************
!!! Subroutine containes 5-dimensional, 3rd-order PCE basis
!!! Try to atomate it by using Matlab UQTK like approach
!SUBROUTINE evaluate_psi_sp5(ndim,term,xi,psi)
!
!    double precision :: psi
!    integer :: term,ndim
!    double precision, dimension(ndim) :: xi
!
!    SELECT CASE (term)
!        CASE (1)                         !! 0th-Order (p=0)
!        psi = 1.0d0
!        CASE (2)                         !! 1st-Order (p=1)
!        psi = xi(1)
!        CASE (3)
!        psi = xi(2)
!        CASE (4)
!        psi = xi(3)
!        CASE (5)
!        psi = xi(4)
!        CASE (6)
!        psi = xi(5)
!        CASE (7)                           !! 2nd-Order (p=2)
!        psi = (xi(1)**2 - 1.0d0)/sqrt(2.0d0)
!        CASE (8)
!        psi = xi(1)*xi(2)
!        CASE (9)
!        psi = xi(1)*xi(3)
!        CASE (10)
!        psi = xi(1)*xi(4)
!        CASE (11)
!        psi = xi(1)*xi(5)
!        CASE (12)
!        psi = (xi(2)**2 - 1.0d0)/sqrt(2.0d0)
!        CASE (13)
!        psi = xi(2)*xi(3)
!        CASE (14)
!        psi = xi(2)*xi(4)
!        CASE (15)
!        psi = xi(2)*xi(5)
!        CASE (16)
!        psi = (xi(3)**2 - 1.0d0)/sqrt(2.0d0)
!        CASE (17)
!        psi = xi(3)*xi(4)
!        CASE (18)
!        psi = xi(3)*xi(5)
!        CASE (19)
!        psi = (xi(4)**2 - 1.0d0)/sqrt(2.0d0)
!        CASE (20)
!        psi = xi(4)*xi(5)
!        CASE (21)
!        psi = (xi(5)**2 - 1.0d0)/sqrt(2.0d0)
!        CASE (22)                          !!3rd-Order (p=3)
!        psi = (xi(1)**3 - 3.0d0*xi(1))/sqrt(6.0d0)
!        CASE (17)
!        psi = (xi(1)**2*xi(2) - xi(2))/sqrt(2.0d0)
!        CASE (18)
!        psi = (xi(1)**2*xi(3) - xi(3))/sqrt(2.0d0)
!        CASE (19)
!        psi = (xi(1)**2*xi(4) - xi(4))/sqrt(2.0d0)
!        CASE (20)
!        psi = (xi(2)**2*xi(1) - xi(1))/sqrt(2.0d0)
!        CASE (21)
!        psi = xi(1)*xi(2)*xi(3)
!        CASE (22)
!        psi = xi(1)*xi(2)*xi(4)
!        CASE (23)
!        psi = (xi(3)**2*xi(1) - xi(1))/sqrt(2.0d0)
!        CASE (24)
!        psi = xi(1)*xi(3)*xi(4)
!        CASE (25)
!        psi = (xi(1)*xi(4)**2 - xi(4))/sqrt(2.0d0)
!        CASE (26)
!        psi = (xi(2)**3 - 3.0d0*xi(2))/sqrt(6.0d0)
!        CASE (27)
!        psi = (xi(2)**2*xi(3) - xi(3))/sqrt(2.0d0)
!        CASE (28)
!        psi = (xi(2)**2*xi(4) - xi(4))/sqrt(2.0d0)
!        CASE (29)
!        psi = (xi(3)**2*xi(2) - xi(2))/sqrt(2.0d0)
!        CASE (30)
!        psi = xi(2)*xi(3)*xi(4)
!        CASE (31)
!        psi = (xi(4)**2*xi(2) - xi(2))/sqrt(2.0d0)
!        CASE (32)
!        psi = (xi(3)**3 - 3.0d0*xi(3))/sqrt(6.0d0)
!        CASE (33)
!        psi = (xi(3)**2*xi(4) - xi(4))/sqrt(2.0d0)
!        CASE (34)
!        psi = (xi(4)**2*xi(3) - xi(3))/sqrt(2.0d0)
!        CASE (35)
!        psi = (xi(4)**3 - 3.0d0*xi(4))/sqrt(6.0d0)
!    END SELECT
!
!END SUBROUTINE evaluate_psi_sp4

SUBROUTINE findFact(number,factorial)

    INTEGER :: number, n, factorial

    factorial = 1
    n = number
    DO WHILE (n/=0)
        factorial = n*factorial
        n = n-1
    END DO

END SUBROUTINE findFact


SUBROUTINE findFactReal(number,factorial)

    INTEGER :: number, n
    DOUBLE PRECISION :: factorial

    factorial = 1
    n = number
    DO WHILE (n/=0)
        factorial = n*factorial
        n = n-1
    END DO

END SUBROUTINE findFactReal



SUBROUTINE npce(nord, ndim, npceout)

    INTEGER :: nord, ndim, npceout
    INTEGER :: n1, d1, d2

    CALL findFact((nord+ndim),n1)
    CALL findFact((nord),d1)
    CALL findFact((ndim),d2)

    npceout = n1/(d1*d2)

END SUBROUTINE npce

END MODULE myCommon
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%** END **%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! MODIFICATIONS FOR LARGE RV CASE : Date 23 July 2015
!
!subroutine construct_Amat(k,np,nb,npceout,ncijk,ijk,cijk,Apce2,Amat2)
!
!integer :: k, np, nb, npceout, ncijk, indexi
!integer, dimension(ncijk,3) :: ijk
!
!double precision, dimension(ncijk) :: cijk
!double precision, dimension(np,np) :: Apce2
!double precision, dimension(np*npceout,np*npceout) :: Amat2
!
!indexi = 1
!do i = 1,npceout
!    do j = 1,npceout    
!        !do k = 1,npcein
!        if ((i .eq. ijk(indexi,1)) .and. (j .eq. ijk(indexi,2)) .and. (k .eq. ijk(indexi,3))) then
!        Amat2(((i-1)*(np-nb)+1):(i*(np-nb)),((j-1)*(np-nb)+1):(j*(np-nb))) = Amat2(((i-1)*(np-nb)+1):(i*(np-nb)), &
!        ((j-1)*(np-nb)+1):(j*(np-nb))) + cijk(indexi)*Apce2(1:(np-nb),1:(np-nb)) !! Aii
!        Amat2(((i-1)*(np-nb)+1):(i*(np-nb)),((np-nb)*npceout)+((j-1)*nb+1):((np-nb)*npceout)+(j*nb)) = &
!        Amat2(((i-1)*(np-nb)+1):(i*(np-nb)),((np-nb)*npceout)+((j-1)*nb+1):((np-nb)*npceout)+(j*nb)) &
!        + cijk(indexi)*Apce2(1:(np-nb),(np-nb)+1:np)    !! Aigamma
!        Amat2(((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+(i*nb),((j-1)*(np-nb)+1):(j*(np-nb))) = &
!        Amat2(((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+(i*nb),((j-1)*(np-nb)+1):(j*(np-nb))) &
!        + cijk(indexi)*Apce2((np-nb)+1:np,1:(np-nb))    !! Agammai
!        Amat2(((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+(i*nb),((np-nb)*npceout)+((j-1)*nb+1): &
!        ((np-nb)*npceout)+(j*nb)) = Amat2(((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+(i*nb), &
!        ((np-nb)*npceout)+((j-1)*nb+1):((np-nb)*npceout)+(j*nb)) &
!        + cijk(indexi)*Apce2((np-nb)+1:np,(np-nb)+1:np) !! Agammagamma
!        indexi = indexi+1
!        end if
!        !end do
!    end do
!end do
!
!
!
!
!end subroutine construct_Amat


!subroutine readnfloat(nfloat)
!integer :: nfloat

!open(unit=1,file='nfloatingDomain.dat',status='old')
!READ(1,*) nfloat
!close(1)

!end subroutine readnfloat


!subroutine findFloatingSubs(nfloat, floatSubs)
!INTEGER :: nfloat
!INTEGER :: floatSubs(nfloat)

!open(unit=1,file='floatingDomians.dat',status='old')
!READ(1,*) floatSubs
!close(1)

!end subroutine findFloatingSubs

!
!!!!!!!!!!!!!!!!!!!!!!! EXTRA USEFUL STUFF !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!SUBROUTINE getRM(pid, nbg, nb, Rmat2)
!!
!!CHARACTER(len=255) :: extension, filename, str1
!!
!!INTEGER :: nbg, nb
!!INTEGER :: i, j, Ri, Rj, pid
!!INTEGER, DIMENSION(:), ALLOCATABLE :: GlobalBN, bnodes12
!!
!!DOUBLE PRECISION, DIMENSION(:,:) :: Rmat2
!!
!!ALLOCATE(GlobalBN(nbg), bnodes12(nb))
!!OPEN(unit=1,file='boundary_nodes.dat',status='old')
!!READ(1,*) GlobalBN
!!CLOSE(1)
!!
!!CALL int2str(str1,pid+1,4)
!!str1 = 'nbnodes' // trim(str1) // '.dat'
!!OPEN(unit=1,file=str1,status='old')
!!READ(1,*) bnodes12
!!CLOSE(1)
!!
!!Rmat2(:,:) = 0.0d0
!!
!!DO i = 1,nb
!!   Ri = bnodes12(i)
!!   DO j = 1,nbg
!!      Rj = GlobalBN(j)
!!      IF (Ri == Rj) THEN
!!         Rmat2(i,j) = 1
!!      END IF
!!   END DO
!!END DO
!!
!!DEALLOCATE(bnodes12, GlobalBN)
!!
!!END SUBROUTINE getRM
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!SUBROUTINE getSerBDM(id,nb,npceout,Dmat)
!!
!!CHARACTER(len=255) :: extension, filename, str1
!!
!!INTEGER :: npg,neg,ntg,nbg,ndom,np,ne,nt,nb,npceout
!!INTEGER :: i, j, k, m, Ri, Rj, Ri2, Rj2, id
!!INTEGER, DIMENSION(:), ALLOCATABLE :: GlobalBN, bnodes12, bnodes13
!!
!!DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: CountR
!!DOUBLE PRECISION, DIMENSION(:,:) :: Dmat
!!
!!
!!open(unit=1,file='meshdim.dat',status='old')
!!read(unit=1,fmt='(5I8)') npg, neg, ntg, nbg, ndom
!!close(1)
!!
!!ALLOCATE(GlobalBN(nbg),CountR(nbg))
!!open(unit=1,file='boundary_nodes.dat',status='old')
!!READ(1,*) GlobalBN
!!close(1)
!!
!!CountR(:) = 0.0d0
!!DO k = 1,ndom
!!CALL int2str(extension,k,4)
!!extension = trim(extension)
!!CALL readmeshdim(np,ne,nt,nb,extension)
!!ALLOCATE(bnodes12(nb))
!!
!!CALL int2str(str1,k,4)
!!str1 = 'nbnodes' // trim(str1) // '.dat'  ! changing from nbnodes to nodes wont work
!!OPEN(unit=1,file=str1,status='old')       ! Need to confirm this change
!!READ(1,*) bnodes12                        ! Date Feb 1, 2014
!!CLOSE(1)
!!
!!DO i = 1,nbg
!!Ri = GlobalBN(i)
!!DO j = 1,nb
!!Rj = bnodes12(j)
!!IF (Ri == Rj) THEN
!!CountR(i) = CountR(i) + 1
!!END IF
!!END DO
!!END DO
!!
!!DEALLOCATE(bnodes12)
!!END DO
!!
!!CALL int2str(extension,id,4)
!!extension = trim(extension)
!!CALL readmeshdim(np,ne,nt,nb,extension)
!!ALLOCATE(bnodes13(nb))
!!
!!CALL int2str(str1,id,4)
!!str1 = 'nbnodes' // trim(str1) // '.dat'
!!OPEN(unit=1,file=str1,status='old')
!!READ(1,*) bnodes13
!!CLOSE(1)
!!
!!Dmat(:,:) = 0.0d0
!!
!!Do i = 1,nb
!!Ri2 = bnodes13(i)
!!DO j = 1,nbg
!!Rj2 = GlobalBN(j)
!!IF (Ri2 == Rj2) THEN
!!Dmat(i,i) = 1/CountR(j)
!!END IF
!!END DO
!!END DO
!!
!!DO m = 2,npceout
!!Dmat((m-1)*nb+1:m*nb,(m-1)*nb+1:m*nb) = Dmat(1:nb,1:nb)
!!END DO
!!
!!DEALLOCATE(bnodes13, GlobalBN, CountR)
!!
!!END SUBROUTINE getSerBDM
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!subroutine getSerub(id,nb,nbg,npceout,Ub_g,Ub)
!!
!!character(len=255) :: str1
!!integer :: nb, nbg, i, j, k, l, id, temp1, npceout
!!double precision, dimension(nbg*npceout) :: Ub_g
!!double precision, dimension(nb*npceout) :: Ub
!!integer, dimension(nb) :: bnodes
!!
!!call int2str(str1,id,4)
!!str1 = 'bnodes' // trim(str1) // '.dat'
!!open(unit=1,file=str1,status='old')
!!do i = 1,nb
!!read(unit=1,fmt=*) bnodes(i)
!!end do
!!close(1)
!!do k = 1,npceout
!!do i = 1,nb
!!Ub((k-1)*nb+i) = Ub_g((k-1)*nbg+bnodes(i))
!!end do
!!end do
!!
!!end subroutine getSerub
!!
!!subroutine getSerubg(id,nb,nbg,npceout,Ub_g,Ub)
!!
!!character(len=255) :: str1
!!integer :: nb, nbg, i, j, k, l, id, temp1, npceout
!!double precision, dimension(nbg*npceout) :: Ub_g
!!double precision, dimension(nb*npceout) :: Ub
!!integer, dimension(nb) :: bnodes
!!
!!call int2str(str1,id,4)
!!str1 = 'bnodes' // trim(str1) // '.dat'
!!open(unit=1,file=str1,status='old')
!!do i = 1,nb
!!read(unit=1,fmt=*) bnodes(i)
!!end do
!!close(1)
!!Ub_g = 0.0d0
!!do k = 1,npceout
!!do i = 1,nb
!!Ub_g((k-1)*nbg+bnodes(i)) = Ub((k-1)*nb+i)
!!end do
!!end do
!!
!!end subroutine getSerubg
