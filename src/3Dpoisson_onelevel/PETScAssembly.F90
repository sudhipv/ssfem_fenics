!! This MODULE is part of HPSFEM package : Ajit Desai : July 2015
!! This is the assembly code addapted from FEniCs/Puffin FEM solver for : Au = f
!! It assembles the left hand side, FEM matrix (A) and right hand side force vector (f)
!! It's updated on 12/2015 to assemble matrix and vector directly into PETSc format
!! Ajit Desai, 7/2015, Mohammad Khalil, 1/2012
!! Contains: 
!! 1. PETScAssembleADVmatrix : PETSc Advection assembly matrix 
!! 2. PETScAssembleDIFmatrix : PETSc Diffusion assembly matrix
!! 3. PETScAssemblevector    : PETSc Source-assembly vector
!
! Record of revisions:
!     Date            Programmer          Description of change
!  ==========         ==========          =====================
!  27/01/2016            AD               Created new module to Assemble subdomian level
!                                         PETSc-Sparse-Block-Matrices(PAsii, PAsig etc)
!!-----------------------------------------------------------------------------------

module PETScAssembly

!! This Module dependends on following Modules 
use common
use variationalform
implicit none

contains


!!***************************************************************************************
!! Subroutine to assemble right hand side force vector (f)
subroutine PETscSubVecAssembly(Fsi, Fsg, points, edges, triangles, np, ne, nt, nb, &
                               ndim, npceout, nomga, amp, dbounds, mIndex, sIndex)

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    Vec              Fsi, Fsg
    PetscInt         one, pi
    PetscErrorCode   ierr

    integer  :: np, ne, nt, i, j, k, ndim,npceout,nomga
    double precision :: detJ, dx, phis, phiedges
    double precision :: v, integral, ds, pe, amp

    double precision, dimension(2,np) :: points
    integer, dimension(3,ne)          :: edges
    integer, dimension(3)             :: edge
    integer, dimension(2)             :: nodese
    integer, dimension(3,nt)          :: triangles
    integer, dimension(3)             :: triangle
    integer, dimension(3)             :: nodest
    double precision, dimension(2,3)  :: midpoints, coordt
    double precision, dimension(2)    :: gausspoints, dphiv, pt, x, dv
    double precision, dimension(2,2)  :: Jmat, Ji, Jtraninv, coorde, dbounds

    integer, dimension(ndim,npceout)  :: mIndex
    integer, dimension(2,ndim)        :: sIndex

    integer :: nb, ni, ii !!jj
    double precision :: intT, intE

    ! initialize vector
    one  = 1
    ni   = np-nb
    intT = 0.0d0
    intE = 0.0d0

    ! quadrature points
    midpoints(1,:) = [0.5d0, 0.5d0, 0.0d0]
    midpoints(2,:) = [0.0d0, 0.5d0, 0.5d0]
    gausspoints = [(1.0d0-1.0d0/sqrt(3.0d0))/2.0d0, (1.0d0+1.0d0/sqrt(3.0d0))/2.0d0]

    ! Iterate over all elements (interior domain)
    do i = 1,nt
        triangle = triangles(:,i)

        ! Get nodes, coordinates, and domain number
        nodest = triangle(1:3)
        coordt = points(:,nodest)

        ! Compute Jacobian of map and area of element
        Jmat(:,:) = 0.0d0
        do j = 1,3
        call dphi(dphiv,j)
        call outer_product(coordt(:,j),dphiv,Ji)
        Jmat = Jmat + Ji
        enddo

        !obtain inv(J')
        call invert(2,transpose(Jmat),Jtraninv)

        detJ = Jmat(1,1)*Jmat(2,2) - Jmat(2,1)*Jmat(1,2)
        dx = 0.5d0*abs(detJ)

        ! Assemble vector
        do j=1,3
        pt = midpoints(:,j)

        x = 0.0d0
        do k = 1,3
        call phi(phis,k,pt)
        x = x + coordt(:,k)*phis
        enddo

        do k = 1,3
        call phi(v,k,pt)
        call dphi(dphiv,k)
        dv = matmul(Jtraninv,dphiv)

        call advdiff(integral,0.0d0,v,[0.0d0,0.0d0],dv,dx,0.0d0,x,3,amp,&
             0.0d0,ndim,npceout,nomga,[0.0d0,0.0d0,0.0d0,0.0d0],[0.0d0,0.0d0,0.0d0,0.0d0],&
             0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],dbounds,0,0,mIndex,sIndex,0)

        intT = integral/3.0d0
        ii = (nodest(k))
        if (ii <= ni ) then
            pi = ii-1
            call VecSetValues(Fsi,one,pi,intT,ADD_VALUES,ierr)
        else if (ii > ni) then
            pi = (ii-1)-ni
            call VecSetValues(Fsg,one,pi,intT,ADD_VALUES,ierr)
        end if

        enddo
        enddo
    enddo

    ! Iterate over all edges (boundary)
    do i = 1,ne
        edge = edges(:,i)

        ! Get nodes, coordinates, and domain number
        nodese = edge(1:2)
        coorde = points(:,nodese)

        ! Compute length of edge
        ds = sqrt(dot_product((coorde(:,1) - coorde(:,2)),(coorde(:,1) - coorde(:,2))))

        ! Assemble vector
        do j = 1,2
        pe = gausspoints(j)

        x = 0.0d0
        do k = 1,2
        call phiedge(phiedges,k,pe)
        x = x + coorde(:,k)*phiedges
        enddo

        do k = 1,2
        call phiedge(v,k,pe)

        call advdiff(integral,0.0d0,v,[0.0d0,0.0d0],[0.0d0,0.0d0],0.0d0,ds,x,3,amp,&
                     0.0d0,ndim,npceout,nomga,[0.0d0,0.0d0,0.0d0,0.0d0],[0.0d0,0.0d0,0.0d0,0.0d0],&
                     0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],dbounds,0,0,mIndex,sIndex,edge(3))

        intE = integral/2.0d0
        ii = (nodese(k))
        if (ii <= ni ) then
            pi = ii-1
            call VecSetValues(Fsi,one,pi,intE,ADD_VALUES,ierr)
        else if (ii > ni) then
            pi = (ii-1)-ni
            call VecSetValues(Fsg,one,pi,intE,ADD_VALUES,ierr)
        end if

        enddo
        enddo
    enddo

end subroutine PETscSubVecAssembly


!!!***************************************************************************************
!!! Subroutine to assemble left hand side FEM matrix (A)
!!! PArr is not assembled because, it is not required for NNC-PCGM (but for FETIDP we need it)
!subroutine PETScTwoLevelAssemADVmatrix(PAmat,PAcc,PAic,PAci,PAir,PAri,PArc,PAcr,points, edges, triangles, &
!np, ne, nt, nb, nci, eq, meanc,omegas,multipliers,sigma,amps,dbounds,casei,term)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!#include <petsc/finclude/petscmat.h>
!
!    Mat              PAmat,PAcc,PAic,PAci
!    Mat              PAir,PAri,PArc,PAcr
!    PetscInt         pi,pj
!    PetscInt         n,one
!    PetscScalar      imatT, imatE
!    PetscErrorCode   ierr
!
!    !! Variable definition & allocation
!    integer, dimension(3,ne) :: edges
!    integer, dimension(3)    :: edge
!    integer, dimension(2)    :: nodese
!    integer, dimension(3,nt) :: triangles
!    integer, dimension(3)    :: triangle
!    integer, dimension(3)    :: nodest
!    integer :: np, ne, nt, nb, eq, i, j, k, l, term, casei
!
!    double precision, dimension(2,np) :: points
!    double precision, dimension(4)    :: omegas, multipliers
!    double precision, dimension(3)    :: amps
!    double precision, dimension(2,3)  :: midpoints, coordt
!    double precision, dimension(2)    :: gausspoints, dphiv, pt, x, dv, du
!    double precision, dimension(2,2)  :: Jmat, Ji, Jtraninv, coorde, dbounds
!
!    double precision :: detJ, dx, phis, phiedges, v, u
!    double precision :: integral, ds, pe, meanc, sigma
!    integer :: ni, ii, jj, nc, nci
!
!!! Defining variables
!n=np
!one=1
!ni = np-nb            !! Total number of interior nodes for each subdomain
!nc = np-nci           !! Total interior+remainin nodes for each subdomain
!
!imatT = 0.0d0
!imatE = 0.0d0
!
!!!initialize matrix
!!!A(:,:) = 0.0d0
!
!!quadrature points
!midpoints(1,:) = [0.5d0, 0.5d0, 0.0d0]
!midpoints(2,:) = [0.0d0, 0.5d0, 0.5d0]
!gausspoints = [(1.0d0-1.0d0/sqrt(3.0d0))/2.0d0, (1.0d0+1.0d0/sqrt(3.0d0))/2.0d0]
!
!! Iterate over all elements (interior domain)
!do i = 1,nt
!! print*, i, ' of ', nt, 'triangles'
!triangle = triangles(:,i)
!
!! Get nodes, coordinates, and domain number
!nodest = triangle(1:3)
!coordt = points(:,nodest)
!
!! Compute Jacobian of map and area of element
!Jmat = 0.0d0
!do j = 1,3
!call dphi(dphiv,j)
!call outer_product(coordt(:,j),dphiv,Ji)
!Jmat = Jmat + Ji
!enddo
!
!! obtain inv(J') used for transforming grad(u) and grad(v)
!! Jtraninv = transpose(Jmat)
!call invert(2,transpose(Jmat),Jtraninv)
!
!detJ = Jmat(1,1)*Jmat(2,2) - Jmat(2,1)*Jmat(1,2)
!dx = 0.5d0*abs(detJ)
!
!! Assemble matrix
!do j=1,3
!pt = midpoints(:,j)
!
!x = 0.0d0
!do k = 1,3
!call phi(phis,k,pt)
!x = x + coordt(:,k)*phis
!enddo
!
!do k = 1,3
!call phi(v,k,pt)
!call dphi(dphiv,k)
!dv = matmul(Jtraninv,dphiv)
!
!do l = 1,3
!call phi(u,l,pt)
!call dphi(dphiv,l)
!du = matmul(Jtraninv,dphiv)
!
!call advdiff(integral,u,v,du,dv,dx,0.0d0,x,eq,0.0d0,&
!meanc,omegas,multipliers,sigma,amps,dbounds,casei,term,0)
!
!!A(nodest(k),nodest(l)) = A(nodest(k),nodest(l)) + integral/3.0d0
!!A(nodest(k),nodest(l)) = integral/3.0d0
!!imatT = A(nodest(k),nodest(l))
!
!imatT = integral/3.0d0
!ii = nodest(k)
!jj = nodest(l)
!if (ii <= nc .AND. jj <= nc) then
!    pi = (ii-1)
!    pj = (jj-1)
!    if (imatT .ne. 0) then
!    call MatSetValues(PAmat,one,pi,one,pj,imatT,ADD_VALUES,ierr)
!    end if
!end if
!
!if (ii <= ni .AND. jj > nc) then
!    pi = ii-1
!    pj = (jj-nc)-1
!    if (imatT .ne. 0) then
!    call MatSetValues(PAic,one,pi,one,pj,imatT,ADD_VALUES,ierr)
!    end if
!else if (ii > nc .AND. jj <= ni) then
!    pi = (ii-nc)-1
!    pj = jj-1
!    if (imatT .ne. 0) then
!    call MatSetValues(PAci,one,pi,one,pj,imatT,ADD_VALUES,ierr)
!    end if
!else if (ii <= ni .AND. jj > ni .AND. jj <= nc) then
!    pi = ii-1
!    pj = (jj-ni)-1
!    if (imatT .ne. 0) then
!    call MatSetValues(PAir,one,pi,one,pj,imatT,ADD_VALUES,ierr)
!    end if
!else if (ii > ni .AND. ii <= nc .AND. jj <= ni) then
!    pi = (ii-ni)-1
!    pj = jj-1
!    if (imatT .ne. 0) then
!    call MatSetValues(PAri,one,pi,one,pj,imatT,ADD_VALUES,ierr)
!    end if
!else if (ii > ni .AND. ii <= nc .AND. jj > nc) then
!    pi = (ii-ni)-1
!    pj = (jj-nc)-1
!    if (imatT .ne. 0) then
!    call MatSetValues(PArc,one,pi,one,pj,imatT,ADD_VALUES,ierr)
!    end if
!else if (ii > nc .AND. jj > ni .AND. jj <= nc) then
!    pi = (ii-nc)-1
!    pj = (jj-ni)-1
!    if (imatT .ne. 0) then
!    call MatSetValues(PAcr,one,pi,one,pj,imatT,ADD_VALUES,ierr)
!    end if
!else if (ii > nc .AND. jj > nc ) then
!    pi = (ii-nc)-1
!    pj = (jj-nc)-1
!    if (imatT .ne. 0) then
!    call MatSetValues(PAcc,one,pi,one,pj,imatT,ADD_VALUES,ierr)
!    end if
!end if
!
!!pi = (nodest(k)-1)
!!pj = (nodest(l)-1)
!!imatT = integral/3.0d0
!!!if (imat .ne. 0) then
!!call MatSetValues(PetMat,one,pi,one,pj,imatT,ADD_VALUES,ierr)
!!!end if
!
!enddo
!enddo
!enddo
!enddo
!
!
!! Iterate over all edges (boundary)
!do i = 1,ne
!edge = edges(:,i)
!
!! Get nodes, coordinates, and domain number
!nodese = edge(1:2)
!coorde = points(:,nodese)
!
!! Compute length of edge
!ds = sqrt(dot_product((coorde(:,1) - coorde(:,2)),(coorde(:,1) - coorde(:,2))))
!
!! Assemble matrix
!do j = 1,2
!pe = gausspoints(j)
!
!x(:) = 0.0d0
!do k = 1,2
!call phiedge(phiedges,k,pe)
!x = x + coorde(:,k)*phiedges
!enddo
!
!do k = 1,2
!call phiedge(v,k,pe)
!
!do l = 1,2
!call phiedge(u,l,pe)
!call advdiff(integral,u,v,[0.0d0,0.0d0],[0.0d0,0.0d0],0.0d0,ds,x,eq,0.0d0,&
!meanc,omegas,multipliers,sigma,amps,dbounds,casei,term,edge(3))
!
!!A(nodese(k),nodese(l)) = A(nodese(k),nodese(l)) + integral/2.0d0
!!A(nodese(k),nodese(l)) = integral/2.0d0
!!imatE = A(nodese(k),nodese(l))
!
!imatE = integral/2.0d0
!ii = nodese(k)
!jj = nodese(l)
!if (ii <= nc .AND. jj <= nc) then
!    pi = (ii-1)
!    pj = (jj-1)
!    if (imatE .ne. 0) then
!    call MatSetValues(PAmat,one,pi,one,pj,imatE,ADD_VALUES,ierr)
!    end if
!end if
!
!if (ii <= ni .AND. jj > nc) then
!    pi = ii-1
!    pj = (jj-nc)-1
!    if (imatE .ne. 0) then
!    call MatSetValues(PAic,one,pi,one,pj,imatE,ADD_VALUES,ierr)
!    end if
!else if (ii > nc .AND. jj <= ni) then
!    pi = (ii-nc)-1
!    pj = jj-1
!    if (imatE .ne. 0) then
!    call MatSetValues(PAci,one,pi,one,pj,imatE,ADD_VALUES,ierr)
!    end if
!else if (ii <= ni .AND. jj > ni .AND. jj <= nc) then
!    pi = ii-1
!    pj = (jj-ni)-1
!    if (imatE .ne. 0) then
!    call MatSetValues(PAir,one,pi,one,pj,imatE,ADD_VALUES,ierr)
!    end if
!else if (ii > ni .AND. ii <= nc .AND. jj <= ni) then
!    pi = (ii-ni)-1
!    pj = jj-1
!    if (imatE .ne. 0) then
!    call MatSetValues(PAri,one,pi,one,pj,imatE,ADD_VALUES,ierr)
!    end if
!else if (ii > ni .AND. ii <= nc .AND. jj > nc) then
!    pi = (ii-ni)-1
!    pj = (jj-nc)-1
!    if (imatE .ne. 0) then
!    call MatSetValues(PArc,one,pi,one,pj,imatE,ADD_VALUES,ierr)
!    end if
!else if (ii > nc .AND. jj > ni .AND. jj <= nc) then
!    pi = (ii-nc)-1
!    pj = (jj-ni)-1
!    if (imatE .ne. 0) then
!    call MatSetValues(PAcr,one,pi,one,pj,imatE,ADD_VALUES,ierr)
!    end if
!else if (ii > nc .AND. jj > nc ) then
!    pi = (ii-nc)-1
!    pj = (jj-nc)-1
!    if (imatE .ne. 0) then
!    call MatSetValues(PAcc,one,pi,one,pj,imatE,ADD_VALUES,ierr)
!    end if
!end if
!
!!pi = (nodese(k)-1)
!!pj = (nodese(l)-1)
!!imatE = integral/2.0d0
!!!if (imat .ne. 0) then
!!call MatSetValues(PetMat,one,pi,one,pj,imatE,ADD_VALUES,ierr)
!!!end if
!
!enddo
!enddo
!enddo
!enddo
!
!
!END SUBROUTINE PETScTwoLevelAssemADVmatrix
!
!
!!!***************************************************************************************
!!! Subroutine to assemble left hand side FEM matrix (A)
!
!subroutine PETScTwoLevelAssemDIFmatrix(PAmat,PAcc,PAic,PAci,PAir,PAri,PArc,PAcr,points, edges, triangles, &
!np, ne, nt, nb, nci,eq, meanc,omegas,multipliers,sigma,amps,dbounds,casei,term)
!
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!#include <petsc/finclude/petscmat.h>
!
!
!Mat              PAmat,PAcc,PAic,PAci
!Mat              PAir,PAri,PArc,PAcr
!PetscInt         pi,pj
!PetscInt         n,one
!PetscScalar      imatT, imatE
!PetscErrorCode   ierr
!
!!! Variable definition & allocation
!integer, dimension(3,ne) :: edges
!integer, dimension(3) :: edge
!integer, dimension(2) :: nodese
!integer, dimension(3,nt) :: triangles
!integer, dimension(3) :: triangle
!integer, dimension(3) :: nodest
!integer :: np, ne, nt, nb, eq, i, j, k, l, term, casei
!
!!double precision, dimension(np,np) :: A
!double precision, dimension(2,np) :: points
!double precision, dimension(4) :: omegas, multipliers
!double precision, dimension(3) :: amps
!double precision, dimension(2,3) :: midpoints, coordt
!double precision, dimension(2) :: gausspoints, dphiv, pt, x, dv, du
!double precision, dimension(2,2) :: Jmat, Ji, Jtraninv, coorde, dbounds
!double precision :: detJ, dx, phis, phiedges, v, u
!double precision :: integral, ds, pe, meanc, sigma
!
!integer :: ni, ii, jj, nci, nc
!
!n=np
!one=1
!
!ni = np-nb
!nc = np-nci
!
!imatT = 0.0d0
!imatE = 0.0d0
!
!!initialize matrix
!!A(:,:) = 0.0d0
!
!!quadrature points
!midpoints(1,:) = [0.5d0, 0.5d0, 0.0d0]
!midpoints(2,:) = [0.0d0, 0.5d0, 0.5d0]
!gausspoints = [(1.0d0-1.0d0/sqrt(3.0d0))/2.0d0, (1.0d0+1.0d0/sqrt(3.0d0))/2.0d0]
!
!! Iterate over all elements (interior domain)
!do i = 1,nt
!! print*, i, ' of ', nt, 'triangles'
!triangle = triangles(:,i)
!
!! Get nodes, coordinates, and domain number
!nodest = triangle(1:3)
!coordt = points(:,nodest)
!
!! Compute Jacobian of map and area of element
!Jmat = 0.0d0
!do j = 1,3
!call dphi(dphiv,j)
!call outer_product(coordt(:,j),dphiv,Ji)
!Jmat = Jmat + Ji
!enddo
!
!!obtain inv(J') used for transforming grad(u) and grad(v)
!!Jtraninv = transpose(Jmat)
!call invert(2,transpose(Jmat),Jtraninv)
!
!detJ = Jmat(1,1)*Jmat(2,2) - Jmat(2,1)*Jmat(1,2)
!dx = 0.5d0*abs(detJ)
!
!! Assemble matrix
!do j=1,3
!pt = midpoints(:,j)
!
!x = 0.0d0
!do k = 1,3
!call phi(phis,k,pt)
!x = x + coordt(:,k)*phis
!enddo
!
!do k = 1,3
!call phi(v,k,pt)
!call dphi(dphiv,k)
!dv = matmul(Jtraninv,dphiv)
!
!do l = 1,3
!call phi(u,l,pt)
!call dphi(dphiv,l)
!du = matmul(Jtraninv,dphiv)
!
!call advdiff(integral,u,v,du,dv,dx,0.0d0,x,eq,0.0d0,&
!meanc,omegas,multipliers,sigma,amps,dbounds,casei,term,0)
!!A(nodest(k),nodest(l)) = A(nodest(k),nodest(l)) + integral/3.0d0
!!A(nodest(k),nodest(l)) = integral/3.0d0
!!imatT = A(nodest(k),nodest(l))
!
!!ii = nodest(k)
!!jj = nodest(l)
!!if (ii <= ni .AND. jj <= ni) then
!!pi = (nodest(k)-1)
!!pj = (nodest(l)-1)
!!imatT = integral/3.0d0
!!if (imatT .ne. 0) then
!!call MatSetValues(PAii,one,pi,one,pj,imatT,ADD_VALUES,ierr)
!!end if
!!end if
!
!imatT = integral/3.0d0
!ii = nodest(k)
!jj = nodest(l)
!if (ii <= nc .AND. jj <= nc) then
!    pi = (ii-1)
!    pj = (jj-1)
!    if (imatT .ne. 0) then
!    call MatSetValues(PAmat,one,pi,one,pj,imatT,ADD_VALUES,ierr)
!    end if
!end if
!
!if (ii <= ni .AND. jj > nc) then
!    pi = ii-1
!    pj = (jj-nc)-1
!    if (imatT .ne. 0) then
!    call MatSetValues(PAic,one,pi,one,pj,imatT,ADD_VALUES,ierr)
!    end if
!else if (ii > nc .AND. jj <= ni) then
!    pi = (ii-nc)-1
!    pj = jj-1
!    if (imatT .ne. 0) then
!    call MatSetValues(PAci,one,pi,one,pj,imatT,ADD_VALUES,ierr)
!    end if
!else if (ii <= ni .AND. jj > ni .AND. jj <= nc) then
!    pi = ii-1
!    pj = (jj-ni)-1
!    if (imatT .ne. 0) then
!    call MatSetValues(PAir,one,pi,one,pj,imatT,ADD_VALUES,ierr)
!    end if
!else if (ii > ni .AND. ii <= nc .AND. jj <= ni) then
!    pi = (ii-ni)-1
!    pj = jj-1
!    if (imatT .ne. 0) then
!    call MatSetValues(PAri,one,pi,one,pj,imatT,ADD_VALUES,ierr)
!    end if
!else if (ii > ni .AND. ii <= nc .AND. jj > nc) then
!    pi = (ii-ni)-1
!    pj = (jj-nc)-1
!    if (imatT .ne. 0) then
!    call MatSetValues(PArc,one,pi,one,pj,imatT,ADD_VALUES,ierr)
!    end if
!else if (ii > nc .AND. jj > ni .AND. jj <= nc) then
!    pi = (ii-nc)-1
!    pj = (jj-ni)-1
!    if (imatT .ne. 0) then
!    call MatSetValues(PAcr,one,pi,one,pj,imatT,ADD_VALUES,ierr)
!    end if
!else if (ii > nc .AND. jj > nc ) then
!    pi = (ii-nc)-1
!    pj = (jj-nc)-1
!    if (imatT .ne. 0) then
!    call MatSetValues(PAcc,one,pi,one,pj,imatT,ADD_VALUES,ierr)
!    end if
!end if
!
!!pi = (nodest(k)-1)
!!pj = (nodest(l)-1)
!!imatT = integral/3.0d0
!!!if (imat .ne. 0) then
!!call MatSetValues(PetMat,one,pi,one,pj,imatT,ADD_VALUES,ierr)
!!!end if
!
!enddo
!enddo
!enddo
!enddo
!
!! Iterate over all edges (boundary)
!
!do i = 1,ne
!edge = edges(:,i)
!
!! Get nodes, coordinates, and domain number
!nodese = edge(1:2)
!coorde = points(:,nodese)
!
!! Compute length of edge
!ds = sqrt(dot_product((coorde(:,1) - coorde(:,2)),(coorde(:,1) - coorde(:,2))))
!
!! Assemble matrix
!do j = 1,2
!pe = gausspoints(j)
!
!x(:) = 0.0d0
!do k = 1,2
!call phiedge(phiedges,k,pe)
!x = x + coorde(:,k)*phiedges
!enddo
!
!do k = 1,2
!call phiedge(v,k,pe)
!
!do l = 1,2
!call phiedge(u,l,pe)
!call advdiff(integral,u,v,[0.0d0,0.0d0],[0.0d0,0.0d0],0.0d0,ds,x,eq,0.0d0,&
!meanc,omegas,multipliers,sigma,amps,dbounds,casei,term,edge(3))
!
!!A(nodese(k),nodese(l)) = A(nodese(k),nodese(l)) + integral/2.0d0
!!A(nodese(k),nodese(l)) = integral/2.0d0
!!imatE = A(nodese(k),nodese(l))
!
!!ii = nodese(k)
!!jj = nodese(l)
!!if (ii <= ni .AND. jj <= ni) then
!!pi = (nodese(k)-1)
!!pj = (nodese(l)-1)
!!imatE = integral/2.0d0
!!if (imatE .ne. 0) then
!!call MatSetValues(PAii,one,pi,one,pj,imatE,ADD_VALUES,ierr)
!!end if
!!end if
!
!
!imatE = integral/2.0d0
!ii = nodese(k)
!jj = nodese(l)
!if (ii <= nc .AND. jj <= nc) then
!    pi = (ii-1)
!    pj = (jj-1)
!    if (imatE .ne. 0) then
!    call MatSetValues(PAmat,one,pi,one,pj,imatE,ADD_VALUES,ierr)
!    end if
!end if
!
!if (ii <= ni .AND. jj > nc) then
!    pi = ii-1
!    pj = (jj-nc)-1
!    if (imatE .ne. 0) then
!    call MatSetValues(PAic,one,pi,one,pj,imatE,ADD_VALUES,ierr)
!    end if
!else if (ii > nc .AND. jj <= ni) then
!    pi = (ii-nc)-1
!    pj = jj-1
!    if (imatE .ne. 0) then
!    call MatSetValues(PAci,one,pi,one,pj,imatE,ADD_VALUES,ierr)
!    end if
!else if (ii <= ni .AND. jj > ni .AND. jj <= nc) then
!    pi = ii-1
!    pj = (jj-ni)-1
!    if (imatE .ne. 0) then
!    call MatSetValues(PAir,one,pi,one,pj,imatE,ADD_VALUES,ierr)
!    end if
!else if (ii > ni .AND. ii <= nc .AND. jj <= ni) then
!    pi = (ii-ni)-1
!    pj = jj-1
!    if (imatE .ne. 0) then
!    call MatSetValues(PAri,one,pi,one,pj,imatE,ADD_VALUES,ierr)
!    end if
!else if (ii > ni .AND. ii <= nc .AND. jj > nc) then
!    pi = (ii-ni)-1
!    pj = (jj-nc)-1
!    if (imatE .ne. 0) then
!    call MatSetValues(PArc,one,pi,one,pj,imatE,ADD_VALUES,ierr)
!    end if
!else if (ii > nc .AND. jj > ni .AND. jj <= nc) then
!    pi = (ii-nc)-1
!    pj = (jj-ni)-1
!    if (imatE .ne. 0) then
!    call MatSetValues(PAcr,one,pi,one,pj,imatE,ADD_VALUES,ierr)
!    end if
!else if (ii > nc .AND. jj > nc ) then
!    pi = (ii-nc)-1
!    pj = (jj-nc)-1
!    if (imatE .ne. 0) then
!    call MatSetValues(PAcc,one,pi,one,pj,imatE,ADD_VALUES,ierr)
!    end if
!end if
!
!
!!pi = (nodese(k)-1)
!!pj = (nodese(l)-1)
!!imatE = integral/2.0d0
!!!if (imat .ne. 0) then
!!call MatSetValues(PetMat,one,pi,one,pj,imatE,ADD_VALUES,ierr)
!!!end if
!
!enddo
!enddo
!enddo
!enddo
!
!
!END SUBROUTINE PETScTwoLevelAssemDIFmatrix
!
!
!
!!!***************************************************************************************
!!! Subroutine to assemble right hand side force vector (f)
!SUBROUTINE PETScTwoLevelAssemVec(PFi, PFg, PFc, points, edges, triangles, np, ne, nt, nb, nci, &
!                                 amp, dbounds,tempi,tempb,tempc)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!!!#include <petsc/finclude/petscmat.h>
!
!Vec              PFi, PFg, PFc
!PetscInt         pi,ii
!PetscInt         ni,one
!PetscScalar      ivecT,ivecE
!PetscErrorCode   ierr
!
!integer, dimension(3,ne) :: edges
!integer, dimension(3) :: edge
!integer, dimension(2) :: nodese
!integer, dimension(3,nt) :: triangles
!integer, dimension(3) :: triangle
!integer, dimension(3) :: nodest
!integer :: np, ne, nt, nb, nci, nc, i, j, k
!
!!double precision, dimension(np) :: b
!double precision, dimension(2,np) :: points
!double precision :: detJ, dx, phis, phiedges, v, integral, ds, pe, amp
!double precision, dimension(2,3) :: midpoints, coordt
!double precision, dimension(2) :: gausspoints, dphiv, pt, x, dv
!double precision, dimension(2,2) :: Jmat, Ji, Jtraninv, coorde, dbounds
!
!integer, dimension(np-nb) :: tempi
!integer, dimension(nb)    :: tempb
!integer, dimension(nci)   :: tempc
!
!
!ni=np-nb
!nc=np-nci
!one=1
!
!!initialize vector
!!b(:) = 0.0d0
!
!!quadrature points
!midpoints(1,:) = [0.5d0, 0.5d0, 0.0d0]
!midpoints(2,:) = [0.0d0, 0.5d0, 0.5d0]
!gausspoints = [(1.0d0-1.0d0/sqrt(3.0d0))/2.0d0, (1.0d0+1.0d0/sqrt(3.0d0))/2.0d0]
!
!! Iterate over all elements (interior domain)
!do i = 1,nt
!triangle = triangles(:,i)
!
!! Get nodes, coordinates, and domain number
!nodest = triangle(1:3)
!coordt = points(:,nodest)
!
!! Compute Jacobian of map and area of element
!Jmat(:,:) = 0.0d0
!do j = 1,3
!call dphi(dphiv,j)
!call outer_product(coordt(:,j),dphiv,Ji)
!Jmat = Jmat + Ji
!enddo
!
!!obtain inv(J')
!call invert(2,transpose(Jmat),Jtraninv)
!
!detJ = Jmat(1,1)*Jmat(2,2) - Jmat(2,1)*Jmat(1,2)
!dx = 0.5d0*abs(detJ)
!
!! Assemble vector
!do j=1,3
!pt = midpoints(:,j)
!
!x = 0.0d0
!do k = 1,3
!call phi(phis,k,pt)
!x = x + coordt(:,k)*phis
!enddo
!
!do k = 1,3
!call phi(v,k,pt)
!call dphi(dphiv,k)
!dv = matmul(Jtraninv,dphiv)
!
!call advdiff(integral,0.0d0,v,[0.0d0,0.0d0],dv,dx,0.0d0,x,3,amp,&
!0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],[0.0d0,0.0d0,0.0d0,0.0d0],&
!0.0d0,[0.0d0,0.0d0,0.0d0],dbounds,0,0,0)
!
!ivecT = integral/3.0d0
!ii = nodest(k)
!if (ii <= ni) then
!    pi = ii-1
!    tempi(ii) = pi
!    call VecSetValues(PFi,one,pi,ivecT,ADD_VALUES,ierr)
!else if (ii > ni) then
!    pi = (ii-ni)-1
!    tempb(ii-ni) = pi
!    call VecSetValues(PFg,one,pi,ivecT,ADD_VALUES,ierr)
!end if
!if (ii > nc) then
!    pi = (ii-nc)-1
!    tempc(ii-nc) = pi
!    call VecSetValues(PFc,one,pi,ivecT,ADD_VALUES,ierr)
!end if
!
!!b(nodest(k)) = b(nodest(k)) + integral/3.0d0
!!ivec = b(nodest(k))
!
!!ivecT = integral/3.0d0
!!ii = (nodest(k)-1)
!!call VecSetValues(PetVec,one,ii,ivecT,ADD_VALUES,ierr)
!
!enddo
!enddo
!enddo
!
!! Iterate over all edges (boundary)
!do i = 1,ne
!edge = edges(:,i)
!
!! Get nodes, coordinates, and domain number
!nodese = edge(1:2)
!coorde = points(:,nodese)
!
!! Compute length of edge
!ds = sqrt(dot_product((coorde(:,1) - coorde(:,2)),(coorde(:,1) - coorde(:,2))))
!
!! Assemble vector
!do j = 1,2
!pe = gausspoints(j)
!
!x = 0.0d0
!do k = 1,2
!call phiedge(phiedges,k,pe)
!x = x + coorde(:,k)*phiedges
!enddo
!
!do k = 1,2
!call phiedge(v,k,pe)
!
!call advdiff(integral,0.0d0,v,[0.0d0,0.0d0],[0.0d0,0.0d0],0.0d0,ds,x,3,amp,&
!0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],[0.0d0,0.0d0,0.0d0,0.0d0],0.0d0,&
![0.0d0,0.0d0,0.0d0],dbounds,0,0,edge(3))
!
!!b(nodese(k)) = b(nodese(k)) + integral/2.0d0
!!ivec2 = b(nodese(k))
!
!ivecE = integral/2.0d0
!ii = nodest(k)
!if (ii <= ni) then
!    pi = ii-1
!    tempi(ii) = pi
!    call VecSetValues(PFi,one,pi,ivecE,ADD_VALUES,ierr)
!else if (ii > ni) then
!    pi = (ii-ni)-1
!    tempb(ii-ni) = pi
!    call VecSetValues(PFg,one,pi,ivecE,ADD_VALUES,ierr)
!end if
!if (ii > nc) then
!    pi = (ii-nc)-1
!    tempc(ii-nc) = pi
!    call VecSetValues(PFc,one,pi,ivecE,ADD_VALUES,ierr)
!end if
!
!!ivecE = integral/2.0d0
!!jj = (nodese(k)-1)
!!call VecSetValues(PetVec,one,jj,ivecE,ADD_VALUES,ierr)
!
!enddo
!enddo
!enddo
!
!END SUBROUTINE PETScTwoLevelAssemVec
!
!!!***************************************************************************************
!
!
!
!
!
!!!***************************************************************************************
!!! Subroutine to assemble left hand side FEM matrix (A)
!
!subroutine PETScOneLevelAssemADVmatrix(PAii,PAgg,PAig,PAgi,points, edges, triangles, np, ne, nt, nb,&
!           eq, meanc,omegas,multipliers,sigma,amps,dbounds,casei,term)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!#include <petsc/finclude/petscmat.h>
!
!Mat              PAii,PAgg,PAig,PAgi
!PetscInt         pi,pj
!PetscInt         n,one
!PetscScalar      imatT, imatE
!PetscErrorCode   ierr
!
!!! Variable definition & allocation
!integer, dimension(3,ne) :: edges
!integer, dimension(3) :: edge
!integer, dimension(2) :: nodese
!integer, dimension(3,nt) :: triangles
!integer, dimension(3) :: triangle
!integer, dimension(3) :: nodest
!integer :: np, ne, nt, nb, eq, i, j, k, l, term, casei
!
!double precision, dimension(2,np) :: points
!double precision, dimension(4) :: omegas, multipliers
!double precision, dimension(3) :: amps
!double precision, dimension(2,3) :: midpoints, coordt
!double precision, dimension(2) :: gausspoints, dphiv, pt, x, dv, du
!double precision, dimension(2,2) :: Jmat, Ji, Jtraninv, coorde, dbounds
!double precision :: detJ, dx, phis, phiedges, v, u
!double precision :: integral, ds, pe, meanc, sigma
!!!double precision, dimension(np,np) :: A
!
!integer :: ni, ii, jj
!
!n=np
!one=1
!ni = np-nb
!imatT = 0.0d0
!imatE = 0.0d0
!
!!!initialize matrix
!!!A(:,:) = 0.0d0
!
!!quadrature points
!    midpoints(1,:) = [0.5d0, 0.5d0, 0.0d0]
!    midpoints(2,:) = [0.0d0, 0.5d0, 0.5d0]
!    gausspoints = [(1.0d0-1.0d0/sqrt(3.0d0))/2.0d0, (1.0d0+1.0d0/sqrt(3.0d0))/2.0d0]
!
!! Iterate over all elements (interior domain)
!do i = 1,nt
!    ! print*, i, ' of ', nt, 'triangles'
!    triangle = triangles(:,i)
!
!    ! Get nodes, coordinates, and domain number
!    nodest = triangle(1:3)
!    coordt = points(:,nodest)
!
!    ! Compute Jacobian of map and area of element
!    Jmat = 0.0d0
!    do j = 1,3
!        call dphi(dphiv,j)
!        call outer_product(coordt(:,j),dphiv,Ji)
!        Jmat = Jmat + Ji
!    enddo
!
!    ! obtain inv(J') used for transforming grad(u) and grad(v)
!    ! Jtraninv = transpose(Jmat)
!    call invert(2,transpose(Jmat),Jtraninv)
!
!    detJ = Jmat(1,1)*Jmat(2,2) - Jmat(2,1)*Jmat(1,2)
!    dx = 0.5d0*abs(detJ)
!
!    ! Assemble matrix
!    do j=1,3
!    pt = midpoints(:,j)
!
!    x = 0.0d0
!    do k = 1,3
!    call phi(phis,k,pt)
!    x = x + coordt(:,k)*phis
!    enddo
!
!    do k = 1,3
!    call phi(v,k,pt)
!    call dphi(dphiv,k)
!    dv = matmul(Jtraninv,dphiv)
!
!    do l = 1,3
!    call phi(u,l,pt)
!    call dphi(dphiv,l)
!    du = matmul(Jtraninv,dphiv)
!
!    call advdiff(integral,u,v,du,dv,dx,0.0d0,x,eq,0.0d0,&
!    meanc,omegas,multipliers,sigma,amps,dbounds,casei,term,0)
!
!    !A(nodest(k),nodest(l)) = A(nodest(k),nodest(l)) + integral/3.0d0
!    !A(nodest(k),nodest(l)) = integral/3.0d0
!    !imatT = A(nodest(k),nodest(l))
!
!imatT = integral/3.0d0
!ii = nodest(k)
!jj = nodest(l)
!if (ii <= ni .AND. jj <= ni) then
!    pi = (nodest(k)-1)
!    pj = (nodest(l)-1)
!    if (imatT .ne. 0) then
!    call MatSetValues(PAii,one,pi,one,pj,imatT,ADD_VALUES,ierr)
!    end if
!else if (ii > ni .AND. jj > ni) then
!    pi = ((nodest(k)-ni)-1)
!    pj = ((nodest(l)-ni)-1)
!    if (imatT .ne. 0) then
!    call MatSetValues(PAgg,one,pi,one,pj,imatT,ADD_VALUES,ierr)
!    end if
!else if (ii <= ni .AND. jj > ni) then
!    pi = (nodest(k)-1)
!    pj = ((nodest(l)-ni)-1)
!    if (imatT .ne. 0) then
!    call MatSetValues(PAig,one,pi,one,pj,imatT,ADD_VALUES,ierr)
!    end if
!else if (ii > ni .AND. jj <= ni) then
!    pi = ((nodest(k)-ni)-1)
!    pj = (nodest(l)-1)
!    if (imatT .ne. 0) then
!    call MatSetValues(PAgi,one,pi,one,pj,imatT,ADD_VALUES,ierr)
!    end if
!end if
!
!    !pi = (nodest(k)-1)
!    !pj = (nodest(l)-1)
!    !imatT = integral/3.0d0
!    !!if (imat .ne. 0) then
!    !call MatSetValues(PetMat,one,pi,one,pj,imatT,ADD_VALUES,ierr)
!    !!end if
!
!    enddo
!    enddo
!    enddo
!enddo
!
!! Iterate over all edges (boundary)
!
!do i = 1,ne
!    edge = edges(:,i)
!
!    ! Get nodes, coordinates, and domain number
!    nodese = edge(1:2)
!    coorde = points(:,nodese)
!
!    ! Compute length of edge
!    ds = sqrt(dot_product((coorde(:,1) - coorde(:,2)),(coorde(:,1) - coorde(:,2))))
!
!    ! Assemble matrix
!    do j = 1,2
!    pe = gausspoints(j)
!
!    x(:) = 0.0d0
!    do k = 1,2
!    call phiedge(phiedges,k,pe)
!    x = x + coorde(:,k)*phiedges
!    enddo
!
!    do k = 1,2
!    call phiedge(v,k,pe)
!
!    do l = 1,2
!    call phiedge(u,l,pe)
!    call advdiff(integral,u,v,[0.0d0,0.0d0],[0.0d0,0.0d0],0.0d0,ds,x,eq,0.0d0,&
!    meanc,omegas,multipliers,sigma,amps,dbounds,casei,term,edge(3))
!
!    !A(nodese(k),nodese(l)) = A(nodese(k),nodese(l)) + integral/2.0d0
!    !A(nodese(k),nodese(l)) = integral/2.0d0
!    !imatE = A(nodese(k),nodese(l))
!
!imatE = integral/2.0d0
!ii = nodese(k)
!jj = nodese(l)
!if (ii <= ni .AND. jj <= ni) then
!    pi = (ii-1)
!    pj = (jj-1)
!    if (imatE .ne. 0) then
!    call MatSetValues(PAii,one,pi,one,pj,imatE,ADD_VALUES,ierr)
!    end if
!else if (ii > ni .AND. jj > ni) then
!    pi = (ii-1)-ni
!    pj = (jj-1)-ni
!    if (imatE .ne. 0) then
!    call MatSetValues(PAgg,one,pi,one,pj,imatE,ADD_VALUES,ierr)
!    end if
!else if (ii <= ni .AND. jj > ni) then
!    pi = (ii-1)
!    pj = (jj-1)-ni
!    if (imatE .ne. 0) then
!    call MatSetValues(PAig,one,pi,one,pj,imatE,ADD_VALUES,ierr)
!    end if
!else if (ii > ni .AND. jj <= ni) then
!    pi = (ii-1)-ni
!    pj = (jj-1)
!    if (imatE .ne. 0) then
!    call MatSetValues(PAgi,one,pi,one,pj,imatE,ADD_VALUES,ierr)
!    end if
!end if
!
!    !pi = (nodese(k)-1)
!    !pj = (nodese(l)-1)
!    !imatE = integral/2.0d0
!    !!if (imat .ne. 0) then
!    !call MatSetValues(PetMat,one,pi,one,pj,imatE,ADD_VALUES,ierr)
!    !!end if
!
!    enddo
!    enddo
!    enddo
!enddo
!
!
!END SUBROUTINE PETScOneLevelAssemADVmatrix
!
!
!!!***************************************************************************************
!!! Subroutine to assemble left hand side FEM matrix (A)
!
!subroutine PETScOneLevelAssemDIFmatrix(PAii,PAgg,PAig,PAgi,points, edges, triangles, np,&
!                  ne, nt, nb, eq, meanc,omegas,multipliers,sigma,amps,dbounds,casei,term)
!
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!#include <petsc/finclude/petscmat.h>
!
!
!Mat              PAii,PAgg,PAig,PAgi
!PetscInt         pi,pj
!PetscInt         n,one
!PetscScalar      imatT, imatE
!PetscErrorCode   ierr
!
!!! Variable definition & allocation
!integer, dimension(3,ne) :: edges
!integer, dimension(3) :: edge
!integer, dimension(2) :: nodese
!integer, dimension(3,nt) :: triangles
!integer, dimension(3) :: triangle
!integer, dimension(3) :: nodest
!integer :: np, ne, nt, nb, eq, i, j, k, l, term, casei
!
!!double precision, dimension(np,np) :: A
!double precision, dimension(2,np) :: points
!double precision, dimension(4) :: omegas, multipliers
!double precision, dimension(3) :: amps
!double precision, dimension(2,3) :: midpoints, coordt
!double precision, dimension(2) :: gausspoints, dphiv, pt, x, dv, du
!double precision, dimension(2,2) :: Jmat, Ji, Jtraninv, coorde, dbounds
!double precision :: detJ, dx, phis, phiedges, v, u
!double precision :: integral, ds, pe, meanc, sigma
!
!integer :: ni, ii, jj
!
!    n=np
!    one=1
!    ni = np-nb
!
!    imatT = 0.0d0
!    imatE = 0.0d0
!!initialize matrix
!!A(:,:) = 0.0d0
!
!!quadrature points
!    midpoints(1,:) = [0.5d0, 0.5d0, 0.0d0]
!    midpoints(2,:) = [0.0d0, 0.5d0, 0.5d0]
!    gausspoints = [(1.0d0-1.0d0/sqrt(3.0d0))/2.0d0, (1.0d0+1.0d0/sqrt(3.0d0))/2.0d0]
!
!! Iterate over all elements (interior domain)
!do i = 1,nt
!    ! print*, i, ' of ', nt, 'triangles'
!    triangle = triangles(:,i)
!
!    ! Get nodes, coordinates, and domain number
!    nodest = triangle(1:3)
!    coordt = points(:,nodest)
!
!    ! Compute Jacobian of map and area of element
!    Jmat = 0.0d0
!    do j = 1,3
!    call dphi(dphiv,j)
!    call outer_product(coordt(:,j),dphiv,Ji)
!    Jmat = Jmat + Ji
!    enddo
!
!    !obtain inv(J') used for transforming grad(u) and grad(v)
!    !Jtraninv = transpose(Jmat)
!    call invert(2,transpose(Jmat),Jtraninv)
!
!    detJ = Jmat(1,1)*Jmat(2,2) - Jmat(2,1)*Jmat(1,2)
!    dx = 0.5d0*abs(detJ)
!
!    ! Assemble matrix
!    do j=1,3
!    pt = midpoints(:,j)
!
!    x = 0.0d0
!    do k = 1,3
!    call phi(phis,k,pt)
!    x = x + coordt(:,k)*phis
!    enddo
!
!    do k = 1,3
!    call phi(v,k,pt)
!    call dphi(dphiv,k)
!    dv = matmul(Jtraninv,dphiv)
!
!    do l = 1,3
!    call phi(u,l,pt)
!    call dphi(dphiv,l)
!    du = matmul(Jtraninv,dphiv)
!
!    call advdiff(integral,u,v,du,dv,dx,0.0d0,x,eq,0.0d0,&
!    meanc,omegas,multipliers,sigma,amps,dbounds,casei,term,0)
!    !A(nodest(k),nodest(l)) = A(nodest(k),nodest(l)) + integral/3.0d0
!    !A(nodest(k),nodest(l)) = integral/3.0d0
!    !imatT = A(nodest(k),nodest(l))
!
!    !ii = nodest(k)
!    !jj = nodest(l)
!    !if (ii <= ni .AND. jj <= ni) then
!    !pi = (nodest(k)-1)
!    !pj = (nodest(l)-1)
!    !imatT = integral/3.0d0
!    !if (imatT .ne. 0) then
!    !call MatSetValues(PAii,one,pi,one,pj,imatT,ADD_VALUES,ierr)
!    !end if
!    !end if
!
!imatT = integral/3.0d0
!ii = nodest(k)
!jj = nodest(l)
!if (ii <= ni .AND. jj <= ni) then
!    pi = (ii-1)
!    pj = (jj-1)
!    !imatT = integral/3.0d0
!    if (imatT .ne. 0) then
!    call MatSetValues(PAii,one,pi,one,pj,imatT,ADD_VALUES,ierr)
!    end if
!else if (ii > ni .AND. jj > ni) then
!    pi = (ii-1)-ni
!    pj = (jj-1)-ni
!    !imatT = integral/3.0d0
!    if (imatT .ne. 0) then
!    call MatSetValues(PAgg,one,pi,one,pj,imatT,ADD_VALUES,ierr)
!    end if
!else if (ii <= ni .AND. jj > ni) then
!    pi = (ii-1)
!    pj = (jj-1)-ni
!    !imatT = integral/3.0d0
!    if (imatT .ne. 0) then
!    call MatSetValues(PAig,one,pi,one,pj,imatT,ADD_VALUES,ierr)
!    end if
!else if (ii > ni .AND. jj <= ni) then
!    pi = (ii-1)-ni
!    pj = (jj-1)
!    !imatT = integral/3.0d0
!    if (imatT .ne. 0) then
!    call MatSetValues(PAgi,one,pi,one,pj,imatT,ADD_VALUES,ierr)
!    end if
!end if
!
!    !pi = (nodest(k)-1)
!    !pj = (nodest(l)-1)
!    !imatT = integral/3.0d0
!    !!if (imat .ne. 0) then
!    !call MatSetValues(PetMat,one,pi,one,pj,imatT,ADD_VALUES,ierr)
!    !!end if
!
!    enddo
!    enddo
!    enddo
!enddo
!
!! Iterate over all edges (boundary)
!
!do i = 1,ne
!    edge = edges(:,i)
!
!    ! Get nodes, coordinates, and domain number
!    nodese = edge(1:2)
!    coorde = points(:,nodese)
!
!    ! Compute length of edge
!    ds = sqrt(dot_product((coorde(:,1) - coorde(:,2)),(coorde(:,1) - coorde(:,2))))
!
!    ! Assemble matrix
!    do j = 1,2
!    pe = gausspoints(j)
!
!    x(:) = 0.0d0
!    do k = 1,2
!    call phiedge(phiedges,k,pe)
!    x = x + coorde(:,k)*phiedges
!    enddo
!
!    do k = 1,2
!    call phiedge(v,k,pe)
!
!    do l = 1,2
!    call phiedge(u,l,pe)
!    call advdiff(integral,u,v,[0.0d0,0.0d0],[0.0d0,0.0d0],0.0d0,ds,x,eq,0.0d0,&
!    meanc,omegas,multipliers,sigma,amps,dbounds,casei,term,edge(3))
!
!    !A(nodese(k),nodese(l)) = A(nodese(k),nodese(l)) + integral/2.0d0
!    !A(nodese(k),nodese(l)) = integral/2.0d0
!    !imatE = A(nodese(k),nodese(l))
!
!    !ii = nodese(k)
!    !jj = nodese(l)
!    !if (ii <= ni .AND. jj <= ni) then
!    !pi = (nodese(k)-1)
!    !pj = (nodese(l)-1)
!    !imatE = integral/2.0d0
!    !if (imatE .ne. 0) then
!    !call MatSetValues(PAii,one,pi,one,pj,imatE,ADD_VALUES,ierr)
!    !end if
!    !end if
!
!imatE = integral/2.0d0
!ii = nodese(k)
!jj = nodese(l)
!if (ii <= ni .AND. jj <= ni) then
!    pi = (ii-1)
!    pj = (jj-1)
!    if (imatE .ne. 0) then
!    call MatSetValues(PAii,one,pi,one,pj,imatE,ADD_VALUES,ierr)
!    end if
!else if (ii > ni .AND. jj > ni) then
!    pi = (ii-1)-ni
!    pj = (jj-1)-ni
!    if (imatE .ne. 0) then
!    call MatSetValues(PAgg,one,pi,one,pj,imatE,ADD_VALUES,ierr)
!    end if
!else if (ii <= ni .AND. jj > ni) then
!    pi = (ii-1)
!    pj = (jj-1)-ni
!    if (imatE .ne. 0) then
!    call MatSetValues(PAig,one,pi,one,pj,imatE,ADD_VALUES,ierr)
!    end if
!else if (ii > ni .AND. jj <= ni) then
!    pi = (ii-1)-ni
!    pj = (jj-1)
!    if (imatE .ne. 0) then
!    call MatSetValues(PAgi,one,pi,one,pj,imatE,ADD_VALUES,ierr)
!    end if
!end if
!
!    !pi = (nodese(k)-1)
!    !pj = (nodese(l)-1)
!    !imatE = integral/2.0d0
!    !!if (imat .ne. 0) then
!    !call MatSetValues(PetMat,one,pi,one,pj,imatE,ADD_VALUES,ierr)
!    !!end if
!
!    enddo
!    enddo
!    enddo
!enddo
!
!
!END SUBROUTINE PETScOneLevelAssemDIFmatrix
!
!


!!!***************************************************************************************
!!! Subroutine to assemble right hand side force vector (f)
!SUBROUTINE PETScOneLevelAssemVec(PFi, PFg, points, edges, triangles, np, ne, nt, nb, amp, dbounds,tempi,tempb)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!!!#include <petsc/finclude/petscmat.h>
!
!    Vec              PFi, PFg
!    PetscInt         pi,ii
!    PetscInt         ni,one
!    PetscScalar      ivecT,ivecE
!    PetscErrorCode   ierr
!
!    integer, dimension(3,ne) :: edges
!    integer, dimension(3) :: edge
!    integer, dimension(2) :: nodese
!    integer, dimension(3,nt) :: triangles
!    integer, dimension(3) :: triangle
!    integer, dimension(3) :: nodest
!    integer :: np, ne, nt, nb, i, j, k
!
!    !double precision, dimension(np) :: b
!    double precision, dimension(2,np) :: points
!    double precision :: detJ, dx, phis, phiedges, v, integral, ds, pe, amp
!    double precision, dimension(2,3) :: midpoints, coordt
!    double precision, dimension(2) :: gausspoints, dphiv, pt, x, dv
!    double precision, dimension(2,2) :: Jmat, Ji, Jtraninv, coorde, dbounds
!
!    integer, dimension(np-nb) :: tempi
!    integer, dimension(nb) :: tempb
!
!
!    ni=np-nb
!    one=1
!
!    !initialize vector
!    !b(:) = 0.0d0
!
!    !quadrature points
!    midpoints(1,:) = [0.5d0, 0.5d0, 0.0d0]
!    midpoints(2,:) = [0.0d0, 0.5d0, 0.5d0]
!    gausspoints = [(1.0d0-1.0d0/sqrt(3.0d0))/2.0d0, (1.0d0+1.0d0/sqrt(3.0d0))/2.0d0]
!
!    ! Iterate over all elements (interior domain)
!    do i = 1,nt
!    triangle = triangles(:,i)
!
!    ! Get nodes, coordinates, and domain number
!    nodest = triangle(1:3)
!    coordt = points(:,nodest)
!
!    ! Compute Jacobian of map and area of element
!    Jmat(:,:) = 0.0d0
!    do j = 1,3
!    call dphi(dphiv,j)
!    call outer_product(coordt(:,j),dphiv,Ji)
!    Jmat = Jmat + Ji
!    enddo
!
!    !obtain inv(J')
!    call invert(2,transpose(Jmat),Jtraninv)
!
!    detJ = Jmat(1,1)*Jmat(2,2) - Jmat(2,1)*Jmat(1,2)
!    dx = 0.5d0*abs(detJ)
!
!    ! Assemble vector
!    do j=1,3
!    pt = midpoints(:,j)
!
!    x = 0.0d0
!    do k = 1,3
!    call phi(phis,k,pt)
!    x = x + coordt(:,k)*phis
!    enddo
!
!    do k = 1,3
!    call phi(v,k,pt)
!    call dphi(dphiv,k)
!    dv = matmul(Jtraninv,dphiv)
!
!    call advdiff(integral,0.0d0,v,[0.0d0,0.0d0],dv,dx,0.0d0,x,3,amp,&
!    0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],[0.0d0,0.0d0,0.0d0,0.0d0],&
!    0.0d0,[0.0d0,0.0d0,0.0d0],dbounds,0,0,0)
!
!    ivecT = integral/3.0d0
!    ii = nodest(k)
!    if (ii <= ni) then
!        pi = ii-1
!        tempi(ii) = pi
!        call VecSetValues(PFi,one,pi,ivecT,ADD_VALUES,ierr)
!    else
!        pi = (ii-1)-ni
!        tempb(ii-ni) = pi
!        call VecSetValues(PFg,one,pi,ivecT,ADD_VALUES,ierr)
!    end if
!
!    !b(nodest(k)) = b(nodest(k)) + integral/3.0d0
!    !ivec = b(nodest(k))
!
!    !ivecT = integral/3.0d0
!    !ii = (nodest(k)-1)
!    !call VecSetValues(PetVec,one,ii,ivecT,ADD_VALUES,ierr)
!
!    enddo
!    enddo
!    enddo
!
!    ! Iterate over all edges (boundary)
!    do i = 1,ne
!    edge = edges(:,i)
!
!    ! Get nodes, coordinates, and domain number
!    nodese = edge(1:2)
!    coorde = points(:,nodese)
!
!    ! Compute length of edge
!    ds = sqrt(dot_product((coorde(:,1) - coorde(:,2)),(coorde(:,1) - coorde(:,2))))
!
!    ! Assemble vector
!    do j = 1,2
!    pe = gausspoints(j)
!
!    x = 0.0d0
!    do k = 1,2
!    call phiedge(phiedges,k,pe)
!    x = x + coorde(:,k)*phiedges
!    enddo
!
!    do k = 1,2
!    call phiedge(v,k,pe)
!
!    call advdiff(integral,0.0d0,v,[0.0d0,0.0d0],[0.0d0,0.0d0],0.0d0,ds,x,3,amp,&
!    0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],[0.0d0,0.0d0,0.0d0,0.0d0],0.0d0,&
!    [0.0d0,0.0d0,0.0d0],dbounds,0,0,edge(3))
!
!    !b(nodese(k)) = b(nodese(k)) + integral/2.0d0
!    !ivec2 = b(nodese(k))
!
!    ivecE = integral/2.0d0
!    ii = nodest(k)
!    if (ii <= ni) then
!        pi = ii-1
!        tempi(ii) = pi
!        call VecSetValues(PFi,one,pi,ivecE,ADD_VALUES,ierr)
!    else
!        pi = (ii-1)-ni
!        tempb(ii-ni) = pi
!        call VecSetValues(PFg,one,pi,ivecE,ADD_VALUES,ierr)
!    end if
!
!    !ivecE = integral/2.0d0
!    !jj = (nodese(k)-1)
!    !call VecSetValues(PetVec,one,jj,ivecE,ADD_VALUES,ierr)
!
!    enddo
!    enddo
!    enddo
!
!END SUBROUTINE PETScOneLevelAssemVec
!
!
!
!!!***************************************************************************************
!!! Subroutine to assemble right hand side force vector (f)
!SUBROUTINE PETScAssemblevector(PetVec, SolVec, points, edges, triangles, np, ne, nt, amp, dbounds)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!!!#include <petsc/finclude/petscmat.h>
!
!Vec              PetVec, SolVec
!PetscInt         ii,jj
!PetscInt         nn,one
!PetscScalar      ivecT,ivecE
!PetscErrorCode   ierr
!
!integer, dimension(3,ne) :: edges
!integer, dimension(3) :: edge
!integer, dimension(2) :: nodese
!integer, dimension(3,nt) :: triangles
!integer, dimension(3) :: triangle
!integer, dimension(3) :: nodest
!integer :: np, ne, nt, i, j, k
!
!!double precision, dimension(np) :: b
!double precision, dimension(2,np) :: points
!double precision :: detJ, dx, phis, phiedges, v, integral, ds, pe, amp
!double precision, dimension(2,3) :: midpoints, coordt
!double precision, dimension(2) :: gausspoints, dphiv, pt, x, dv
!double precision, dimension(2,2) :: Jmat, Ji, Jtraninv, coorde, dbounds
!
!
!nn=np
!one=1
!
!!initialize vector
!!b(:) = 0.0d0
!
!!quadrature points
!    midpoints(1,:) = [0.5d0, 0.5d0, 0.0d0]
!    midpoints(2,:) = [0.0d0, 0.5d0, 0.5d0]
!    gausspoints = [(1.0d0-1.0d0/sqrt(3.0d0))/2.0d0, (1.0d0+1.0d0/sqrt(3.0d0))/2.0d0]
!
!! Iterate over all elements (interior domain)
!do i = 1,nt
!    triangle = triangles(:,i)
!
!    ! Get nodes, coordinates, and domain number
!    nodest = triangle(1:3)
!    coordt = points(:,nodest)
!
!    ! Compute Jacobian of map and area of element
!    Jmat(:,:) = 0.0d0
!    do j = 1,3
!    call dphi(dphiv,j)
!    call outer_product(coordt(:,j),dphiv,Ji)
!    Jmat = Jmat + Ji
!    enddo
!
!    !obtain inv(J')
!    call invert(2,transpose(Jmat),Jtraninv)
!
!    detJ = Jmat(1,1)*Jmat(2,2) - Jmat(2,1)*Jmat(1,2)
!    dx = 0.5d0*abs(detJ)
!
!    ! Assemble vector
!    do j=1,3
!    pt = midpoints(:,j)
!
!    x = 0.0d0
!    do k = 1,3
!    call phi(phis,k,pt)
!    x = x + coordt(:,k)*phis
!    enddo
!
!    do k = 1,3
!    call phi(v,k,pt)
!    call dphi(dphiv,k)
!    dv = matmul(Jtraninv,dphiv)
!
!    call advdiff(integral,0.0d0,v,[0.0d0,0.0d0],dv,dx,0.0d0,x,3,amp,&
!    0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],[0.0d0,0.0d0,0.0d0,0.0d0],&
!    0.0d0,[0.0d0,0.0d0,0.0d0],dbounds,0,0,0)
!
!    !b(nodest(k)) = b(nodest(k)) + integral/3.0d0
!    !ivec = b(nodest(k))
!
!    ivecT = integral/3.0d0
!    ii = (nodest(k)-1)
!    call VecSetValues(PetVec,one,ii,ivecT,ADD_VALUES,ierr)
!
!    enddo
!    enddo
!enddo
!
!! Iterate over all edges (boundary)
!do i = 1,ne
!    edge = edges(:,i)
!
!    ! Get nodes, coordinates, and domain number
!    nodese = edge(1:2)
!    coorde = points(:,nodese)
!
!    ! Compute length of edge
!    ds = sqrt(dot_product((coorde(:,1) - coorde(:,2)),(coorde(:,1) - coorde(:,2))))
!
!    ! Assemble vector
!    do j = 1,2
!    pe = gausspoints(j)
!
!    x = 0.0d0
!    do k = 1,2
!    call phiedge(phiedges,k,pe)
!    x = x + coorde(:,k)*phiedges
!    enddo
!
!    do k = 1,2
!    call phiedge(v,k,pe)
!
!    call advdiff(integral,0.0d0,v,[0.0d0,0.0d0],[0.0d0,0.0d0],0.0d0,ds,x,3,amp,&
!    0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],[0.0d0,0.0d0,0.0d0,0.0d0],0.0d0,&
!    [0.0d0,0.0d0,0.0d0],dbounds,0,0,edge(3))
!
!    !b(nodese(k)) = b(nodese(k)) + integral/2.0d0
!    !ivec2 = b(nodese(k))
!
!    ivecE = integral/2.0d0
!    jj = (nodese(k)-1)
!    call VecSetValues(PetVec,one,jj,ivecE,ADD_VALUES,ierr)
!
!    enddo
!    enddo
!enddo
!
!END SUBROUTINE PETScAssemblevector
!
!
!!!***************************************************************************************
!!! Subroutine to assemble left hand side FEM matrix (A)
!
!subroutine PETScAssembleADVmatrix(PetMat, points, edges, triangles, np, ne, nt, eq, meanc,omegas, &
!multipliers,sigma,amps,dbounds,casei,term)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!#include <petsc/finclude/petscmat.h>
!
!Mat              PetMat
!PetscInt         ii,jj
!PetscInt         n,one
!PetscScalar      imatT, imatE
!PetscErrorCode   ierr
!
!!! Variable definition & allocation
!integer, dimension(3,ne) :: edges
!integer, dimension(3) :: edge
!integer, dimension(2) :: nodese
!integer, dimension(3,nt) :: triangles
!integer, dimension(3) :: triangle
!integer, dimension(3) :: nodest
!integer :: np, ne, nt, eq, i, j, k, l, term, casei
!
!double precision, dimension(2,np) :: points
!double precision, dimension(4) :: omegas, multipliers
!double precision, dimension(3) :: amps
!double precision, dimension(2,3) :: midpoints, coordt
!double precision, dimension(2) :: gausspoints, dphiv, pt, x, dv, du
!double precision, dimension(2,2) :: Jmat, Ji, Jtraninv, coorde, dbounds
!double precision :: detJ, dx, phis, phiedges, v, u
!double precision :: integral, ds, pe, meanc, sigma
!!!double precision, dimension(np,np) :: A
!
!
!n=np
!one=1
!
!!!initialize matrix
!!!A(:,:) = 0.0d0
!
!!quadrature points
!midpoints(1,:) = [0.5d0, 0.5d0, 0.0d0]
!midpoints(2,:) = [0.0d0, 0.5d0, 0.5d0]
!gausspoints = [(1.0d0-1.0d0/sqrt(3.0d0))/2.0d0, (1.0d0+1.0d0/sqrt(3.0d0))/2.0d0]
!
!! Iterate over all elements (interior domain)
!do i = 1,nt
!! print*, i, ' of ', nt, 'triangles'
!triangle = triangles(:,i)
!
!! Get nodes, coordinates, and domain number
!nodest = triangle(1:3)
!coordt = points(:,nodest)
!
!! Compute Jacobian of map and area of element
!Jmat = 0.0d0
!do j = 1,3
!call dphi(dphiv,j)
!call outer_product(coordt(:,j),dphiv,Ji)
!Jmat = Jmat + Ji
!enddo
!
!! obtain inv(J') used for transforming grad(u) and grad(v)
!! Jtraninv = transpose(Jmat)
!call invert(2,transpose(Jmat),Jtraninv)
!
!detJ = Jmat(1,1)*Jmat(2,2) - Jmat(2,1)*Jmat(1,2)
!dx = 0.5d0*abs(detJ)
!
!! Assemble matrix
!do j=1,3
!pt = midpoints(:,j)
!
!x = 0.0d0
!do k = 1,3
!call phi(phis,k,pt)
!x = x + coordt(:,k)*phis
!enddo
!
!do k = 1,3
!call phi(v,k,pt)
!call dphi(dphiv,k)
!dv = matmul(Jtraninv,dphiv)
!
!do l = 1,3
!call phi(u,l,pt)
!call dphi(dphiv,l)
!du = matmul(Jtraninv,dphiv)
!
!call advdiff(integral,u,v,du,dv,dx,0.0d0,x,eq,0.0d0,&
!meanc,omegas,multipliers,sigma,amps,dbounds,casei,term,0)
!
!!A(nodest(k),nodest(l)) = A(nodest(k),nodest(l)) + integral/3.0d0
!!A(nodest(k),nodest(l)) = integral/3.0d0
!!imatT = A(nodest(k),nodest(l))
!
!ii = (nodest(k)-1)
!jj = (nodest(l)-1)
!imatT = integral/3.0d0
!!if (imat .ne. 0) then
!call MatSetValues(PetMat,one,ii,one,jj,imatT,ADD_VALUES,ierr)
!!end if
!
!enddo
!enddo
!enddo
!enddo
!
!! Iterate over all edges (boundary)
!
!do i = 1,ne
!edge = edges(:,i)
!
!! Get nodes, coordinates, and domain number
!nodese = edge(1:2)
!coorde = points(:,nodese)
!
!! Compute length of edge
!ds = sqrt(dot_product((coorde(:,1) - coorde(:,2)),(coorde(:,1) - coorde(:,2))))
!
!! Assemble matrix
!do j = 1,2
!pe = gausspoints(j)
!
!x(:) = 0.0d0
!do k = 1,2
!call phiedge(phiedges,k,pe)
!x = x + coorde(:,k)*phiedges
!enddo
!
!do k = 1,2
!call phiedge(v,k,pe)
!
!do l = 1,2
!call phiedge(u,l,pe)
!call advdiff(integral,u,v,[0.0d0,0.0d0],[0.0d0,0.0d0],0.0d0,ds,x,eq,0.0d0,&
!meanc,omegas,multipliers,sigma,amps,dbounds,casei,term,edge(3))
!
!!A(nodese(k),nodese(l)) = A(nodese(k),nodese(l)) + integral/2.0d0
!!A(nodese(k),nodese(l)) = integral/2.0d0
!!imatE = A(nodese(k),nodese(l))
!
!ii = (nodese(k)-1)
!jj = (nodese(l)-1)
!imatE = integral/2.0d0
!!if (imat .ne. 0) then
!call MatSetValues(PetMat,one,ii,one,jj,imatE,ADD_VALUES,ierr)
!!end if
!
!enddo
!enddo
!enddo
!enddo
!
!
!END SUBROUTINE PETScAssembleADVmatrix
!
!
!!!***************************************************************************************
!!! Subroutine to assemble left hand side FEM matrix (A)
!
!subroutine PETScAssembleDIFmatrix(PetMat, points, edges, triangles, np, ne, nt, eq, meanc,omegas, &
!multipliers,sigma,amps,dbounds,casei,term)
!
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!#include <petsc/finclude/petscmat.h>
!
!
!Mat              PetMat
!PetscInt         ii,jj
!PetscInt         n,one
!PetscScalar      imatT, imatE
!PetscErrorCode   ierr
!
!!! Variable definition & allocation
!integer, dimension(3,ne) :: edges
!integer, dimension(3) :: edge
!integer, dimension(2) :: nodese
!integer, dimension(3,nt) :: triangles
!integer, dimension(3) :: triangle
!integer, dimension(3) :: nodest
!integer :: np, ne, nt, eq, i, j, k, l, term, casei
!
!!double precision, dimension(np,np) :: A
!double precision, dimension(2,np) :: points
!double precision, dimension(4) :: omegas, multipliers
!double precision, dimension(3) :: amps
!double precision, dimension(2,3) :: midpoints, coordt
!double precision, dimension(2) :: gausspoints, dphiv, pt, x, dv, du
!double precision, dimension(2,2) :: Jmat, Ji, Jtraninv, coorde, dbounds
!double precision :: detJ, dx, phis, phiedges, v, u
!double precision :: integral, ds, pe, meanc, sigma
!
!
!n=np
!one=1
!
!!initialize matrix
!!A(:,:) = 0.0d0
!
!!quadrature points
!midpoints(1,:) = [0.5d0, 0.5d0, 0.0d0]
!midpoints(2,:) = [0.0d0, 0.5d0, 0.5d0]
!gausspoints = [(1.0d0-1.0d0/sqrt(3.0d0))/2.0d0, (1.0d0+1.0d0/sqrt(3.0d0))/2.0d0]
!
!! Iterate over all elements (interior domain)
!do i = 1,nt
!! print*, i, ' of ', nt, 'triangles'
!triangle = triangles(:,i)
!
!! Get nodes, coordinates, and domain number
!nodest = triangle(1:3)
!coordt = points(:,nodest)
!
!! Compute Jacobian of map and area of element
!Jmat = 0.0d0
!do j = 1,3
!call dphi(dphiv,j)
!call outer_product(coordt(:,j),dphiv,Ji)
!Jmat = Jmat + Ji
!enddo
!
!!obtain inv(J') used for transforming grad(u) and grad(v)
!!Jtraninv = transpose(Jmat)
!call invert(2,transpose(Jmat),Jtraninv)
!
!detJ = Jmat(1,1)*Jmat(2,2) - Jmat(2,1)*Jmat(1,2)
!dx = 0.5d0*abs(detJ)
!
!! Assemble matrix
!do j=1,3
!pt = midpoints(:,j)
!
!x = 0.0d0
!do k = 1,3
!call phi(phis,k,pt)
!x = x + coordt(:,k)*phis
!enddo
!
!do k = 1,3
!call phi(v,k,pt)
!call dphi(dphiv,k)
!dv = matmul(Jtraninv,dphiv)
!
!do l = 1,3
!call phi(u,l,pt)
!call dphi(dphiv,l)
!du = matmul(Jtraninv,dphiv)
!
!call advdiff(integral,u,v,du,dv,dx,0.0d0,x,eq,0.0d0,&
!meanc,omegas,multipliers,sigma,amps,dbounds,casei,term,0)
!!A(nodest(k),nodest(l)) = A(nodest(k),nodest(l)) + integral/3.0d0
!!A(nodest(k),nodest(l)) = integral/3.0d0
!!imatT = A(nodest(k),nodest(l))
!
!ii = (nodest(k)-1)
!jj = (nodest(l)-1)
!imatT = integral/3.0d0
!!if (imat .ne. 0) then
!call MatSetValues(PetMat,one,ii,one,jj,imatT,ADD_VALUES,ierr)
!!end if
!
!enddo
!enddo
!enddo
!enddo
!
!! Iterate over all edges (boundary)
!
!do i = 1,ne
!edge = edges(:,i)
!
!! Get nodes, coordinates, and domain number
!nodese = edge(1:2)
!coorde = points(:,nodese)
!
!! Compute length of edge
!ds = sqrt(dot_product((coorde(:,1) - coorde(:,2)),(coorde(:,1) - coorde(:,2))))
!
!! Assemble matrix
!do j = 1,2
!pe = gausspoints(j)
!
!x(:) = 0.0d0
!do k = 1,2
!call phiedge(phiedges,k,pe)
!x = x + coorde(:,k)*phiedges
!enddo
!
!do k = 1,2
!call phiedge(v,k,pe)
!
!do l = 1,2
!call phiedge(u,l,pe)
!call advdiff(integral,u,v,[0.0d0,0.0d0],[0.0d0,0.0d0],0.0d0,ds,x,eq,0.0d0,&
!meanc,omegas,multipliers,sigma,amps,dbounds,casei,term,edge(3))
!
!!A(nodese(k),nodese(l)) = A(nodese(k),nodese(l)) + integral/2.0d0
!!A(nodese(k),nodese(l)) = integral/2.0d0
!!imatE = A(nodese(k),nodese(l))
!
!ii = (nodese(k)-1)
!jj = (nodese(l)-1)
!imatE = integral/2.0d0
!!if (imat .ne. 0) then
!call MatSetValues(PetMat,one,ii,one,jj,imatE,ADD_VALUES,ierr)
!!end if
!
!enddo
!enddo
!enddo
!enddo
!
!
!END SUBROUTINE PETScAssembleDIFmatrix


!!***************************************************************************************
!! Subroutine for Basis functions and gradiend of the basis functions
   subroutine phi(y,indx,x)!Basis functions on the reference element

      double precision, dimension(2) :: x
      double precision :: y
      integer :: indx

      select case (indx)
      case (1)
         y = 1 - x(1) - x(2)
      case (2)
         y = x(1)
      case (3)
         y = x(2)
      end select

   end subroutine phi

   subroutine phiedge(y,indx,x)!Basis functions on the reference edge

      double precision :: y, x
      integer :: indx

      select case (indx)
      case (1)
         y = 1 - x
      case (2)
         y = x
      end select

   end subroutine phiedge

   subroutine dphi(y,indx)!Gradients of basis functions on the reference element

      double precision, dimension(2) :: y
      integer :: indx

      select case (indx)
      case (1)
         y = [-1.0d0, -1.0d0]
      case (2)
         y = [1.0d0, 0.0d0]
      case (3)
         y = [0.0d0, 1.0d0]
      end select

   end subroutine dphi


end module PETScAssembly
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%** END **%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
