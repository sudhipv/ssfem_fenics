!! Assembly Module: By Ajit Desai, 8/2015, Mohammad Khalil, 1/2012
!! purpose: This is the assembly code addapted from FEniCs/Puffin FEM solver for : Au = f
!! It assembles the left hand side FEM matrix (A) and right hand side FEM vector (f)
!
! input:
!      points: nodal coordinates (x , y for 2D)
!      edges : 1D nodal connectivity (physical boundary)
!      triangles: 2D nodel connectivity (physical surface)
!      np    : number of points 
!      ne    : number of edges
!      nt    : number of triangles
!      ndim  : random dimension of KLE/PCE expansion used for input represenations 
!      eq    : equation type : diffusion or convection-diffusion 
!      meanc : 
!      omegas: w's of KLE expansion terms : Rogers Book
!      multipliers : √λ / √(a1 = (sin(2*w_i a1)/2*w_i a_1)): KLE expansion terms : Rogers Book
!      sigma : σ for underlying Guassian process 
!      amps  : amp : 
!      dbounds : physical domain limits : x=[0, 1], y=[0,1] : for unit square 
!      casei : to select the case: 1: deterministi, 2: Stochastic MC 3: Stochastic PCE for Lognormal 
!      term  : to select terms for PCE expansion for Lognormal Process  
!
! output:
!      Amat/A: FEM assembled matrix  
!      b/f   : FEM assembled vector 
! NOTE : for theory refer Puffin/FEniCs codes & documentaion "http://fenicsproject.org"
!
!!---------------------------------------------------------------------------------------
module assembly

!! This Module dependends on following Modules 
use common
use myCommon
use variationalform

implicit none

contains


!!***************************************************************************************
!! Subroutine to assemble left hand side FEM matrix (A)
subroutine SubAssembleMatrix(Adii,Adgg,Adgi,Adri,Adci,Adcr,Adcc,Adrr,points,edges,triangles, &
                             np,ne,nt,nb,nci,ndim,npceout,nomga,eq,meanc,omegas,&
                             multipliers,sigma,amps,dbounds,mIndex,sIndex,casei,term) !mIndex

    integer, dimension(3,ne) :: edges
    integer, dimension(3)    :: edge
    integer, dimension(2)    :: nodese
    integer, dimension(3,nt) :: triangles
    integer, dimension(3)    :: triangle
    integer, dimension(3)    :: nodest
    integer :: np, ne, nt, eq, i, j, k, l, term
    integer :: casei,ndim,npceout,nomga

    !!double precision, dimension(np,np) :: A
    double precision, dimension(2,np) :: points
    double precision, dimension(nomga):: omegas, multipliers
    double precision, dimension(ndim) :: amps
    double precision, dimension(2,3)  :: midpoints, coordt
    double precision, dimension(2)    :: gausspoints, dphiv, pt, x, dv, du
    double precision, dimension(2,2)  :: Jmat, Ji, Jtraninv, coorde, dbounds
    double precision :: meanc, detJ, dx, phis, phiedges, v, u, integral, ds, pe, sigma

    !! NEW : Feb 2016
    integer :: nb, ni, nci, nc, ii, jj
    double precision :: intT, intE
    double precision, dimension((np-nb),(np-nb))  :: Adii
    double precision, dimension(nb,nb)            :: Adgg
    !double precision, dimension((np-nb),nb)       :: Adig
    double precision, dimension(nb,(np-nb))       :: Adgi
    !double precision, dimension((np-nb),(nb-nci)) :: Adir
    double precision, dimension((nb-nci),(np-nb)) :: Adri
    !double precision, dimension((np-nb),nci)      :: Adic
    double precision, dimension(nci,(np-nb))      :: Adci
    !double precision, dimension((nb-nci),nci)     :: Adrc
    double precision, dimension(nci,(nb-nci))     :: Adcr
    double precision, dimension(nci,nci)          :: Adcc
    double precision, dimension((nb-nci),(nb-nci)):: Adrr

!! UnderStudy: July 7, 2016
integer, dimension(ndim,npceout) :: mIndex
integer, dimension(2,ndim) :: sIndex



    !! initialize matrix
    !! A(:,:) = 0.0d0

    !! NEW : Feb 2016
    ni = np-nb
    nc = np-nci
    intT = 0.0d0
    intE = 0.0d0

    !! quadrature points
    midpoints(1,:) = [0.5d0, 0.5d0, 0.0d0]
    midpoints(2,:) = [0.0d0, 0.5d0, 0.5d0]
    gausspoints = [(1.0d0-1.0d0/sqrt(3.0d0))/2.0d0, (1.0d0+1.0d0/sqrt(3.0d0))/2.0d0]

    !! Iterate over all elements (interior domain)
    do i = 1,nt

        ! print*, i, ' of ', nt, 'triangles'
        triangle = triangles(:,i)

        ! Get nodes, coordinates, and domain number
        nodest = triangle(1:3)
        coordt = points(:,nodest)

        ! Compute Jacobian of map and area of element
        Jmat = 0.0d0
        do j = 1,3
        call dphi(dphiv,j)
        call outer_product(coordt(:,j),dphiv,Ji)
        Jmat = Jmat + Ji
        enddo

        !obtain inv(J') used for transforming grad(u) and grad(v)
        !Jtraninv = transpose(Jmat)
        call invert(2,transpose(Jmat),Jtraninv)

        detJ = Jmat(1,1)*Jmat(2,2) - Jmat(2,1)*Jmat(1,2)
        dx = 0.5d0*abs(detJ)

        ! Assemble matrix
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

        do l = 1,3
        call phi(u,l,pt)
        call dphi(dphiv,l)
        du = matmul(Jtraninv,dphiv)

        call advdiff(integral,u,v,du,dv,dx,0.0d0,x,eq,0.0d0,meanc,ndim,npceout,nomga, &
                   omegas,multipliers,sigma,amps,dbounds,casei,term,mIndex,sIndex,0) !mIndex

        !!A(nodest(k),nodest(l)) = A(nodest(k),nodest(l)) + integral/3.0d0

        !! NEW : Feb 2016
        intT = integral/3.0d0
        ii = nodest(k)
        jj = nodest(l)
        if (ii <= ni .AND. jj <= ni) then
            Adii(ii,jj) = Adii(ii,jj) + intT
        else if (ii > ni .AND. jj > ni) then
            Adgg((ii-ni),(jj-ni)) = Adgg((ii-ni),(jj-ni)) + intT
!        else if (ii <= ni .AND. jj > ni) then
!            Adig(ii,(jj-ni)) = Adig(ii,(jj-ni)) + intT
        else if (ii > ni .AND. jj <= ni) then
            Adgi((ii-ni),jj) = Adgi((ii-ni),jj) + intT
        end if


        if (ii > ni .AND. ii <= nc .AND. jj <= ni) then
            Adri((ii-ni),jj) = Adri((ii-ni),jj) + intT
!        else if (ii <= ni .AND. jj > ni .AND. jj <= nc) then
!            Adir(ii,(jj-ni)) = Adir(ii,(jj-ni)) + intT
!        else if (ii <= ni .AND. jj > nc) then
!            Adic(ii,(jj-nc)) = Adic(ii,(jj-nc)) + intT
        else if (ii > nc .AND. jj <= ni) then
            Adci((ii-nc),jj) = Adci((ii-nc),jj) + intT
!        else if (ii > ni .AND. ii <= nc .AND. jj > nc) then
!            Adrc((ii-ni),(jj-nc)) = Adrc((ii-ni),(jj-nc)) + intT
        else if (ii > nc .AND. jj > ni .AND. jj <= nc) then
            Adcr((ii-nc),(jj-ni)) = Adcr((ii-nc),(jj-ni)) + intT
        else if (ii > nc .AND. jj > nc ) then
            Adcc((ii-nc),(jj-nc)) = Adcc((ii-nc),(jj-nc)) + intT
        else if (ii > ni .AND. ii <= nc .AND. jj > ni .AND. jj <= nc) then
            Adrr((ii-ni),(jj-ni)) = Adrr((ii-ni),(jj-ni)) + intT
        end if

        enddo
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

        ! Assemble matrix
        do j = 1,2
        pe = gausspoints(j)

        x(:) = 0.0d0
        do k = 1,2
        call phiedge(phiedges,k,pe)
        x = x + coorde(:,k)*phiedges
        enddo

        do k = 1,2
        call phiedge(v,k,pe)

        do l = 1,2
        call phiedge(u,l,pe)
        call advdiff(integral,u,v,[0.0d0,0.0d0],[0.0d0,0.0d0],0.0d0,ds,x,eq,0.0d0,meanc,ndim,npceout,nomga, &
                    omegas,multipliers,sigma,amps,dbounds,casei,term,mIndex,sIndex,edge(3)) !mIndex

        !!A(nodese(k),nodese(l)) = A(nodese(k),nodese(l)) + integral/2.0d0

        !! NEW : Feb 2016
        intE = integral/2.0d0
        ii = nodese(k)
        jj = nodese(l)
        if (ii <= ni .AND. jj <= ni) then
            Adii(ii,jj) = Adii(ii,jj) + intE
        else if (ii > ni .AND. jj > ni) then
            Adgg((ii-ni),(jj-ni)) = Adgg((ii-ni),(jj-ni)) + intE
!        else if (ii <= ni .AND. jj > ni) then
!            Adig(ii,(jj-ni)) = Adig(ii,(jj-ni)) + intE
        else if (ii > ni .AND. jj <= ni) then
            Adgi((ii-ni),jj) = Adgi((ii-ni),jj) + intE
        end if


        if (ii > ni .AND. ii <= nc .AND. jj <= ni) then
            Adri((ii-ni),jj) = Adri((ii-ni),jj) + intE
!        else if (ii <= ni .AND. jj > ni .AND. jj <= nc) then
!            Adir(ii,(jj-ni)) = Adir(ii,(jj-ni)) + intE
!        else if (ii <= ni .AND. jj > nc) then
!            Adic(ii,(jj-nc)) = Adic(ii,(jj-nc)) + intE
        else if (ii > nc .AND. jj <= ni) then
            Adci((ii-nc),jj) = Adci((ii-nc),jj) + intE
!        else if (ii > ni .AND. ii <= nc .AND. jj > nc) then
!            Adrc((ii-ni),(jj-nc)) = Adrc((ii-ni),(jj-nc)) + intE
        else if (ii > nc .AND. jj > ni .AND. jj <= nc) then
            Adcr((ii-nc),(jj-ni)) = Adcr((ii-nc),(jj-ni)) + intE
        else if (ii > nc .AND. jj > nc ) then
            Adcc((ii-nc),(jj-nc)) = Adcc((ii-nc),(jj-nc)) + intE
        else if (ii > ni .AND. ii <= nc .AND. jj > ni .AND. jj <= nc) then
            Adrr((ii-ni),(jj-ni)) = Adrr((ii-ni),(jj-ni)) + intE
        end if

        enddo
        enddo
        enddo
    enddo

end subroutine SubAssembleMatrix


!!!***************************************************************************************
!!! Subroutine to assemble right hand side force vector (f)
!subroutine SubAssembleVector(Fsi, Fsg, points, edges, triangles, np, ne, nt, nb, npceout, amp, dbounds)
!
!    !!double precision, dimension(np)   :: b
!    double precision, dimension(2,np) :: points
!    integer, dimension(3,ne) :: edges
!    integer, dimension(3)    :: edge
!    integer, dimension(2)    :: nodese
!    integer, dimension(3,nt) :: triangles
!    integer, dimension(3)    :: triangle
!    integer, dimension(3)    :: nodest
!    integer :: np, ne, nt, i, j, k, ndim
!    double precision :: detJ, dx, phis, phiedges, v, integral, ds, pe, amp
!    double precision, dimension(2,3) :: midpoints, coordt
!    double precision, dimension(2)   :: gausspoints, dphiv, pt, x, dv
!    double precision, dimension(2,2) :: Jmat, Ji, Jtraninv, coorde, dbounds
!
!    !! NEW : Feb 2016
!    integer :: nb, ni, ii, jj, npceout
!    double precision :: intT, intE
!    double precision, dimension((np-nb)*npceout)  :: Fsi
!    double precision, dimension(nb*npceout)       :: Fsg
!
!    !initialize vector
!    !!b(:) = 0.0d0
!    intT = 0.0d0
!    intE = 0.0d0
!    ni = np-nb
!
!    !quadrature points
!    midpoints(1,:) = [0.5d0, 0.5d0, 0.0d0]
!    midpoints(2,:) = [0.0d0, 0.5d0, 0.5d0]
!    gausspoints = [(1.0d0-1.0d0/sqrt(3.0d0))/2.0d0, (1.0d0+1.0d0/sqrt(3.0d0))/2.0d0]
!
!    ! Iterate over all elements (interior domain)
!    do i = 1,nt
!        triangle = triangles(:,i)
!
!        ! Get nodes, coordinates, and domain number
!        nodest = triangle(1:3)
!        coordt = points(:,nodest)
!
!        ! Compute Jacobian of map and area of element
!        Jmat(:,:) = 0.0d0
!        do j = 1,3
!        call dphi(dphiv,j)
!        call outer_product(coordt(:,j),dphiv,Ji)
!        Jmat = Jmat + Ji
!        enddo
!
!        !obtain inv(J')
!        call invert(2,transpose(Jmat),Jtraninv)
!
!        detJ = Jmat(1,1)*Jmat(2,2) - Jmat(2,1)*Jmat(1,2)
!        dx = 0.5d0*abs(detJ)
!
!        ! Assemble vector
!        do j=1,3
!        pt = midpoints(:,j)
!
!        x = 0.0d0
!        do k = 1,3
!        call phi(phis,k,pt)
!        x = x + coordt(:,k)*phis
!        enddo
!
!        do k = 1,3
!        call phi(v,k,pt)
!        call dphi(dphiv,k)
!        dv = matmul(Jtraninv,dphiv)
!
!        call advdiff(integral,0.0d0,v,[0.0d0,0.0d0],dv,dx,0.0d0,x,3,amp,&
!        0.0d0,ndim,[0.0d0,0.0d0,0.0d0,0.0d0],[0.0d0,0.0d0,0.0d0,0.0d0],0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],dbounds,0,0,0)
!
!        !!b(nodest(k)) = b(nodest(k)) + integral/3.0d0
!
!        !! NEW : Feb 2016
!        intT = integral/3.0d0
!        ii = (nodest(k))
!        if (ii <= ni ) then
!        Fsi(ii) = Fsi(ii) + intT
!        else if (ii > ni) then
!        Fsg(ii-ni) = Fsg(ii-ni) + intT
!        end if
!
!        enddo
!        enddo
!    enddo
!
!    ! Iterate over all edges (boundary)
!    do i = 1,ne
!        edge = edges(:,i)
!
!        ! Get nodes, coordinates, and domain number
!        nodese = edge(1:2)
!        coorde = points(:,nodese)
!
!        ! Compute length of edge
!        ds = sqrt(dot_product((coorde(:,1) - coorde(:,2)),(coorde(:,1) - coorde(:,2))))
!
!        ! Assemble vector
!        do j = 1,2
!        pe = gausspoints(j)
!
!        x = 0.0d0
!        do k = 1,2
!        call phiedge(phiedges,k,pe)
!        x = x + coorde(:,k)*phiedges
!        enddo
!
!        do k = 1,2
!        call phiedge(v,k,pe)
!
!        call advdiff(integral,0.0d0,v,[0.0d0,0.0d0],[0.0d0,0.0d0],0.0d0,ds,x,3,amp,&
!        0.0d0,ndim,[0.0d0,0.0d0,0.0d0,0.0d0],[0.0d0,0.0d0,0.0d0,0.0d0],0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],dbounds,0,0,edge(3))
!
!        !!b(nodese(k)) = b(nodese(k)) + integral/2.0d0
!
!        !! NEW : Feb 2016
!        intE = integral/2.0d0
!        ii = (nodese(k))
!        if (ii <= ni ) then
!        Fsi(ii) = Fsi(ii) + intE
!        else if (ii > ni) then
!        Fsg(ii-ni) = Fsg(ii-ni) + intE
!        end if
!
!        enddo
!        enddo
!    enddo
!
!end subroutine SubAssembleVector
!
!
!
!!!***************************************************************************************
!!! Subroutine to assemble left hand side FEM matrix (A)
!subroutine assemblematrix(A,points,edges,triangles,np,ne,nt, ndim, eq, &
!                          meanc,omegas,multipliers,sigma,amps,dbounds,casei,term)
!
!    integer, dimension(3,ne) :: edges
!    integer, dimension(3) :: edge
!    integer, dimension(2) :: nodese
!    integer, dimension(3,nt) :: triangles
!    integer, dimension(3) :: triangle
!    integer, dimension(3) :: nodest
!    integer :: np, ne, nt, eq, i, j, k, l, term,casei,ndim
!
!    double precision, dimension(np,np) :: A
!    double precision, dimension(2,np) :: points
!    double precision, dimension(ndim) :: omegas, multipliers
!    double precision, dimension(ndim) :: amps
!    double precision :: meanc, detJ, dx, phis, phiedges, v, u, integral, ds, pe, sigma
!    double precision, dimension(2,3) :: midpoints, coordt
!    double precision, dimension(2) :: gausspoints, dphiv, pt, x, dv, du
!    double precision, dimension(2,2) :: Jmat, Ji, Jtraninv, coorde, dbounds
!
!
!      !initialize matrix
!      A(:,:) = 0.0d0
!
!      !quadrature points
!      midpoints(1,:) = [0.5d0, 0.5d0, 0.0d0]
!      midpoints(2,:) = [0.0d0, 0.5d0, 0.5d0]
!      gausspoints = [(1.0d0-1.0d0/sqrt(3.0d0))/2.0d0, (1.0d0+1.0d0/sqrt(3.0d0))/2.0d0]
!
!      ! Iterate over all elements (interior domain)
!      do i = 1,nt
!      !     print*, i, ' of ', nt, 'triangles'
!         triangle = triangles(:,i)
!
!         ! Get nodes, coordinates, and domain number
!         nodest = triangle(1:3)
!         coordt = points(:,nodest)
!
!         ! Compute Jacobian of map and area of element
!         Jmat = 0.0d0
!         do j = 1,3
!            call dphi(dphiv,j)
!            call outer_product(coordt(:,j),dphiv,Ji)
!            Jmat = Jmat + Ji
!         enddo
!
!         !obtain inv(J') used for transforming grad(u) and grad(v)
!         !Jtraninv = transpose(Jmat)
!         call invert(2,transpose(Jmat),Jtraninv)
!
!         detJ = Jmat(1,1)*Jmat(2,2) - Jmat(2,1)*Jmat(1,2)
!         dx = 0.5d0*abs(detJ)
!
!         ! Assemble matrix
!         do j=1,3
!            pt = midpoints(:,j)
!
!            x = 0.0d0
!            do k = 1,3
!               call phi(phis,k,pt)
!               x = x + coordt(:,k)*phis
!            enddo
!
!            do k = 1,3
!               call phi(v,k,pt)
!               call dphi(dphiv,k)
!               dv = matmul(Jtraninv,dphiv)
!
!               do l = 1,3
!                  call phi(u,l,pt)
!                  call dphi(dphiv,l)
!                  du = matmul(Jtraninv,dphiv)
!
!                  call advdiff(integral,u,v,du,dv,dx,0.0d0,x,eq,0.0d0,&
!                     meanc,ndim,omegas,multipliers,sigma,amps,dbounds,casei,term,0)
!                  A(nodest(k),nodest(l)) = A(nodest(k),nodest(l)) + integral/3.0d0
!
!               enddo
!            enddo
!         enddo
!      enddo
!
!      ! Iterate over all edges (boundary)
!
!      do i = 1,ne
!         edge = edges(:,i)
!
!         ! Get nodes, coordinates, and domain number
!         nodese = edge(1:2)
!         coorde = points(:,nodese)
!
!         ! Compute length of edge
!         ds = sqrt(dot_product((coorde(:,1) - coorde(:,2)),(coorde(:,1) - coorde(:,2))))
!
!         ! Assemble matrix
!         do j = 1,2
!            pe = gausspoints(j)
!
!            x(:) = 0.0d0
!            do k = 1,2
!               call phiedge(phiedges,k,pe)
!               x = x + coorde(:,k)*phiedges
!            enddo
!         
!            do k = 1,2
!               call phiedge(v,k,pe)
!
!               do l = 1,2
!                  call phiedge(u,l,pe)
!                  call advdiff(integral,u,v,[0.0d0,0.0d0],[0.0d0,0.0d0],0.0d0,ds,x,eq,0.0d0,&
!                     meanc,ndim,omegas,multipliers,sigma,amps,dbounds,casei,term,edge(3))
!                  A(nodese(k),nodese(l)) = A(nodese(k),nodese(l)) + integral/2.0d0
!
!               enddo
!            enddo
!         enddo
!      enddo
!
!   end subroutine assemblematrix
!   
!
!!!***************************************************************************************
!!! Subroutine to assemble right hand side force vector (f)
!   subroutine assemblevector(b, points, edges, triangles, np, ne, nt, amp, dbounds)
!
!      double precision, dimension(np) :: b
!      double precision, dimension(2,np) :: points
!      integer, dimension(3,ne) :: edges
!      integer, dimension(3) :: edge
!      integer, dimension(2) :: nodese
!      integer, dimension(3,nt) :: triangles
!      integer, dimension(3) :: triangle
!      integer, dimension(3) :: nodest
!      integer :: np, ne, nt, i, j, k, ndim
!      double precision :: detJ, dx, phis, phiedges, v, integral, ds, pe, amp
!      double precision, dimension(2,3) :: midpoints, coordt
!      double precision, dimension(2) :: gausspoints, dphiv, pt, x, dv
!      double precision, dimension(2,2) :: Jmat, Ji, Jtraninv, coorde, dbounds
!
!      !initialize vector
!      b(:) = 0.0d0
!
!      !quadrature points
!      midpoints(1,:) = [0.5d0, 0.5d0, 0.0d0]
!      midpoints(2,:) = [0.0d0, 0.5d0, 0.5d0]
!      gausspoints = [(1.0d0-1.0d0/sqrt(3.0d0))/2.0d0, (1.0d0+1.0d0/sqrt(3.0d0))/2.0d0]
!
!      ! Iterate over all elements (interior domain)
!      do i = 1,nt
!         triangle = triangles(:,i)
!
!         ! Get nodes, coordinates, and domain number
!         nodest = triangle(1:3)
!         coordt = points(:,nodest)
!
!         ! Compute Jacobian of map and area of element
!         Jmat(:,:) = 0.0d0
!         do j = 1,3
!            call dphi(dphiv,j)
!            call outer_product(coordt(:,j),dphiv,Ji)
!            Jmat = Jmat + Ji
!         enddo
!
!         !obtain inv(J')
!         call invert(2,transpose(Jmat),Jtraninv)
!
!         detJ = Jmat(1,1)*Jmat(2,2) - Jmat(2,1)*Jmat(1,2)
!         dx = 0.5d0*abs(detJ)
!
!         ! Assemble vector
!         do j=1,3
!            pt = midpoints(:,j)
!
!            x = 0.0d0
!            do k = 1,3
!               call phi(phis,k,pt)
!               x = x + coordt(:,k)*phis
!            enddo
!
!            do k = 1,3
!               call phi(v,k,pt)
!               call dphi(dphiv,k)
!               dv = matmul(Jtraninv,dphiv)
!
!               call advdiff(integral,0.0d0,v,[0.0d0,0.0d0],dv,dx,0.0d0,x,3,amp,&
!               0.0d0,ndim,[0.0d0,0.0d0,0.0d0,0.0d0],[0.0d0,0.0d0,0.0d0,0.0d0],0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],dbounds,0,0,0)
!               b(nodest(k)) = b(nodest(k)) + integral/3.0d0
!            enddo
!         enddo
!      enddo
!
!      ! Iterate over all edges (boundary)
!      do i = 1,ne
!         edge = edges(:,i)
!
!         ! Get nodes, coordinates, and domain number
!         nodese = edge(1:2)
!         coorde = points(:,nodese)
!
!         ! Compute length of edge
!         ds = sqrt(dot_product((coorde(:,1) - coorde(:,2)),(coorde(:,1) - coorde(:,2))))
!
!         ! Assemble vector
!         do j = 1,2
!            pe = gausspoints(j)
!
!            x = 0.0d0
!            do k = 1,2
!               call phiedge(phiedges,k,pe)
!               x = x + coorde(:,k)*phiedges
!            enddo
!         
!            do k = 1,2
!               call phiedge(v,k,pe)
!
!               call advdiff(integral,0.0d0,v,[0.0d0,0.0d0],[0.0d0,0.0d0],0.0d0,ds,x,3,amp,&
!               0.0d0,ndim,[0.0d0,0.0d0,0.0d0,0.0d0],[0.0d0,0.0d0,0.0d0,0.0d0],0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],dbounds,0,0,edge(3))
!               b(nodese(k)) = b(nodese(k)) + integral/2.0d0
!            enddo
!         enddo
!      enddo
!
!   end subroutine assemblevector

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

end module assembly
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%