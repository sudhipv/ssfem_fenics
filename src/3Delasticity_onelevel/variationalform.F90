!! Variational Formulation Module : By Ajit Desai, 8/2015, Mohammad Khalil, 1/2012
!! purpose: is to define Stochast advection-diffusion differential equation in weak form 
!
!!---------------------------------------------------------------------------------------
module variationalform

use myCommon

implicit none

contains

!!***************************************************************************************
!! Subroutine for definin advection-diffusion differential equaiton in weak form
   subroutine advdiff(integral,u,v,du,dv,dx,ds,x,eq,amp,meanc,ndim,npceout,nomga, &
               omegas,multipliers,sigma,amps,dbounds,casei,term,mIndex,sIndex,ident)

    double precision :: integral, u, v, dx, ds, amp, meanc, gs, gds, fs, gns, cs, sigma
    integer :: eq, term, casei, ident, ndim, npceout, nomga
    double precision, dimension(nomga):: omegas, multipliers
    double precision, dimension(ndim) :: amps
    double precision, dimension(2)    :: du, dv, x, bv
    double precision, dimension(2,2)  :: dbounds

    integer, dimension(ndim,npceout)  :: mIndex
    integer, dimension(2,ndim)        :: sIndex

    SELECT CASE (eq)
        CASE (1)
            call b(bv,x)
             bv(:) = 0.0d0
            call g(gs,x,dbounds,ident)
            integral = dot_product(bv,du)*v*dx + gs*u*v*ds
        CASE (2)
            call c(cs,x,ndim,npceout,nomga,omegas,multipliers,sigma,amps,casei,meanc,dbounds,term,mIndex,sIndex)
            integral = cs*dot_product(du,dv)*dx
        CASE (3)
            call g(gs,x,dbounds,ident)
            call gd(gds,x,ident)
            call f(fs,x,amp)
            call gn(gns,x,dbounds,ident)
            integral = fs*v*dx + (gs*gds - gns)*v*ds
    END SELECT

   end subroutine advdiff

!!***************************************************************************************
!! Advection/convection term 
   subroutine b(y,x) 

      double precision, dimension(2) :: x,y

      y(1) = 0.0d0                           !!y(1) = 1.0d0
      y(2) = 0.0d0

   end subroutine b

!!***************************************************************************************
!! Conductivit or Penalty factor 
   subroutine g(y,x,dbounds,ident) 

      integer :: ident
      double precision, dimension(2) :: x
      double precision :: y
      double precision, dimension(2,2) :: dbounds

     !! Specify the boundaries on which condition should be imposed 
     if ((ident .eq. 1) .or. (ident .eq. 2) .or. (ident .eq. 3) .or. (ident .eq. 4)) then
         y = 1.0d50
      else
         y = 0.0d0
      endif

   end subroutine g

!!***************************************************************************************
!! Dirichlet boundary condition
   subroutine gd(y,x,ident) 
    
      integer :: ident
      double precision, dimension(2) :: x
      double precision :: y
      !! if (ident .ge. 5) then
      !!     y = 0.0d0
      !!  else
         y = 0.0d0
      !!  end if

   end subroutine gd

!!***************************************************************************************
!! Neumann boundary condition
   subroutine gn(y,x,dbounds,ident) 

      integer :: ident
      double precision, dimension(2) :: x
      double precision :: y !! eps
      double precision, dimension(2,2) :: dbounds

      !! if (ident .eq. 2) then
      !!    y = -0.6d0*exp(-100.0d0*(x(2)-1.0d0)**2)
      !! else
      !!   y = 0.0d0
      !! endif
         y = 0.0d0

   end subroutine gn

!!***************************************************************************************
!! Right-hand side, source term
   subroutine f(y, x, amp) 

      double precision, dimension(2) :: x, centre
      double precision :: y, amp

      centre(1) = 0.5d0
      centre(2) = 0.5d0
      
      !! point force at the center of the plate
      !! y = amp*exp(-2000.0d0*dot_product((x-centre),(x-centre)))
      y = 1.0d0

   end subroutine f

!!***************************************************************************************
!! Random Diffusion coefficient
   subroutine c(y, x, ndim,npceout,nomga,omegas,multipliers,sigma,xi,casei,meanc,dbounds,term,mIndex,sIndex)
            !!c(cs,x, omegas, multipliers, sigma,amps,casei, meanc, dbounds, term)

    integer :: casei, term, ndim, npceout, nomga
    double precision, dimension(2)    :: x, xdash,offsets
    double precision, dimension(nomga):: omegas, multipliers
    double precision, dimension(ndim) :: xi, g
    double precision, dimension(2,2)  :: dbounds
    double precision, dimension(20)   :: lterms
    double precision :: sigma, y, meanc
    !!integer, dimension(2,ndim) :: lmdindex

    integer :: i, j, index
    double precision :: newY, yy, factorial

    integer :: Xindex, Yindex
    double precision :: gg1, gg2

    integer, dimension(ndim,npceout)  :: mIndex
    integer, dimension(2,ndim)        :: sIndex

    SELECT CASE (casei)
    !!---------------------------------------------------------------------------------------
    !! Deterministic and random variable cases
    !! CASE (1)
    !! y = meanc

    !!---------------------------------------------------------------------------------------
    !! MC with trunctaed PCE of lognormal process
    CASE (1)
    offsets(1) = sum(dbounds(1,:))/2.0d0
    offsets(2) = sum(dbounds(2,:))/2.0d0
    xdash = x-offsets

!!print*, "offsets = ", xdash
!!print*, "-------------------------------------------------------"

    ! g(1) = sigma*multipliers(1)*multipliers(1)*(cos(omegas(1)*xdash(1))*cos(omegas(1)*xdash(2)))
    ! g(2) = sigma*multipliers(1)*multipliers(2)*(1.0d0/sqrt(2.0d0))*(cos(omegas(1)*xdash(1))*&
    !    sin(omegas(2)*xdash(2))+sin(omegas(2)*xdash(1))*cos(omegas(1)*xdash(2)))
    ! g(3) = sigma*multipliers(1)*multipliers(3)*(1.0d0/sqrt(2.0d0))*(cos(omegas(1)*xdash(1))*&
    !    cos(omegas(3)*xdash(2))+cos(omegas(3)*xdash(1))*cos(omegas(1)*xdash(2)))

    g(1) = sigma*multipliers(1)*multipliers(1)*(cos(omegas(1)*xdash(1))*cos(omegas(1)*xdash(2)))
    g(2) = sigma*multipliers(1)*multipliers(2)*(cos(omegas(1)*xdash(1))*sin(omegas(2)*xdash(2)))
    g(3) = sigma*multipliers(1)*multipliers(2)*(sin(omegas(2)*xdash(1))*cos(omegas(1)*xdash(2)))

    lterms(1) = 1.0d0
    lterms(2) = xi(1)*g(1)
    lterms(3) = xi(2)*g(2)
    lterms(4) = xi(3)*g(3)
    lterms(5) = (xi(1)**2 - 1.0d0)*g(1)**2/2.0d0
    lterms(6) = (xi(2)**2 - 1.0d0)*g(2)**2/2.0d0
    lterms(7) = (xi(3)**2 - 1.0d0)*g(3)**2/2.0d0
    lterms(8) = xi(1)*xi(2)*g(1)*g(2)
    lterms(9) = xi(1)*xi(3)*g(1)*g(3)
    lterms(10) = xi(2)*xi(3)*g(2)*g(3)
    lterms(11) = (xi(1)**3 - 3.0d0*xi(1))*g(1)**3/6.0d0
    lterms(12) = (xi(2)**3 - 3.0d0*xi(2))*g(2)**3/6.0d0
    lterms(13) = (xi(3)**3 - 3.0d0*xi(3))*g(3)**3/6.0d0
    lterms(14) = (xi(1)*xi(1)*xi(2) - xi(2))*g(1)**2*g(2)/2.0d0
    lterms(15) = (xi(1)*xi(1)*xi(3) - xi(3))*g(1)**2*g(3)/2.0d0
    lterms(16) = (xi(2)*xi(2)*xi(1) - xi(1))*g(2)**2*g(1)/2.0d0
    lterms(17) = (xi(2)*xi(2)*xi(3) - xi(3))*g(2)**2*g(3)/2.0d0
    lterms(18) = (xi(3)*xi(3)*xi(1) - xi(1))*g(3)**2*g(1)/2.0d0
    lterms(19) = (xi(3)*xi(3)*xi(2) - xi(2))*g(3)**2*g(2)/2.0d0
    lterms(20) = xi(1)*xi(2)*xi(3)*g(1)*g(2)*g(3)
    y = meanc*sum(lterms)
    !! ! for exact lognormal process
    !! y = meanc*exp(-0.5d0*sum(g**2))*exp(xi(1)*g(1) + xi(2)*g(2) + xi(3)*g(3)) 


    !!---------------------------------------------------------------------------------------
    !! Trunctaed PCE of lognormal process
    !! UnderStudy: July 14, 2016 : General case L=n, O=m

    CASE (2) !! UnderStudy: July 14, 2016 :
    offsets(1) = sum(dbounds(1,:))/2.0d0     !! It can be calculated upfront of defined cost
    offsets(2) = sum(dbounds(2,:))/2.0d0     !! make changes to make it offset(1),(2) = 0.5
    xdash = x-offsets

    !! obtain KLE terms
    do i = 1,ndim
        Xindex = sIndex(1,i)
        if (MOD(Xindex,2) .eq. 0) then
            gg1 = multipliers(Xindex)*(sin(omegas(Xindex)*xdash(1)))
        else
            gg1 = multipliers(Xindex)*(cos(omegas(Xindex)*xdash(1)))
        end if

        Yindex = sIndex(2,i)
        if (MOD(Yindex,2) .eq. 0) then
            gg2 = multipliers(Yindex)*(sin(omegas(Yindex)*xdash(2)))
        else
            gg2 = multipliers(Yindex)*(cos(omegas(Yindex)*xdash(2)))
        end if
    g(i) = sigma*gg1*gg2
    end do

    !! obtain PCE terms
    newY = 1.0d0
    i = term
    do j = 1,ndim
    index = mIndex(j,i)
        if (index .eq. 0) then
            yy = 1.0d0
        elseif (index .eq. 1) then
            yy = g(j)
        else
            Call findFactReal(index,factorial)
            yy = (g(j)**(index))/sqrt(factorial)
        end if
    newY = newY*yy
    end do
    y = meanc*newY


    !!---------------------------------------------------------------------------------------
    !! Trunctaed PCE of lognormal process
     CASE (3) !! L=3, O=3
    offsets(1) = sum(dbounds(1,:))/2.0d0
    offsets(2) = sum(dbounds(2,:))/2.0d0
    xdash = x-offsets

    g(1) = sigma*multipliers(1)*multipliers(1)*(cos(omegas(1)*xdash(1))*cos(omegas(1)*xdash(2)))
    g(2) = sigma*multipliers(1)*multipliers(2)*(cos(omegas(1)*xdash(1))*sin(omegas(2)*xdash(2)))
    g(3) = sigma*multipliers(2)*multipliers(1)*(sin(omegas(2)*xdash(1))*cos(omegas(1)*xdash(2)))

    SELECT CASE (term)
       CASE (1)                           !! 0th-Order (p=0)
          y = 1.0d0
       CASE (2)                           !! 1st-Order (p=1)
          y = g(1)
       CASE (3)
          y = g(2)
       CASE (4)
          y = g(3)
       CASE (5)                           !! 2nd-Order (p=2)
          y = sqrt(2.0d0)*g(1)**2/2.0d0
       CASE (6)
          y = g(1)*g(2)
       CASE (7)
          y = g(1)*g(3)
       CASE (8)
          y = sqrt(2.0d0)*g(2)**2/2.0d0
       CASE (9)
          y = g(2)*g(3)
       CASE (10)
          y = sqrt(2.0d0)*g(3)**2/2.0d0
       CASE (11)                          !! 3rd-Order (p=3)
          y = sqrt(6.0d0)*g(1)**3/6.0d0
       CASE (12)
          y = sqrt(2.0d0)*g(1)**2*g(2)/2.0d0
       CASE (13)
          y = sqrt(2.0d0)*g(1)**2*g(3)/2.0d0
       CASE (14)
          y = sqrt(2.0d0)*g(2)**2*g(1)/2.0d0
       CASE (15)
          y = g(1)*g(2)*g(3)
       CASE (16)
          y = sqrt(2.0d0)*g(3)**2*g(1)/2.0d0
       CASE (17)
          y = sqrt(6.0d0)*g(2)**3/6.0d0
       CASE (18)
          y = sqrt(2.0d0)*g(2)**2*g(3)/2.0d0
       CASE (19)
          y = sqrt(2.0d0)*g(3)**2*g(2)/2.0d0
       CASE (20)
          y = sqrt(6.0d0)*g(3)**3/6.0d0
    END SELECT
    y = meanc*y


    !!---------------------------------------------------------------------------------------
    !! Trunctaed PCE of lognormal process
    CASE (4) !! L=4, O=3
    offsets(1) = sum(dbounds(1,:))/2.0d0
    offsets(2) = sum(dbounds(2,:))/2.0d0
    xdash = x-offsets

    g(1) = sigma*multipliers(1)*multipliers(1)*(cos(omegas(1)*xdash(1))*cos(omegas(1)*xdash(2)))
    g(2) = sigma*multipliers(1)*multipliers(2)*(cos(omegas(1)*xdash(1))*sin(omegas(2)*xdash(2)))
    g(3) = sigma*multipliers(2)*multipliers(1)*(sin(omegas(2)*xdash(1))*cos(omegas(1)*xdash(2)))
    g(4) = sigma*multipliers(1)*multipliers(3)*(sin(omegas(1)*xdash(1))*sin(omegas(3)*xdash(2)))

    SELECT CASE (term)
       CASE (1)                           !! 0th-Order (p=0)
          y = 1.0d0
       CASE (2)                           !! 1st-Order (p=1)
          y = g(1)
       CASE (3)
          y = g(2)
       CASE (4)
          y = g(3)
       CASE (5)
          y = g(4)
       CASE (6)                           !! 2nd-Order (p=2)
          y = sqrt(2.0d0)*g(1)**2/2.0d0
       CASE (7)
          y = g(1)*g(2)   
       CASE (8)
          y = g(1)*g(3)
       CASE (9)
          y = g(1)*g(4)   
       CASE (10)   
          y = sqrt(2.0d0)*g(2)**2/2.0d0 
       CASE (11)
          y = g(2)*g(3)
       CASE (12) 
          y = g(2)*g(4)
       CASE (13)
          y = sqrt(2.0d0)*g(3)**2/2.0d0   
       CASE (14)
          y = g(3)*g(4)
       CASE (15)
          y = sqrt(2.0d0)*g(4)**2/2.0d0 
       CASE (16)                          !!3rd-Order (p=3)
          y = sqrt(6.0d0)*g(1)**3/6.0d0
       CASE (17)
          y = sqrt(2.0d0)*g(1)**2*g(2)/2.0d0
       CASE (18)
          y = sqrt(2.0d0)*g(1)**2*g(3)/2.0d0
       CASE (19)
          y = sqrt(2.0d0)*g(1)**2*g(4)/2.0d0
       CASE (20)
          y = sqrt(2.0d0)*g(2)**2*g(1)/2.0d0
       CASE (21)
          y = g(1)*g(2)*g(3)
       CASE (22)
          y = g(1)*g(2)*g(4)
       CASE (23)
          y = sqrt(2.0d0)*g(3)**2*g(1)/2.0d0
       CASE (24)
          y = g(1)*g(3)*g(4)
       CASE (25)
          y = sqrt(2.0d0)*g(1)*g(4)**2/2.0d0
       CASE (26)
          y = sqrt(6.0d0)*g(2)**3/6.0d0
       CASE (27)
          y = sqrt(2.0d0)*g(2)**2*g(3)/2.0d0
       CASE (28)
          y = sqrt(2.0d0)*g(2)**2*g(4)/2.0d0
       CASE (29)
          y = sqrt(2.0d0)*g(3)**2*g(2)/2.0d0
       CASE (30)
          y = g(2)*g(3)*g(4)
       CASE (31)
          y = sqrt(2.0d0)*g(4)**2*g(2)/2.0d0
       CASE (32)
          y = sqrt(6.0d0)*g(3)**3/6.0d0
       CASE (33)
          y = sqrt(2.0d0)*g(3)**2*g(4)/2.0d0
       CASE (34)
          y = sqrt(2.0d0)*g(4)**2*g(3)/2.0d0
       CASE (35)
          y = sqrt(6.0d0)*g(4)**3/6.0d0
    END SELECT
    y = meanc*y

    !!---------------------------------------------------------------------------------------
    !! Trunctaed PCE of lognormal process
    CASE (5) !! L=5, O=3
    offsets(1) = sum(dbounds(1,:))/2.0d0
    offsets(2) = sum(dbounds(2,:))/2.0d0
    xdash = x-offsets

    g(1) = sigma*multipliers(1)*multipliers(1)*(cos(omegas(1)*xdash(1))*cos(omegas(1)*xdash(2)))
    g(2) = sigma*multipliers(1)*multipliers(2)*(cos(omegas(1)*xdash(1))*sin(omegas(2)*xdash(2)))
    g(3) = sigma*multipliers(2)*multipliers(1)*(sin(omegas(2)*xdash(1))*cos(omegas(1)*xdash(2)))
    g(4) = sigma*multipliers(2)*multipliers(2)*(cos(omegas(2)*xdash(1))*cos(omegas(2)*xdash(2)))
    g(5) = sigma*multipliers(1)*multipliers(3)*(cos(omegas(1)*xdash(1))*cos(omegas(3)*xdash(2)))

    SELECT CASE (term)
       CASE (1)                           !! 0th-Order (p=0)
          y = 1.0d0
       CASE (2)                           !! 1st-Order (p=1)
          y = g(1)
       CASE (3)
          y = g(2)
       CASE (4)
          y = g(3)
       CASE (5)
          y = g(4)
       CASE (6)
          y = g(5)
       CASE (7)                           !! 2nd-Order (p=2)
          y = sqrt(2.0d0)*g(1)**2/2.0d0
       CASE (8)
          y = g(1)*g(2)
       CASE (9)
          y = g(1)*g(3)
       CASE (10)
          y = g(1)*g(4)
       CASE (11)
          y = g(1)*g(5)
       CASE (12)
          y = sqrt(2.0d0)*g(2)**2/2.0d0
       CASE (13)
          y = g(2)*g(3)
       CASE (14)
          y = g(2)*g(4)
       CASE (15)
          y = g(2)*g(5)
       CASE (16)
          y = sqrt(2.0d0)*g(3)**2/2.0d0
       CASE (17)
          y = g(3)*g(4)
       CASE (18)
          y = g(3)*g(5)
       CASE (19)
          y = sqrt(2.0d0)*g(4)**2/2.0d0
       CASE (20)
          y = g(4)*g(5)
       CASE (21)
          y = sqrt(2.0d0)*g(5)**2/2.0d0
       CASE (22)                          !!3rd-Order (p=3)
          y = sqrt(6.0d0)*g(1)**3/6.0d0
       CASE (23)
          y = sqrt(2.0d0)*g(1)**2*g(2)/2.0d0
       CASE (24)
          y = sqrt(2.0d0)*g(1)**2*g(3)/2.0d0
       CASE (25)
          y = sqrt(2.0d0)*g(1)**2*g(4)/2.0d0
       CASE (26)
          y = sqrt(2.0d0)*g(1)**2*g(5)/2.0d0
       CASE (27)
          y = sqrt(2.0d0)*g(1)*g(2)**2/2.0d0
       CASE (28)
          y = g(1)*g(2)*g(3)
       CASE (29)
          y = g(1)*g(2)*g(4)
       CASE (30)
          y = g(1)*g(2)*g(5)
       CASE (31)
          y = sqrt(2.0d0)*g(1)*g(3)**2/2.0d0
       CASE (32)
          y = g(1)*g(3)*g(4)
       CASE (33)
          y = g(1)*g(3)*g(5)
       CASE (34)
          y = sqrt(2.0d0)*g(1)*g(4)**2/2.0d0
       CASE (35)
          y = g(1)*g(4)*g(5)
       CASE (36)
          y = sqrt(2.0d0)*g(1)*g(5)**2/2.0d0
       CASE (37)
          y = sqrt(6.0d0)*g(2)**3/6.0d0
       CASE (38)
          y = sqrt(2.0d0)*g(2)**2*g(3)/2.0d0
       CASE (39)
          y = sqrt(2.0d0)*g(2)**2*g(4)/2.0d0
       CASE (40)
          y = sqrt(2.0d0)*g(2)**2*g(5)/2.0d0
       CASE (41)
          y = sqrt(2.0d0)*g(2)*g(3)**2/2.0d0
       CASE (42)
          y = g(2)*g(3)*g(4)
       CASE (43)
          y = g(2)*g(3)*g(5)
       CASE (44)
          y = sqrt(2.0d0)*g(2)*g(4)**2/2.0d0
       CASE (45)
          y = g(2)*g(4)*g(5)
       CASE (46)
          y = sqrt(2.0d0)*g(2)*g(5)**2/2.0d0
       CASE (47)
          y = sqrt(6.0d0)*g(3)**3/6.0d0
       CASE (48)
          y = sqrt(2.0d0)*g(3)**2*g(4)/2.0d0
       CASE (49)
          y = sqrt(2.0d0)*g(3)**2*g(5)/2.0d0
       CASE (50)
          y = sqrt(2.0d0)*g(3)*g(4)**2/2.0d0
       CASE (51)
          y = g(3)*g(4)*g(5)
       CASE (52)
          y = sqrt(2.0d0)*g(3)*g(5)**2/2.0d0
       CASE (53)
          y = sqrt(6.0d0)*g(4)**3/6.0d0
       CASE (54)
          y = sqrt(2.0d0)*g(4)**2*g(5)/2.0d0
       CASE (55)
          y = sqrt(2.0d0)*g(4)*g(5)**2/2.0d0
       CASE (56)
          y = sqrt(6.0d0)*g(5)**3/6.0d0
    END SELECT
    y = meanc*y


    END SELECT

    end subroutine c


end module variationalform
