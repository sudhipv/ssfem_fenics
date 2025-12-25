!! Pre-processing code-1: By Ajit Desai, 8/2015, Mohammad Khalil, 1/2012
!! purpose: is to extract "GLOBAL" mesh data from "gmsh.msh" file 
!
!input:
!    gmsh.msh  : mesh file generated usig GMSH-ver2.8.5 or later
!
!output:
!    points.dat    : global nodal co-ordinates
!    edges.dat     : global 1D-connectivity 
!    triangles.dat : global 2D-connectivity 
!    dbounds.dat   : physical domain limits 
!    meshdim.at    : global mesh dimension (Np,Ne,Nt,Nb)
!    boundary_nodes.dat: global boundary nodes
!!
!! TEMP : For high nDim there is issue with "npceout" : Need to fix
!!---------------------------------------------------------------------------------------
PROGRAM KLE_PCE_Data

implicit none

    character(len=255) :: filename, str1, str2, str3
    integer, parameter :: nOrd = 2
    integer, parameter :: nDim = 50

    !!Temp: Assign PCE value = number of rows in multiIndex
    integer, parameter :: nPCE = 1326    !!Temp


    !!integer            :: nPCE         !!Temp
    integer            :: i

    integer, dimension(:,:), allocatable :: mIndex, SortIndex

    print*, '-----------------------'
    print*, 'Order     =', nOrd
    print*, '-----------------------'
    print*, 'Dimen     =', nDim
    print*, '-----------------------'

    ! CALL npceout(nOrd, nDim, npce) !!Order of PCE expansion      !!Temp
    print*, 'nPCE      =', nPCE
    print*, '-----------------------'

    !allocate(mIndex(nDim,nPCE))
    call int2str(str1,nord,4)
    call int2str(str2,ndim,4)

    filename = 'cp multiIndex' // trim(str1) // trim(str2) // '.dat multiIndex.dat'
    CALL system(filename)

    filename = 'cp cijk' // trim(str1) // trim(str2) // ' cijk'
    CALL system(filename)


    open(unit=1,file='pcedata.dat',status='replace')
    write(unit=1,fmt='(3I8)') nOrd, nDim, nPCE
close(1)

    allocate(SortIndex(2,nDim))
    open(unit=2,file='SortIndices.dat',status='old')
    read(2,*) SortIndex
close(2)

    open(unit=3,file='sortIndex.dat',status='replace')
    do i = 1,nDim
        write(unit=3,fmt='(2I8)') SortIndex(1,i), SortIndex(2,i)
    end do
close(3)


    CALL system("mv multiIndex.dat ../")
    CALL system("mv pcedata.dat ../")
    CALL system("mv sortIndex.dat ../")
    CALL system("mv cijk ../")


END PROGRAM KLE_PCE_Data

!!Temp
!!!**************************************************************************
!SUBROUTINE findFact(number,factorial)
!
!    INTEGER :: number, n, factorial
!
!    factorial = 1
!    n = number
!    DO WHILE (n/=0)
!        factorial = n*factorial
!        n = n-1
!    END DO
!
!END SUBROUTINE findFact
!
!
!SUBROUTINE npceout(nord, ndim, npce)
!
!    INTEGER :: nord, ndim, npce
!    INTEGER :: n1, d1, d2
!
!    CALL findFact((nord+ndim),n1)
!    CALL findFact((nord),d1)
!    CALL findFact((ndim),d2)
!
!    npce = n1/(d1*d2)
!
!END SUBROUTINE npceout

!!*********************************************************
!SUBROUTINE findFactReal(number,factorial)
!
!    INTEGER :: number, n
!    DOUBLE PRECISION :: factorial
!
!    factorial = 1
!    n = number
!    DO WHILE (n/=0)
!        factorial = n*factorial
!        n = n-1
!    END DO
!
!END SUBROUTINE findFactReal


!SUBROUTINE multiIndex(ndim,npce,mIndex)
!
!    integer :: ndim, npce
!    integer, dimension(3,20) :: mIndex
!
!    character(len=255) :: filename,str1
!    call int2str(str1,ndim,4)
!    filename = 'multiIndex' // trim(str1) // '.dat'
!
!    !open(unit=1,file='multiIndex.dat',status='old')
!    open(unit=1,file=filename,status='old')
!    read(1,*) mIndex
!    close(1)
!
!END SUBROUTINE multiIndex


!!**************************************************************************
!! Subroutine to convert integer to string to read "*00*.dat" files
SUBROUTINE int2str(str1,int1,n0pads)

    character(len=255) :: str1, str2
    integer :: int1, intdigs, i, curdig, curint
    integer :: n0pads

    intdigs = floor(log10(dble(int1))) + 1
    curint = int1
    str1 = ''
    do i = 1,intdigs
    curdig = mod(curint,10)
    curint = floor(dble(curint)/10.0d0)
    str1 = achar(48 + curdig) // str1
    enddo
    do i = intdigs+1,n0pads
    str1 = '0' // str1
    enddo

END SUBROUTINE int2str
