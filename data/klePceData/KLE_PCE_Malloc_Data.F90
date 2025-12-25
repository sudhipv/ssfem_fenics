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
!!---------------------------------------------------------------------------------------
PROGRAM KLE_PCE_Data

implicit none

    character(len=255) :: filename, str1, str2, str3
    integer, parameter :: nOrd = 3
    integer, parameter :: nDim = 3
    integer            :: nPCE

    CALL system("cp ~/currentRW/DDMpackage/SparseDDM/Stochastic/OneLevPCGM/VerMallocS/gmsh.msh ../")
    CALL system("cp ~/currentRW/DDMpackage/SparseDDM/Stochastic/OneLevPCGM/VerMallocS/multiIndex.dat ../")
    CALL system("cp ~/currentRW/DDMpackage/SparseDDM/Stochastic/OneLevPCGM/VerMallocS/sortIndex.dat ../")
    CALL system("cp ~/currentRW/DDMpackage/SparseDDM/Stochastic/OneLevPCGM/VerMallocS/pcedata.dat ../")
    CALL system("cp ~/currentRW/DDMpackage/SparseDDM/Stochastic/OneLevPCGM/VerMallocS/cijk ../")
    CALL system("cp ~/currentRW/DDMpackage/SparseDDM/Stochastic/OneLevPCGM/VerMallocS/nnzi0* ../")
    CALL system("cp ~/currentRW/DDMpackage/SparseDDM/Stochastic/OneLevPCGM/VerMallocS/nnzb0* ../")
    CALL system("cp ~/currentRW/DDMpackage/SparseDDM/Stochastic/OneLevPCGM/VerMallocS/nnzbi0* ../")

END PROGRAM KLE_PCE_Data

!!**************************************************************************
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
!SUBROUTINE npceout2(nord, ndim, npce)
!
!INTEGER :: nord, ndim, npce
!INTEGER :: aa, bb, cc
!
!    aa = nord+ndim
!    bb = aa
!    IF (nord/=1) THEN
!        DO i=1,(nord-1)
!           bb = bb*(aa-i)
!        END DO
!    END IF
!    CALL findFact((nord),cc)
!    npce = bb/cc
!
!END SUBROUTINE npceout2
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
