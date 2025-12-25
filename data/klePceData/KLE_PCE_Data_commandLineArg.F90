!! KLE/PCE Pre-processing code: By Ajit Desai, 8/2015
!! purpose: is to extract "KLE/PCE" data for selected case
!
!input:
!    KLEord  : Order of input KLE: currently hardcoded to 2, can be changed here
!    nOrd    : Order of output-PCE: currently hardcoded to 3, can be changed here
!    nDim    : Dimension (number of RVs) for both input and output PCE
!            : Takes in command-line input after code compilation
!
!output:
!    pcedata.dat    : create pcedata file with "nOrd, nDim, nPCE, npcein"
!    multiIndex.dat : multiIndices for selected case
!    nZijk.dat      : nonZero block structure for selected case
!    cijk           : nonZero tripple products and the respective indices
!                   : for the selected case
!!
!!---------------------------------------------------------------------------------------
PROGRAM KLE_PCE_Data

implicit none

    character(len=255) :: filename, str1, str2, str3
    integer, parameter :: KLEord = 2
    integer, parameter :: nOrd = 3
    integer            :: nDim

    integer            :: nPCE, npcein
    integer            :: i

    integer, dimension(:,:), allocatable :: mIndex, SortIndex

    print*, '=============================='
    print*, 'Generate & Import KLE/PCE Data'

    print*, 'Enter PCE Dimension (RVs): '
    read(*,*), nDim


if (nDim == 1) THEN
    print*, "Special case: for nDim=1, nOrd=1"
    call system("cp pcedata00010001.dat pcedata.dat")
    call system("cp multiIndex00010001.dat multiIndex")
    call system("cp nZijk00010001.dat nZijk.dat")
    call system("cp cijk00010001 cijk")
else
    print*, '-----------------------'
    print*, 'Order     =', nOrd
    print*, '-----------------------'
    print*, 'Dimen     =', nDim
    print*, '-----------------------'

    !!CALL npceout(nOrd, nDim, npce) !! for small nOrd and small nDim
    CALL npceout2(nOrd, nDim, npce)  !! for small nOrd and large nDim
    CALL npceout2(KLEord, nDim, npcein)  !! for small nOrd and large nDim

    print*, 'nPCEin    =', npcein
    print*, '-----------------------'
    print*, 'nPCEout   =', nPCE
    print*, '-----------------------'

    !!allocate(mIndex(nDim,nPCE))
    call int2str(str1,nord,4)
    call int2str(str2,ndim,4)

    filename = 'cp multiIndex' // trim(str1) // trim(str2) // '.dat multiIndex.dat'
    CALL system(filename)

    filename = 'cp cijk' // trim(str1) // trim(str2) // ' cijk'
    CALL system(filename)

    filename = 'cp nZijk' // trim(str1) // trim(str2) // '.dat nZijk.dat'
    CALL system(filename)

    open(unit=1,file='pcedata.dat',status='replace')
    write(unit=1,fmt='(4I8)') nOrd, nDim, nPCE, npcein
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

    print*, '=============================='

end if

END PROGRAM KLE_PCE_Data

!CALL system("mv multiIndex.dat ../")
!CALL system("mv pcedata.dat ../")
!!**************************************************************************
SUBROUTINE findFact(number,factorial)

    INTEGER :: number, n, factorial

    factorial = 1
    n = number
    DO WHILE (n/=0)
        factorial = n*factorial
        n = n-1
    END DO

END SUBROUTINE findFact


SUBROUTINE npceout2(nord, ndim, npce)

INTEGER :: nord, ndim, npce
INTEGER :: aa, bb, cc

    aa = nord+ndim
    bb = aa
    IF (nord/=1) THEN
        DO i=1,(nord-1)
           bb = bb*(aa-i)
        END DO
    END IF
    CALL findFact((nord),cc)
    npce = bb/cc

END SUBROUTINE npceout2


SUBROUTINE npceout(nord, ndim, npce)

    INTEGER :: nord, ndim, npce
    INTEGER :: n1, d1, d2

    CALL findFact((nord+ndim),n1)
    CALL findFact((nord),d1)
    CALL findFact((ndim),d2)

    npce = n1/(d1*d2)

END SUBROUTINE npceout

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
