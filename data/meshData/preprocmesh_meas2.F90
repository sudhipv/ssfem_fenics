!! Require for PCKF
program preprocmesh_meas2

implicit none

integer, parameter :: n = 255

integer :: nmeas,i,j,k,nmeasi,curnode,np,ne,nt,nb,nd
integer, allocatable, dimension(:) :: measnodes,measnodesi
character(len=n) :: filename
character(len=n) :: str1


filename = 'nmeas.dat'
open(unit=1,file=filename,status='old')
read(1,*) nmeas
close(1)

allocate(measnodes(nmeas),measnodesi(nmeas))
filename = 'measnodes.dat'
open(unit=1,file=filename,status='old')
do i = 1,nmeas
   read(1,*) measnodes(i)
end do
close(1)

open(unit=1,file='meshdim.dat',status='old')
read(unit=1,fmt='(6I8)') np, ne, nt, nb, nd
close(1)

! !!!!!!!!!!!! Measurements on Boundary !!!!!!!!!!!!!!!!!!!
! filename = 'boundary_nodes.dat'
! open(unit=1,file=filename,status='old')
! filename = 'measnodes_boundary.dat'
! open(unit=2,file=filename,status='replace')
! filename = 'measnodes_boundary_map.dat'
! open(unit=3,file=filename,status='replace')
! nmeasi = 0
! do i = 1,nb
!    read(1,*) curnode
! !    print*,curnode
!    do j = 1,nmeas
!       if (curnode .eq. measnodes(j)) then
!          nmeasi = nmeasi + 1
!          measnodes(j) = 0
!          write(unit=2,fmt='(I8)') i
!          write(unit=3,fmt='(I8)') j
!       end if
!    end do
! end do
! close(1)
! close(2)
! close(3)
! filename = 'nmeas_boundary.dat'
! open(unit=1,file=filename,status='replace')
! write(unit=1,fmt='(I8)') nmeasi
! close(1)

!!!!!!!!!!!! Interior measurements !!!!!!!!!!!!!!!!!!!
do k = 1,nd
   call int2str(str1,k,4)
   str1 = 'meshdim' // trim(str1) // '.dat'
   open(unit=1,file=str1,status='old')
   read(unit=1,fmt=*) np, ne, nt, nb
   close(1)

   call int2str(str1,k,4)
   str1 = 'nodes' // trim(str1) // '.dat'
   open(unit=1,file=str1,status='old')
   call int2str(str1,k,4)
   str1 = 'measnodes' // trim(str1) // '.dat'
   open(unit=2,file=str1,status='replace')
   call int2str(str1,k,4)
   str1 = 'measnodes' // trim(str1) // '_map.dat'
   open(unit=3,file=str1,status='replace')
   nmeasi = 0
   do i = 1,np
      read(1,*) curnode
      do j = 1,nmeas
         if (curnode .eq. measnodes(j)) then
            nmeasi = nmeasi + 1
            measnodes(j) = 0
            write(unit=2,fmt='(I8)') i
            write(unit=3,fmt='(I8)') j
         end if
      end do
   end do
   close(1)
   close(2)
   close(3)
   call int2str(str1,k,4)
   str1 = 'nmeas' // trim(str1) // '.dat'
   open(unit=1,file=str1,status='replace')
   write(unit=1,fmt='(I8)') nmeasi
   close(1)
end do

end program preprocmesh_meas2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine int2str(str1,int1,n0pads)

  character(len=255) :: str1, str2
  integer :: int1, intdigs, i, curdig, curint,n0pads

  intdigs = floor(log10(dble(int1))) + 1
  curint = int1
  str1 = ''
  do i = 1,intdigs
    curdig = mod(curint,10)
    curint = floor(dble(curint)/10.0)
    str1 = achar(48 + curdig) // str1
  enddo
  do i = intdigs+1,n0pads
    str1 = '0' // str1
  enddo

end subroutine int2str
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!