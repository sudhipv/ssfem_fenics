!! Require for PCKF 
program preprocmesh_meas

implicit none

integer, parameter :: nmeasx = 8,nmeasy = 6   !!nmeasx = 8,nmeasy = 6 for 48 measurements
integer, parameter :: nmeas = nmeasx*nmeasy

integer, parameter :: n = 255
character(len=n) :: filename
character(len=n) :: str, str_o
character(len=4) :: str1, str2
integer :: np, ne, nt, nb, i, j, net, nedgeinfo, temp4, l, cur_node, cur_dom, found, ndom, nb_sus
double precision, dimension(4) :: temp1
double precision :: curd
double precision, dimension(nmeasx) :: locsx
double precision, dimension(nmeasy) :: locsy
double precision, dimension(nmeas,2) :: locs, locs_act
double precision, dimension(nmeas) :: mind
integer, dimension(nmeas) :: measnodes
integer, dimension(12) :: temp2
double precision, dimension(2) :: temp3
integer, dimension(:,:), allocatable :: nodes
double precision, dimension(2,2) :: dbounds


filename = 'meshdim.dat'
open(unit=1,file=filename,status='old')
read(1,*) np
close(1)

open(unit=1,file='dbounds.dat',status='old')
read(1,*) dbounds(1,:)
read(1,*) dbounds(2,:)
close(1)

mind(:) = (dbounds(1,2)-dbounds(1,1))**2 + (dbounds(2,2)-dbounds(2,1))**2

do i = 1,nmeasx
   locsx(i) = dbounds(1,1) + dble(i)*(dbounds(1,2)-dbounds(1,1))/dble(nmeasx+1)
end do
! locsx(:) = [0.5d0]
do i = 1,nmeasy
   locsy(i) = dbounds(2,1) + dble(i)*(dbounds(2,2)-dbounds(2,1))/dble(nmeasy+1)
end do

!print*, locsx
!print*, locsy
! locsy(:) = [0.5d0]

! locsx(1) = 5.0
! locsy(:) = (dbounds(2,2)+dbounds(2,1))/2.0

temp4 = 1
do i = 1,nmeasx
   do j = 1,nmeasy
      locs(temp4,1) = locsx(i)
      locs(temp4,2) = locsy(j)
      temp4 = temp4+1
   end do
end do
! locs(1,1) = 1.0
! locs(2,1) = 2.0
! locs(3,1) = 3.0
! locs(4,1) = 4.0
! locs(:,2) = (dbounds(2,2)+dbounds(2,1))/2.0

filename = 'points.dat'
open(unit=1,file=filename,status='old')
do i = 1,np
   read(unit=1,fmt=*) temp3
   do j = 1,nmeas
      curd = (locs(j,1)-temp3(1))**2 + (locs(j,2)-temp3(2))**2
      if (curd .lt. mind(j)) then
         mind(j) = curd
         measnodes(j) = i
         locs_act(j,:) = temp3
      end if
   end do
end do
close(1)

open(unit=1,file='measnodes.dat',status='replace')
open(unit=2,file='measlocs.dat',status='replace')
do i = 1,nmeas
   write(unit=1,fmt='(I8)') measnodes(i)
   write(unit=2,fmt='(2ES17.8E3)') locs_act(i,:)
end do
close(1)
close(2)

open(unit=1,file='nmeas.dat',status='replace')
write(unit=1,fmt='(I8)') nmeas
close(1)

end program preprocmesh_meas
