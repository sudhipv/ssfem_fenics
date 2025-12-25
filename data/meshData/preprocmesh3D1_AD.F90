!! 3D Pre-processing code-1: Global decomposer
! Copyright (C) 2017 Ajit Desai, PhD Candidate, ajit.ndesai@gmail.com
!! purpose: is to extract "GLOBAL" mesh data from "gmsh.msh" file 
!
!input:
!    3D-GMSH.msh     : mesh file generated usig GMSH-ver2.8.5 or later
!
!outputs:
!    points.dat      : global nodal co-ordinates
!    edges.dat       : global 1D-connectivity
!    triangles.dat   : global 2D-connectivity
!    tetrahedrons.dat: global 3D-connectivity
!    dbounds.dat     : physical domain limits
!    meshdim.at      : global mesh dimension (Np,Ne,Nt,Nh,Nb,Nd)
!    boundary_nodes.dat: global boundary nodes
!!
!!---------------------------------------------------------------------------------------
PROGRAM preprocmesh3D1_AD

implicit none

integer, parameter :: n = 255
character(len=n) :: filename = 'foo3D.msh'
character(len=n) :: str
character(len=4) :: str1, str2
integer :: len_str, k, np, ne, nt, nh, nb, i, j, net, nedgeinfo
integer :: temp4, l, cur_node, cur_dom, found, ndom, nb_sus
double precision, dimension(4) :: temp1
double precision, dimension(3) :: xlims,ylims,zlims
integer, dimension(12) :: temp2
integer, dimension(3) :: temp3
integer, dimension(:,:), allocatable :: nodes
double precision :: mesh3D

print*, '----------------------------------------'
print*, '3D Local Mesh Decomposition Initiated...'
print*, '----------------------------------------'
str1 = '$Nod'
str2 = '$End'

!! Opend input "GMSH-file"
open(unit=1,file=filename,status='old')

do while(.true.)
  read (1,'(a)') str
  len_str = len_trim(str)
  k = index(str(1:len_str),str1)
  if(k.ne.0) then
    read (1,*) np
    exit
  end if
end do

read (1,*) temp1
open(unit=2,file='points.dat',status='replace')
write(unit=2,fmt='(3ES17.8E3)') temp1(2:4)

xlims(:) = temp1(2)
ylims(:) = temp1(3)
zlims(:) = temp1(4)

do i = 2,np
    read (1,*) temp1
    write(unit=2,fmt='(3ES17.8E3)') temp1(2:4)
    xlims(1) = min(xlims(1),temp1(2))
    xlims(2) = max(xlims(2),temp1(2))
    ylims(1) = min(ylims(1),temp1(3))
    ylims(2) = max(ylims(2),temp1(3))
    zlims(1) = min(zlims(1),temp1(4))
    zlims(2) = max(zlims(2),temp1(4))
end do
close(2)

open(unit=8,file='dbounds.dat',status='replace')
write(unit=8,fmt='(3ES17.8E3)') xlims
write(unit=8,fmt='(3ES17.8E3)') ylims
write(unit=8,fmt='(3ES17.8E3)') zlims
close(8)

if (maxval(zlims) .eq. 0) then
    print*, "**** NOT A 3D MESH ****"
    print*, '----------------------------------'
    stop
end if


read (1,'(a)') str
read (1,'(a)') str
read (1,*) net
read (1,'(a)') str
read (str,*) temp3

nedgeinfo = 5+temp3(3)
read (str,*) temp2(1:nedgeinfo)
open(unit=3,file='edges.dat',status='replace')
write(unit=3,fmt='(4I8)') temp2(nedgeinfo-1), temp2(nedgeinfo), temp2(4), temp2(nedgeinfo-2)
ne = 1
do while(.true.)
    read (1,'(a)') str
    read (str,*) temp2(1:nedgeinfo)
    if (temp2(2) .eq. 1) then
        write(unit=3,fmt='(4I8)') temp2(nedgeinfo-1), temp2(nedgeinfo), temp2(4), temp2(nedgeinfo-2)
    else
        exit
    end if
    ne = ne+1
end do
close(3)

nedgeinfo = 10
read (str,*) temp2(1:nedgeinfo)
open(unit=4,file='triangles.dat',status='replace')
write(unit=4,fmt='(4I8)') temp2(nedgeinfo-2), temp2(nedgeinfo-1), temp2(nedgeinfo), temp2(7)
nt = 1
do while(.true.)
    read (1,'(a)') str
    read (str,*) temp2(1:nedgeinfo)
    if (temp2(2) .eq. 2) then
    write(unit=4,fmt='(4I8)') temp2(nedgeinfo-2), temp2(nedgeinfo-1), temp2(nedgeinfo), temp2(7)
    else
    exit
end if
nt = nt+1
end do
close(4)

nh = net - ne - nt
read (str,*) temp3
nedgeinfo = 7+temp3(3)
read (str,*) temp2(1:nedgeinfo)
open(unit=6,file='tetrahedrons_on_boundary.dat',status='replace')
open(unit=9,file='tetrahedrons.dat',status='replace')
write(unit=9,fmt='(5I8)') temp2(nedgeinfo-3), temp2(nedgeinfo-2), temp2(nedgeinfo-1), temp2(nedgeinfo), temp2(7)
nb_sus = 0
if (temp3(3) .gt. 4) then
    write(unit=6,fmt='(2I8)') temp2(nedgeinfo-3), temp2(7)
    write(unit=6,fmt='(2I8)') temp2(nedgeinfo-2), temp2(7)
    write(unit=6,fmt='(2I8)') temp2(nedgeinfo-1), temp2(7)
    write(unit=6,fmt='(2I8)') temp2(nedgeinfo), temp2(7)
    nb_sus = nb_sus+4
end if
do i = 2,nh
    read (1,'(a)') str
    read (str,*) temp3
    nedgeinfo = 7+temp3(3)
    read (str,*) temp2(1:nedgeinfo)
    write(unit=9,fmt='(5I8)') temp2(nedgeinfo-3), temp2(nedgeinfo-2), temp2(nedgeinfo-1), temp2(nedgeinfo), temp2(7)
    if (temp3(3) .gt. 4) then
        write(unit=6,fmt='(2I8)') temp2(nedgeinfo-3), temp2(7)
        write(unit=6,fmt='(2I8)') temp2(nedgeinfo-2), temp2(7)
        write(unit=6,fmt='(2I8)') temp2(nedgeinfo-1), temp2(7)
        write(unit=6,fmt='(2I8)') temp2(nedgeinfo), temp2(7)
        nb_sus = nb_sus+4
    end if
end do
close(6)
close(9)
close(1)


open(unit=6,file='tetrahedrons_on_boundary.dat',status='old')
allocate(nodes(nb_sus,3))
nodes(:,:) = 0
j = 0
ndom = 0
nb = 0
do i = 1,nb_sus
    read(6,*) cur_node, cur_dom
    ndom = max(ndom, cur_dom)
    found = 0
    do k = 1,j
        if (cur_node .eq. nodes(k,1)) then
            found = 1
            if (cur_dom .ne. nodes(k,2)) then
                nb = nb+1
                nodes(k,3) = 1
            end if
            exit
        end if
    end do
    if (found .eq. 0) then
        j = j+1
        nodes(j,1) = cur_node
        nodes(j,2) = cur_dom
    end if
end do

open(unit=7,file='boundary_nodes.dat',status='replace')
nb = 0
do i = 1,j
    if (nodes(i,3) .eq. 1) then
        nb = nb+1
        write(unit=7,fmt='(I8)') nodes(i,1)
    end if
end do

open(unit=5,file='meshdim.dat',status='replace')
write(unit=5,fmt='(7I8)') np, ne, nt, nh, nb, ndom, nb_sus
close(5)
close(7)

END PROGRAM preprocmesh3D1_AD
