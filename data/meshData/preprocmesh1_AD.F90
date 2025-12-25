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
PROGRAM preprocmesh1_AD

implicit none

integer, parameter :: n = 255
character(len=n) :: filename = 'foo.msh'
character(len=n) :: str, str_o
character(len=4) :: str1, str2
integer :: len_str, k, np, ne, nt, nb, i, j, net, nedgeinfo
integer :: temp4, l, cur_node, cur_dom, found, ndom, nb_sus
double precision, dimension(4) :: temp1
double precision, dimension(2) :: xlims,ylims
integer, dimension(12) :: temp2
integer, dimension(3) :: temp3
integer, dimension(:,:), allocatable :: nodes

print*, '-------------------------------'
print*, 'Mesh Decomposition Initiated...'
print*, '-------------------------------'
str1 = '$Nod'
str2 = '$End'

open(unit=1,file=filename,status='old')
open(unit=2,file='points.dat',status='replace')
open(unit=3,file='edges.dat',status='replace')
open(unit=4,file='triangles.dat',status='replace')
open(unit=5,file='meshdim.dat',status='replace')
open(unit=6,file='triangles_on_boundary.dat',status='replace')
open(unit=7,file='boundary_nodes.dat',status='replace')
open(unit=8,file='dbounds.dat',status='replace')

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
write(unit=2,fmt='(2ES17.8E3)') temp1(2:3)
xlims(:) = temp1(2)
ylims(:) = temp1(3)
do i = 2,np
  read (1,*) temp1
  write(unit=2,fmt='(2ES17.8E3)') temp1(2:3)
  xlims(1) = min(xlims(1),temp1(2))
  xlims(2) = max(xlims(2),temp1(2))
  ylims(1) = min(ylims(1),temp1(3))
  ylims(2) = max(ylims(2),temp1(3))
end do
write(unit=8,fmt='(2ES17.8E3)') xlims
write(unit=8,fmt='(2ES17.8E3)') ylims

read (1,'(a)') str
read (1,'(a)') str
read (1,*) net

read (1,'(a)') str
read (str,*) temp3
nedgeinfo = 5+temp3(3)
read (str,*) temp2(1:nedgeinfo)
write(unit=3,fmt='(4I8)') temp2(nedgeinfo-1), temp2(nedgeinfo), temp2(5), temp2(nedgeinfo-2)
ne = 1
do while(.true.)
    read (1,'(a)') str
    read (str,*) temp2(1:nedgeinfo)
    if (temp2(2) .eq. 1) then
        write(unit=3,fmt='(4I8)') temp2(nedgeinfo-1), temp2(nedgeinfo), temp2(5), temp2(nedgeinfo-2)
    else
        exit
    end if
    ne = ne+1
end do

nt = net - ne
read (str,*) temp3
nedgeinfo = 6+temp3(3)
read (str,*) temp2(1:nedgeinfo)
write(unit=4,fmt='(4I8)') temp2(nedgeinfo-2), temp2(nedgeinfo-1), temp2(nedgeinfo), temp2(7)
nb_sus = 0
if (temp3(3) .gt. 4) then
    write(unit=6,fmt='(2I8)') temp2(nedgeinfo-2), temp2(7)
    write(unit=6,fmt='(2I8)') temp2(nedgeinfo-1), temp2(7)
    write(unit=6,fmt='(2I8)') temp2(nedgeinfo), temp2(7)
    nb_sus = nb_sus+3
end if
do i = 2,nt
    read (1,'(a)') str
    read (str,*) temp3
    nedgeinfo = 6+temp3(3)
    read (str,*) temp2(1:nedgeinfo)
    write(unit=4,fmt='(4I8)') temp2(nedgeinfo-2), temp2(nedgeinfo-1), temp2(nedgeinfo), temp2(7)
    if (temp3(3) .gt. 4) then
        write(unit=6,fmt='(2I8)') temp2(nedgeinfo-2), temp2(7)
        write(unit=6,fmt='(2I8)') temp2(nedgeinfo-1), temp2(7)
        write(unit=6,fmt='(2I8)') temp2(nedgeinfo), temp2(7)
        nb_sus = nb_sus+3
    end if
end do
close(6)

open(unit=6,file='triangles_on_boundary.dat',status='old')
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

nb = 0
do i = 1,j
    if (nodes(i,3) .eq. 1) then
        nb = nb+1
        write(unit=7,fmt='(I8)') nodes(i,1)
    end if
end do

write(unit=5,fmt='(6I8)') np, ne, nt, nb, ndom, nb_sus

close(1)
close(2)
close(3)
close(4)
close(5)
close(6)
close(7)
close(8)

END PROGRAM preprocmesh1_AD
