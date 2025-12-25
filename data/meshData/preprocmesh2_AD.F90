!! Pre-processing code-2: By Ajit Desai, 8/2015, Mohammad Khalil, 1/2012
!! purpose: is to extract "LOCAL" mesh data from "gmsh.msh" file 
!
!input:
!    gmsh.msh      : mesh file generated usig GMSH-ver2.8.5 or later
!    global-data   : global data generated uinsg preprocmesh1_AD.f90
!
!output:
!    nodes00*.dat  : local nodes  
!    bnodes00*.dat : local boundary nodes
!    cnodes00*.dat : local corner nodes
!    rnodes00*.dat : local remaining nodes
!    dimcrn00*.dat : local dimension for corner & remaining 
!    points00*.dat : local nodal co-ordinates
!    edges00*.dat  : local 1D-connectivity 
!    meshdim00*.dat: local mesh dimension (np,ne,nr,nb) 
!    meshdim.at    : global mesh dimension (Np,Ne,Nt,Nb)
!    bnodes00*.dat : new-rearranged-renumbered boundary nodes
!    boundary_nodes.dat: global boundary nodes
!    boundary_nodes.dat: global corner nodes
!    triangles00*.dat  : local 2D-connectivity 
!!
!!---------------------------------------------------------------------------------------

PROGRAM preprocmesh2_AD

! use common
implicit none

!!First-Level-Data
integer, parameter :: n = 255
character(len=n) :: filename = 'gmsh.msh'
character(len=n) :: str, str_o
character(len=n) :: str1, str2, extension
integer(kind=8) :: len_str, k, np, ne, nt, nb
integer(kind=8) :: i, j, net, nedgeinfo, temp4, l
integer(kind=8) :: cur_node, cur_dom, found, ndom
integer(kind=8) :: nbsus, nbicur, npicur, Ri, Rj
integer(kind=8), dimension(3) :: temp3, trianglei
integer(kind=8), dimension(4) :: triangle, edge
integer(kind=8), dimension(12) :: temp2
double precision, dimension(4) :: temp1
integer, dimension(:), allocatable :: nti, npi, nei, bnodes
integer, dimension(:), allocatable :: nbi, nodes, node_renum
double precision, dimension(:,:), allocatable :: points

!! Second-Level-Data :Coarse-Mesh
integer(kind=8) :: cn, nc1, nci, nri, nc2, nbgc, nbgr
integer(kind=8), dimension(:), allocatable :: edges2,countit,corner,cnodes,rnodes
integer(kind=8), dimension(:), allocatable :: cornerEdg2, cornerTemp, remaining
integer(kind=8), dimension(:,:), allocatable :: nodes2

!!---------------------------------------------------------------------------------------
open(unit=1,file='meshdim.dat',status='old')
read(unit=1,fmt='(6I8)') np, ne, nt, nb, ndom, nbsus
close(1)

open(unit=2,file='num_partition.dat',status='replace')
write(unit=2,fmt='(I8)') ndom
close(2)

allocate(npi(ndom),nti(ndom),nei(ndom),nbi(ndom),bnodes(nb),points(2,np),node_renum(np), &
& nodes(floor(2.0d0*dble(np)/dble(ndom) + 10.0d0)))
npi(:) = 0
nti(:) = 0
nei(:) = 0
nbi(:) = 0

!!!!!!!!!!!!!  BOUNDARY  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
open(unit=1,file='boundary_nodes.dat',status='old')
do i = 1,nb
    read(1,*) bnodes(i)
end do
close(1)

open(unit=1,file='triangles_on_boundary.dat',status='old')
do i=1,ndom
    call int2str(str1,i,4)
    str1 = 'bnodestemp' // trim(str1) // '.dat'
    open(unit=10+i,file=str1,status='replace')
end do
do i = 1,nbsus
    read(1,*) cur_node, cur_dom
    do j = 1,nb
        if (cur_node .eq. bnodes(j)) then
            write(unit=10+cur_dom,fmt='(I8)') cur_node
            nbi(cur_dom) = nbi(cur_dom)+1
            exit
        end if
    end do
end do
! close(unit=1,status='delete')
close(1)

do i=1,ndom
    close(unit=10+i)
end do

do j=1,ndom
    call int2str(str1,j,4)
    str1 = 'bnodestemp' // trim(str1) // '.dat'
    open(unit=1,file=str1,status='old')

    nodes(:) = 0
    nbicur = 0
    do i = 1,nbi(j)
        read(1,*) cur_node
        found = 0
        do k = 1,nbicur
            if (cur_node .eq. nodes(k)) then
                found = 1
                exit
            end if
        end do
        if (found .eq. 0) then
            nbicur=nbicur+1
            if (nbicur .gt. size(nodes)) then
                print*, 'not enough allocation for nodes in preprocmesh2.f'
                call exit(0)
            end if
            nodes(nbicur) = cur_node
        end if
    end do
    nbi(j) = nbicur
    close(unit=1,status='delete')
    
    call int2str(str1,j,4)
    str1 = 'bnodes' // trim(str1) // '.dat'
    open(unit=1,file=str1,status='replace')
    do i = 1,nbi(j)
        write(unit=1,fmt='(I8)') nodes(i)
    end do
    close(1)
end do


!!!!!!!!!!!!!  EDGES  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
open(unit=1,file='edges.dat',status='old')
do i=1,ndom
    call int2str(str1,i,4)
    str1 = 'edges' // trim(str1) // '.dat'
    open(unit=10+i,file=str1,status='replace')
end do

do i=1,ne
    read(unit=1,fmt=*) edge
    nei(edge(4)) = nei(edge(4))+1
    write(unit=10+edge(4),fmt='(3I8)') edge(1:3)
end do

close(1)
do i=1,ndom
    close(unit=10+i)
end do


!!!!!!!!!!!!!  TRIANGLES  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
open(unit=1,file='triangles.dat',status='old')
do i=1,ndom
    call int2str(str1,i,4)
    str1 = 'triangles' // trim(str1) // '.dat'
    open(unit=10+i,file=str1,status='replace')
end do

do i=1,nt
    read(unit=1,fmt='(4I8)') triangle
    nti(triangle(4)) = nti(triangle(4))+1
    write(unit=10+triangle(4),fmt='(3I8)') triangle(1:3)
end do

close(1)
do i=1,ndom
    close(unit=10+i)
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! MODIFIED FOR CORNER & REMAINING NODES BY AJIT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PRINT*, 'Finding Corner & Remaining Nodes...'
!! The array to assign count for each node (the particular node shared between how many sub-domains).
 allocate(countit(nb))
 countit = 0.0d0
do i=1,ndom
  
    allocate (nodes2(nbi(i),2))
    nodes2(:,:) = 0.0d0
    call int2str(str1,i,4)
    str1 = 'bnodes' // trim(str1) // '.dat'
    open(unit=2,file=str1,status='old')
    do j = 1,nbi(i)
	read(unit=2,fmt='(2I8)') nodes2(j,1)
    end do
    close(2)
    
    do j = 1,nb      !Global boundary_nodes
      do k = 1,nbi(i)   !local boundary_nodes
	  if (bnodes(j) .eq. nodes2(k,1)) then
	  countit(j) = countit(j) + 1          
	  end if
      end do
    end do   
    deallocate(nodes2)
end do
open(unit=2,file='countit.dat',status='replace')
write(unit=2,fmt='(I8)') countit
close(2)

!!----------------------------------------------------------------------------------------
!!  To Extract the corner nodes: Nodes exist in more than two sub-domain
open(unit=3,file='corner_nodesWithoutEdges.dat',status='replace')
nc1 = 0.0d0
do i=1,nb   
  if (countit(i) .gt. 2.0) then
    cn = bnodes(i)
    ! WRITE(*,*) cn
    nc1 = nc1 + 1
    write(unit=3,fmt='(I8)') cn
  end if
end do
close(3)

allocate(corner(nc1))
open(unit=2,file='corner_nodesWithoutEdges.dat',status='old')
do i = 1,nc1
    read(unit=2,fmt='(I8)') corner(i)
end do
close(2)


!!----------------------------------------------------------------------------------------
!! To Separate the corner-nodes(+plus the nodes at the ends of interface edges), 
!! the remaining-nodes & their dimensions for each sub-domain.

allocate(edges2(ne))
open(unit=2,file='edges.dat',status='old')
do i = 1,ne
    read(unit=2,fmt='(4I8)') edges2(i)
end do
close(2)


!! To Find The Global Interface Corner nodes 'corner_nodes.dat'
!! Added on Feb 13,2014: Need to check throughtly for large number of partions
allocate(cornerEdg2(nb))
cornerEdg2(:) = 0.0d0
DO i = 1,nb
   Ri = bnodes(i)
   DO j = 1,ne
      Rj = edges2(j)
      IF (Ri == Rj) THEN
         cornerEdg2(i) = 1
      END IF
   END DO
END DO

open(unit=7,file='corner_nodes.dat',status='replace')
open(unit=8,file='nbgc.dat',status='replace')

nc2 = 0
do i = 1,nb
   if (cornerEdg2(i) .eq. 1) then
   nc2 = nc2+1
   write(unit=7,fmt='(I8)') bnodes(i)
   end if
end do
close(7)

open(unit=7,file='corner_nodes.dat', position="append", status='old')
write(unit=7,fmt='(I8)') corner
nbgc = nc2 + nc1
write(unit=8,fmt='(I8)') nbgc
close(7)
close(8)

deallocate(cornerEdg2)
!!!!!!!!!!!!!!!!!!!
! To find the Global Remaining-Interface nodes 'remaining_nodes.dat'
! Added on March 26, 2014:
allocate(cornerTemp(nbgc))
open(unit=2,file='corner_nodes.dat',status='old')
do i = 1,nbgc
    read(unit=2,fmt='(I8)') cornerTemp(i)
end do
close(2)
allocate(cornerEdg2(nb))
cornerEdg2(:) = 0.0d0
DO i = 1,nb
   Ri = bnodes(i)
   DO j = 1,nbgc
      Rj = cornerTemp(j)
      IF (Ri == Rj) THEN
         cornerEdg2(i) = 1
      END IF
   END DO
END DO

open(unit=7,file='remaining_nodes.dat',status='replace')
open(unit=8,file='nbgr.dat',status='replace')

nbgr = 0
do i = 1,nb
   if (cornerEdg2(i) .eq. 0) then
   nbgr = nbgr+1
   write(unit=7,fmt='(I8)') bnodes(i)
   end if
end do
write(unit=8,fmt='(I8)') nbgr
close(7)
close(8)
deallocate(cornerEdg2,cornerTemp)
!!!!!!!!!!!!!!!!!!!

do i=1,ndom

    allocate (nodes2(nbi(i),2))
    nodes2(:,:) = 0
    call int2str(str1,i,4)
    str1 = 'bnodes' // trim(str1) // '.dat'
    open(unit=2,file=str1,status='old')
    do j = 1,nbi(i)
	read(2,*) nodes2(j,1)
    end do
    close(2)

    do j = 1,nbi(i)
      do k = 1,ne
	if (nodes2(j,1) .eq. edges2(k)) then
	nodes2(j,2) = 1
	end if
      end do
      nri = nri + 1
    end do

    do j = 1,nbi(i)
      do k = 1,nc1
	if (nodes2(j,1) .eq. corner(k)) then
	nodes2(j,2) = 1
	end if
      end do
    end do
  
    call int2str(str1,i,4)
    str1 = 'rnodes' // trim(str1) // '.dat' 
    open(unit=2,file=str1,status='replace')

    call int2str(str1,i,4)
    str1 = 'cnodes' // trim(str1) // '.dat' 
    open(unit=3,file=str1,status='replace')
    
    call int2str(str1,i,4)
    str1 = 'dimcrn' // trim(str1) // '.dat' 
    open(unit=4,file=str1,status='replace')
	
    nci = 0.0d0
    nri = 0.0d0
	
    do j = 1,nbi(i)
      if (nodes2(j,2) .ne. 1) then
      write(unit=2,fmt='(I8)') nodes2(j,1)
      nri = nri + 1
      else 
      write(unit=3,fmt='(I8)') nodes2(j,1)
      nci = nci + 1
      end if
    end do   
    close(2) 
    close(3)

    write(unit=4,fmt='(2I8)') nci, nri
    close(4)
    deallocate(nodes2)
end do

print*, 'Re-arranging Initiated...'
do i=1,ndom
    call int2str(str1,i,4)
    str1 = 'nbnodes' // trim(str1) // '.dat' 
    open(unit=4,file=str1,status='replace')

    call int2str(str1,i,4)
    str1 = 'dimcrn' // trim(str1) // '.dat' 
    open(unit=1,file=str1,status='old')
    read(1,fmt='(2I8)') nci, nri
    close(1)
    
    allocate(rnodes(nri),cnodes(nci))
    
    call int2str(str1,i,4)
    str1 = 'rnodes' // trim(str1) // '.dat' 
    open(unit=2,file=str1,status='old')
    
    do j = 1,nri
       read(2,*) rnodes(j)
       write(unit=4,fmt='(I8)') rnodes(j)
    end do
    close(2)
    close(4)
   
    call int2str(str1,i,4)
    str1 = 'cnodes' // trim(str1) // '.dat' 
    open(unit=3,file=str1,status='old')
    
    call int2str(str1,i,4)
    str1 = 'nbnodes' // trim(str1) // '.dat' 
    open(unit=4,file=str1 ,position="append", status='old')

    do j = 1,nci
       read(3,*) cnodes(j)
       write(unit=4,fmt='(I8)') cnodes(j)
    end do
    close(3)
    close(4)
    deallocate(rnodes, cnodes)   
end do  

deallocate(edges2, countit, corner)
!!!!!!!!!!! UNTIL HERE MODIFIED BY AJIT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!  POINTS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=1,ndom
    call int2str(str1,i,4)
    str1 = 'nbnodes' // trim(str1) // '.dat'
    open(unit=1,file=str1,status='old')
    do j = 1,nbi(i)
        read(1,*) nodes(j)
    end do
    close(1)

    call int2str(str1,i,4)
    str1 = 'triangles' // trim(str1) // '.dat'
    open(unit=1,file=str1,status='old')
    call int2str(str1,i,4)
    str1 = 'nodes' // trim(str1) // '.dat'
    open(unit=2,file=str1,status='replace')

    npicur = nbi(i)
    do j = 1,nti(i)
        read(unit=1,fmt='(3I8)') trianglei
        do l = 1,3
        found = 0
            do k = 1,npicur
                if (nodes(k) .eq. trianglei(l)) then
                    found = 1
                    exit
                end if
            end do
            if (found .eq. 0) then
                npicur = npicur+1
                if (npicur .gt. size(nodes)) then
                    print*, 'not enough allocation for nodes in preprocmesh2.f'
                    call exit(0)
                end if
                nodes(npicur) = trianglei(l)
                write(unit=2,fmt='(I8)') trianglei(l)
            end if
        end do
    end do
    do j = 1,nbi(i)
        write(unit=2,fmt='(I8)') nodes(j)
    end do

    close(1)
    close(2)
    npi(i) = npicur
end do

open(unit=1,file='points.dat',status='old')
do i=1,np
    read(1,*) points(:,i)
end do
close(1)

do i=1,ndom
    call int2str(str1,i,4)
    str1 = 'nodes' // trim(str1) // '.dat'
    open(unit=1,file=str1,status='old')

    call int2str(str1,i,4)
    str1 = 'points' // trim(str1) // '.dat'
    open(unit=2,file=str1,status='replace')

    do j = 1,npi(i)
        read(unit=1,fmt=*) temp4
        write(unit=2,fmt='(2ES17.8E3)') points(:,temp4)
    end do
    
    close(1)
    close(2)
end do

!!!!!!!!!!!!!  MESH DIMENSIONS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=1,ndom
    call int2str(str1,i,4)
    str1 = 'meshdim' // trim(str1) // '.dat'
    open(unit=1,file=str1,status='replace')
    write(unit=1,fmt='(4I8)') npi(i), nei(i), nti(i), nbi(i)
    close(1)
end do

!!!!!!!!!!!!!  RENUMBER EDGES AND TRIANGLES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
print*, 'Re-numbering Initiated...'
do i=1,ndom
    call int2str(str1,i,4)
    str1 = 'nodes' // trim(str1) // '.dat'
    open(unit=1,file=str1,status='old')
    do j = 1,npi(i)
        read(unit=1,fmt=*) temp4
        node_renum(temp4) = j
    end do
    close(1)

    call int2str(str1,i,4)
    str1 = 'edges' // trim(str1) // '.dat'
    open(unit=1,file=str1,status='old')
    open(unit=2,file="temp.dat",status='replace')
    do j = 1,nei(i)
        read(unit=1,fmt=*) temp3
        write(unit=2,fmt='(3I8)') node_renum(temp3(1)), node_renum(temp3(2)), temp3(3)
    end do
    close(1)
    close(2)

    call int2str(str1,i,4)
    str1 = 'mv temp.dat edges' // trim(str1) // '.dat'
    CALL system(str1)

    call int2str(str1,i,4)
    str1 = 'triangles' // trim(str1) // '.dat'
    open(unit=1,file=str1,status='old')
    open(unit=2,file="temp.dat",status='replace')
    do j = 1,nti(i)
        read(unit=1,fmt=*) temp3
        write(unit=2,fmt='(3I8)') node_renum(temp3(1)), node_renum(temp3(2)), node_renum(temp3(3))
    end do
    close(1)
    close(2)

    call int2str(str1,i,4)
    str1 = 'mv temp.dat triangles' // trim(str1) // '.dat'
    CALL system(str1)
end do


!!!!!!!!!!!!!  RENUMBER BOUNDARY NODES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call int2str(str1,i,4)
str1 = 'boundary_nodes.dat'
open(unit=1,file=str1,status='old')
do j = 1,nb
    read(unit=1,fmt=*) temp4
    node_renum(temp4) = j
end do
close(1)

do i=1,ndom
    call int2str(str1,i,4)
    str1 = 'nbnodes' // trim(str1) // '.dat'
    open(unit=1,file=str1,status='old')
    open(unit=2,file="temp.dat",status='replace')
    do j = 1,nbi(i)
        read(unit=1,fmt=*) temp4
        write(unit=2,fmt='(3I8)') node_renum(temp4)
    end do
    close(1)
    close(2)

    call int2str(str1,i,4)
    str1 = 'mv temp.dat bnodes' // trim(str1) // '.dat'
    CALL system(str1)
end do

deallocate(npi,nti,nei,nbi,bnodes,points,nodes,node_renum)

call getBrSCount(nbgr)   !!! Added May 3, 2014


END PROGRAM preprocmesh2_AD


!!**************************************************************************
!! Subroutine used to extract 'Br_s' matrix
!! it finds "how many domains the interface nodes is shared with" 
SUBROUTINE getBrSCount(nbgr)

    CHARACTER(len=255) :: str1 ,extension   !!filename

    INTEGER(kind=8) :: nri,nbgr,nci, npg, neg, ntg, nbg, ndom
    INTEGER(kind=8) :: i, j, Ri, Rj, pid, k
    INTEGER(kind=8), DIMENSION(:), ALLOCATABLE :: rnode12, GlobalRN, FinalCount
    INTEGER(kind=8), DIMENSION(nbgr, 3) :: BrSmatCount

    open(unit=1,file='meshdim.dat',status='old')
    read(unit=1,fmt='(5I8)') npg, neg, ntg, nbg, ndom
    close(1)

    ALLOCATE(GlobalRN(nbgr),FinalCount(nbgr))
    open(unit=2,file='remaining_nodes.dat',status='old')
    read(2,*) GlobalRN
    close(2)

    BrSmatCount(:,:) = 0
    FinalCount(:) = 0.0d0

    DO k = 1,ndom
        CALL int2str(extension,k,4)
        extension = trim(extension)
        call readcrdim(nci,nri,extension)     !! Check use of this subroutine
        ALLOCATE(rnode12(nri))

        call int2str(str1,k,4)
        str1 = 'rnodes' // trim(str1) // '.dat'
        open(unit=2,file=str1,status='old')
        READ(2,*) rnode12
        close(2)

        Do i = 1,nbgr
            Ri = GlobalRN(i)
            BrSmatCount(i:i,1) = Ri
            DO j = 1,nri
                Rj = rnode12(j)
                IF ((Ri == Rj) .AND. (BrSmatCount(i,3)== 0)) THEN
                   BrSmatCount(i,2) = k
                   BrSmatCount(i,3) = 1
                   FinalCount(i) = k
                END IF
            END DO
        END DO

        DEALLOCATE(rnode12)
        END DO

    open(unit=3,file='final_count.dat',status='replace')
    write(unit=3,fmt='(I8)') FinalCount
    close(3)

    DEALLOCATE(FinalCount, GlobalRN)

END SUBROUTINE getBrSCount

!!**************************************************************************
!! Subroutine to read 2nd level local dimensions 
!! nri : number of remaining nodes for each partition
!! nci : number of corner nodes for each partition
SUBROUTINE readcrdim(nci,nri,extension)

    integer(kind=8) :: nci, nri
    character(len=255) :: filename,extension

    filename = 'dimcrn' // trim(extension) // '.dat'
    open(unit=1,file=filename,status='old')
    read(1,*) nci, nri
    close(1)

END SUBROUTINE readcrdim

!!**************************************************************************
!! Subroutine to convert integer to string to read "*00*.dat" files 
SUBROUTINE int2str(str1,int1,n0pads)

    character(len=255) :: str1, str2
    integer(kind=8) :: int1, intdigs, i, curdig, curint
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
