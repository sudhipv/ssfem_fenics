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

PROGRAM preprocmesh3D2_AD

! use common
implicit none

!!First-Level-Data

integer, parameter :: wb=1  !! 0 for vertex grid, 1 for wire-basket

integer, parameter :: n = 255
character(len=n) :: str1, extension
integer(kind=8) :: np, ne, nt, nh, nb
integer(kind=8) :: i, j, k, net, temp4, l
integer(kind=8) :: cur_node, cur_dom, found, ndom
integer(kind=8) :: nbsus, nbicur, npicur, Ri, Rj
integer(kind=8), dimension(2) :: edges
integer(kind=8), dimension(3) :: temp3
integer(kind=8), dimension(4) :: triangle, temp4h !edge
integer(kind=8), dimension(5) :: tetrahedrons
integer(kind=8), dimension(4) :: tetrahedronsi
integer(kind=8), dimension(12) :: temp2
integer, dimension(:), allocatable :: nti, npi, nei, nhi, nbi
integer, dimension(:), allocatable :: bnodes, nodes, node_renum
double precision, dimension(:,:), allocatable :: points

!! Second-Level-Data :Coarse-Mesh
integer(kind=8) :: nci, nri, nbgr !cn, nc1, nc2, nbgc
integer, dimension(:,:), allocatable :: nodes2
integer(kind=8), dimension(:), allocatable :: countit,cnodes,rnodes,gCnodes
integer(kind=8), dimension(:), allocatable :: interface_edges, interface_triangles

INTEGER(kind=8) :: outcome, elm, nelm, ii, jj, kk, count
INTEGER(kind=8) :: nCnodes, nRnodes

!!---------------------------------------------------------------------------------------
open(unit=1,file='meshdim.dat',status='old')
read(unit=1,fmt='(7I8)') np, ne, nt, nh, nb, ndom, nbsus
close(1)

allocate(npi(ndom),nti(ndom),nei(ndom),nhi(ndom),nbi(ndom))
! Changed to 5 by sudhi since gives error with some parameters because of less allocated memory 
allocate(bnodes(nb),points(3,np),node_renum(np),nodes(floor(5.0d0*dble(np)/dble(ndom) + 10.0d0)))
allocate(interface_edges(nb),interface_triangles(nb))

npi(:) = 0
nti(:) = 0
nei(:) = 0
nhi(:) = 0
nbi(:) = 0

print*, '-----------------------------------------------------------------'
print*, '3D Local Mesh Decomposition, Renumbering & Renaming Initiated ...'
if (wb == 1) then
    print*, 'Using wire-basket based coarse grid'
    print*, '-----------------------------------------------------------------'
else
    print*, 'Using vertex based coarse grid'
    print*, '-----------------------------------------------------------------'
end if

open(unit=2,file='num_partition.dat',status='replace')
write(unit=2,fmt='(I8)') ndom
close(2)

!!!!!!!!!!!!!  BOUNDARY  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
open(unit=1,file='boundary_nodes.dat',status='old')
do i = 1,nb
    read(1,*) bnodes(i)
end do
close(1)

open(unit=1,file='tetrahedrons_on_boundary.dat',status='old')
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
                print*, 'not enough allocation for nodes in boundary nodes loop :preprocmesh3D2_AD.f90'
                call exit(0)
            end if
            nodes(nbicur) = cur_node
        end if
    end do
    nbi(j) = nbicur
    close(unit=1,status='delete')

    call int2str(str1,j,4)
    str1 = 'bnodes' // trim(str1) // '.dat'
    open(unit=2,file=str1,status='replace')
    do i = 1,nbi(j)
        write(unit=2,fmt='(I8)') nodes(i)
    end do
    close(2)
    close(1)
end do


!!!!!!!!!!!!!  EDGES  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!open(unit=1,file='edges.dat',status='old')
do i=1,ndom
    call int2str(str1,i,4)
    str1 = 'edges' // trim(str1) // '.dat'
    open(unit=10+i,file=str1,status='replace')
    write(unit=10+i,fmt='(I8)') 0
end do

!do i=1,ne
!    read(unit=1,fmt=*) edge
!    !! There is no edge dividation between subdomain
!    !nei(edge(4)) = nei(edge(4))+1
!    !nei(edge(4)) = nei(edge(4))+1
!    write(unit=10+edge(4),fmt='(3I8)') edge(1:3)
!end do
!close(1)
do i=1,ndom
    close(unit=10+i)
end do

!!!  Edge-boundary nodes (type of corner nodes)
interface_edges(:) = 0.0
open(unit=1,file='edges.dat',status='old')
count = 0
do i=1,ne
    read(unit=1,fmt=*) edges
    ii = edges(1)
    jj = edges(2)
    do j=1,nb
        if ((ii .eq. bnodes(j)) .or. (jj .eq. bnodes(j))) then
            call ismemberof(bnodes(j), nb, interface_edges, outcome)
            if (outcome .eq. 0) then
                count = count + 1
                interface_edges(count) = bnodes(j)
            end if
        end if
    end do
end do
close(1)

open(unit=3,file='edges_boundary_nodes.dat',status='replace')
write(unit=3,fmt=*) interface_edges
close(3)


!!!!!!!!!!!!!  TRIANGLES  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
open(unit=1,file='triangles.dat',status='old')
do i=1,ndom
    call int2str(str1,i,4)
    str1 = 'triangles' // trim(str1) // '.dat'
    open(unit=10+i,file=str1,status='replace')
end do

!!!  Triangles-boundary nodes (type of corner nodes)
interface_triangles(:) = 0.0
count = 0

do i=1,nt
    read(unit=1,fmt='(4I8)') triangle
    ii = triangle(1)
    jj = triangle(2)
    kk = triangle(3)
    do j=1,nb
        if ((ii .eq. bnodes(j)) .or. (jj .eq. bnodes(j)) .or. (kk .eq. bnodes(j))) then
        call ismemberof(bnodes(j), nb, interface_triangles, outcome)
            if (outcome .eq. 0) then
            count = count + 1
            interface_triangles(count) = bnodes(j)
            end if
        end if
    end do

    !! For each subdomain
    nti(triangle(4)) = nti(triangle(4))+1
    write(unit=10+triangle(4),fmt='(3I8)') triangle(1:3)
end do

close(1)
do i=1,ndom
    close(unit=10+i)
end do

open(unit=3,file='triangles_boundary_nodes.dat',status='replace')
write(unit=3,fmt=*) interface_triangles
close(3)


!!!!!!!!!!!!!  TETRAHEDRONS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
open(unit=1,file='tetrahedrons.dat',status='old')
open(unit=3,file='tetrahedrons4c.dat',status='replace')  !we need only first 4 columns
do i=1,ndom
    call int2str(str1,i,4)
    str1 = 'tetrahedrons' // trim(str1) // '.dat'
    open(unit=10+i,file=str1,status='replace')
end do

do i=1,nh
    read(unit=1,fmt='(5I8)') tetrahedrons
    write(unit=3,fmt='(4I8)') tetrahedrons(1:4)
    nhi(tetrahedrons(5)) = nhi(tetrahedrons(5))+1
    write(unit=10+tetrahedrons(5),fmt='(4I8)') tetrahedrons(1:4)
end do

close(3)
close(1)
do i=1,ndom
    close(unit=10+i)
end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! MODIFIED FOR CORNER & REMAINING NODES BY AJIT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!print*, 'Finding Corner & Remaining Nodes...'

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
!!!!!!!!!!!!!!!!!!!! NEW LOGIC TO GET CORNER & REMAINING NODES !!!!!!!!!!!!!!!!!
!! To Extracting the GLOBAL corner nodes: Node exist in more than two sub-domain
nCnodes = 0
nRnodes = 0
open(unit=7,file='corner_nodes.dat',status='replace')
open(unit=9,file='remaining_nodes.dat',status='replace')
do i=1,nb
    if (countit(i) .gt. 2.0) then
        nCnodes = nCnodes + 1
        write(unit=7,fmt='(I8)') bnodes(i)
    else
        if (wb == 1) then
            call ismemberof(bnodes(i), nb, interface_triangles, outcome) !wire-basket
        else
            call ismemberof(bnodes(i), nb, interface_edges, outcome)  !vertex-only
        end if

        if (outcome .eq. 0) then
            nRnodes = nRnodes + 1
            write(unit=9,fmt='(I8)') bnodes(i)
        else
            nCnodes = nCnodes + 1
            write(unit=7,fmt='(I8)') bnodes(i)
        end if
    end if
end do
close(7)
close(9)
open(unit=6,file='nbgc.dat',status='replace')
write(unit=6,fmt='(I8)') nCnodes
close(6)
open(unit=8,file='nbgr.dat',status='replace')
write(unit=8,fmt='(I8)') nRnodes
close(8)

deallocate(interface_edges, interface_triangles)


!!----------------------------------------------------------------------------------------
!! To Extracting the LOCAL-corner nodes: Node exist in more than two sub-domain
allocate(gCnodes(nCnodes))
open(unit=1,file='corner_nodes.dat',status='old')
do i = 1,nCnodes
    read(1,*) gCnodes(i)
end do
close(1)

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

    call int2str(str1,i,4)
    str1 = 'rnodes' // trim(str1) // '.dat'
    open(unit=2,file=str1,status='replace')

    call int2str(str1,i,4)
    str1 = 'cnodes' // trim(str1) // '.dat'
    open(unit=3,file=str1,status='replace')

    call int2str(str1,i,4)
    str1 = 'dimcrn' // trim(str1) // '.dat'
    open(unit=4,file=str1,status='replace')

    nci = 0
    nri = 0
    do j = 1,nbi(i)
        call ismemberof(nodes2(j,1), nCnodes, gCnodes, outcome)
        if (outcome .eq. 0) then
            nri = nri + 1
            write(unit=2,fmt='(I8)') nodes2(j,1)
        else
            nci = nci + 1
            write(unit=3,fmt='(I8)') nodes2(j,1)
        end if
    end do
    close(2)
    close(3)

    write(unit=4,fmt='(2I8)') nci, nri
    close(4)
    deallocate(nodes2)
end do

deallocate(gCnodes, countit)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!print*, 'Re-arranging Initiated...'

do i=1,ndom

    call int2str(str1,i,4)
    str1 = 'dimcrn' // trim(str1) // '.dat'
    open(unit=1,file=str1,status='old')
    read(1,fmt='(2I8)') nci, nri
    close(1)

    allocate(rnodes(nri),cnodes(nci))

    call int2str(str1,i,4)
    str1 = 'nbnodes' // trim(str1) // '.dat'
    open(unit=4,file=str1,status='replace')

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
    open(unit=4,file=str1,position="append", status='old')

    do j = 1,nci
       read(3,*) cnodes(j)
       write(unit=4,fmt='(I8)') cnodes(j)
    end do
    close(3)
    close(4)
    deallocate(rnodes, cnodes)
end do

!!!!!!!!!! UNTIL HERE MODIFIED BY AJIT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!  POINTS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=1,ndom
    call int2str(str1,i,4)
    str1 = 'nbnodes' // trim(str1) // '.dat'
    open(unit=1,file=str1,status='old')
    do j = 1,nbi(i)
        read(1,*) nodes(j)
    end do
    close(1)

    ! print*, size(nodes)
    ! print*, nodes
    call int2str(str1,i,4)
    str1 = 'tetrahedrons' // trim(str1) // '.dat'
    open(unit=1,file=str1,status='old')
    call int2str(str1,i,4)
    str1 = 'nodes' // trim(str1) // '.dat'
    open(unit=2,file=str1,status='replace')

    npicur = nbi(i)
    ! print*, "npicur is "
    ! print*, npicur
    do j = 1,nhi(i)
        read(unit=1,fmt='(4I8)') tetrahedronsi
        do l = 1,4
        found = 0
            do k = 1,npicur
                if (nodes(k) .eq. tetrahedronsi(l)) then
                    found = 1
                    exit
                end if
            end do
            if (found .eq. 0) then
                npicur = npicur+1
                if (npicur .gt. size(nodes)) then
                    print*, 'not enough allocation for nodes in tetrahedrons loop : preprocmesh3D2_AD.f90'
                    call exit(0)
                end if
                nodes(npicur) = tetrahedronsi(l)
                write(unit=2,fmt='(I8)') tetrahedronsi(l)
            end if
        end do
    end do
    do j = 1,nbi(i)
        write(unit=2,fmt='(I8)') nodes(j)
    end do

    close(1)
    close(2)
    npi(i) = npicur
    !print*, "final npicur"
    !print*, npicur
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
        write(unit=2,fmt='(3ES17.8E3)') points(:,temp4)
    end do

    close(1)
    close(2)
end do

!!!!!!!!!!!!!  MESH DIMENSIONS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=1,ndom
    call int2str(str1,i,4)
    str1 = 'meshdim' // trim(str1) // '.dat'
    open(unit=1,file=str1,status='replace')
    write(unit=1,fmt='(5I8)') npi(i), nei(i), nti(i), nhi(i), nbi(i)
    close(1)
end do

!!!!!!!!!!!!!  RENUMBER EDGES, TRIANGLES & TETRADEDRONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!print*, 'Re-numbering Initiated...'

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


    call int2str(str1,i,4)
    str1 = 'tetrahedrons' // trim(str1) // '.dat'
    open(unit=1,file=str1,status='old')
    open(unit=2,file="temp.dat",status='replace')
    do j = 1,nhi(i)
        read(unit=1,fmt=*) temp4h
        write(unit=2,fmt='(4I8)') node_renum(temp4h(1)), node_renum(temp4h(2)), node_renum(temp4h(3)), node_renum(temp4h(4))
    end do
    close(1)
    close(2)

    call int2str(str1,i,4)
    str1 = 'mv temp.dat tetrahedrons' // trim(str1) // '.dat'

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

deallocate(npi,nti,nei,nhi,nbi,bnodes,points,nodes,node_renum)

call getBrSCount(nRnodes)   !!! Added May 3, 2014

END PROGRAM preprocmesh3D2_AD


!!**************************************************************************
SUBROUTINE ismemberof(elm, nelm, set, outcome)

    INTEGER(kind=8) :: outcome, nelm, i, setElm
    INTEGER(kind=4) :: elm
    INTEGER(kind=8), DIMENSION(nelm) :: set

    outcome = 0
    DO i = 1,nelm
        setElm = set(i)
        IF (setElm == elm) then
            outcome = outcome+1
        END IF
    END DO

END SUBROUTINE ismemberof


!!**************************************************************************
!! Subroutine used to extract 'Br_s' matrix
!! it finds "how many domains the interface nodes is shared with"
SUBROUTINE getBrSCount(nbgr)

    CHARACTER(len=255) :: str1 ,extension   !!filename

    INTEGER(kind=8) :: nri,nbgr,nci, npg, neg, ntg, nbg, ndom, nhg
    INTEGER(kind=8) :: i, j, Ri, Rj, pid, k
    INTEGER(kind=8), DIMENSION(:), ALLOCATABLE :: rnode12, GlobalRN, FinalCount
    INTEGER(kind=8), DIMENSION(nbgr, 3) :: BrSmatCount

    open(unit=1,file='meshdim.dat',status='old')
    read(unit=1,fmt='(6I8)') npg, neg, ntg, nhg, nbg, ndom
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

    character(len=255) :: str1
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
