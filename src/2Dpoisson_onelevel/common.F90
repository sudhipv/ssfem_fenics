!! Common Module: By Ajit Desai, 8/2015, Mohammad Khalil, 1/2012
!! purpose: it contains various subrountines used in this package 
!!
!!-----------------------------------------------------------------------------------

module common

implicit none

contains

!!*********************************************************
!! Subroutine 
subroutine construct_coln_boundary(Ub_g,U_g)

  character(len=255) :: str1
  integer :: nbg, i, temp1
  double precision, dimension(:) :: Ub_g, U_g

  nbg = size(Ub_g)

  str1 = '../../data/meshData/boundary_nodes.dat'    !!**
  open(unit=1,file=str1,status='old')
  do i = 1,nbg
    read(unit=1,fmt=*) temp1
    U_g(temp1) = Ub_g(i)
  end do
  close(1)

end subroutine construct_coln_boundary

!!*********************************************************
!! Subroutine 
subroutine construct_coln(pid,Ui,U_gi)

  character(len=255) :: str1
  integer :: npi, i, pid, temp1
  double precision, dimension(:) :: Ui, U_gi

  npi = size(Ui)

  call int2str(str1,pid+1,4)
  str1 = '../../data/meshData/nodes' // trim(str1) // '.dat'   !!**
  open(unit=1,file=str1,status='old')
  U_gi(:) = 0.0d0
  do i = 1,npi
    read(unit=1,fmt=*) temp1
    U_gi(temp1) = Ui(i)
  end do
  close(1)

end subroutine construct_coln

!!*********************************************************
!! Subroutine 
subroutine int2str(str1,int1,n0pads)

  character(len=255) :: str1
  integer :: int1, intdigs, i, curdig, curint,n0pads

  intdigs = floor(log10(dble(int1))) + 1
  curint = int1
  str1 = ''
  do i = 1,intdigs
    curdig = mod(curint,10)
    curint = floor(dble(curint)/10.0d0)
    str1 = trim(achar(48 + curdig) // str1)
  enddo
  do i = intdigs+1,n0pads
    str1 = trim('0' // str1)
  enddo

end subroutine int2str

!!*********************************************************
!! Subroutine 
subroutine writevtk(filename,p,t,u,np,nt)

integer :: np, nt, i
double precision, dimension(np) :: u
double precision, dimension(2,np) :: p
integer, dimension(3,nt) :: t
character(len=255) :: filename, str1, str2

open(unit=1001,file=filename)
write(1001,'(a)') '# vtk DataFile Version 1.0'
write(1001,'(a)') 'LDFEMLSS'
write(1001,'(a)') 'ASCII'
write(1001,'(a)') ' '
write(1001,'(a)') 'DATASET POLYDATA'
call int2str(str1,np,0)
write(1001,'(a)') 'POINTS     ' // trim(str1) // '  double'
do i = 1,np
    write(1001,*) p(1,i), p(2,i), 0.0d0
enddo
write(1001,'(a)') ' '
call int2str(str1,nt,0)
call int2str(str2,4*nt,0)
write(1001,'(a)') 'POLYGONS   ' // trim(str1) // '    ' // trim(str2)
do i = 1,nt
    write(1001,*) 3 , t(1,i)-1 ,t(2,i)-1, t(3,i)-1
enddo
write(1001,'(a)') ' '
call int2str(str1,np,0)
write(1001,'(a)') 'POINT_DATA    ' // trim(str1)
write(1001,'(a)') 'SCALARS u double'
write(1001,'(a)') 'LOOKUP_TABLE default'
do i = 1,np
    write(1001,*) u(i)
enddo
close(1001)

end subroutine writevtk

!!*********************************************************
!! Subroutine 
subroutine readmeshdata(p,e,t,np,ne,nt,extension)          !!**

integer :: np, ne, nt
double precision, dimension(2,np) :: p
integer, dimension(3,ne) :: e
integer, dimension(3,nt) :: t
integer :: i
character(len=255) :: extension, filename

filename = '../../data/meshData/points' // trim(extension) // '.dat'     !!**
open(unit=1,file=filename,status='old')
do i = 1,np
    read(unit=1,fmt=*) p(:,i)
end do
close(1)

filename = '../../data/meshData/edges' // trim(extension) // '.dat'      !!**
open(unit=1,file=filename,status='old')
do i = 1,ne
    read(unit=1,fmt=*) e(:,i)
end do
close(1)

filename = '../../data/meshData/triangles' // trim(extension) // '.dat'   !!**
open(unit=1,file=filename,status='old')
do i = 1,nt
    read(unit=1,fmt=*) t(:,i)
end do
close(1)

end subroutine readmeshdata

!!*********************************************************
!! Subroutine 
subroutine readmeshdim(np,ne,nt,nb,extension)

integer :: np, ne, nt, nb
character(len=255) :: filename,extension

filename = '../../data/meshData/meshdim' // trim(extension) // '.dat'  !!**
open(unit=1,file=filename,status='old')
read(1,*) np, ne, nt, nb
close(1)

end subroutine readmeshdim

!!*********************************************************
!! Subroutine
subroutine readGlobMeshDim(nParts)

    integer :: nParts, aa, bb, cc, dd

    open(unit=1,file='../../data/meshData/meshdim.dat',status='old')    !!**
    read(unit=1,fmt='(6I8)') aa, bb, cc, dd, nParts
    close(1)

end subroutine readGlobMeshDim

!!*********************************************************
!! Subroutine
subroutine readMeshPart(nParts)

    integer :: nParts
    open(unit=1,file='../../data/meshData/num_partition.dat',status='old')    !!**
    read(1,*) nParts
    close(1)

end subroutine readMeshPart


!!*********************************************************
!! Subroutine
subroutine readPceData(nord, ndim, npceout, npcein)

integer :: nord, ndim, npceout, npcein

open(unit=2,file='../../data/klePceData/pcedata.dat',status='old')      !!**
read(unit=2,fmt='(4I8)') nord, ndim, npceout, npcein
close(2)

end subroutine readPceData


!!*********************************************************
!! Subroutine
subroutine readIndices(nDim,npceout,mIndex,sIndex)

    integer  :: i, nDim, npceout
    integer, dimension(2,ndim)  :: sIndex
    integer, dimension(nDim,npceout)  :: mIndex


    open(unit=1,file='../../data/klePceData/multiIndex.dat',status='old')  !!**
    read(1,*) mIndex
    close(1)

    open(unit=2,file='../../data/klePceData/sortIndex.dat',status='old')   !!**
    do i = 1,nDim
    read(unit=2,fmt='(2I8)') sIndex(1,i), sIndex(2,i)
    end do
    close(2)

end subroutine readIndices


!!*********************************************************
!! Subroutine 
subroutine readprocess(nomga,omegas,multipliers)

integer :: nomga
double precision, dimension(nomga) :: omegas, multipliers

open(unit=1,file='../../data/klePceData/omegas',status='old')       !!**
read(1,*) omegas
close(1)
open(unit=1,file='../../data/klePceData/multipliers',status='old')  !!**
read(1,*) multipliers
close(1)

end subroutine readprocess

!!*********************************************************
!! Subroutine 
subroutine readdbounds(dbounds)

double precision, dimension(2,2) :: dbounds

open(unit=1,file='../../data/meshData/dbounds.dat',status='old')  !!**
read(1,*) dbounds(1,:)
read(1,*) dbounds(2,:)
close(1)

end subroutine readdbounds

!!*********************************************************
!! Subroutine 
subroutine readncijk(ncijk)

integer :: ncijk

open(unit=1,file="../../data/klePceData/cijk",status='old')    !!**
read(unit=1,fmt=*) ncijk
close(1)

end subroutine readncijk

!!*********************************************************
!! Subroutine 
subroutine readcijk(ncijk,ijk,cijk)

integer :: ncijk, i
double precision, dimension(ncijk) :: cijk
integer, dimension(ncijk,3) :: ijk

open(unit=1,file="../../data/klePceData/cijk",status='old')      !!**
read(unit=1,fmt=*) i
do i = 1,ncijk
   read(unit=1,fmt=*) ijk(i,:), cijk(i)
end do
close(1)

end subroutine readcijk

!!*********************************************************
!! Subroutine 
subroutine invert(n,A,Ainv)

double precision :: A(n,n),Ainv(n,n)
double precision, dimension(n) :: work
integer, dimension(n) :: ipvt
integer :: info, n

Ainv = A
call dgetrf(n,n,Ainv,n,ipvt,info)

call dgetri(n,Ainv,n,ipvt,work,n,info)

end subroutine invert

!!*********************************************************
!! Subroutine 
subroutine solve(n,nrhs,A,u,f)

double precision :: A(n,n), Atemp(n,n)
double precision :: u(n,nrhs),f(n,nrhs)
integer, dimension(n) :: ipvt
integer :: info, n, nrhs

Atemp = A
u = f
call dgesv(n,nrhs,Atemp,n,ipvt,u,n,info)

end subroutine solve

!!*********************************************************
!! Subroutine 
subroutine minmax(in,minin,maxin)

double precision, dimension(:) :: in
double precision :: minin,maxin
integer :: length, i

length = size(in)

minin = in(1)
maxin = in(1)
do i = 2,length
    if (minin > in(i)) then
        minin = in(i)
    else if (maxin < in(i)) then
        maxin = in(i)
    end if
end do

end subroutine minmax


!!*********************************************************
!! Subroutine 
subroutine covariance(in,nr,nc,out)

integer :: nr, nc
double precision, dimension(:,:) :: in
double precision, dimension(:,:) :: out
double precision :: meani, meanj, sum
integer :: i, j

do i=1, nr
do j = 1,nr
sum = dot_product(in(i,:),in(j,:))
call mean(in(i,:),meani)
call mean(in(j,:),meanj)
out(i,j) = (sum - dble(nc)*meani*meanj)/dble(nc-1)
end do
end do

end subroutine covariance

!!*********************************************************
!! Subroutine 
subroutine outer_product(a,b,out)

double precision, dimension(:) :: a,b
double precision, dimension(:,:) :: out
integer :: length_i, length_j, i, j

length_i = size(a)
length_j = size(b)

do i = 1,length_i
do j = 1,length_j
out(i,j) = a(i)*b(j)
end do
end do

end subroutine outer_product

!!*********************************************************
!! Subroutine 
subroutine diag(in,out)

double precision, dimension(:) :: in
double precision, dimension(:,:) :: out
integer :: length, i,j

length = size(in)

out = 0.0d0

do i = 1,length
do j = 1,length
out(i,j) = 0.0d0
end do
out(i,i) = in(i)
end do

end subroutine diag

!!*********************************************************
!! Subroutine 
subroutine variance(in,out)

double precision, dimension(:) :: in
double precision :: out, meani
integer :: length, i

length = size(in)

out = 0.0d0

do i = 1,length
out = out + in(i)**2
end do

call mean(in,meani)

out = (out - dble(length)*meani**2)/dble(length-1)

end subroutine variance

!!!*********************************************************
!!! Subroutine 
!subroutine variance_matrix(in,out)
!
!double precision, dimension(:,:) :: in
!double precision :: out, meani
!integer :: i, j
!integer :: length_i, length_j
!
!length_i = size(in,1)
!length_j = size(in,2)
!
!out = 0.0d0
!
!do i = 1,length_i
!    do j = 1,length_j
!        out = out + in(i,j)**2
!    end do
!end do
!
!call mean_matrix(in,meani)
!
!out = (out - dble(length_i*length_j)*meani**2)/dble(length_i*length_j - 1)
!
!end subroutine variance_matrix

!!*********************************************************
!! Subroutine 
subroutine mean(in,out)

double precision, dimension(:) :: in
double precision :: out
integer :: length, i

length = size(in)

out = 0.0d0

do i = 1,length
out = out + in(i)
end do

out = out/dble(length)

end subroutine mean

!!!*********************************************************
!!! Subroutine 
!subroutine mean_matrix(in,out)
!
!double precision, dimension(:,:) :: in
!double precision :: out
!integer :: i, j
!integer :: length_i, length_j
!
!length_i = size(in,1)
!length_j = size(in,2)
!
!out = 0.0d0
!
!do i = 1,length_i
!do j = 1,length_j
!out = out + in(i,j)
!end do
!end do
!
!out = out/dble(length_i*length_j)
!
!end subroutine mean_matrix
!
!!!*********************************************************
!!! Subroutine 
!subroutine linspace_d(a,b,c,out)
!
!double precision :: a, b, delta
!integer :: c, i
!double precision, dimension(c) :: out
!
!delta = (b - a)/dble(c-1)
!
!do i=1,c
!out(i) = a + (b-a)*(dble(i-1)/(c-1))
!end do
!
!end subroutine linspace_d

!!*********************************************************
!! Subroutine 
subroutine randperm(out)

double precision :: sample
integer :: i, index, left
integer, dimension(:) :: out
integer :: temp, length

length = size(out)

do i=1,length
out(i) = i
end do

left = length
do i=1,length
call random_number(sample)
index = length - floor(sample*dble(left))

temp = out(i)
out(i) = out(index)
out(index) = temp

left = left - 1
end do

end subroutine randperm

!!*********************************************************
!! Subroutine 
subroutine myerfinv(inout)

double precision, dimension(:) :: inout
double precision, dimension(4) :: a, b, c
double precision, dimension(2) :: d
double precision :: y0, temp, z
integer :: i

a = [ 0.886226899d0, -1.645349621d0,  0.914624893d0, -0.140543331d0];
b = [-2.118377725d0,  1.442710462d0, -0.329097515d0,  0.012229801d0];
c = [-1.970840454d0, -1.624906493d0,  3.429567803d0,  1.641345311d0];
d = [ 3.543889200d0,  1.637067800d0];

y0 = 0.7d0

do i = 1,size(inout)
temp = inout(i)
if (abs(temp) <= y0) then
z = temp*temp
inout(i) = temp*(((a(4)*z+a(3))*z+a(2))*z+a(1))/((((b(4)*z+b(3))*z+b(2))*z+b(1))*z+1)
else if (y0 < temp) then
z = sqrt(-log((1-temp)/2))
inout(i) = (((c(4)*z+c(3))*z+c(2))*z+c(1))/((d(2)*z+d(1))*z+1)
else
z = sqrt(-log((1+temp)/2))
inout(i) = -(((c(4)*z+c(3))*z+c(2))*z+c(1))/((d(2)*z+d(1))*z+1)
end if
end do

end subroutine myerfinv

!!*********************************************************
!! Subroutine 
subroutine norm_gen(out,n)

   double precision, dimension(:) :: out
   integer :: n

   integer, dimension(n) :: perm
!   double precision, dimension(n) :: rands

   integer :: length

   length = n

   call randperm(perm)
!    call random_number(rands)
!    out = 2.0d0*(dble(perm) - rands)/length - 1.0d0
   out = 2.0d0*(dble(perm) - 0.5d0)/length - 1.0d0

   call myerfinv(out)

   out = sqrt(2.0d0)*out
      
end subroutine norm_gen

!!*********************************************************
!! Subroutine 
subroutine getub(pid,nb,nbg,npceout,Ub_g,Ub)

   character(len=255) :: str1
   integer :: nb, nbg, i, k, pid, npceout
   double precision, dimension(nbg*npceout) :: Ub_g
   double precision, dimension(nb*npceout) :: Ub
   integer, dimension(nb) :: bnodes

   call int2str(str1,pid+1,4)
   str1 = '../../data/meshData/bnodes' // trim(str1) // '.dat'    !!**
   open(unit=1,file=str1,status='old')
   do i = 1,nb
      read(unit=1,fmt=*) bnodes(i)
   end do
   close(1)
   do k = 1,npceout
      do i = 1,nb
         Ub((k-1)*nb+i) = Ub_g((k-1)*nbg+bnodes(i))
      end do
   end do

end subroutine getub

!!*********************************************************
!! Subroutine 
subroutine getubg(pid,nb,nbg,npceout,Ub_g,Ub)

   character(len=255) :: str1
   integer :: nb, nbg, i, k, pid, npceout
   double precision, dimension(nbg*npceout) :: Ub_g
   double precision, dimension(nb*npceout) :: Ub
   integer, dimension(nb) :: bnodes

   call int2str(str1,pid+1,4)
   str1 = '../../data/meshData/bnodes' // trim(str1) // '.dat'   !!**
   open(unit=1,file=str1,status='old')
   do i = 1,nb
      read(unit=1,fmt=*) bnodes(i)
   end do
   close(1)
   Ub_g = 0.0d0
   do k = 1,npceout
      do i = 1,nb
         Ub_g((k-1)*nbg+bnodes(i)) = Ub((k-1)*nb+i) 
      end do
   end do

end subroutine getubg

!!*********************************************************
!! Subroutine 
subroutine multiply(A,b,c,n1,n2,n3)

integer :: n1,n2,n3
double precision :: A(n1,n2),b(n2,n3),c(n1,n3)

call dgemm('N','N',n1,n3,n2,1.0d0,A,n1,b,n2,0.0d0,c,n1)

end subroutine multiply

!!!*********************************************************
!!! Subroutine
!subroutine pcgm(pid,n,A,U,F,maxiter,tol,set,Minv)
!
!   integer :: n,maxiter,set
!   double precision :: tol
!   double precision :: A(n,n),F(n),U(n)
!
!   integer :: i, pid
!   double precision :: rho_next,rho_curr, alpha, beta, err, norm1
!   double precision :: r(n),p(n),q(n),z(n),Minv(n,n)
!
!   if (set .eq. 1) then
!      U(:) = 0.0d0
!   end if
!
!   call multiply(A,U,r,n,n,1)
!   r = F - r
!
!   call multiply(Minv,r,z,n,n,1)
!
!   p = z
!   rho_next = dot_product(r,z)
!
!   do i = 1,maxiter
!      call multiply(A,p,q,n,n,1)
!      rho_curr = rho_next
!      alpha = rho_curr/dot_product(q,p)
!      err = alpha*alpha*dot_product(p,p)/dot_product(U,U)
!      U = U + alpha*p
!      r = r - alpha*q
!      !!if (pid .eq. 0) print*, 'iteration # ', i, ', relative error of ', err
!      if (err .lt. tol) exit
!      call multiply(Minv,r,z,n,n,1)
!      rho_next = dot_product(r,z)
!      beta = rho_next/rho_curr
!      p = z + beta*p
!   end do
!
!end subroutine pcgm
!
!!!*********************************************************
!!! Subroutine
!subroutine create_r(pid,nb,nbg,npceout,Rmat)
!
!  character(len=255) :: str1
!  integer :: nb, nbg, i, pid, temp1, npceout
!  double precision, dimension(:,:) :: Rmat
!
!  call int2str(str1,pid+1,4)
!  str1 = 'bnodes' // trim(str1) // '.dat'
!  open(unit=1,file=str1,status='old')
!  Rmat(:,:) = 0.0d0
!  do i = 1,nb
!    read(unit=1,fmt=*) temp1
!    Rmat(i,temp1) = 1
!  end do
!  close(1)
!  do i = 2,npceout
!    Rmat((i-1)*nb+1:i*nb,(i-1)*nbg+1:i*nbg) = Rmat(1:nb,1:nbg)
!  end do
!
!end subroutine create_r
!
!!!*********************************************************
!!! Subroutine
!subroutine histcount(outpos,outn,in,nhist,hist_min,hist_max)
!
!integer :: nhist, length
!double precision, dimension(:) :: in
!double precision, dimension(:) :: outpos
!integer, dimension(:) :: outn
!double precision, dimension(nhist+1) :: edges
!double precision :: hist_min, hist_max
!integer :: i, index1
!double precision :: current
!
!length = size(in)
!
!call linspace_d(hist_min,hist_max,nhist+1,edges)
!
!call linspace_d(hist_min,hist_max,nhist,outpos)
!
!outn = 0
!do i=1,length
!index1 = 1
!current = in(i)
!do while ((current > edges(index1+1)) .and. index1 < nhist)
!index1 = index1 + 1
!end do
!outn(index1) = outn(index1) + 1
!end do
!
!end subroutine histcount
!
!!!*********************************************************
!!! Subroutine
!subroutine hist(outpos,outn,in,nhist)
!
!integer :: nhist, length
!double precision, dimension(:) :: in
!double precision, dimension(:) :: outpos, outn
!double precision, dimension(nhist+1) :: edges
!double precision :: minin, maxin
!integer :: i, index1
!double precision :: current
!
!length = size(in)
!
!call minmax(in,minin,maxin)
!
!call linspace_d(minin,maxin,nhist+1,edges)
!
!call linspace_d(minin,maxin,nhist,outpos)
!
!outn = 0.0d0
!do i=1,length
!index1 = 1
!current = in(i)
!do while ((current > edges(index1+1)) .and. index1 < nhist)
!index1 = index1 + 1
!end do
!outn(index1) = outn(index1) + 1.0d0
!end do
!
!outn = outn / (dble(length)*(maxin - minin)/dble(nhist - 1))
!
!end subroutine hist
!
!!!*********************************************************
!!! Subroutine
!subroutine rms_error_state(est,exact,rms)
!
!double precision, dimension(:) :: est, exact
!double precision :: rms
!
!rms = sqrt(dot_product(exact-est,exact-est)/dot_product(exact,exact))
!
!end subroutine rms_error_state
!
!!!*********************************************************
!!! Subroutine
!subroutine rms_error_param(est,exact,rms)
!
!double precision, dimension(:) :: est
!double precision :: rms, exact
!
!rms = sqrt(dot_product(exact-est,exact-est)/(dble(size(est))*(exact**2)))
!
!end subroutine rms_error_param
!
!
!!!*********************************************************
subroutine readnmeas(nmeas,extension)

integer :: nmeas
character(len=255) :: filename,extension

filename = '../../data/meshData/nmeas' // trim(extension) // '.dat'
open(unit=1,file=filename,status='old')
read(1,*) nmeas
close(1)

end subroutine readnmeas

!!!*********************************************************
subroutine readmeasnodes(nmeas,measnodes,extension)

integer, dimension(nmeas) :: measnodes
integer :: nmeas, i
character(len=255) :: filename,extension

filename = '../../data/meshData/measnodes' // trim(extension) // '.dat'
open(unit=1,file=filename,status='old')
do i = 1,nmeas
   read(1,*) measnodes(i)
end do
close(1)

end subroutine readmeasnodes

!!!*********************************************************
subroutine create_df_enkf(nmeas,nmc,np,measnodes,Aenkf,Df)

double precision, dimension(nmeas,nmc) :: Df
double precision, dimension(np,nmc) :: Aenkf

integer, dimension(nmeas) :: measnodes
integer :: nmeas,i,nmc,np

do i = 1,nmc
   Df(:,i) = Aenkf(measnodes,i)
end do

end subroutine create_df_enkf


subroutine create_d(nmeas,np,measnodes,u,d,gamma)

double precision, dimension(np) :: u
double precision, dimension(nmeas) :: d
double precision, dimension(nmeas) :: error


integer, dimension(nmeas) :: measnodes
integer :: nmeas,np
double precision :: gamma

! call norm_gen(error,nmeas)
call random_number(error)

d = u(measnodes) + sqrt(gamma)*error

end subroutine create_d

subroutine read_d(nmeas,extension,d)

integer :: nmeas,i
double precision, dimension(nmeas) :: d
character(len=255) :: str1, extension

str1 = '../../data/meshData/meas' // trim(extension) // '.dat'
open(unit=1,file=str1,status='old')
do i = 1,nmeas
   read(unit=1,fmt=*) d(i)
end do
close(1)

end subroutine read_d


subroutine create_d_pce_rv(nmeas,np,extension,npceout,measnodes,Af,d,gamma,xi)

double precision, dimension(np,npceout+nmeas) :: Af
double precision, dimension(nmeas) :: d
double precision, dimension(nmeas) :: error

integer, dimension(nmeas) :: measnodes
integer :: nmeas,np,i,npceout
double precision :: gamma,xi

character(len=255) :: str1,extension

call random_number(error)
d = Af(measnodes,1) + Af(measnodes,2)*xi + Af(measnodes,3)*(xi**2-1.0d0)/sqrt(2.0d0) + &
   Af(measnodes,4)*(xi**3-3.0d0*xi)/sqrt(6.0d0) + sqrt(gamma)*error

str1 = '../../data/meshData/meas' // trim(extension) // '.dat'
open(unit=1,file=str1,status='replace')
do i = 1,nmeas
   write(unit=1,fmt='(1ES17.8E3)') d(i)
end do
close(1)

end subroutine create_d_pce_rv


subroutine create_d_pce_sp(nmeas,np,extension,npceout,measnodes,Af,d,gamma,xi)

double precision, dimension(np,npceout+nmeas) :: Af
double precision, dimension(nmeas) :: d
double precision, dimension(nmeas) :: error

integer, dimension(nmeas) :: measnodes
integer :: nmeas,np,i,npceout
double precision :: gamma,psi
double precision, dimension(3) :: xi

character(len=255) :: str1,extension

d(:) = 0.0d0
do i = 1,npceout
   call evaluate_psi_sp(i,xi,psi)
   d = d + Af(measnodes,i)*psi
end do

call random_number(error)
d = d + sqrt(gamma)*error

str1 = '../../data/meshData/meas' // trim(extension) // '.dat'
open(unit=1,file=str1,status='replace')
do i = 1,nmeas
   write(unit=1,fmt='(1ES17.8E3)') d(i)
end do
close(1)

end subroutine create_d_pce_sp

!!!*********************************************************
subroutine evaluate_psi_sp(term,xi,psi)

double precision :: psi
integer :: term
double precision, dimension(3) :: xi

SELECT CASE (term)
   CASE (1)
      psi = 1.0d0
   CASE (2)
      psi = xi(1)
   CASE (3)
      psi = xi(2)
   CASE (4)
      psi = xi(3)
   CASE (5)
      psi = (xi(1)**2 - 1.0d0)/sqrt(2.0d0)
   CASE (6)
      psi = (xi(2)**2 - 1.0d0)/sqrt(2.0d0)
   CASE (7)
      psi = (xi(3)**2 - 1.0d0)/sqrt(2.0d0)
   CASE (8)
      psi = xi(1)*xi(2)
   CASE (9)
      psi = xi(1)*xi(3)
   CASE (10)
      psi = xi(2)*xi(3)
   CASE (11)
      psi = (xi(1)**3 - 3.0d0*xi(1))/sqrt(6.0d0)
   CASE (12)
      psi = (xi(2)**3 - 3.0d0*xi(2))/sqrt(6.0d0)
   CASE (13)
      psi = (xi(3)**3 - 3.0d0*xi(3))/sqrt(6.0d0)
   CASE (14)
      psi = (xi(1)*xi(1)*xi(2) - xi(2))/sqrt(2.0d0)
   CASE (15)
      psi = (xi(1)*xi(1)*xi(3) - xi(3))/sqrt(2.0d0)
   CASE (16)
      psi = (xi(2)*xi(2)*xi(1) - xi(1))/sqrt(2.0d0)
   CASE (17)
      psi = (xi(2)*xi(2)*xi(3) - xi(3))/sqrt(2.0d0)
   CASE (18)
      psi = (xi(3)*xi(3)*xi(1) - xi(1))/sqrt(2.0d0)
   CASE (19)
      psi = (xi(3)*xi(3)*xi(2) - xi(2))/sqrt(2.0d0)
   CASE (20)
      psi = xi(1)*xi(2)*xi(3)
END SELECT

end subroutine evaluate_psi_sp

!!!*********************************************************
!subroutine solve_woodbury(n,n1,n2,A,u,f)
!
!double precision :: A(n,n), Atemp(n,n), S(n1,n1), out1(n2,n1)
!double precision :: u(n),f(n), u1(n1), u2(n2), f1(n1), f2(n2)
!integer, dimension(n) :: ipvt
!integer :: info, n, nrhs, n1, n2
!
!call solve(n2,n1,A((n1+1):n,(n1+1):n),out1,A((n1+1):n,1:n1))
!
!S = A(1:n1,1:n1) - matmul(A(1:n1,(n1+1):n),out1)
!
!call solve(n2,1,A((n1+1):n,(n1+1):n),u2,f((n1+1):n))
!
!f1 = f(1:n1) - matmul(A(1:n1,(n1+1):n),u2)
!
!call solve(n1,1,S,u1,f1)
!
!f2 = f((n1+1):n) - matmul(A((n1+1):n,1:n1),u1)
!call solve(n2,1,A((n1+1):n,(n1+1):n),u2,f2)
!
!u(1:n1) = u1
!u((n1+1):n2) = u2
!
!end subroutine solve_woodbury

!!!*********************************************************
!subroutine create_realization_rv(np,npceout,nmeas,Af,uvec,xi)
!
!double precision, dimension(np,npceout+nmeas) :: Af
!double precision, dimension(np) :: uvec
!double precision :: xi
!
!integer :: np,npceout,nmeas
!
!uvec = Af(:,1) + Af(:,2)*xi + Af(:,3)*(xi**2-1.0d0)/sqrt(2.0d0) + &
!   Af(:,4)*(xi**3-3.0d0*xi)/sqrt(6.0d0)
!
!end subroutine create_realization_rv


!subroutine create_realization_rv_i(np,nb,npceout,nmeas,Af,uveci,xi)
!
!double precision, dimension(np,npceout+nmeas) :: Af
!double precision, dimension(np-nb) :: uveci
!double precision :: xi
!
!integer :: np,npceout,nmeas,nb
!
!uveci = Af(1:(np-nb),1) + Af(1:(np-nb),2)*xi + Af(1:(np-nb),3)*(xi**2-1.0d0)/sqrt(2.0d0) + &
!   Af(1:(np-nb),4)*(xi**3-3.0d0*xi)/sqrt(6.0d0)
!
!end subroutine create_realization_rv_i
!
!subroutine create_realization_rv_b(nbg,npceout,nmeas,Afb,uvecb,xi)
!
!double precision, dimension(nbg,npceout+nmeas) :: Afb
!double precision, dimension(nbg) :: uvecb
!double precision :: xi
!
!integer :: nbg,npceout,nmeas
!
!uvecb = Afb(:,1) + Afb(:,2)*xi + Afb(:,3)*(xi**2-1.0d0)/sqrt(2.0d0) + &
!   Afb(:,4)*(xi**3-3.0d0*xi)/sqrt(6.0d0)
!
!end subroutine create_realization_rv_b

!!!*********************************************************
subroutine create_realization_sp(np,npceout,nmeas,Af,uvec,xi)

double precision, dimension(np,npceout+nmeas) :: Af
double precision, dimension(np) :: uvec
double precision, dimension(3) :: xi
double precision :: psi

integer :: np,npceout,nmeas,i

uvec(:) = 0.0d0
do i = 1,npceout
   call evaluate_psi_sp(i,xi,psi)
   uvec = uvec + Af(:,i)*psi
end do

end subroutine create_realization_sp


!!!*********************************************************
subroutine create_realization_sp_i(np,nb,npceout,nmeas,Af,uveci,xi)

double precision, dimension(np,npceout+nmeas) :: Af
double precision, dimension(np-nb) :: uveci
double precision, dimension(3) :: xi
double precision :: psi

integer :: np,npceout,nmeas,nb,i

uveci(:) = 0.0d0
do i = 1,npceout
   call evaluate_psi_sp(i,xi,psi)
   uveci = uveci + Af(1:(np-nb),i)*psi
end do

end subroutine create_realization_sp_i


!!!*********************************************************
subroutine create_realization_sp_b(nbg,npceout,nmeas,Afb,uvecb,xi)

double precision, dimension(nbg,npceout+nmeas) :: Afb
double precision, dimension(nbg) :: uvecb
double precision, dimension(3) :: xi
double precision :: psi

integer :: nbg,npceout,nmeas,i

uvecb(:) = 0.0d0
do i = 1,npceout
   call evaluate_psi_sp(i,xi,psi)
   uvecb = uvecb + Afb(:,i)*psi
end do

end subroutine create_realization_sp_b


!!!*********************************************************
subroutine create_dp_enkf(nmeas,nmc,d,Dp,gamma)

double precision, dimension(nmeas,nmc) :: Dp
double precision, dimension(nmeas) :: d
double precision, dimension(nmc) :: error

integer :: nmeas,nmc,i
double precision :: gamma

do i = 1,nmeas
   call norm_gen(error,nmc)
!    call random_number(error)
   Dp(i,:) = d(i) + sqrt(gamma)*error
end do

end subroutine create_dp_enkf

!!!!*********************************************************
!subroutine create_cmat_enkf(nmeas,nmc,Df,Dp,Cmat)
!
!double precision, dimension(nmc,nmc) :: Cmat
!double precision, dimension(nmeas,nmc) :: Dp, Df, Dzero, temp
!double precision, dimension(nmeas,nmeas) :: Gamma, P, inv
!double precision, dimension(nmeas) :: meand
!
!integer :: nmeas,nmc,i
!
!meand =  sum(Dp,2)/dble(nmc)
!do i = 1,nmc
!   Dzero(:,i) = Dp(:,i) - meand
!end do
!call dgemm('N','T',nmeas,nmeas,nmc,1.0d0,Dzero,nmeas,Dzero,nmeas,0.0d0,Gamma,nmeas)
!
!meand =  sum(Df,2)/dble(nmc)
!do i = 1,nmc
!   Dzero(:,i) = Df(:,i) - meand
!end do
!call dgemm('N','T',nmeas,nmeas,nmc,1.0d0,Dzero,nmeas,Dzero,nmeas,0.0d0,P,nmeas)
!
!call invert(nmeas,(P+Gamma),inv)
!
!call dgemm('N','N',nmeas,nmc,nmeas,1.0d0,inv,nmeas,(Dp-Df),nmeas,0.0d0,temp,nmeas)
!call dgemm('T','N',nmc,nmc,nmeas,1.0d0,Dzero,nmeas,temp,nmeas,0.0d0,Cmat,nmc)
!
!do i = 1,nmc
!   Cmat(i,i) = Cmat(i,i) + 1.0d0
!end do
!
!end subroutine create_cmat_enkf


!!!*********************************************************
subroutine create_Af_pckf(nmeas,npceout,np,Upce,Af)

double precision, dimension(np,npceout+nmeas) :: Af
double precision, dimension(np*npceout) :: Upce

integer :: npceout,i,np,nmeas

do i = 1,npceout
   Af(:,i) = Upce((i-1)*np+1:i*np)
end do
Af(:,(npceout+1):(npceout+nmeas)) = 0.0d0

end subroutine create_Af_pckf

subroutine create_dp_pckf(nmeas,npceout,d,Dp,gamma)

double precision, dimension(nmeas,npceout+nmeas) :: Dp
double precision, dimension(nmeas) :: d

integer :: nmeas,npceout,i
double precision :: gamma

Dp(:,:) = 0.0d0
Dp(:,1) = d

do i = 1,nmeas
   Dp(i,npceout+i) = sqrt(gamma)
end do

end subroutine create_dp_pckf


!!!*********************************************************
subroutine create_df_pckf(nmeas,npceout,np,measnodes,Apckf,Df)

double precision, dimension(nmeas,npceout+nmeas) :: Df
double precision, dimension(np,npceout+nmeas) :: Apckf

integer, dimension(nmeas) :: measnodes
integer :: nmeas,npceout,np

Df = Apckf(measnodes,:)

end subroutine create_df_pckf

subroutine create_Dfg_pckf(pid,nmeas,nmeasg,np,npceout,Af,RDf)

double precision, dimension(nmeasg,npceout+nmeasg) :: RDf
double precision, dimension(np,npceout+nmeasg) :: Af

integer, dimension(nmeas) :: measnodes, measnodes_map
integer :: pid,nmeas,nmeasg,i,npceout,np

character(len=255) :: str1

call int2str(str1,pid+1,4)
str1 = '../../data/meshData/measnodes' // trim(str1) // '.dat'
open(unit=1,file=str1,status='old')
do i = 1,nmeas
   read(unit=1,fmt=*) measnodes(i)
end do
close(1)

call int2str(str1,pid+1,4)
str1 = '../../data/meshData/measnodes' // trim(str1) // '_map.dat'
open(unit=1,file=str1,status='old')
do i = 1,nmeas
   read(unit=1,fmt=*) measnodes_map(i)
end do
close(1)

RDf(:,:) = 0.0d0
do i = 1,nmeas
   RDf(measnodes_map(i),:) = Af(measnodes(i),:)
end do

end subroutine create_Dfg_pckf


!!!*********************************************************
subroutine create_Dpg_pckf(pid,nmeas,nmeasg,np,npceout,d,RDp,gamma)

double precision, dimension(nmeasg,npceout+nmeasg) :: RDp
double precision, dimension(nmeas) :: d

integer, dimension(nmeas) :: measnodes_map
integer :: pid,nmeas,nmeasg,i,npceout,np
double precision :: gamma

character(len=255) :: str1

call int2str(str1,pid+1,4)
str1 = '../../data/meshData/measnodes' // trim(str1) // '_map.dat'
open(unit=1,file=str1,status='old')
do i = 1,nmeas
   read(unit=1,fmt=*) measnodes_map(i)
end do
close(1)

RDp(:,:) = 0.0d0

do i = 1,nmeas
   RDp(measnodes_map(i),1) = d(i)
   RDp(measnodes_map(i),npceout+measnodes_map(i)) = sqrt(gamma)
end do

end subroutine create_Dpg_pckf

!!!*********************************************************
subroutine create_cmat_pckf(nmeas,npceout,Df,Dp,Cmat)

double precision, dimension(npceout+nmeas,npceout+nmeas) :: Cmat
double precision, dimension(nmeas,npceout+nmeas) :: Dp, Df, temp
double precision, dimension(nmeas,npceout+nmeas) :: Dzero
double precision, dimension(nmeas,nmeas) :: Gamma, P, inv
!!double precision, dimension(nmeas) :: meand

integer :: nmeas,npceout,i

Dzero(:,1) = 0.0d0
Dzero(:,2:npceout+nmeas) = Dp(:,2:npceout+nmeas)
call dgemm('N','T',nmeas,nmeas,npceout+nmeas,1.0d0,Dzero,nmeas,Dzero,nmeas,0.0d0,Gamma,nmeas)

Dzero(:,2:npceout+nmeas) = Df(:,2:npceout+nmeas)
call dgemm('N','T',nmeas,nmeas,npceout+nmeas,1.0d0,Dzero,nmeas,Dzero,nmeas,0.0d0,P,nmeas)

call invert(nmeas,(P+Gamma),inv)

call dgemm('N','N',nmeas,npceout+nmeas,nmeas,1.0d0,inv,nmeas,(Dp-Df),nmeas,0.0d0,temp,nmeas)
call dgemm('T','N',npceout+nmeas,npceout+nmeas,nmeas,1.0d0,Dzero,nmeas,temp,nmeas,0.0d0,Cmat,npceout+nmeas)

do i = 1,npceout+nmeas
   Cmat(i,i) = Cmat(i,i) + 1.0d0
end do

end subroutine create_cmat_pckf

!!!*********************************************************
!subroutine create_Dfg_enkf(pid,nmeas,nmeasg,nmc,np,Af,RDf)
!
!double precision, dimension(nmeasg,nmc) :: RDf
!double precision, dimension(np,nmc) :: Af
!
!integer, dimension(nmeas) :: measnodes, measnodes_map
!integer :: pid,nmeas,nmeasg,i,nmc,np
!
!character(len=255) :: str1
!
!call int2str(str1,pid+1,4)
!str1 = 'measnodes' // trim(str1) // '.dat'
!open(unit=1,file=str1,status='old')
!do i = 1,nmeas
!   read(unit=1,fmt=*) measnodes(i)
!end do
!close(1)
!
!call int2str(str1,pid+1,4)
!str1 = 'measnodes' // trim(str1) // '_map.dat'
!open(unit=1,file=str1,status='old')
!do i = 1,nmeas
!   read(unit=1,fmt=*) measnodes_map(i)
!end do
!close(1)
!
!RDf(:,:) = 0.0d0
!do i = 1,nmeas
!   RDf(measnodes_map(i),:) = Af(measnodes(i),:)
!end do
!
!end subroutine create_Dfg_enkf
!
!
!subroutine create_Dpg_enkf(pid,nmeas,nmeasg,nmc,np,d,RDp,gamma)
!
!double precision, dimension(nmeasg,nmc) :: RDp
!double precision, dimension(nmeas) :: d
!double precision, dimension(nmc) :: error
!
!integer, dimension(nmeas) :: measnodes_map
!integer :: pid,nmeas,nmeasg,i,nmc,np
!double precision :: gamma
!
!character(len=255) :: str1
!
!call int2str(str1,pid+1,4)
!str1 = 'measnodes' // trim(str1) // '_map.dat'
!open(unit=1,file=str1,status='old')
!do i = 1,nmeas
!   read(unit=1,fmt=*) measnodes_map(i)
!end do
!close(1)
!
!RDp(:,:) = 0.0d0
!
!do i = 1,nmeas
!   call norm_gen(error,nmc)
!!    call random_number(error)
!   RDp(measnodes_map(i),:) = d(i) + sqrt(gamma)*error
!end do
!
!end subroutine create_Dpg_enkf


!!!*********************************************************
subroutine create_Af_pckfddm(nmeasg,npceout,np,nb,Ui,Ub,Af)

double precision, dimension(np,npceout+nmeasg) :: Af
double precision, dimension((np-nb)*npceout) :: Ui
double precision, dimension(nb*npceout) :: Ub

integer :: npceout,i,np,nmeasg,nb

do i = 1,npceout
   Af(1:(np-nb),i) = Ui((i-1)*(np-nb)+1:i*(np-nb))
   Af((np-nb+1):np,i) = Ub((i-1)*nb+1:i*nb)
end do
Af(:,(npceout+1):(npceout+nmeasg)) = 0.0d0

end subroutine create_Af_pckfddm


!!!*********************************************************
subroutine create_Afb_pckfddm(nmeasg,npceout,nbg,Ub_g,Afb)

double precision, dimension(nbg,npceout+nmeasg) :: Afb
double precision, dimension(nbg*npceout) :: Ub_g

integer :: npceout,i,nmeasg,nbg

do i = 1,npceout
   Afb(:,i) = Ub_g((i-1)*nbg+1:i*nbg)
end do
Afb(:,(npceout+1):(npceout+nmeasg)) = 0.0d0

end subroutine create_Afb_pckfddm


!subroutine create_logw(d,Af,gamma,nmc,np,nmeas,measnodes,w)
!
!integer :: nmc,np,nmeas,i,j
!double precision :: gamma
!double precision, dimension(np,nmc) :: Af
!double precision, dimension(nmc) :: w
!double precision, dimension(nmeas) :: d
!integer, dimension(nmeas) :: measnodes
!
!w(:) = 0.0d0
!do i = 1,nmc
!   do j = 1,nmeas
!      w(i) = w(i) + (-1.0d0/(2.0d0*gamma))*(d(j)-Af(measnodes(j),i))**2
!   end do
!end do
!
!end subroutine create_logw


end module common
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%** END **%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
