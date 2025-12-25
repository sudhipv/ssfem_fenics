!!Solvers Module: by Ajit Desai, 7/2015, Mohammad Khalil, 1/2012 
!!contains: various solver for stochastic-DDM based PCGM with different Preconditioners 
!input:
!      DDM-Data  : Descretized blocks of matrices & vectors 
!      Mesh-Data : Global & local moes data
!      MPI-Data  : For MPI routines  
!      PCE-Data  : For stochastic formulation 
!output:
!      Ui, Ub, Urs, Ucs : Local solution vectors 
!      Ub_g, Uc         : Global solution vectors 
!
!Record of revisions:
!     Date        Programmer     Description of change
!  ==========     ==========     =====================
!  23/05/2014        MK          Stochastic-Lumped-PCGM solver(L)
!  12/06/2014        AD          Stochastic-Weighted-Lumped-PCGM Solver added (WL)
!  16/06/2014        AD          Stochastic-OneLevel-Naumann-Naumann-PCGM Solver added(NN)
!  23/07/2014        AD          Stochastic-TwoLevel-Naumann-Naumann-PCGM Solver added(NNC)
!  01/08/2014        AD          Both one/two Level Solver are optimized & cleaned
!  07/08/2015        AD          Coarse Problem: LUMPED-PCGM subroutine added
!  12/08/2015        AD          FETI-DP with Dirichlet-PCGM Solver is added
!!
!!----------------------------------------------------------------------------------------------
MODULE solvers

use common
use myCommon

implicit none

contains

!!!----------------------------------------------------------------------------------------
!!!The Direct Sub-structuring Method : Algorithm:1 from the report 
!!!----------------------------------------------------------------------------------------
SUBROUTINE solve_direct(pid,nb,np,nbg,npceout,Aii,Aig,Agi,Agg,Fi,Fg,Ui,Ub,Ub_g)

   use mpi

   integer :: pid,nb,np,npceout,nbg,ierr
   double precision :: Aii((np-nb)*npceout,(np-nb)*npceout),Aig((np-nb)*npceout,nb*npceout)
   double precision :: Agi(nb*npceout,(np-nb)*npceout),Agg(nb*npceout,nb*npceout)
   double precision :: Fi((np-nb)*npceout),Fg(nb*npceout),Ui((np-nb)*npceout),Ub(nb*npceout)
   double precision :: out1((np-nb)*npceout,nb*npceout),out2((np-nb)*npceout),RGi(nbg*npceout)
   double precision :: Si(nb*npceout,nb*npceout), Gi(nb*npceout), Rmat(nb*npceout,nbg*npceout)
   double precision :: RSiR(nbg*npceout,nbg*npceout),S(nbg*npceout,nbg*npceout)
   double precision :: Ub_g(nbg*npceout),calG(nbg*npceout)

!!!-------------------------------------------------------------------------------------------
   if (pid .eq. 0) print*, '-------------------------------------'
   if (pid .eq. 0) print*, 'Initializing Stochastic-Direct-DDM...'
   if (pid .eq. 0) print*, '-------------------------------------'
   
   if (pid .eq. 0) print*, 'solving once, matrix of dimension ',(np-nb)*npceout
   call solve((np-nb)*npceout,nb*npceout,Aii,out1,Aig)
   Si = Agg - matmul(Agi,out1)
   
   if (pid .eq. 0) print*, 'solving twice'
   call solve((np-nb)*npceout,1,Aii,out2,Fi)
   Gi = Fg - matmul(Agi,out2)

   call create_r(pid,nb,nbg,npceout,Rmat)
   RSiR = matmul(transpose(Rmat),matmul(Si,Rmat))
   RGi = matmul(transpose(Rmat),Gi)
   
   CALL MPI_REDUCE(RSiR,S,(nbg*npceout)*(nbg*npceout),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
   CALL MPI_REDUCE(RGi,calG,nbg*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

   if (pid .eq. 0) then
         print*, 'Su=g'
      call solve(nbg*npceout,1,S,Ub_g,calG)
      !!call solve_woodbury(nbg*npceout,nbg*npceout/2,nbg*npceout-(nbg*npceout/2),S,Ub_g,calG)
   end if

   CALL MPI_BCAST(Ub_g,nbg*npceout,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   Ub = matmul(Rmat,Ub_g)   
   
   if (pid .eq. 0) print*, 'solving thrice'
   call solve((np-nb)*npceout,1,Aii,Ui,Fi-matmul(Aig,Ub))

END SUBROUTINE solve_direct

!!!---------------------------------------------------------------------------------------------
!!!Lumped-Preconditioned PCGM Solver : Algorithm:2 from the report
!!!---------------------------------------------------------------------------------------------
SUBROUTINE solve_lpcgm(pid,nb,np,nbg,npceout,Aii,Aig,Agi,Agg,Fi,Fg,Ui,Ub,Ub_g,maxiter,tol,direct)

   use mpi

   integer :: pid,nb,np,npceout,nbg,maxiter,direct,ierr,i
   double precision :: Aii((np-nb)*npceout,(np-nb)*npceout),Aig((np-nb)*npceout,nb*npceout)
   double precision :: Agi(nb*npceout,(np-nb)*npceout),Agg(nb*npceout,nb*npceout)
   double precision :: Fi((np-nb)*npceout),Fg(nb*npceout),Ui((np-nb)*npceout)
   double precision :: Ub_g(nbg*npceout),Ub(nb*npceout),tol

   double precision :: rho_next, rho_curr, alpha, beta, err
   double precision :: out2((np-nb)*npceout), RGi(nbg*npceout),Pb(nb*npceout)
   double precision :: Gi(nb*npceout), Rmat(nb*npceout,nbg*npceout),Zb(nb*npceout)
   double precision :: rb_g(nbg*npceout),rb(nb*npceout),Minv(nb*npceout,nb*npceout)
   double precision :: RZb(nbg*npceout),Zb_g(nbg*npceout),Pb_g(nbg*npceout)
   double precision :: Qb(nb*npceout),RQb(nbg*npceout),Qb_g(nbg*npceout),out3((np-nb)*npceout)
   double precision :: MAii((np-nb)*npceout,(np-nb)*npceout), temp((np-nb),(np-nb))

   Ub_g(:) = 0.0d0
   Minv = Agg

!!!-------------------------------------------------------------------------------------------
!! Cholesky Decomposition of As to solve x = inv(As)*b : direct = 1
   if (direct .eq. 1) then
      if (pid .eq. 0) print*, 'Obtaining Cholesky Decomposition of Aii (',(np-nb)*npceout,')'
      call dpotrf('L',(np-nb)*npceout,Aii,(np-nb)*npceout,ierr)
      if (ierr .ne. 0) print*, 'Error in factorizing Aii! Processor # ', pid
   else
      MAii(:,:) = 0.0d0
      call invert((np-nb),Aii(1:(np-nb),1:(np-nb)),temp)
      do i = 1,npceout
         MAii((i-1)*(np-nb)+1:i*(np-nb),(i-1)*(np-nb)+1:i*(np-nb)) = temp
      end do
   end if
   
   if (pid .eq. 0) print*, 'Obtaining Cholesky Decomposition of Minv(',nb*npceout,')'
   call dpotrf('L',nb*npceout,Minv,nb*npceout,ierr)
   if (ierr .ne. 0) print*, 'Error in factorizing Minv! Processor # ', pid

   if (pid .eq. 0) print*, '--------------------------------------'
   if (pid .eq. 0) print*, 'Initializing Stochastic-Lumped-PCGM...'
   if (pid .eq. 0) print*, '--------------------------------------'
   
   if (direct .eq. 1) then
      out2 = Fi
      call dpotrs('L',(np-nb)*npceout,1,Aii,(np-nb)*npceout,out2,(np-nb)*npceout,ierr)
   else
      call pcgm(pid,(np-nb)*npceout,Aii,out2,Fi,10000,1.0d-10,1,MAii)
   end if
   call multiply(Agi,out2,Gi,nb*npceout,(np-nb)*npceout,1)
   Gi = Fg - Gi
   
!!!-------------------------------------------------------------------------------------------
!! Parallel-Preconditioning effect
   call getubg(pid,nb,nbg,npceout,RGi,Gi)  
   call MPI_ALLREDUCE(RGi,rb_g,nbg*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr) 
   !!call create_r(pid,nb,nbg,npceout,Rmat)    !!for fair time comparison uncoment this line
   !!RGi = matmul(transpose(Rmat),Gi)
   call getub(pid,nb,nbg,npceout,rb_g,rb) 
   !!rb = matmul(Rmat,rb_g)                   !!for fair time comparison uncoment this line
   Zb = rb
   call dpotrs('L',nb*npceout,1,Minv,nb*npceout,Zb,nb*npceout,ierr)
   call getubg(pid,nb,nbg,npceout,RZb,Zb) 
   !!RZb = matmul(transpose(Rmat),Zb)          !!for fair time comparison uncoment this line
   CALL MPI_REDUCE(RZb,Zb_g,nbg*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

!!!---------------------------------------------------------------------------------------------
   if (pid .eq. 0) Pb_g = Zb_g
   if (pid .eq. 0) rho_next = dot_product(rb_g,Zb_g)

!!!---------------------------------------------------------------------------------------------
   out2(:) = 0.0d0
   do i = 1,maxiter
      CALL MPI_BCAST(Pb_g,nbg*npceout,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call getub(pid,nb,nbg,npceout,Pb_g,Pb)
      call multiply(Aig,Pb,out2,(np-nb)*npceout,nb*npceout,1)
      if (direct .eq. 1) then
         Ui = out2
         call dpotrs('L',(np-nb)*npceout,1,Aii,(np-nb)*npceout,Ui,(np-nb)*npceout,ierr)
      else
         call pcgm(pid,(np-nb)*npceout,Aii,Ui,out2,10000,1.0d-10,1,MAii)
      end if

      call multiply(Agi,Ui,out2,nb*npceout,(np-nb)*npceout,1)
      call multiply(Agg,Pb,Qb,nb*npceout,nb*npceout,1)
      Qb = Qb - out2

      call getubg(pid,nb,nbg,npceout,RQb,Qb)
      CALL MPI_REDUCE(RQb,Qb_g,nbg*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if (pid .eq. 0) then
         rho_curr = rho_next
         alpha = rho_curr/dot_product(Qb_g,Pb_g)
         err = alpha*alpha*dot_product(Pb_g,Pb_g)/dot_product(Ub_g,Ub_g)
         Ub_g = Ub_g + alpha*Pb_g
         rb_g = rb_g - alpha*Qb_g

         print*, 'iteration # ', i, ', relative error of ', err
         !!print*, 'and sum of Ub of ', sum(Ub_g)
      end if
     
      CALL MPI_BCAST(err,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      if (err .lt. tol) exit
      
!!!-------------------------------------------------------------------------------------------
!! Parallel-Preconditioning effect 
      CALL MPI_BCAST(rb_g,nbg*npceout,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call getub(pid,nb,nbg,npceout,rb_g,rb) 
    
      Zb = rb
      call dpotrs('L',nb*npceout,1,Minv,nb*npceout,Zb,nb*npceout,ierr)
      call getubg(pid,nb,nbg,npceout,RZb,Zb) 
    
      CALL MPI_REDUCE(RZb,Zb_g,nbg*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!!!---------------------------------------------------------------------------------------------

      if (pid .eq. 0) then
         rho_next = dot_product(rb_g,Zb_g)
         beta = rho_next/rho_curr
         Pb_g = Zb_g + beta*Pb_g
      end if
   end do

   CALL MPI_BCAST(Ub_g,nbg*npceout,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   call getub(pid,nb,nbg,npceout,Ub_g,Ub)
   
   call multiply(Aig,Ub,out2,(np-nb)*npceout,nb*npceout,1)
   if (direct .eq. 1) then
      Ui = Fi-out2
      call dpotrs('L',(np-nb)*npceout,1,Aii,(np-nb)*npceout,Ui,(np-nb)*npceout,ierr)
   else
      call pcgm(pid,(np-nb)*npceout,Aii,Ui,Fi-out2,10000,1.0d-10,1,MAii)
   end if

END SUBROUTINE solve_lpcgm


!!!---------------------------------------------------------------------------------------------
!!! Weighted-Lumped-Preconditioned PCGM Solver : Algorithm:3 from the report 
!!!---------------------------------------------------------------------------------------------
SUBROUTINE solve_wlpcgm(pid,nb,np,nbg,npceout,Aii,Aig,Agi,Agg,Fi,Fg,Ui,Ub,Ub_g,maxiter,tol,direct)

    use mpi

    integer :: pid,nb,np,npceout,nbg,maxiter,direct,ierr,i
    double precision :: Aii((np-nb)*npceout,(np-nb)*npceout),Aig((np-nb)*npceout,nb*npceout)
    double precision :: Agi(nb*npceout,(np-nb)*npceout),Agg(nb*npceout,nb*npceout)
    double precision :: Fi((np-nb)*npceout),Fg(nb*npceout),Ui((np-nb)*npceout),Ub(nb*npceout)
    double precision :: Ub_g(nbg*npceout), tol

    double precision :: rho_next, rho_curr, alpha, beta, err
    double precision :: Gi(nb*npceout), Rmat(nb*npceout,nbg*npceout)
    double precision :: out2((np-nb)*npceout),RGi(nbg*npceout),Zb(nb*npceout)
    double precision :: rb_g(nbg*npceout),rb(nb*npceout),Minv(nb*npceout,nb*npceout)
    double precision :: RZb(nbg*npceout),Zb_g(nbg*npceout),Pb_g(nbg*npceout)
    double precision :: Qb(nb*npceout),RQb(nbg*npceout),Qb_g(nbg*npceout),Pb(nb*npceout)
    double precision :: MAii((np-nb)*npceout,(np-nb)*npceout), temp((np-nb),(np-nb))
    double precision :: Dmat(nb*npceout,nb*npceout),yb(nb*npceout),Utau(nb*npceout)

    Ub_g(:) = 0.0d0
    Minv = Agg

!!!-------------------------------------------------------------------------------------------
!! Cholesky Decomposition of As to solve x = inv(As)*b : direct = 1
    if (direct .eq. 1) then
       if (pid .eq. 0) print*, 'Obtaining Cholesky Decomposition of Aii (',(np-nb)*npceout,')'
       call dpotrf('L',(np-nb)*npceout,Aii,(np-nb)*npceout,ierr)
       if (ierr .ne. 0) print*, 'Error in factorizing Aii! Processor # ', pid
    else
       MAii(:,:) = 0.0d0
       call invert((np-nb),Aii(1:(np-nb),1:(np-nb)),temp)
       do i = 1,npceout
          MAii((i-1)*(np-nb)+1:i*(np-nb),(i-1)*(np-nb)+1:i*(np-nb)) = temp
       end do
    end if
    
    if (pid .eq. 0) print*, 'Obtaining Cholesky Decomposition of Minv(',nb*npceout,')'
    call dpotrf('L',nb*npceout,Minv,nb*npceout,ierr)
    if (ierr .ne. 0) print*, 'Error in factorizing Minv! Processor # ', pid

    if (pid .eq. 0) print*, '-----------------------------------------------'
    if (pid .eq. 0) print*, 'Initializing Stochastic-Weighted-Lumped-PCGM...'
    if (pid .eq. 0) print*, '-----------------------------------------------'
    
    if (direct .eq. 1) then
        out2 = Fi
        call dpotrs('L',(np-nb)*npceout,1,Aii,(np-nb)*npceout,out2,(np-nb)*npceout,ierr)
    else
        call pcgm(pid,(np-nb)*npceout,Aii,out2,Fi,10000,1.0d-10,1,MAii)
    end if
    call multiply(Agi,out2,Gi,nb*npceout,(np-nb)*npceout,1)
    Gi = Fg - Gi

!!!-------------------------------------------------------------------------------------------
!! Parallel-Preconditioning effect
    call getubg(pid,nb,nbg,npceout,RGi,Gi)  
    call MPI_ALLREDUCE(RGi,rb_g,nbg*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    !!CALL MPI_REDUCE(RGi,rb_g,nbg*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    !!CALL MPI_BCAST(rb_g,nbg*npceout,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    
    call getub(pid,nb,nbg,npceout,rb_g,rb) 
    call getBDM(pid,nb,npceout,Dmat)
    
    yb = matmul(Dmat,rb)
    call dpotrs('L',nb*npceout,1,Minv,nb*npceout,yb,nb*npceout,ierr)
    Zb = matmul(Dmat,yb)
    
    call getubg(pid,nb,nbg,npceout,RZb,Zb) 
    call MPI_REDUCE(RZb,Zb_g,nbg*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!!!---------------------------------------------------------------------------------------------

    if (pid .eq. 0) Pb_g = Zb_g
    if (pid .eq. 0) rho_next = dot_product(rb_g,Zb_g)

!!!---------------------------------------------------------------------------------------------
    out2(:) = 0.0d0
    do i = 1,maxiter
    
!! Parallel Matrix-Vector product
        CALL MPI_BCAST(Pb_g,nbg*npceout,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call getub(pid,nb,nbg,npceout,Pb_g,Pb)
        
        call multiply(Aig,Pb,out2,(np-nb)*npceout,nb*npceout,1)
        if (direct .eq. 1) then
           Ui = out2
           call dpotrs('L',(np-nb)*npceout,1,Aii,(np-nb)*npceout,Ui,(np-nb)*npceout,ierr)
        else
           call pcgm(pid,(np-nb)*npceout,Aii,Ui,out2,10000,1.0d-10,1,MAii)
        end if

        call multiply(Agi,Ui,out2,nb*npceout,(np-nb)*npceout,1)
        call multiply(Agg,Pb,Qb,nb*npceout,nb*npceout,1)
        Qb = Qb - out2

        call getubg(pid,nb,nbg,npceout,RQb,Qb)
        CALL MPI_REDUCE(RQb,Qb_g,nbg*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        
        if (pid .eq. 0) then
            rho_curr = rho_next
            alpha = rho_curr/dot_product(Qb_g,Pb_g)
            err = alpha*alpha*dot_product(Pb_g,Pb_g)/dot_product(Ub_g,Ub_g)
            Ub_g = Ub_g + alpha*Pb_g
            rb_g = rb_g - alpha*Qb_g
            print*, 'iteration # ', i, ', relative error of ', err
            !!print*,  ' sum of Ub of ', sum(Ub_g)
        end if
        CALL MPI_BCAST(err,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        if (err .lt. tol) exit
        
!!!-------------------------------------------------------------------------------------------
!! Parallel-Preconditioning effect
        call MPI_BCAST(rb_g,nbg*npceout,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call getub(pid,nb,nbg,npceout,rb_g,rb) 

        yb = matmul(Dmat,rb)
        call dpotrs('L',nb*npceout,1,Minv,nb*npceout,yb,nb*npceout,ierr)
        Zb = matmul(Dmat,yb)

        call getubg(pid,nb,nbg,npceout,RZb,Zb) 
        call MPI_REDUCE(RZb,Zb_g,nbg*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!!!---------------------------------------------------------------------------------------------

        if (pid .eq. 0) then
        rho_next = dot_product(rb_g,Zb_g)
        beta = rho_next/rho_curr
        Pb_g = Zb_g + beta*Pb_g
        end if
    end do

!!!-------------------------------------------------------------------------------------------
    CALL MPI_BCAST(Ub_g,nbg*npceout,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call getub(pid,nb,nbg,npceout,Ub_g,Ub)
    
    call multiply(Aig,Ub,out2,(np-nb)*npceout,nb*npceout,1)
    if (direct .eq. 1) then
        Ui = Fi -out2
        call dpotrs('L',(np-nb)*npceout,1,Aii,(np-nb)*npceout,Ui,(np-nb)*npceout,ierr)
    else
        call pcgm(pid,(np-nb)*npceout,Aii,Ui,Fi-out2,10000,1.0d-10,1,MAii)
    end if

END SUBROUTINE solve_wlpcgm

!!!----------------------------------------------------------------------------------------------
!!! One-Level Neumann-Neumann Preconditioned PCGM solver : Algorithm:5 from the report
!!!----------------------------------------------------------------------------------------------
SUBROUTINE solve_nnpcgm(pid,nb,np,nbg,npceout,Aii,Aig,Agi,Agg,Fi,Fg,Ui,Ub,Ub_g,maxiter,tol,direct)

    use mpi

    integer :: pid,nb,np,nbg,ierr,i,npceout,maxiter,direct,j,k
    double precision :: Aii((np-nb)*npceout,(np-nb)*npceout),Aig((np-nb)*npceout,nb*npceout)
    double precision :: Agi(nb*npceout,(np-nb)*npceout),Agg(nb*npceout,nb*npceout)
    double precision :: Fi((np-nb)*npceout),Fg(nb*npceout),Ui((np-nb)*npceout),Ub(nb*npceout)
    double precision :: Ub_g(nbg*npceout), tol

    double precision :: rho_next, rho_curr, alpha, beta, err
    double precision :: out2((np-nb)*npceout),RGi(nbg*npceout),Gi(nb*npceout)
    double precision :: rb_g(nbg*npceout),rb(nb*npceout),Minv(nb*npceout,nb*npceout),Zb(nb*npceout)
    double precision :: RZb(nbg*npceout),Zb_g(nbg*npceout),Pb_g(nbg*npceout),Pb(nb*npceout)
    double precision :: Qb(nb*npceout),RQb(nbg*npceout),Qb_g(nbg*npceout)
    double precision :: MAii((np-nb)*npceout,(np-nb)*npceout), temp((np-nb),(np-nb))
    double precision :: Utau(nb*npceout), yb(nb*npceout), Dmat(nb*npceout,nb*npceout)
    double precision :: Amat(np*npceout,np*npceout), ybc(np*npceout), out12(np*npceout)

!!-----------------------------------------------------------------------------------------------
!!!!! Uncomment this section for Execution Time Calculations !!!!!
!!-----------------------------------------------------------------------------------------------
    !integer :: c1,c2,cr
    !double precision :: time1
    !double precision :: time2
    !call system_clock(count_rate=cr)
    !call cpu_time(time1)
    !call system_clock(c1)

!!-----------------------------------------------------------------------------------------------
!! To create Block Matrix Amat, used for Dirichlet & Naumann Solve
    Amat(:,:) = 0.0d0
    Amat(1:(np-nb)*npceout,1:(np-nb)*npceout) = Aii
    Amat(1:(np-nb)*npceout,((np-nb)*npceout+1):np*npceout) = Aig
    Amat(((np-nb)*npceout+1):np*npceout,1:(np-nb)*npceout) = Agi
    Amat(((np-nb)*npceout+1):np*npceout,((np-nb)*npceout+1):np*npceout) = Agg

!!-----------------------------------------------------------------------------------------------
    Ub_g(:) = 0.0d0
    Minv = Agg

    if (direct .eq. 1) then
      if (pid .eq. 0) print*, 'Obtaining Cholesky Decomposition of Aii (',(np-nb)*npceout,')'
      call dpotrf('L',(np-nb)*npceout,Aii,(np-nb)*npceout,ierr)
      if (ierr .ne. 0) print*, 'Error in factorizing Aii! Processor # ', pid
    else
      MAii(:,:) = 0.0d0
      call invert((np-nb),Aii(1:(np-nb),1:(np-nb)),temp)
      do i = 1,npceout
         MAii((i-1)*(np-nb)+1:i*(np-nb),(i-1)*(np-nb)+1:i*(np-nb)) = temp
      end do
    end if
    
    if (pid .eq. 0) print*, 'Obtaining Cholesky Decomposition of Minv(',nb*npceout,')'
    call dpotrf('L',nb*npceout,Minv,nb*npceout,ierr)
    if (ierr .ne. 0) print*, 'Error in factorizing Minv! Processor # ', pid

    if (pid .eq. 0) print*, 'Obtaining Cholesky Decomposition of Amat(',np*npceout,')'
    call dpotrf('L',np*npceout,Amat,np*npceout,ierr)
    if (ierr .ne. 0) print*, 'Error in factorizing Amat! Processor # ', pid

    !if (pid .eq. 0) then
    !call cpu_time(time2)
    !call system_clock(c2)
    !print*, 'Time taken for Cholesky Decomposition Amat', dble(c2-c1)/dble(cr),&
    !'seconds (', time2-time1, 'cpu time)'
    !end if

!!!----------------------------------------------------------------------------------------
    if (pid .eq. 0) print*, '---------------------------------------------------------'
    if (pid .eq. 0) print*, 'Initializing Stochastic-One-Level-Neumann-Neumann-PCGM...'
    if (pid .eq. 0) print*, '---------------------------------------------------------'
    
    if (direct .eq. 1) then
      out2 = Fi
      call dpotrs('L',(np-nb)*npceout,1,Aii,(np-nb)*npceout,out2,(np-nb)*npceout,ierr)
    else
      call pcgm(pid,(np-nb)*npceout,Aii,out2,Fi,10000,1.0d-10,1,MAii)
    end if
    call multiply(Agi,out2,Gi,nb*npceout,(np-nb)*npceout,1)
    Gi = Fg - Gi
    call getubg(pid,nb,nbg,npceout,RGi,Gi)
    CALL MPI_REDUCE(RGi,rb_g,nbg*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    
!!!----------------------------------------------------------------------------------------
!!! OneLevel Newmann-Newmann Precondtioner
    CALL MPI_BCAST(rb_g,nbg*npceout,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)      
    call getub(pid,nb,nbg,npceout,rb_g,rb)
    call getBDM(pid, nb, npceout, Dmat)
    yb = matmul(Dmat,rb)                  !Zb = rb

    ybc(1:(np-nb)*npceout) = 0.0d0
    ybc(((np-nb)*npceout+1):np*npceout) = yb
    out12 = ybc
    call dpotrs('L',np*npceout,1,Amat,np*npceout,out12,np*npceout,ierr)
    Utau = out12(((np-nb)*npceout+1):np*npceout)
    Zb = matmul(Dmat,Utau)

    call getubg(pid,nb,nbg,npceout,RZb,Zb)
    CALL MPI_REDUCE(RZb,Zb_g,nbg*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!!----------------------------------------------------------------------------------------

    if (pid .eq. 0) Pb_g = Zb_g
    if (pid .eq. 0) rho_next = dot_product(rb_g,Zb_g)
    
!!-----------------------------------------------------------------------------------------------
    out2(:) = 0.0d0
    do i = 1,maxiter
      CALL MPI_BCAST(Pb_g,nbg*npceout,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call getub(pid,nb,nbg,npceout,Pb_g,Pb)
      call multiply(Aig,Pb,out2,(np-nb)*npceout,nb*npceout,1)
      if (direct .eq. 1) then
         Ui = out2
         call dpotrs('L',(np-nb)*npceout,1,Aii,(np-nb)*npceout,Ui,(np-nb)*npceout,ierr)
      else
         call pcgm(pid,(np-nb)*npceout,Aii,Ui,out2,10000,1.0d-10,1,MAii)
      end if
    
      call multiply(Agi,Ui,out2,nb*npceout,(np-nb)*npceout,1)
      call multiply(Agg,Pb,Qb,nb*npceout,nb*npceout,1)
      Qb = Qb - out2
    
      call getubg(pid,nb,nbg,npceout,RQb,Qb)
      CALL MPI_REDUCE(RQb,Qb_g,nbg*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if (pid .eq. 0) then
         rho_curr = rho_next
         alpha = rho_curr/dot_product(Qb_g,Pb_g)
         err = alpha*alpha*dot_product(Pb_g,Pb_g)/dot_product(Ub_g,Ub_g)
         Ub_g = Ub_g + alpha*Pb_g
         rb_g = rb_g - alpha*Qb_g
         print*, 'iteration # ', i, ', relative error of ', err
         !!print*, 'and sum of Ub of ', sum(Ub_g)
      end if
      CALL MPI_BCAST(err,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      if (err .lt. tol) exit
      
!!!----------------------------------------------------------------------------------------
!! OneLevel Newmann-Newmann Precondtioner
      CALL MPI_BCAST(rb_g,nbg*npceout,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)      
      call getub(pid,nb,nbg,npceout,rb_g,rb)
      yb = matmul(Dmat,rb)                  !Zb = rb

      !ybc(1:(np-nb)*npceout) = 0.0d0
      ybc(((np-nb)*npceout+1):np*npceout) = yb
      out12 = ybc
      call dpotrs('L',np*npceout,1,Amat,np*npceout,out12,np*npceout,ierr)
      Utau = out12(((np-nb)*npceout+1):np*npceout)
      Zb = matmul(Dmat,Utau)

      call getubg(pid,nb,nbg,npceout,RZb,Zb)
      CALL MPI_REDUCE(RZb,Zb_g,nbg*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

!!!----------------------------------------------------------------------------------------
      if (pid .eq. 0) then
         rho_next = dot_product(rb_g,Zb_g)
         beta = rho_next/rho_curr
         Pb_g = Zb_g + beta*Pb_g
      end if
    end do
    
    CALL MPI_BCAST(Ub_g,nbg*npceout,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    
    call getub(pid,nb,nbg,npceout,Ub_g,Ub)
    call multiply(Aig,Ub,out2,(np-nb)*npceout,nb*npceout,1)
    
    if (direct .eq. 1) then
      Ui = Fi-out2
      call dpotrs('L',(np-nb)*npceout,1,Aii,(np-nb)*npceout,Ui,(np-nb)*npceout,ierr)
    else
      call pcgm(pid,(np-nb)*npceout,Aii,Ui,Fi-out2,10000,1.0d-10,1,MAii)
    end if

END SUBROUTINE solve_nnpcgm


!!!----------------------------------------------------------------------------------------------
!!!!!!!!!!!!!! Two-Level Neumann-Neumann Preconditioned PCGM Solver !!!!!!!!!!!!!!
!!!----------------------------------------------------------------------------------------------
SUBROUTINE solve_nncpcgm(pid,nb,np,nbg,npceout,nri,nci,nbgc,Aii,Aig,Agi,Agg,Air,Aic,&
                         Ari,Aci,Arr,Acc,Arc,Acr,Fi,Fg,Ui,Ub,Ub_g,maxiter,tol)
    use mpi

    integer :: pid,nb,np,nbg,ierr,i,j,npceout,nci,nri,nbgc, maxiter !! integer :: direct = 1

    double precision :: rho_next,rho_curr, alpha, beta, err,residual,tol
    double precision :: Aii((np-nb)*npceout,(np-nb)*npceout),Aig((np-nb)*npceout,nb*npceout)
    double precision :: Agi(nb*npceout,(np-nb)*npceout),Agg(nb*npceout,nb*npceout)
    double precision :: Fi((np-nb)*npceout),Fg(nb*npceout),out2((np-nb)*npceout)
    double precision :: Ub_g(nbg*npceout),Ui((np-nb)*npceout),Ub(nb*npceout)
    double precision :: Gi(nb*npceout),yb(nb*npceout),RGi(nbg*npceout)

    double precision :: rb_g(nbg*npceout),rb(nb*npceout),Minv(nb*npceout,nb*npceout)
    double precision :: RZb(nbg*npceout),Zb_g(nbg*npceout),Pb_g(nbg*npceout),Pb(nb*npceout)
    double precision :: Qb(nb*npceout),RQb(nbg*npceout),Qb_g(nbg*npceout)
    double precision :: MAii((np-nb)*npceout,(np-nb)*npceout),temp((np-nb),(np-nb))

    double precision :: Air((np-nb)*npceout,(nb-nci)*npceout),Aic((np-nb)*npceout,nci*npceout)
    double precision :: Arr((nb-nci)*npceout,(nb-nci)*npceout),Acc(nci*npceout,nci*npceout)
    double precision :: Ari((nb-nci)*npceout,(np-nb)*npceout),Aci(nci*npceout,(np-nb)*npceout)
    double precision :: Arc((nb-nci)*npceout,nci*npceout),Acr(nci*npceout,(nb-nci)*npceout)

    double precision :: Dmat(nb*npceout,nb*npceout), RcSmat(nci*npceout,nb*npceout)
    double precision :: RrSmat(nri*npceout,nb*npceout), BcSmat(nci*npceout,nbgc*npceout)
    double precision :: BcSmatT(nbgc*npceout,nci*npceout),FrS(nri*npceout), FcS(nci*npceout)
    double precision :: RcSmatT(nb*npceout,nci*npceout), RrSmatT(nb*npceout,nri*npceout)

    double precision :: v1(nri*npceout), DcBc(nbgc*npceout), dcS(nci*npceout)
    double precision :: Dc(nbgc*npceout), Zc(nbgc*npceout), ZcS(nci*npceout)
    double precision :: v2(nri*npceout), ZrS(nri*npceout), Zr2(nb*npceout)
    double precision :: Zc2(nb*npceout), Zs(nb*npceout), ZsD(nb*npceout)
    double precision :: Amat((np-nci)*npceout,(np-nci)*npceout),ybc((np-nci)*npceout)

    double precision :: u11((np-nb)*npceout),u12((np-nb)*npceout),u13(nri*npceout)
    double precision :: u14(nri*npceout), u23(nci*npceout), u24(nci*npceout)
    double precision :: Minc(nci*npceout,nci*npceout),out12((np-nci)*npceout)
    double precision :: tolc = 1.0d-4

!!-----------------------------------------------------------------------------------------------
!!  Uncomment this section for Execution Time Calculations !!!!!
!!-----------------------------------------------------------------------------------------------
    !integer :: c1,c2,cr
    !double precision :: time1
    !double precision :: time2
    !call system_clock(count_rate=cr)
    
!!-----------------------------------------------------------------------------------------------
!! To create Block Matrix Amat, used for Dirichlet & Naumann Solve
    Amat(1:(np-nb)*npceout,1:(np-nb)*npceout) = Aii
    Amat(1:(np-nb)*npceout,((np-nb)*npceout+1):(np-nci)*npceout) = Air
    Amat(((np-nb)*npceout+1):(np-nci)*npceout,1:(np-nb)*npceout) = Ari
    Amat(((np-nb)*npceout+1):(np-nci)*npceout,((np-nb)*npceout+1):(np-nci)*npceout) = Arr

!!-----------------------------------------------------------------------------------------------
!! Initiate Array
    Ub_g(:) = 0.0d0   !Minv = Agg
    Minc = Acc
    
!!-----------------------------------------------------------------------------------------------
!! Cholesky Decomposition of As to solve x = inv(As)*b
    if (pid .eq. 0) print*, 'Obtaining Cholesky Decomposition of Aii (',(np-nb)*npceout,')'
    call dpotrf('L',(np-nb)*npceout,Aii,(np-nb)*npceout,ierr)
    if (ierr .ne. 0) print*, 'Error in factorizing Aii! Processor # ', pid

    if (pid .eq. 0) print*, 'Obtaining Cholesky Decomposition of Amat(',(np-nci)*npceout,')'
    call dpotrf('L',(np-nci)*npceout,Amat,(np-nci)*npceout,ierr)
    if (ierr .ne. 0) print*, 'Error in factorizing Amat! Processor # ', pid

    if (pid .eq. 0) print*, 'Obtaining Cholesky Decomposition of Acc(',nci*npceout,')'
    call dpotrf('L',nci*npceout,Minc,nci*npceout,ierr)
    if (ierr .ne. 0) print*, 'Error in factorizing Acc! Processor # ', pid

!!----------------------------------------------------------------------------------------------
    if (pid .eq. 0) print*, '---------------------------------------------------------'
    if (pid .eq. 0) print*, 'Initializing Stochastic-Two-Level-Neumann-Neumann-PCGM...'
    if (pid .eq. 0) print*, '---------------------------------------------------------'
    
    !! STEP 3 : Solve
    
    out2 = Fi
    call dpotrs('L',(np-nb)*npceout,1,Aii,(np-nb)*npceout,out2,(np-nb)*npceout,ierr)
    
    !! STEP 4 : Compute
    call multiply(Agi,out2,Gi,nb*npceout,(np-nb)*npceout,1)
    Gi = Fg - Gi
    
    !! STEP 5 : Gather 
    call getubg(pid,nb,nbg,npceout,RGi,Gi)
    CALL MPI_ALLREDUCE(RGi,rb_g,nbg*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    !!CALL MPI_REDUCE(RGi,rb_g,nbg*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

!!-----------------------------------------------------------------------------------------------
!!  Uncomment this section for Residual Calculation !!!!!
!!-----------------------------------------------------------------------------------------------
    !if (pid .eq. 0) then
    !residual = SQRT(dot_product(rb_g,rb_g))
    !open (unit=7,file="residual.dat",action="write",status="replace")
    !print*,residual
    !write (7,*), residual
    !close (7)
    !end if

!!-----------------------------------------------------------------------------------------------
!!! Two-Level Newmann-Newmann Precondtioner: Using : Algorithm 7 from the report
    !! CALL MPI_BCAST(rb_g,nbg*npceout,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    !! SUB-STEP 3 : Scatter    
    call getub(pid,nb,nbg,npceout,rb_g,rb)                  ! -> rb
    
    !! SUB-STEP 4 : Weighting
    call getBDM(pid, nb, npceout, Dmat)                     ! -> Dmat
    call multiply(Dmat,rb,yb,nb*npceout,nb*npceout,1)       ! yb = matmul(Dmat,rb)
    
    !! SUB-STEP 5 & 6 : Compute 
    call getRrS(pid, nb, nri, npceout, RrSmat)
    call multiply(RrSmat,yb,FrS,nri*npceout,nb*npceout,1)   !FrS = matmul(RrSmat,yb) !error here
    call getRcS(pid, nb, nci, npceout, RcSmat)
    call multiply(RcSmat,yb,FcS,nci*npceout,nb*npceout,1)   !FcS = matmul(RcSmat,yb)

!!-----------------------------------------------------------------------------------------------
    !! SUB-STEP 7 : Solve
    ybc(:) = 0.0d0
    ybc(((np-nb)*npceout+1):(np-nci)*npceout) = FrS
    out12 = ybc
    call dpotrs('L',(np-nci)*npceout,1,Amat,(np-nci)*npceout,out12,(np-nci)*npceout,ierr)
    v1 = out12(((np-nb)*npceout+1):(np-nci)*npceout)
    
    !! SUB-STEP 8 : Compute 
    !!!------------------- Section3/Parallel MatrixVector Product ---------------------!!!
    !!! Scr*v1 = [Acr - (Aci *inv(Aii)*Air)]*v1
    call multiply(Air,v1,u11,(np-nb)*npceout,nri*npceout,1)
    u12 = u11
    call dpotrs('L',(np-nb)*npceout,1,Aii,(np-nb)*npceout,u12,(np-nb)*npceout,ierr)
    call multiply(Aci,u12,u23,nci*npceout,(np-nb)*npceout,1)
    call multiply(Acr,v1,u24,nci*npceout,nri*npceout,1)
    dcS = u24 - u23
    
    !! SUB-STEP 9 : Update 
    dcS = FcS - dcS
    
    !! SUB-STEP 10 : Gather     
    call getBcS(pid, nci, nbgc, npceout, BcSmat)
    BcSmatT = TRANSPOSE(BcSmat)
    call multiply(BcSmatT,dcS,DcBc,nbgc*npceout,nci*npceout,1)
    CALL MPI_ALLREDUCE(DcBc,Dc,nbgc*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    !!CALL MPI_REDUCE(DcBc,Dc,nbgc*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

!!-----------------------------------------------------------------------------------------------
!! SUB-STEP 12 : Coarse Problem: LUMPED-PCGM : Algorithm 9 from the report
!!-----------------------------------------------------------------------------------------------
    call solve_coarsepcgm(pid,nb,np,nci,nri,nbg,nbgc,npceout,BcSmat,BcSmatT, & 
                          Aii,Acc,Aic,Air,Aci,Ari,Arc,Acr,Amat,Minc,Dc,Zc,tolc)

!!-----------------------------------------------------------------------------------------------
    !! SUB-STEP 14 : Scatter
    CALL MPI_BCAST(Zc,nbgc*npceout,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call multiply(BcSmat,Zc,ZcS,nci*npceout,nbgc*npceout,1)           !! ZcS = matmul(BcSmat,Zc)
    
    !! SUB-STEP 15 : Compute 
    !!!-------------------Section1/Parallel MatrixVector Product ---------------------!!!
    !!! Src*ZcS = [Arc - (Ari *inv(Aii)*Air)]*ZcS
    call multiply(Aic,ZcS,u11,(np-nb)*npceout,nci*npceout,1)
    u12 = u11
    call dpotrs('L',(np-nb)*npceout,1,Aii,(np-nb)*npceout,u12,(np-nb)*npceout,ierr)
    call multiply(Ari,u12,u13,nri*npceout,(np-nb)*npceout,1)
    call multiply(Arc,ZcS,u14,nri*npceout,nci*npceout,1)
    v2 = u14 - u13
    
    !! SUB-STEP 16 : Update    
    v2 = FrS - v2

!!-----------------------------------------------------------------------------------------------
!! Solving Naumann-Naumann Problem: solve(nri*npceout,1,Srr,ZrS,v2)
!!-----------------------------------------------------------------------------------------------
    !! SUB-STEP 17 : Solve
    ybc(((np-nb)*npceout+1):(np-nci)*npceout) = v2
    out12 = ybc
    call dpotrs('L',(np-nci)*npceout,1,Amat,(np-nci)*npceout,out12,(np-nci)*npceout,ierr)
    ZrS = out12(((np-nb)*npceout+1):(np-nci)*npceout)

!!-----------------------------------------------------------------------------------------------
    !! SUB-STEP 18 & 19 : Compute
    RcSmatT = TRANSPOSE(RcSmat)
    RrSmatT = TRANSPOSE(RrSmat)
    call multiply(RrSmatT,ZrS,Zr2,nb*npceout,nri*npceout,1)
    call multiply(RcSmatT,ZcS,Zc2,nb*npceout,nci*npceout,1)
    
    !! SUB-STEP 20 : Update
    Zs = Zr2 + Zc2
    
    !! SUB-STEP 21 : Weighting     
    ZsD = matmul(Dmat,Zs)   
    
    !! SUB-STEP : Gather
    call getubg(pid,nb,nbg,npceout,RZb,ZsD)
    CALL MPI_REDUCE(RZb,Zb_g,nbg*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

!!-----------------------------------------------------------------------------------------------
!!-----------------------------------------------------------------------------------------------
    !! STEP 8 : Compute
    if (pid .eq. 0) Pb_g = Zb_g
    if (pid .eq. 0) rho_next = dot_product(rb_g,Zb_g)

!!-----------------------------------------------------------------------------------------------
    !! STEP 10 : For loop : Main Iterations Start
    out2(:) = 0.0d0                        
    do i = 1,maxiter
    
        !! STEP 12 : Scatter
        CALL MPI_BCAST(Pb_g,nbg*npceout,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call getub(pid,nb,nbg,npceout,Pb_g,Pb)
        
        !! STEP 13 : Solve 
        call multiply(Aig,Pb,out2,(np-nb)*npceout,nb*npceout,1)
        Ui = out2
        call dpotrs('L',(np-nb)*npceout,1,Aii,(np-nb)*npceout,Ui,(np-nb)*npceout,ierr)
        
        !! STEP 14 : Compute
        call multiply(Agi,Ui,out2,nb*npceout,(np-nb)*npceout,1)
        call multiply(Agg,Pb,Qb,nb*npceout,nb*npceout,1)
        Qb = Qb - out2
        
        !! STEP 15 : Gather 
        call getubg(pid,nb,nbg,npceout,RQb,Qb)
        CALL MPI_REDUCE(RQb,Qb_g,nbg*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        
        !! STEP 16, 17, 18 & 19 : Compute & Update
        if (pid .eq. 0) then
        rho_curr = rho_next
        alpha = rho_curr/dot_product(Qb_g,Pb_g)
        err = alpha*alpha*dot_product(Pb_g,Pb_g)/dot_product(Ub_g,Ub_g)
        Ub_g = Ub_g + alpha*Pb_g
        rb_g = rb_g - alpha*Qb_g

!!-----------------------------------------------------------------------------------------------
!!  Uncomment this section for Residual Calculation !!!!!
!!-----------------------------------------------------------------------------------------------
        !residual = SQRT(dot_product(rb_g,rb_g))
        !print*,residual
        !open(unit=7,file="residual.dat", position="append", status='old')
        !write (7,*), residual
        !close (7)

!!-----------------------------------------------------------------------------------------------
        !! STEP 20 : Exit
        print*, '----------------'
        print*, 'Main Iteration #',i, ', relative error of ',err
        print*, '----------------'
        !!print*, 'and sum of Ub of ', sum(Ub_g)
        end if
        
        CALL MPI_BCAST(err,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        if (err .lt. tol) exit

!!-----------------------------------------------------------------------------------------------
!! Two-Level Newmann-Newmann Precondtioner Using : Algorithm 7 from the report
!!-----------------------------------------------------------------------------------------------
        !! SUB-STEP 3 : Scatter
        CALL MPI_BCAST(rb_g,nbg*npceout,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call getub(pid,nb,nbg,npceout,rb_g,rb)                  ! -> rb
        
        !! SUB-STEP 4 : Weighting 
        call multiply(Dmat,rb,yb,nb*npceout,nb*npceout,1)       ! yb = matmul(Dmat,rb)

        !! SUB-STEP 5 & 6 : Compute
        call multiply(RrSmat,yb,FrS,nri*npceout,nb*npceout,1)   !FrS = matmul(RrSmat,yb) !
        call multiply(RcSmat,yb,FcS,nci*npceout,nb*npceout,1)   !FcS = matmul(RcSmat,yb)

!!-----------------------------------------------------------------------------------------------
        !! SUB-STEP 7 : Solve 
        !call solve(nri*npceout,1,Srr,v1,FrS)                    ! -> v1
        !! ybc(:) = 0.0d0
        ybc(((np-nb)*npceout+1):(np-nci)*npceout) = FrS
        out12 = ybc
        call dpotrs('L',(np-nci)*npceout,1,Amat,(np-nci)*npceout,out12,(np-nci)*npceout,ierr)
        v1 = out12(((np-nb)*npceout+1):(np-nci)*npceout)

        !! SUB-STEP 8 : Compute 
        !!!------------------- Section3/Parallel MatrixVector Product ---------------------!!!
        !!! Scr*v1 = [Acr - (Aci *inv(Aii)*Air)]*v1
        call multiply(Air,v1,u11,(np-nb)*npceout,nri*npceout,1)
        u12 = u11
        call dpotrs('L',(np-nb)*npceout,1,Aii,(np-nb)*npceout,u12,(np-nb)*npceout,ierr)
        call multiply(Aci,u12,u23,nci*npceout,(np-nb)*npceout,1)
        call multiply(Acr,v1,u24,nci*npceout,nri*npceout,1)
        
        !! SUB-STEP 9 : Update 
        dcS = u24 - u23

        !! SUB-STEP 10 : Gather
        dcS = FcS - dcS
        call multiply(BcSmatT,dcS,DcBc,nbgc*npceout,nci*npceout,1)     
        CALL MPI_ALLREDUCE(DcBc,Dc,nbgc*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        !!CALL MPI_REDUCE(DcBc,Dc,nbgc*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

!!-----------------------------------------------------------------------------------------------
!!  SUB-STEP 12 : Coarse Problem: LUMPED-PCGM : Algorithm 9 from the report
!!-----------------------------------------------------------------------------------------------
        call solve_coarsepcgm(pid,nb,np,nci,nri,nbg,nbgc,npceout,BcSmat,BcSmatT, & 
                     Aii,Acc,Aic,Air,Aci,Ari,Arc,Acr,Amat,Minc,Dc,Zc,tolc)

!!-----------------------------------------------------------------------------------------------
        !! SUB-STEP 14 : Scatter
        CALL MPI_BCAST(Zc,nbgc*npceout,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call multiply(BcSmat,Zc,ZcS,nci*npceout,nbgc*npceout,1)           

        !! SUB-STEP 15 : Compute
        !!!-------------------Section1/Parallel MatrixVector Product ---------------------!!!
        !!! Src*ZcS = [Arc - (Ari *inv(Aii)*Air)]*ZcS
        call multiply(Aic,ZcS,u11,(np-nb)*npceout,nci*npceout,1)
        u12 = u11
        call dpotrs('L',(np-nb)*npceout,1,Aii,(np-nb)*npceout,u12,(np-nb)*npceout,ierr)
        call multiply(Ari,u12,u13,nri*npceout,(np-nb)*npceout,1)
        call multiply(Arc,ZcS,u14,nri*npceout,nci*npceout,1)
        v2 = u14 - u13
        
        !! SUB-STEP 16 : Update
        v2 = FrS - v2

!!-----------------------------------------------------------------------------------------------
        !! SUB-STEP 17 : Solve
        ybc(((np-nb)*npceout+1):(np-nci)*npceout) = v2
        out12 = ybc
        call dpotrs('L',(np-nci)*npceout,1,Amat,(np-nci)*npceout,out12,(np-nci)*npceout,ierr)
        ZrS = out12(((np-nb)*npceout+1):(np-nci)*npceout)

!!-----------------------------------------------------------------------------------------------
        !! SUB-STEP 18 & 19 : Compute
        call multiply(RrSmatT,ZrS,Zr2,nb*npceout,nri*npceout,1)
        call multiply(RcSmatT,ZcS,Zc2,nb*npceout,nci*npceout,1)
        
        !! SUB-STEP 20 : Update
        Zs = Zr2 + Zc2
        
        !! SUB-STEP 21 : Weighting 
        ZsD = matmul(Dmat,Zs)    
        
        !! SUB-STEP 22 : Gather
        call getubg(pid,nb,nbg,npceout,RZb,ZsD)
        CALL MPI_REDUCE(RZb,Zb_g,nbg*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

!!-----------------------------------------------------------------------------------------------
        !! STEP 22, 23 & 24 : Compute & Update
        if (pid .eq. 0) then
            rho_next = dot_product(rb_g,Zb_g)
            beta = rho_next/rho_curr
            Pb_g = Zb_g + beta*Pb_g
        end if
    end do

!! Main PCGM -iterations ends here
!!-----------------------------------------------------------------------------------------------
!!-----------------------------------------------------------------------------------------------
    !! STEP 28 : Scatter
    CALL MPI_BCAST(Ub_g,nbg*npceout,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call getub(pid,nb,nbg,npceout,Ub_g,Ub)
    
    !! STEP 29 : Compute & Solve 
    call multiply(Aig,Ub,out2,(np-nb)*npceout,nb*npceout,1)
    Ui = Fi-out2
    call dpotrs('L',(np-nb)*npceout,1,Aii,(np-nb)*npceout,Ui,(np-nb)*npceout,ierr)


END SUBROUTINE solve_nncpcgm



!!!!---------------------------------------------------------------------------------------------
!!!!!!!!!!!!!!! The FETI-DP Procedure with Dirichalet-PCGM Solver !!!!!!!!!!!!!!
!!!!---------------------------------------------------------------------------------------------
SUBROUTINE solve_fetidp(pid,nb,np,nbg,npceout,nri,nci,nbgc,nbgr,Aii,Agg,Air,Aic,Ari,Aci,&
                        Arr,Acc,Arc,Acr,Fi,Fr,Fc,UiS,UrS,UcS,Uc)

    use mpi

    !! Variable Definitions
    integer :: maxiter = 100
    integer :: i, j, ierr, nb, np, nbg, pid
    integer :: nci, nri, nbgc, nbgr, npceout
    double precision :: alpha, beta, err, rho_next, rho_curr, residual
    double precision :: tol = 1.0d-11     !! Main PCGM
    double precision :: tolc = 1.0d-11    !! Coarse PCGM

    double precision :: Aii((np-nb)*npceout,(np-nb)*npceout),Agg(nb*npceout,nb*npceout)
    double precision :: Air((np-nb)*npceout,(nb-nci)*npceout), Aic((np-nb)*npceout,nci*npceout) 
    double precision :: Arr((nb-nci)*npceout,(nb-nci)*npceout), Acc(nci*npceout,nci*npceout)
    double precision :: Ari((nb-nci)*npceout,(np-nb)*npceout), Aci(nci*npceout,(np-nb)*npceout)
    double precision :: Arc((nb-nci)*npceout,nci*npceout), Acr(nci*npceout,(nb-nci)*npceout)
    double precision :: Amat((np-nci)*npceout,(np-nci)*npceout),Minc(nci*npceout,nci*npceout)

    double precision :: BcSmat(nci*npceout,nbgc*npceout), BcSmatT(nbgc*npceout,nci*npceout) 
    double precision :: BrSmat(nbgr*npceout,nri*npceout), BrSmatT(nri*npceout,nbgr*npceout)
    double precision :: Dc(nbgc*npceout), Dr(nbgr*npceout), DrS2(nbgr*npceout)
    double precision :: Dr2(nbgr*npceout),DrSmat(nri*npceout,nri*npceout) 
    double precision :: DcBc(nbgc*npceout),DrBr(nbgr*npceout),Dcc(nbgc*npceout),dcS(nci*npceout)

    double precision :: Fi((np-nb)*npceout), Fc(nci*npceout), Fr(nri*npceout)
    double precision :: Gr(nri*npceout),Gc(nci*npceout), rhsR2(nri*npceout)
    double precision :: out2((np-nb)*npceout),out12((np-nci)*npceout),out3c((np-nci)*npceout)
    double precision :: Pb_g(nbgr*npceout), Pb(nri*npceout)
    double precision :: Qbr(nri*npceout), Qbr2(nri*npceout),Qb_g(nbgr*npceout)
    double precision :: RZbr(nbgr*npceout), rb_g(nbgr*npceout)
    double precision :: rb(nri*npceout), rbs(nri*npceout), rhsC(nci*npceout)
    double precision :: RQbr(nbgr*npceout), rhsR(nri*npceout)

    double precision :: Ub_g(nbgr*npceout),Uc(nbgc*npceout)
    double precision :: Ubot(nri*npceout),Utop((np-nb)*npceout),UiUr((np-nci)*npceout)
    double precision :: u11((np-nb)*npceout),u12((np-nb)*npceout),u13(nri*npceout)
    double precision :: u14(nri*npceout), u15(nri*npceout), ZcS(nci*npceout)
    double precision :: u22((np-nci)*npceout),u23(nci*npceout),u24(nci*npceout)
    double precision :: Vs1(nri*npceout), Vs2(nci*npceout), V2S(nbgc*npceout)
    double precision :: V2(nbgc*npceout),V3(nbgc*npceout),Vs3(nci*npceout),Vs4(nri*npceout)
    double precision :: Vs5(nri*npceout), Vc2(nbgc*npceout), ybc((np-nci)*npceout)
    double precision :: Zb_g(nbgr*npceout), ZsDr(nri*npceout), Zc(nbgc*npceout)
    double precision :: UrS(nri*npceout), UiS((np-nb)*npceout), UcS(nci*npceout)

!!-----------------------------------------------------------------------------------------------
!! To create Block Matrix Amat, used for Dirichlet & Naumann Solve
    Amat(1:(np-nb)*npceout,1:(np-nb)*npceout) = Aii
    Amat(1:(np-nb)*npceout,((np-nb)*npceout+1):(np-nci)*npceout) = Air
    Amat(((np-nb)*npceout+1):(np-nci)*npceout,1:(np-nb)*npceout) = Ari
    Amat(((np-nb)*npceout+1):(np-nci)*npceout,((np-nb)*npceout+1):(np-nci)*npceout) = Arr
    !Aii =  Amat(1:(np-nb),1:(np-nb))
    !Air  = Amat(1:(np-nb),(np-nb)+1:(np-nci))
    !Ari  = Amat((np-nb)+1:(np-nci),1:(np-nb))
    !Arr  = Amat((np-nb)+1:(np-nci),(np-nb)+1:(np-nci))

!! Initiate Array    !!Minv = Agg
    Ub_g(:) = 0.0d0
    Minc = Acc

!!-----------------------------------------------------------------------------------------------
!! Cholesky Decomposition of As to solve x = inv(As)*b
    if (pid .eq. 0) print*, 'Obtaining Cholesky Decomposition of Aii (',(np-nb)*npceout,')'
    call dpotrf('L',(np-nb)*npceout,Aii,(np-nb)*npceout,ierr)
    if (ierr .ne. 0) print*, 'Error in factorizing Aii! Processor # ', pid


    if (pid .eq. 0) print*, 'Obtaining Cholesky Decomposition of Amat(',(np-nci)*npceout,')'
    call dpotrf('L',(np-nci)*npceout,Amat,(np-nci)*npceout,ierr)
    if (ierr .ne. 0) print*, 'Error in factorizing Amat! Processor # ', pid


    if (pid .eq. 0) print*, 'Obtaining Cholesky Decomposition of Acc(',nci*npceout,')'
    call dpotrf('L',(nb-nri)*npceout,Minc,(nb-nri)*npceout,ierr)
    if (ierr .ne. 0) print*, 'Error in factorizing Acc! Processor # ', pid

!!----------------------------------------------------------------------------------------------
    if (pid .eq. 0) print*, '----------------------------------------'
    if (pid .eq. 0) print*, 'Initializing Stochastic-FETIDP-Solver...'
    if (pid .eq. 0) print*, '----------------------------------------'
    
    call getBrS(pid, nri, nbgr, npceout, BrSmatT)
    BrSmat = TRANSPOSE(BrSmatT)
    call getBcS(pid, nci, nbgc, npceout, BcSmat)
    BcSmatT = TRANSPOSE(BcSmat)

!!!----------------------------------------------------------------------------------------------
!!! Construct RHS Part-1 : Dc = SUM_1^ns[BcS'(Gc - Scr*inv(Srr)*Gr)]
!!!----------------------------------------------------------------------------------------------
    !! STEP 3 : Solve
    out2 = Fi
    call dpotrs('L',(np-nb)*npceout,1,Aii,(np-nb)*npceout,out2,(np-nb)*npceout,ierr)
    call multiply(Ari,out2,Gr,nri*npceout,(np-nb)*npceout,1)

    !! STEP 4 & 7 : Compute
    Gr = Fr - Gr
    call multiply(Aci,out2,Gc,nci*npceout,(np-nb)*npceout,1)
    Gc = Fc - Gc

    !! STEP 5 : Solve
    !!!-------------------inv(Srr)*rhs---------------------!!!
    ybc(:) = 0.0d0
    ybc(((np-nb)*npceout+1):(np-nci)*npceout) = Gr
    out12 = ybc
    call dpotrs('L',(np-nci)*npceout,1,Amat,(np-nci)*npceout,out12,(np-nci)*npceout,ierr)
    rhsR = out12(((np-nb)*npceout+1):(np-nci)*npceout)

    !! STEP 6 : Solve
    !!!------------------- Scr*rhs---------------------!!!
    !!! Scr*v1 = [Acr - (Aci *inv(Aii)*Air)]*v1
    call multiply(Air,rhsR,u11,(np-nb)*npceout,nri*npceout,1)
    u12 = u11
    call dpotrs('L',(np-nb)*npceout,1,Aii,(np-nb)*npceout,u12,(np-nb)*npceout,ierr)
    call multiply(Aci,u12,u23,nci*npceout,(np-nb)*npceout,1)
    call multiply(Acr,rhsR,u24,nci*npceout,nri*npceout,1)
    rhsC = u24 - u23

    dcS = Gc - rhsC
    call multiply(BcSmatT,dcS,DcBc,nbgc*npceout,nci*npceout,1)
    call multiply(BrSmat,rhsR,DrBr,nbgr*npceout,nri*npceout,1)

    CALL MPI_ALLREDUCE(DcBc,Dc,nbgc*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    CALL MPI_REDUCE(DrBr,Dr,nbgr*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

    if (pid .eq. 0) then
    Dcc = Dc
    end if

!!!----------------------------------------------------------------------------------------------
!!! STEP 11: Solve RHS Coarse Problem: LUMPED-PCGM : Algorithm 9 from the report : Fcc Zc = Dc
!!!----------------------------------------------------------------------------------------------
    if (pid .eq. 0) print*, 'RHS CoarseProblem Initiated...'
    call solve_coarsepcgm(pid,nb,np,nci,nri,nbg,nbgc,npceout,BcSmat,BcSmatT,Aii,Acc, &
                          Aic,Air,Aci,Ari,Arc,Acr,Amat,Minc,Dc,Zc,tolc)

!!!----------------------------------------------------------------------------------------------
!!! Construct RHS Part-2 : Frc*Zc =  Frc = SUM_1^ns[BrS(inv(Srr)*Src)BcS]*Zc
!!!----------------------------------------------------------------------------------------------
    CALL MPI_BCAST(Zc,nbgc*npceout,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)  !! This line was missing
    call multiply(BcSmat,Zc,ZcS,nci*npceout,nbgc*npceout,1)

    !!!------------------- Src*rhs ---------------------!!!
    !!! Src*v1 = [Arc - (Ari *inv(Aii)*Aic)]*v1
    call multiply(Aic,ZcS,u11,(np-nb)*npceout,nci*npceout,1)
    u12 = u11
    call dpotrs('L',(np-nb)*npceout,1,Aii,(np-nb)*npceout,u12,(np-nb)*npceout,ierr)
    call multiply(Ari,u12,u13,nri*npceout,(np-nb)*npceout,1)
    call multiply(Arc,ZcS,u14,nri*npceout,nci*npceout,1)
    rhsR = u14 - u13

    !!!-------------------inv(Srr)*rhsR---------------------!!!
    u22(:) = 0.0d0
    u22(((np-nb)*npceout+1):(np-nci)*npceout)=rhsR          !u1
    out3c = u22
    call dpotrs('L',(np-nci)*npceout,1,Amat,(np-nci)*npceout,out3c,(np-nci)*npceout,ierr)
    rhsR2 = out3c(((np-nb)*npceout+1):(np-nci)*npceout)     !u2

    call multiply(BrSmat,rhsR2,DrS2,nbgr*npceout,nri*npceout,1)
    CALL MPI_REDUCE(DrS2,Dr2,nbgr*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

    if (pid .eq. 0)  then                 
    rb_g = Dr - Dr2
    end if

!!!----------------------------------------------------------------------------------------------
!!! STEP 20: Preconditioning : Z0 = inv(M_d) R0  =  M^-1 = SUM_1^ns[BrS*DrS (Srr) DrS BrS^T]
!!!----------------------------------------------------------------------------------------------
    CALL MPI_BCAST(rb_g,nbgr*npceout,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call getDrS(pid, nbgr, nri, npceout, DrSmat)                  ! -> Dmat
    call multiply(BrSmatT,rb_g,rb,nri*npceout,nbgr*npceout,1)
    call multiply(DrSmat,rb,rbs,nri*npceout,nri*npceout,1)

    !!!------------------- Srr*rhs ---------------------!!!
    !!! Srr*v1 = [Arr - (Ari *inv(Aii)*Air)]*v1
    call multiply(Air,rbs,u11,(np-nb)*npceout,nri*npceout,1)
    u12 = u11
    call dpotrs('L',(np-nb)*npceout,1,Aii,(np-nb)*npceout,u12,(np-nb)*npceout,ierr)
    call multiply(Ari,u12,u13,nri*npceout,(np-nb)*npceout,1)
    call multiply(Arr,rbs,u14,nri*npceout,nri*npceout,1)
    rhsR = u14 - u13

    call multiply(DrSmat,rhsR,ZsDr,nri*npceout,nri*npceout,1)
    call multiply(BrSmat,ZsDr,RZbr,nbgr*npceout,nri*npceout,1)
    CALL MPI_REDUCE(RZbr,Zb_g,nbgr*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

!!----------------------------------------------------------------------------------------------
    if (pid .eq. 0) Pb_g = Zb_g
    if (pid .eq. 0) rho_next = dot_product(rb_g,Zb_g)

!!!----------------------------------------------------------------------------------------------
!!! STEP 29 : Main Iteration Start
!!!----------------------------------------------------------------------------------------------
    do i = 1,maxiter
        CALL MPI_BCAST(Pb_g,nbgr*npceout,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call multiply(BrSmatT,Pb_g,Pb,nri*npceout,nbgr*npceout,1)


        !!!-------------------inv(Srr)*rhsR---------------------!!! 
        u22(:) = 0.0d0
        u22(((np-nb)*npceout+1):(np-nci)*npceout)=Pb
        out3c = u22
        call dpotrs('L',(np-nci)*npceout,1,Amat,(np-nci)*npceout,out3c,(np-nci)*npceout,ierr)
        Vs1 = out3c(((np-nb)*npceout+1):(np-nci)*npceout)     !u2

        !!!------------------- Scr*rhs ---------------------!!!
        call multiply(Air,Vs1,u11,(np-nb)*npceout,nri*npceout,1)
        u12 = u11
        call dpotrs('L',(np-nb)*npceout,1,Aii,(np-nb)*npceout,u12,(np-nb)*npceout,ierr)
        call multiply(Aci,u12,u23,nci*npceout,(np-nb)*npceout,1)
        call multiply(Acr,Vs1,u24,nci*npceout,nri*npceout,1)
        Vs2 = u24 - u23

        call multiply(BcSmatT,Vs2,V2S,nbgc*npceout,nci*npceout,1)
        CALL MPI_ALLREDUCE(V2S,V2,nbgc*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

!!!----------------------------------------------------------------------------------------------
!!! STEP 36 : Solve LHS Coarse Problem: LUMPED-PCGM : Algorithm 9 from the report : Fcc V3 = V2
!!!----------------------------------------------------------------------------------------------
        if (pid .eq. 0) print*, 'LHS CoarseProblem Initiated...'
        call solve_coarsepcgm(pid,nb,np,nci,nri,nbg,nbgc,npceout,BcSmat,BcSmatT,Aii,Acc, &
                              Aic,Air,Aci,Ari,Arc,Acr,Amat,Minc,V2,V3,tolc)

!!----------------------------------------------------------------------------------------------

        CALL MPI_BCAST(V3,nbgc*npceout,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        CALL multiply(BcSmat,V3,Vs3,nci*npceout,nbgc*npceout,1)

        !!!------------------- Src*v1 ---------------------!!!
        call multiply(Aic,Vs3,u11,(np-nb)*npceout,nci*npceout,1)
        !      if (direct .eq. 1) then
        u12 = u11
        call dpotrs('L',(np-nb)*npceout,1,Aii,(np-nb)*npceout,u12,(np-nb)*npceout,ierr)
        !       else
        !          call pcgm(pid,(np-nb)*npceout,Aii,Ui,out2,10000,1.0d-10,1,MAii)
        !       end if
        call multiply(Ari,u12,u13,nri*npceout,(np-nb)*npceout,1)
        call multiply(Arc,Vs3,u14,nri*npceout,nci*npceout,1)
        Vs4 = u14 - u13

        Vs5 = Pb + Vs4

        !!!-------------------inv(Srr)*rhsR---------------------!!! 
        u22(:) = 0.0d0
        u22(((np-nb)*npceout+1):(np-nci)*npceout)=Vs5          !u1
        out3c = u22
        call dpotrs('L',(np-nci)*npceout,1,Amat,(np-nci)*npceout,out3c,(np-nci)*npceout,ierr)
        Qbr2 = out3c(((np-nb)*npceout+1):(np-nci)*npceout)     !u2
        call multiply(BrSmat,Qbr2,RQbr,nbgr*npceout,nri*npceout,1)
        CALL MPI_REDUCE(RQbr,Qb_g,nbgr*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

        if (pid .eq. 0) then
        rho_curr = rho_next
        alpha = rho_curr/dot_product(Qb_g,Pb_g)  
        err = alpha*alpha*dot_product(Pb_g,Pb_g)/dot_product(Ub_g,Ub_g)
        Ub_g = Ub_g + alpha*Pb_g
        rb_g = rb_g - alpha*Qb_g
        
        print*, '----------------'
        print*, 'Main Iteration #',i, ', relative error of ',err
        print*, '----------------'
        !!print*, 'and sum of Ub of ', sum(Ub_g)
        end if

        CALL MPI_BCAST(err,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        if (err .lt. tol) exit

!!!----------------------------------------------------------------------------------------------
!!! STEP 49: Preconditioning : Zi = inv(M_d) Ri  =  M^-1 = SUM_1^ns[BrS DrS (Srr) DrS BrS^T]
!!!----------------------------------------------------------------------------------------------
        CALL MPI_BCAST(rb_g,nbgr*npceout,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call multiply(BrSmatT,rb_g,rb,nri*npceout,nbgr*npceout,1)
        call multiply(DrSmat,rb,rbs,nri*npceout,nri*npceout,1)

        !!!-------------------Srr/Parallel MatrixVector Product ---------------------!!!
        !!! Srr*v1 = [Arr - (Ari *inv(Aii)*Air)]*v1
        call multiply(Air,rbs,u11,(np-nb)*npceout,nri*npceout,1)
        !      if (direct .eq. 1) then
        u12 = u11
        call dpotrs('L',(np-nb)*npceout,1,Aii,(np-nb)*npceout,u12,(np-nb)*npceout,ierr)
        !       else
        !          call pcgm(pid,(np-nb)*npceout,Aii,Ui,out2,10000,1.0d-10,1,MAii)
        !       end if
        call multiply(Ari,u12,u13,nri*npceout,(np-nb)*npceout,1)
        call multiply(Arr,rbs,u14,nri*npceout,nri*npceout,1)
        rhsR = u14 - u13

        call multiply(DrSmat,rhsR,ZsDr,nri*npceout,nri*npceout,1)
        call multiply(BrSmat,ZsDr,RZbr,nbgr*npceout,nri*npceout,1)
        CALL MPI_REDUCE(RZbr,Zb_g,nbgr*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

!!!----------------------------------------------------------------------------------------------
        if (pid .eq. 0) then
        rho_next = dot_product(rb_g,Zb_g)
        beta = rho_next/rho_curr
        Pb_g = Zb_g + beta*Pb_g
        end if
    end do

!!!-----------------------------------------------------------------------------------------------
!! Uc = inv(Fcc)*(Dc + Fcr*Ub_g) = inv(Fcc)*Dcc
!!!-----------------------------------------------------------------------------------------------
    CALL MPI_BCAST(Ub_g,nbgr*npceout,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call multiply(BrSmatT,Ub_g,Pb,nri*npceout,nbgr*npceout,1)

    !! call solve(nri*1,1,Srr,Vs1,Pb)
    !!!-------------------Section2/inv(Srr)*rhsR---------------------!!! 
    u22(:) = 0.0d0
    u22(((np-nb)*npceout+1):(np-nci)*npceout)=Pb          !u1
    out3c = u22
    call dpotrs('L',(np-nci)*npceout,1,Amat,(np-nci)*npceout,out3c,(np-nci)*npceout,ierr)
    Vs1 = out3c(((np-nb)*npceout+1):(np-nci)*npceout)     !u2

    !!!------------------- Scr*rhs ---------------------!!!
    call multiply(Air,Vs1,u11,(np-nb)*npceout,nri*npceout,1)
    u12 = u11
    call dpotrs('L',(np-nb)*npceout,1,Aii,(np-nb)*npceout,u12,(np-nb)*npceout,ierr)
    call multiply(Aci,u12,u23,nci*npceout,(np-nb)*npceout,1)
    call multiply(Acr,Vs1,u24,nci*npceout,nri*npceout,1)
    Vs2 = u24 - u23

    call multiply(BcSmatT,Vs2,V2S,nbgc*npceout,nci*npceout,1)
    CALL MPI_REDUCE(V2S,Vc2,nbgc*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

    if (pid .eq. 0) then
    Dcc = Dcc + Vc2
    endif

!!!-----------------------------------------------------------------------------------------------
!!! Coarse Problem To Find Uc : LUMPED-PCGM : Algorithm 9 from the report : Fcc Uc = Dcc
!!!-----------------------------------------------------------------------------------------------
    if (pid .eq. 0) print*, 'PCGM for CoarseProblem to find Uc...'
    CALL MPI_BCAST(Dcc,nbgc*npceout,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call solve_coarsepcgm(pid,nb,np,nci,nri,nbg,nbgc,npceout,BcSmat,BcSmatT,Aii,Acc, &
                          Aic,Air,Aci,Ari,Arc,Acr,Amat,Minc,Dcc,Uc,tolc)

!!-----------------------------------------------------------------------------------------------
    CALL MPI_BCAST(Uc,nbgc*npceout,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
    call multiply(BcSmat,Uc,UcS,nci*npceout,nbgc*npceout,1)
    call multiply(Aic,UcS,u11,(np-nb)*npceout,nci*npceout,1)
    Utop = Fi - u11

    call multiply(Arc,UcS,u15,nri*npceout,nci*npceout,1)
    Ubot = Fr - Pb - u15

!!!-------------------inv(Srr)*rhsR---------------------!!! 
    u22(:) = 0.0d0
    u22(1:((np-nb)*npceout))= Utop          
    u22(((np-nb)*npceout+1):(np-nci)*npceout)= Ubot          !u1
    UiUr = u22
    call dpotrs('L',(np-nci)*npceout,1,Amat,(np-nci)*npceout,UiUr,(np-nci)*npceout,ierr)
    UiS = UiUr(1:((np-nb)*npceout))     
    UrS = UiUr(((np-nb)*npceout+1):(np-nci)*npceout)


END SUBROUTINE solve_fetidp
!!-----------------------------------------------------------------------------------------------

!!!!---------------------------------------------------------------------------------------------
!!!!!!!!!!!!!!! Coarse Problem: Lumped-PCGM Solver !!!!!!!!!!!!!!
!!!!---------------------------------------------------------------------------------------------
SUBROUTINE solve_coarsepcgm(pid,nb,np,nci,nri,nbg,nbgc,npceout,BcSmat,BcSmatT,Aii,Acc,Aic,Air,Aci,Ari,Arc,Acr,Amat,Minc,Dc,Zc,tolc)

    use mpi

    integer :: maxiter = 100                                        
    integer :: pid,nb,np,nbg,ierr,j,npceout,nci,nri,nbgc
    double precision :: tolc 

    double precision :: Aii((np-nb)*npceout,(np-nb)*npceout), Acc(nci*npceout,nci*npceout)
    double precision :: Air((np-nb)*npceout,(nb-nci)*npceout), Aic((np-nb)*npceout,nci*npceout)
    double precision :: Ari((nb-nci)*npceout,(np-nb)*npceout), Aci(nci*npceout,(np-nb)*npceout)
    double precision :: Arc((nb-nci)*npceout,nci*npceout), Acr(nci*npceout,(nb-nci)*npceout)
    double precision :: BcSmat(nci*npceout,nbgc*npceout), BcSmatT(nbgc*npceout,nci*npceout)
    double precision :: Amat((np-nci)*npceout,(np-nci)*npceout)
    double precision :: Dc(nbgc*npceout), Zc(nbgc*npceout)

    double precision :: u11((np-nb)*npceout),u12((np-nb)*npceout),u13(nri*npceout)
    double precision :: u14(nri*npceout),DcR1(nri*npceout),DcR2(nri*npceout),u22((np-nci)*npceout)
    double precision :: u23(nci*npceout),u24(nci*npceout),DcR3(nci*npceout),DcR4(nci*npceout)
    double precision :: Minc(nci*npceout,nci*npceout),Pb_gc(nbgc*npceout),out3c((np-nci)*npceout)
    double precision :: RZbc(nbgc*npceout),DcR(nci*npceout),out13(nci*npceout),Zb_gc(nbgc*npceout)
    double precision :: Qbr(nci*npceout),RQbc(nbgc*npceout),Qb_gc(nbgc*npceout)
    double precision :: rho_nextc,rho_currc, alphac, betac, errc

    
!!STEP 10 : Preconditioning 
    !if (pid .eq. 0) print*, 'PCGM for CourseProblem...'
!!STEP 11 : Scatter
    call multiply(BcSmat,Dc,DcR,nci*npceout,nbgc*npceout,1)
!!STEP 12 : Solve
    out13 = DcR
    call dpotrs('L',nci*npceout,1,Minc,nci*npceout,out13,nci*npceout,ierr)
!!STEP 14 : Gather    
    call multiply(BcSmatT,out13,RZbc,nbgc*npceout,nci*npceout,1)
    CALL MPI_REDUCE(RZbc,Zb_gc,nbgc*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

!!STEP 16 & 17 : Compute 
    if (pid .eq. 0) Pb_gc = Zb_gc
    if (pid .eq. 0) rho_nextc = dot_product(Dc,Zb_gc)
    
!!STEP 17 : Main Iteration 
    Zc(:) = 0.0d0
    do j = 1,maxiter
    
!!STEP 19 : Scatter 
        CALL MPI_BCAST(Pb_gc,nbgc*npceout,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call multiply(BcSmat,Pb_gc,DcR,nci*npceout,nbgc*npceout,1)

        !!!-------------------Section1/Parallel MatrixVector Product ---------------------!!!
        !!! Src*v1 = [Arc - (Ari *inv(Aii)*Aic)]*DcR
        call multiply(Aic,DcR,u11,(np-nb)*npceout,nci*npceout,1)
        u12 = u11
        call dpotrs('L',(np-nb)*npceout,1,Aii,(np-nb)*npceout,u12,(np-nb)*npceout,ierr)
        call multiply(Ari,u12,u13,nri*npceout,(np-nb)*npceout,1)
        call multiply(Arc,DcR,u14,nri*npceout,nci*npceout,1)
        DcR1 = u14 - u13
!!STEP 20 : Compute
        !!!-------------------Section2/Parallel MatrixVector Product---------------------!!!
        !!! inv(Srr)*DcR2
        u22(:) = 0.0d0
        u22(((np-nb)*npceout+1):(np-nci)*npceout)=DcR1          
        out3c = u22
        call dpotrs('L',(np-nci)*npceout,1,Amat,(np-nci)*npceout,out3c,(np-nci)*npceout,ierr)
        DcR2 = out3c(((np-nb)*npceout+1):(np-nci)*npceout)      
!!STEP 21 : Solve
        !!!-------------------Section3/Parallel MatrixVector Product ---------------------!!!
        !!! Scr*DcR2 = [Acr - (Aci *inv(Aii)*Air)]*DcR2
        call multiply(Air,DcR2,u11,(np-nb)*npceout,nri*npceout,1)
        u12 = u11
        call dpotrs('L',(np-nb)*npceout,1,Aii,(np-nb)*npceout,u12,(np-nb)*npceout,ierr)
        call multiply(Aci,u12,u23,nci*npceout,(np-nb)*npceout,1)
        call multiply(Acr,DcR2,u24,nci*npceout,nri*npceout,1)
        DcR3 = u24 - u23
!!STEP 22 : Compute
        !!!-------------------Section4/Parallel MatrixVector Product ---------------------!!!
        !!! Scc*DcR = [Acc - (Aci *inv(Aii)*Aic)]*DcR
        call multiply(Aic,DcR,u11,(np-nb)*npceout,nci*npceout,1)
        u12 = u11
        call dpotrs('L',(np-nb)*npceout,1,Aii,(np-nb)*npceout,u12,(np-nb)*npceout,ierr)
        call multiply(Aci,u12,u23,nci*npceout,(np-nb)*npceout,1)
        call multiply(Acc,DcR,u24,nci*npceout,nci*npceout,1)
        DcR4 = u24 - u23
!!STEP 23 : Compute
        !!----------------------------------------------------------------------------------------
        
!!STEP 24 : Update 
        Qbr = DcR4 - DcR3
!!STEP 25 : Gather 
        call multiply(BcSmatT,Qbr,RQbc,nbgc*npceout,nci*npceout,1)
        CALL MPI_REDUCE(RQbc,Qb_gc,nbgc*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

!!STEP 27 to 30 : Compute & Update 
        if (pid .eq. 0) then
            rho_currc = rho_nextc
            alphac = rho_currc/dot_product(Qb_gc,Pb_gc)
            errc = alphac*alphac*dot_product(Pb_gc,Pb_gc)/dot_product(Zc,Zc)
            Zc = Zc + alphac*Pb_gc
            Dc = Dc - alphac*Qb_gc
            !print*, 'iteration # ',j,'relative error of', errc
        end if
!!STEP 31 : Exit         
        CALL MPI_BCAST(errc,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        if (errc .lt. tolc) exit

!!STEP 32: Preconditioing 
        CALL MPI_BCAST(Dc,nbgc*npceout,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!!STEP 33 : Scatter        
        call multiply(BcSmat,Dc,DcR,nci*npceout,nbgc*npceout,1)
!!STEP 34 : Solve        
        out13 = DcR
        call dpotrs('L',nci*npceout,1,Minc,nci*npceout,out13,nci*npceout,ierr)
!!STEP 35 : Gather        
        call multiply(BcSmatT,out13,RZbc,nbgc*npceout,nci*npceout,1)
        CALL MPI_REDUCE(RZbc,Zb_gc,nbgc*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

!!STEP 37 to 39 : Compute & Update 
        if (pid .eq. 0) then
        rho_nextc = dot_product(Dc,Zb_gc)
        betac = rho_nextc/rho_currc
        Pb_gc = Zb_gc + betac*Pb_gc
        end if
    end do
    
    
END SUBROUTINE solve_coarsepcgm


!!!---------------------------------------------------------------------------------------------
!!!!!!!!!!!!!! Lumped-Preconditioned BICG-STAB Solver !!!!!!!!!!!!!!
!!!---------------------------------------------------------------------------------------------
SUBROUTINE solve_bicgstab(pid,nb,np,nbg,npceout,Aii,Aig,Agi,Agg,Fi,Fg,Ui,Ub,Ub_g,maxiter,tol,direct)

    use mpi

    integer :: pid,nb,np,npceout,nbg,maxiter,direct
    double precision :: tol
    double precision :: Aii((np-nb)*npceout,(np-nb)*npceout),Aig((np-nb)*npceout,nb*npceout)
    double precision :: Agi(nb*npceout,(np-nb)*npceout),Agg(nb*npceout,nb*npceout)
    double precision :: Fi((np-nb)*npceout),Fg(nb*npceout),Ui((np-nb)*npceout),Ub(nb*npceout)
    double precision :: Ub_g(nbg*npceout)

    integer :: ierr,i
    double precision :: rho_prev,rho_prev2,omega,rho_curr, alpha, beta, err
    double precision :: out1((np-nb)*npceout,nb*npceout),out2((np-nb)*npceout),RGi(nbg*npceout)
    double precision :: Si(nb*npceout,nb*npceout), Gi(nb*npceout), Rmat(nb*npceout,nbg*npceout)
    double precision :: RSiR(nbg*npceout,nbg*npceout),S(nbg*npceout,nbg*npceout),calG(nbg*npceout)
    double precision :: rb_g(nbg*npceout),rb(nb*npceout),Minv(nb*npceout,nb*npceout),Zb(nb*npceout)
    double precision :: RZb(nbg*npceout),Zb_g(nbg*npceout),Pb_g(nbg*npceout),Pb(nb*npceout)
    double precision :: Qb(nb*npceout),RQb(nbg*npceout),Qb_g(nbg*npceout),out3((np-nb)*npceout)
    double precision :: MAii((np-nb)*npceout,(np-nb)*npceout), temp((np-nb),(np-nb))
    double precision :: Aiiinv((np-nb)*npceout,(np-nb)*npceout),Minvinv(nb*npceout,nb*npceout)
    double precision :: rtilde(nbg*npceout), vb_g(nbg*npceout), Ptilde(nbg*npceout)
    double precision :: sb_g(nbg*npceout), stilde(nbg*npceout), tb_g(nbg*npceout)
    double precision :: incr(nbg*npceout), Mt(nbg*npceout), Ms(nbg*npceout)

    !!!---------------------------------------------------------------------------------------------
    Ub_g(:) = 0.0d0
    Minv = Agg

    if (pid .eq. 0) print*, 'inverting Aii, matrix of dimension ',(np-nb)*npceout
    call invert((np-nb)*npceout,Aii,Aiiinv)

    if (pid .eq. 0) print*, 'inverting M, matrix of dimension ',nb*npceout
    call invert(nb*npceout,Minv,Minvinv)

    if (pid .eq. 0) print*, 'Initializing BICGSTAB...'

    ! Get b, right-hand side
    call multiply(Aiiinv,Fi,out2,(np-nb)*npceout,(np-nb)*npceout,1)
    call multiply(Agi,out2,Gi,nb*npceout,(np-nb)*npceout,1)
    Gi = Fg - Gi
    call getubg(pid,nb,nbg,npceout,RGi,Gi)
    CALL MPI_REDUCE(RGi,rb_g,nbg*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    if (pid .eq. 0) rtilde = rb_g

    do i = 1,maxiter
        if (pid .eq. 0) then
           rho_prev = dot_product(rtilde,rb_g)
            if (i .eq. 1) then
               Pb_g = rb_g
            else
               beta = (rho_prev/rho_prev2)*(alpha/omega)
               Pb_g = rb_g + beta*(Pb_g-omega*vb_g)
            end if
        end if

        !solve M ptilde = p
        CALL MPI_BCAST(Pb_g,nbg*npceout,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call getub(pid,nb,nbg,npceout,Pb_g,rb)
        call multiply(Minvinv,rb,Zb,nb*npceout,nb*npceout,1)
        call getubg(pid,nb,nbg,npceout,RZb,Zb)
        CALL MPI_REDUCE(RZb,Ptilde,nbg*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

        ! v = A ptilde
        CALL MPI_BCAST(Ptilde,nbg*npceout,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call getub(pid,nb,nbg,npceout,Ptilde,Pb)
        call multiply(Aig,Pb,out2,(np-nb)*npceout,nb*npceout,1)
        call multiply(Aiiinv,out2,Ui,(np-nb)*npceout,(np-nb)*npceout,1)
        call multiply(Agi,Ui,out2,nb*npceout,(np-nb)*npceout,1)
        call multiply(Agg,Pb,Qb,nb*npceout,nb*npceout,1)
        Qb = Qb - out2
        call getubg(pid,nb,nbg,npceout,RQb,Qb)
        CALL MPI_REDUCE(RQb,vb_g,nbg*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

        if (pid .eq. 0) then
           alpha = rho_prev/dot_product(rtilde,vb_g)
           sb_g = rb_g - alpha*vb_g
        end if

        !solve M stilde = s
        CALL MPI_BCAST(sb_g,nbg*npceout,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call getub(pid,nb,nbg,npceout,sb_g,rb)
        call multiply(Minvinv,rb,Zb,nb*npceout,nb*npceout,1)
        call getubg(pid,nb,nbg,npceout,RZb,Zb)
        CALL MPI_REDUCE(RZb,stilde,nbg*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

        ! t = A stilde
        CALL MPI_BCAST(stilde,nbg*npceout,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call getub(pid,nb,nbg,npceout,stilde,Pb)
        call multiply(Aig,Pb,out2,(np-nb)*npceout,nb*npceout,1)
        call multiply(Aiiinv,out2,Ui,(np-nb)*npceout,(np-nb)*npceout,1)
        call multiply(Agi,Ui,out2,nb*npceout,(np-nb)*npceout,1)
        call multiply(Agg,Pb,Qb,nb*npceout,nb*npceout,1)
        Qb = Qb - out2
        call getubg(pid,nb,nbg,npceout,RQb,Qb)
        CALL MPI_REDUCE(RQb,tb_g,nbg*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

!!!---------------------------------------------------------------------------------------------
        !!!solve M Ms = s
        !!CALL MPI_BCAST(sb_g,nbg*npceout,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        !!call getub(pid,nb,nbg,npceout,sb_g,rb)
        !!call multiply(Minvinv,rb,Zb,nb*npceout,nb*npceout,1)
        !!call getubg(pid,nb,nbg,npceout,RZb,Zb)
        !!CALL MPI_REDUCE(RZb,Ms,nbg*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        !!
        !!!solve M Mt = t
        !!CALL MPI_BCAST(tb_g,nbg*npceout,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        !!call getub(pid,nb,nbg,npceout,tb_g,rb)
        !!call multiply(Minvinv,rb,Zb,nb*npceout,nb*npceout,1)
        !!call getubg(pid,nb,nbg,npceout,RZb,Zb)
        !!CALL MPI_REDUCE(RZb,Mt,nbg*npceout,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

!!!---------------------------------------------------------------------------------------------
        if (pid .eq. 0) then
            omega = dot_product(tb_g,sb_g)/dot_product(tb_g,tb_g)
            !!omega = dot_product(Mt,Ms)/dot_product(Mt,Mt)
            incr = alpha*Ptilde + omega*stilde
            err = dot_product(incr,incr)/dot_product(Ub_g,Ub_g)
            Ub_g = Ub_g + incr
            rb_g = sb_g - omega*tb_g
            rho_prev2 = rho_prev
            print*, 'iteration # ', i, ', relative error of ', err, 'and sum of Ub of ', sum(Ub_g)
        end if

        CALL MPI_BCAST(err,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        if (err .lt. tol) exit
    end do

    CALL MPI_BCAST(Ub_g,nbg*npceout,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

    call getub(pid,nb,nbg,npceout,Ub_g,Ub)
    call multiply(Aig,Ub,out2,(np-nb)*npceout,nb*npceout,1)
    call multiply(Aiiinv,Fi-out2,Ui,(np-nb)*npceout,(np-nb)*npceout,1)

END SUBROUTINE solve_bicgstab


END MODULE solvers

!!-----------------------------------------------------------------------------------------------
!!!!! Use #1 & #2 for Execution Time Calculations !!!!!
!!-----------------------------------------------------------------------------------------------
!!#1
!call cpu_time(time1)
!call system_clock(c1)

!!#2
!if (pid .eq. 0) then
!call cpu_time(time2)
!call system_clock(c2)
!print*, 'Time taken for....', dble(c2-c1)/dble(cr), 'seconds (', time2-time1, 'cpu time)'
!end if
!!-----------------------------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!END!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!