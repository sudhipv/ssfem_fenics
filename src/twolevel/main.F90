!! Main Code: By Ajit Desai, Augest/2015, PETSc Version : March/2016
!! purpose: is to define inputs, load FEM/DDM data and call DDM solvers
!
!input:
!  PCE:KLE-Data:
!      ndim  : random dimension of KLE/PCE expansion used for input represenations
!      nord  : order of basis polynomials used for PCE expansion 
!      npcein: number of PCE terms : npceout
!      Apce  : stochastic-stiffness matrix 
!      Fpce  : stochastic-source vector 
!      Cijk  : multidimentional moments : tripple product 
!      sigma : σ for underlying Guassian process : Strenght of uncertainty
!      omegas: w's of KLE expansion terms : Rogers Book : Scale of uncertainty
!      multipliers  : √λ / √(a1 = (sin(2*w_i a1)/2*w_i a_1)): KLE expansion terms : Rogers Book
!      
!  MESH-Data: 
!      p e t : points: elements: triangles: for local or each subdomain
!      np ne nt     : number of points, elements & triangles for local or each subdomain 
!      pg eg tg     : points: elements: triangles: for global or whole domain 
!      npg neg ntg  : number of points, elements & triangles for global or whole domain 
!      nb nci nri   : number of nodes on boundary, corner & remaining for local or each subdomain 
!      nbg nbgc nbgr: number of nodes on boundary, corner & remaining for global or whole domain 
!      dbounds : physical domain limits : x=[0, 1], y=[0,1] : for unit square 
!
!  MPI-Data: 
!      pid   : Processor ID
!      nproc : number of processors 
!      ierr  : MPI errors 
!
!  FEM-Data: 
!      Amat  : FEM assembly matrix 
!      Aadv  : advection matrix 
!      Adif  : diffusion matrix 
!      fvec  : FEM assembly vector
!
!  DDM-Data:
!      Aii   : interior-interior sub-block of Amat
!      Aig   : interior-boundary sub-block of Amat
!      Air   : interior-remainig sub-block of Amat
!      Aic   : interior-corner sub-block of Amat
!      Fi Fg : interior & boundary sub-block of Fvec
!
! output:
!      Ui Ur Uc : local interior, remaining & corner nodes solution vector 
!      Ub Ub_g  : local & global boundary nodes solution vector
!!
!!---------------------------------------------------------------------------------------
PROGRAM main

    use common
    use myCommon
    use assembly
    use PETScSolvers

    implicit none

include 'mpif.h'

!!------------------------------------------------------
!! Two-Level-NNC-PCGM with PETSc Sparse Iterative Solver
!!------------------------------------------------------

!!MPI-Data:
    integer :: ierr, pid

!!Exit-Criteria
    integer, parameter :: maxiter = 1000
    double precision, parameter :: tol = 1.0d-6

!!PCE-Data:
    !! Select casep
    !! 2: for automatic 
    !! 3: for 3dim, 3ord
    !! 4: for 4dim, 3ord
    !! 5: for 5dim, 3ord
    integer, parameter :: casep = 2
    integer, parameter :: nomga = 9
    integer            :: nord,ndim,npceout,npcein

!!Mallocs-Data:
    !! To precalculate "mallocsCals = 1" : First time run
    !! To use calculated "mallocsCals = 0" : Repeated runs
    integer, parameter :: mallocsCals = 1

!!KLE-Data-Inputs
    integer          :: i, ncijk  !,indexi,j, k, 
    double precision :: amp,const_diff,sigma,gamma
    double precision, dimension(2,2)              :: dbounds
    !double precision, dimension(2)                   :: centre
    integer, allocatable, dimension(:,:)          :: ijk
    double precision, allocatable, dimension(:)   :: cijk
    double precision, allocatable, dimension(:)   :: omegas, multipliers, xi
    integer, allocatable, dimension(:,:)          :: mIndex, sIndex

!!MESH-Data-Inputs:
    integer    :: npg, neg, ntg, nbg, nbgc, nbgr, nParts
    integer    :: np, ne, nt, nb, nci, nri
    character(len=255) :: label, extension
    integer, allocatable, dimension(:,:)          :: e, eg
    integer, allocatable, dimension(:,:)          :: t, tg
    double precision, allocatable, dimension(:,:) :: p, pg

!!FEM & DDM-Data-Outputs:
    double precision, allocatable, dimension(:)   :: Ui, Ub, Ub_g, U_gi, U_g

!!PCKF-Data-Inputs/Outputs                                                       !!** PCKF ***
    integer :: nmeas, nmeasg
    integer, allocatable, dimension(:) :: measnodes  !!,measnodesg
    double precision, allocatable, dimension(:)   :: d
    double precision, allocatable, dimension(:,:) :: Af, Afb
    double precision, allocatable, dimension(:)   :: uveci, uvecb, var_g

!!Time Calculation Data:
    integer :: c1,c2,cr
    double precision :: time1
    double precision :: time2
    !integer :: seed_size
    !integer, dimension (:), allocatable :: seed

!!-----------------------------------------------------------------------------------------------
!! MPI initiation 
    CALL MPI_INIT(ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,pid,ierr)

    !! CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
    !! call random_seed (seed_size)
    !! allocate(seed(seed_size))
    !! do i=1, seed_size
    !!     seed(i) = pid*100
    !! end do
    !! call random_seed(put = seed)
    !! deallocate(seed)

!!-----------------------------------------------------------------------------------------------
!!!!!!!!!!> 1: Pre-Processing: Memory-Allocation / Assembly / Distribution of Amat
!!-----------------------------------------------------------------------------------------------
!!*** Pre-Processing Time T1 ***!!!
    call system_clock(count_rate=cr)
    call cpu_time(time1)
    call system_clock(c1)

    if (pid .eq. 0) then
        call readdbounds(dbounds)
    end if
    CALL MPI_BCAST(dbounds,4,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

!!-----------------------------------------------------------------------------------------------
!! PCE/KLE predefined parameters 
    amp   = 45.23d0
    sigma = 0.3d0       !!making sigma very small  0.0003d0
    gamma = 0.05d0**2   !!0.01d0**2
    const_diff = 0.1d0  !!0.009d0! 1e-1 for coarse mesh


!!-----------------------------------------------------------------------------------------------
!! PCE data
    if (pid .eq. 0) then
        !!print*, 'casep     =', casep
        call readGlobMeshDim(nParts)
        print*, '--------------------------------------'
        print*, 'nParts    =',nParts
        call readPceData(nord, ndim, npceout, npcein)
        print*, 'Order     =', nord
        print*, 'Dimen     =', ndim
        print*, 'nPCEin    =', npcein
        print*, 'nPCE      =', npceout
    end if

    CALL MPI_BCAST(nord,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(ndim,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(npceout,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(npcein,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)


!!-----------------------------------------------------------------------------------------------
!! PCE/KLE Indices
    allocate(omegas(nomga),multipliers(nomga),xi(ndim),mIndex(ndim,npceout),sIndex(2,ndim))
    if (pid .eq. 0) then
        call readIndices(ndim,npceout,mIndex,sIndex)
    end if

    CALL MPI_BCAST(mIndex,ndim*npceout,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(sIndex,2*ndim,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

!!-----------------------------------------------------------------------------------------------
!! Stochastic Parameters
    if (pid .eq. 0) then
        open(unit=1,file='../../data/meshData/params.dat')
        write(1,*) amp, const_diff, sigma, gamma
        close(1)
    end if

!!-----------------------------------------------------------------------------------------------
!! KLE expansion turms & Tripple Product (Cijk =  <Psi_i Psi_j Psi_k> )
    if (pid .eq. 0) then
        call readncijk(ncijk)
        print*, 'nCijk     =',ncijk
        call readprocess(nomga,omegas,multipliers)
    end if
    CALL MPI_BCAST(omegas,nomga,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(multipliers,nomga,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(ncijk,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

!!-----------------------------------------------------------------------------------------------
!! Read mesh data : local & global
    if (pid .eq. 0) then
        extension = ''
        call readmeshdim(npg,neg,ntg,nbg,extension)
        print*, 'nNodes    =',npg
        call readgcrdim(nbgc,nbgr)
        allocate(pg(2,npg),eg(3,neg),tg(3,ntg)) !!,measnodesg(nmeasg))
        call readmeshdata(pg,eg,tg,npg,neg,ntg,extension)
        call readnmeas(nmeasg,extension)                                          !!** PCKF ***
        !!call readmeasnodes(nmeasg,measnodesg,extension)
    end if
    CALL MPI_BCAST(nbg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(npg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(nbgr,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(nbgc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(nmeasg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)                    !!** PCKF ***

    call int2str(extension,pid+1,4)
    extension = trim(extension)
    call readmeshdim(np,ne,nt,nb,extension)
    call readcrdim(nci,nri,extension)              
    call readnmeas(nmeas,extension)                                               !!** PCKF ***

!!-----------------------------------------------------------------------------------------------
!! Dynamic array allocation
    allocate(p(2,np),e(3,ne),t(3,nt),Ub(nb*npceout),ijk(ncijk,3),cijk(ncijk),&
             Ui((np-nb)*npceout),Ub_g(nbg*npceout),U_g(npg),U_gi(npg))
    allocate(Af(np,npceout+nmeasg),uveci(np-nb),uvecb(nbg),d(nmeas), &
             Afb(nbg,npceout+nmeasg),var_g(npg),measnodes(nmeas))                 !!** PCKF ***

    !! allocate(Ur(nri*npceout),Uc(nci*npceout),UiUr((np-nci)*npceout),&
    !! FETIuveci(np-nci),FETIuvecb(nbgc),FETIAfb(nbgc,npceout+nmeasg),Uc_g(nbgc*npceout))

!!-----------------------------------------------------------------------------------------------
    if (pid .eq. 0) then
        call readcijk(ncijk,ijk,cijk)
    end if
    CALL MPI_BCAST(ijk,ncijk*3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(cijk,ncijk,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

    call readmeshdata(p,e,t,np,ne,nt,extension)
    call readmeasnodes(nmeas,measnodes,extension)                                 !!** PCKF ***

!!!-----------------------------------------------------------------------------------------------
!!!!!!!!!> 2: Calling Solver Only for Mallocs Calculations
!!!-----------------------------------------------------------------------------------------------
!! If you want to precalculate mallocs "mallocsCals = 1'      : First time run : Necessory
!! If have calculated mallocs then just use "mallocsCals = 0" : Repeated runs
!! mallocsCals = 1

!    if (mallocsCals .eq. 1) then
!
!    call PETSc_stonncpcgm(pid,p,e,t,np,ne,nt,nb,nbg,ndim,npcein,npceout,nri,nci,nbgc,nomga,ncijk,&
!                    ijk,cijk,mIndex,sIndex,casep,dbounds,const_diff,omegas,multipliers,sigma,amp,&
!                    Ui,Ub,Ub_g,mallocsCals,maxiter,tol)
!    end if

!!!-----------------------------------------------------------------------------------------------
!!*** Pre-Processing Time T2 ***!!!
    if (pid .eq. 0) then
        call cpu_time(time2)
        call system_clock(c2)
        if (pid .eq. 0) print*, '--------------------------------------------------------------'
        print*, 'Time taken for pre-processing', dble(c2-c1)/dble(cr),'seconds'
        print*, 'CPU-Time for pre-processing', time2-time1,'seconds'
        if (pid .eq. 0) print*, '--------------------------------------------------------------'
    end if

!!!-----------------------------------------------------------------------------------------------
!!*** Processing Time T1 ***!!! This includes assembly
    call cpu_time(time1)
    call system_clock(c1)

!!-----------------------------------------------------------------------------------------------
!!!!!!!!> 3: Main Process: Calling Solver
!!-----------------------------------------------------------------------------------------------

    call PETSc_stonncpcgm(pid,p,e,t,np,ne,nt,nb,nbg,ndim,npcein,npceout,nri,nci,nbgc,nomga,ncijk,&
                    ijk,cijk,mIndex,sIndex,casep,dbounds,const_diff,omegas,multipliers,sigma,amp,&
                    Ui,Ub,Ub_g,mallocsCals,maxiter,tol)


!!-----------------------------------------------------------------------------------------------
!!!-----------------------------------------------------------------------------------------------
!!*** Processing Time T2 ***!!!
    if (pid .eq. 0) then
        call cpu_time(time2)
        call system_clock(c2)
        if (pid .eq. 0) print*, '---------------------------------------------------------'
        print*, 'Time taken by solver', dble(c2-c1)/dble(cr),'seconds'
        print*, 'CPU-Time taken by olver',time2-time1,'seconds'
        if (pid .eq. 0) print*, '---------------------------------------------------------'
    end if


!!-----------------------------------------------------------------------------------------------
!!!!!!!!!!> 4: Post-Processing: Creating Final Solution/Writing Output Files
!!-----------------------------------------------------------------------------------------------
!!!*** Post-Processing Time T1 ***!!!
!   call cpu_time(time1)
!   call system_clock(c1)
!!-----------------------------------------------------------------------------------------------
!!! FETIDP-Post-Processing
!!-----------------------------------------------------------------------------------------------
!   IF (SolverID .eq. 6) THEN
!       call create_Ub_FETIddm(pid,nb,npceout,nri,nci,Ur,Uc,Ub)
!       call create_Af_FETIddm(nmeasg,npceout,np,nb,Ui,Ub,Af)
!       if (pid .eq. 0) then
!           call create_Afb_FETIddm(nmeasg,npceout,nbgc,Uc_g,FETIAfb)
!       end if
!       xi(:) = 0.0d0
!       xi(1) = -1.0d0
!       call create_d_pce_sp4(nmeas,np,ndim,extension,npceout,measnodes,Af,d,gamma,xi)

!!!***** Actual Solution *****!!!
        !call create_realization_FETIsp_i(np,nci,npceout,nmeas,Af,FETIuveci,xi)
        !!! call create_realization_sp_i4(np,nb,ndim,npceout,nmeas,Af,uveci,xi)
        !call construct_coln(pid,FETIuveci,U_gi)
        !call MPI_REDUCE(U_gi,U_g,npg,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        !
        !if (pid .eq. 0) then
        !   call create_realization_FETIsp_b(nbgc,npceout,nmeas,Afb,FETIuvecb,xi)
        !   !! call create_realization_sp_b4(nbg,ndim,npceout,nmeas,Afb,uvecb,xi)
        !   call construct_FETIcoln_boundary(FETIuvecb,U_g)
        !   label='out_actual.vtk'                     !!Error Here
        !   call writevtk(label,pg,tg,U_g,npg,ntg)
        !   open(unit=3,file='out_actual.dat',status='replace')
        !   write(unit=3,fmt=*) U_g
        !   close(3)
        !end 

!!!!***** Mean *****!!!
!       call construct_coln(pid,Af(1:(np-nci),1),U_gi)
!       CALL MPI_REDUCE(U_gi,U_g,npg,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!       if (pid .eq. 0) then
!           call construct_FETIcoln_boundary(FETIAfb(:,1),U_g)
!       label='out_mean_prior_pce.vtk'
!       call writevtk(label,pg,tg,U_g,npg,ntg)
!       print*, maxval(U_g)
!       end if
!       
!!!!***** SD *****!!!    
!       var_g(:) = 0.0d0
!       do i = 2,npceout
!           call construct_coln(pid,Af(1:(np-nci),i),U_gi)
!           call MPI_REDUCE(U_gi,U_g,npg,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!           if (pid .eq. 0) then
!           call construct_FETIcoln_boundary(Uc_g((i-1)*nbgc+1:i*nbgc),U_g)
!           var_g = var_g + U_g**2
!           end if
!       end do
!       if (pid .eq. 0) then
!           var_g = sqrt(var_g)
!           label='out_sd_prior_pce.vtk'
!           call writevtk(label,pg,tg,var_g,npg,ntg)
!           print*, maxval(var_g)
!       end if
!   
!!!!***** Higher Order PCE Coefficients *****!!!
!        do i = 2,npceout
!            call construct_coln(pid,Af(1:(np-nci),i),U_gi)
!            CALL MPI_REDUCE(U_gi,U_g,npg,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!            if (pid .eq. 0) then
!            call construct_FETIcoln_boundary(Uc_g((i-1)*nbgc+1:i*nbgc),U_g)
!            call int2str(extension,i,2)
!            label= 'out_pce_prior_' // trim(extension) // '.vtk'
!            call writevtk(label,pg,tg,U_g,npg,ntg)
!            end if
!        end do
!!-----------------------------------------------------------------------------------------------
!    ELSE


!-----------------------------------------------------------------------------------------------
!! PCGM-Post-Processing
!!-----------------------------------------------------------------------------------------------
        call create_Af_pckfddm(nmeasg,npceout,np,nb,Ui,Ub,Af)
        if (pid .eq. 0) then
          call create_Afb_pckfddm(nmeasg,npceout,nbg,Ub_g,Afb)
        end if

        xi(:) = 0.0d0
        xi(1) = -1.0d0

!!!***** Mean *****!!!
       call construct_coln(pid,Af(1:(np-nb),1),U_gi)
       CALL MPI_REDUCE(U_gi,U_g,npg,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
       if (pid .eq. 0) then
          call construct_coln_boundary(Afb(:,1),U_g)
          label='../../data/vtkOutputs/out_mean_prior_pce.vtk'
          call writevtk(label,pg,tg,U_g,npg,ntg)
          print*, maxval(U_g)
       end if

!!!***** SD *****!!!
       var_g(:) = 0.0d0
       do i = 2,npceout
          call construct_coln(pid,Ui((i-1)*(np-nb)+1:i*(np-nb)),U_gi)
          call MPI_REDUCE(U_gi,U_g,npg,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
          if (pid .eq. 0) then
             call construct_coln_boundary(Ub_g((i-1)*nbg+1:i*nbg),U_g)
             var_g = var_g + U_g**2
          end if
       end do
       if (pid .eq. 0) then
          var_g = sqrt(var_g)
          label='../../data/vtkOutputs/out_sd_prior_pce.vtk'
          call writevtk(label,pg,tg,var_g,npg,ntg)
          print*, maxval(var_g)
       end if




!!!***** Higher Order PCE Coefficients *****!!! 
!       do i = 2,npceout
!          call construct_coln(pid,Af(1:(np-nb),i),U_gi)
!          CALL MPI_REDUCE(U_gi,U_g,npg,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!          if (pid .eq. 0) then
!             call construct_coln_boundary(Afb(:,i),U_g)
!             call int2str(extension,i,2)
!             label= 'out_pce_prior_' // trim(extension) // '.vtk'
!             call writevtk(label,pg,tg,U_g,npg,ntg)
!          end if
!       end do

!CALL MPI_BCAST(ndim,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!!!***** Acrual Solution *****!!!
!!! 3-Dimensional RV
!if (ndim .eq. 3)
!    call create_d_pce_sp(nmeas,np,ndim,extension,npceout,measnodes,Af,d,gamma,xi)
!    call create_realization_sp_i(np,nb,npceout,nmeas,Af,uveci,xi)
!    call construct_coln(pid,uveci,U_gi)
!    call MPI_REDUCE(U_gi,U_g,npg,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!
!    if (pid .eq. 0) then
!        call create_realization_sp_b(nbg,npceout,nmeas,Afb,uvecb,xi)
!        call construct_coln_boundary(uvecb,U_g)
!        label='out_actual.vtk'
!        call writevtk(label,pg,tg,U_g,npg,ntg)
!    end if
!! 4-Dimensional RV
!else !! if (ndim .eq. 4)
!           call create_d_pce_sp4(nmeas,np,ndim,extension,npceout,measnodes,Af,d,gamma,xi)
!           call create_realization_sp_i4(np,nb,ndim,npceout,nmeas,Af,uveci,xi)
!           call construct_coln(pid,uveci,U_gi)
!           call MPI_REDUCE(U_gi,U_g,npg,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!
!           if (pid .eq. 0) then
!              call create_realization_sp_b4(nbg,ndim,npceout,nmeas,Afb,uvecb,xi)
!              call construct_coln_boundary(uvecb,U_g)
!              label='out_actual.vtk'
!              call writevtk(label,pg,tg,U_g,npg,ntg)
!           end if
!end if
!
!       call construct_coln(pid,Af(1:(np-nb),1),U_gi)
!       CALL MPI_REDUCE(U_gi,U_g,npg,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!       if (pid .eq. 0) then
!          uvecb(:) = -1.0d0
!          call construct_coln_boundary(uvecb,U_g)
!          label='out_mean_prior_pce_decom.vtk'
!          call writevtk(label,pg,tg,U_g,npg,ntg)
!       end if
!       var_g(:) = 0.0d0
!       do i = 2,npceout
!          call construct_coln(pid,Ui((i-1)*(np-nb)+1:i*(np-nb)),U_gi)
!          call MPI_REDUCE(U_gi,U_g,npg,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!          if (pid .eq. 0) then
!             var_g = var_g + U_g**2
!          end if
!       end do
!       if (pid .eq. 0) then
!          var_g = sqrt(var_g)
!          uvecb(:) = -1.0d0
!          call construct_coln_boundary(uvecb,var_g)
!          label='out_sd_prior_pce_decom.vtk'
!          call writevtk(label,pg,tg,var_g,npg,ntg)
!       end if
!
!   END IF
!
!!-----------------------------------------------------------------------------------------------
!!!!!!!!!!! PCKF Update
!!-----------------------------------------------------------------------------------------------
!!FromHere
!!call cpu_time(time1)
!!call system_clock(c1)
!!
!!   call create_Dfg_pckf(pid,nmeas,nmeasg,np,npceout,Af,RDf)
!!   call MPI_REDUCE(RDf,Df,nmeasg*(npceout+nmeasg),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!!
!!   call create_Dpg_pckf(pid,nmeas,nmeasg,np,npceout,d,RDp,gamma)
!    call MPI_REDUCE(RDp,Dp,nmeasg*(npceout+nmeasg),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!!
!!   if (pid .eq. 0) then
!!      call create_cmat_pckf(nmeasg,npceout,Df,Dp,Cmat)
!!   end if
!!   CALL MPI_BCAST(Cmat,(npceout+nmeasg)**2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!!
!!   call multiply(Af,Cmat,Aa,np,(npceout+nmeasg),(npceout+nmeasg))
!!   if (pid .eq. 0) then
!!      call multiply(Afb,Cmat,Aab,nbg,(npceout+nmeasg),(npceout+nmeasg))
!!   end if
!!
!!   if (pid .eq. 0) then
!!        call cpu_time(time2)
!!        call system_clock(c2)
!!        print*, 'Time taken for Update', dble(c2-c1)/dble(cr), 'seconds (', time2-time1, 'cpu time)'
!!        open(unit=7,file="results.txt", position="append", status='old')
!!        write (7,*) 'Time taken for Update', dble(c2-c1)/dble(cr), 'seconds (', time2-time1, 'cpu time)'
!!        close (7)
!!   end if

!!   do i = 2,npceout+nmeasg
!!      call construct_coln(pid,Aa(1:(np-nb),i),U_gi)
!!      CALL MPI_REDUCE(U_gi,U_g,npg,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!!      if (pid .eq. 0) then
!!         call construct_coln_boundary(Aab(:,i),U_g)
!!         call int2str(extension,i,2)
!!         label= 'out_pce_posterior_' // trim(extension) // '.vtk'
!!         call writevtk(label,pg,tg,U_g,npg,ntg)
!!      end if
!!   end do
!!
!!   call construct_coln(pid,Aa(1:(np-nb),1),U_gi)
!!   CALL MPI_REDUCE(U_gi,U_g,npg,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!!   if (pid .eq. 0) then
!!     call construct_coln_boundary(Aab(:,1),U_g)
!!   end if
!!   if (pid .eq. 0) then
!!     label='out_mean_posterior_pckf.vtk'
!!     call writevtk(label,pg,tg,U_g,npg,ntg)
!!   end if
!!
!!   var_g(:) = 0.0d0
!!  do i = 2,npceout+nmeasg
!!     call construct_coln(pid,Aa(1:(np-nb),i),U_gi)
!!     call MPI_REDUCE(U_gi,U_g,npg,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!!   if (pid .eq. 0) then
!!      call construct_coln_boundary(Aab(:,i),U_g)
!!      var_g = var_g + U_g**2
!!   end if
!!  end do
!!  if (pid .eq. 0) then
!!    var_g = sqrt(var_g)
!!    label='out_sd_posterior_pckf.vtk'
!!    call writevtk(label,pg,tg,var_g,npg,ntg)
!!    print*, maxval(var_g)
!!  end if
!
!!  do j = 1,nmeasg
!!      cross_g(:) = 0.0d0
!!      do i = 2,npceout
!!         call construct_coln(pid,Ui((i-1)*(np-nb)+1:i*(np-nb)),U_gi)
!!         call MPI_REDUCE(U_gi,U_g,npg,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!!         if (pid .eq. 0) then
!!            call construct_coln_boundary(Ub_g((i-1)*nbg+1:i*nbg),U_g)
!!            cross_g = cross_g + U_g*Df(j,i)
!!         end if
!!      end do
!!      if (pid .eq. 0) then
!!         call int2str(extension,j,2)
!!         label= 'out_cross_corr_' // trim(extension) // '.vtk'
!!         call writevtk(label,pg,tg,cross_g,npg,ntg)
!!      end if
!!  end do
!!
!!-----------------------------------------------------------------------------------------------
!!!*** Post-Processing Time T1 ***!!!
!!  if (pid .eq. 0) then
!!      call cpu_time(time2)
!!      call system_clock(c2)
!!      print*, 'Time taken by solver', dble(c2-c1)/dble(cr),'seconds'
!!      print*, 'CPU-Time taken by olver',time2-time1,'seconds'
!!  end if
!!
!!!!UntillHere
!
!!-----------------------------------------------------------------------------------------------
!! Dynamic array de-allocation 
   if (pid .eq. 0) then
      deallocate(pg,eg,tg)
   end if

    deallocate(p,e,t,Ub,Ui,Ub_g,U_g,U_gi,ijk,cijk)
    deallocate(omegas,multipliers,xi,mIndex,sIndex)
    deallocate(Af,Afb,measnodes,d,var_g,uveci,uvecb)                               !!** PCKF ***

   CALL MPI_FINALIZE(ierr)

   
END PROGRAM main
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%** END **%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!-----------------------------------------------------------------------------------------------
!!!!!!!!!!!! END !!!!!!!!!!
!!!-----------------------------------------------------------------------------------------------
!!!!!!!! Extra IMP Stuffs !!!!!!!
!!!!!-----------------------------------------------------------------------------------------------
!!!! To assemble vector : Method-4
!!    Fsi(:) = 0.0d0
!!    Fsg(:) = 0.0d0
!!    call SubAssembleVector(Fsi, Fsg, p, e, t, np, ne, nt, nb, npceout, amp, dbounds)
!!
!!!!-----------------------------------------------------------------------------------------------
!!!! To assemble matrix
!!!! Stochastic Matices
!!    Asii(:,:) = 0.0d0
!!    Asig(:,:) = 0.0d0
!!    Asgi(:,:) = 0.0d0
!!    Asgg(:,:) = 0.0d0
!!    Asir(:,:) = 0.0d0
!!    Asri(:,:) = 0.0d0
!!    Asic(:,:) = 0.0d0
!!    Asci(:,:) = 0.0d0
!!    Asrc(:,:) = 0.0d0
!!    Ascr(:,:) = 0.0d0
!!    Ascc(:,:) = 0.0d0
!!    Asrr(:,:) = 0.0d0
!!
!!!!-----------------------------------------------------------------------------------------------
!!!! Method : 3
!!    indexi = 1
!!    do k = 1,npcein
!!        !! Deterministic Matices
!!        Adii(:,:) = 0.0d0
!!        Adgg(:,:) = 0.0d0
!!        Adig(:,:) = 0.0d0
!!        Adgi(:,:) = 0.0d0
!!        Adir(:,:) = 0.0d0
!!        Adri(:,:) = 0.0d0
!!        Adic(:,:) = 0.0d0
!!        Adci(:,:) = 0.0d0
!!        Adrc(:,:) = 0.0d0
!!        Adcr(:,:) = 0.0d0
!!        Adcc(:,:) = 0.0d0
!!        Adrr(:,:) = 0.0d0
!!
!!        if (k .eq. 1) then
!!            !! Det-Advection Matrix
!!            call SubAssembleMatrix(Adii,Adgg,Adig,Adgi,Adir,Adri,Adic,Adci,Adrc,Adcr,Adcc,Adrr,p,e,t,np,ne,nt,nb,nci,ndim,&
!!                     1,0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],[0.0d0,0.0d0,0.0d0,0.0d0],0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],dbounds,0,0)
!!            !! Det-Diffusion Matrix
!!            call SubAssembleMatrix(Adii,Adgg,Adig,Adgi,Adir,Adri,Adic,Adci,Adrc,Adcr,Adcc,Adrr,p,e,t,np,ne,nt,nb,nci,ndim,2,&
!!                                   const_diff,omegas,multipliers,sigma,[0.0d0,0.0d0,0.0d0,0.0d0],dbounds,3,1)
!!        else
!!            !! Sto-Diffusion Matrices
!!            call SubAssembleMatrix(Adii,Adgg,Adig,Adgi,Adir,Adri,Adic,Adci,Adrc,Adcr,Adcc,Adrr,p,e,t,np,ne,nt,nb,nci,ndim,2,&
!!                                   const_diff,omegas,multipliers,sigma,[0.0d0,0.0d0,0.0d0,0.0d0],dbounds,3,k)
!!        end if
!!
!!        do i = 1,npceout
!!        do j = 1,npceout
!!
!            if ((k .eq. ijk(indexi,1)) .and. (i .eq. ijk(indexi,2)) .and. (j .eq. ijk(indexi,3))) then
!!                Asii(((i-1)*(np-nb)+1):(i*(np-nb)),((j-1)*(np-nb)+1):(j*(np-nb))) =&
!!                Asii(((i-1)*(np-nb)+1):(i*(np-nb)),((j-1)*(np-nb)+1):(j*(np-nb))) &
!!                + cijk(indexi)*Adii !!Aii !!Apce(1:(np-nb),1:(np-nb),k) !!Adif(1:(np-nb),1:(np-nb))
!!
!!                Asig(((i-1)*(np-nb)+1):(i*(np-nb)),((j-1)*nb+1):(j*nb)) = &
!!                Asig(((i-1)*(np-nb)+1):(i*(np-nb)),((j-1)*nb+1):(j*nb)) &
!!                + cijk(indexi)*Adig !!Aig !!Apce(1:(np-nb),(np-nb)+1:np,k) !!Adif(1:(np-nb),(np-nb)+1:np)
!!
!!                Asgi(((i-1)*nb+1):(i*nb),((j-1)*(np-nb)+1):(j*(np-nb))) = &
!!                Asgi(((i-1)*nb+1):(i*nb),((j-1)*(np-nb)+1):(j*(np-nb))) &
!!                + cijk(indexi)*Adgi !!Agi !!Apce((np-nb)+1:np,1:(np-nb),k) !!Adif((np-nb)+1:np,1:(np-nb))
!!
!!                Asgg(((i-1)*nb+1):(i*nb),((j-1)*nb+1):(j*nb)) =&
!!                Asgg(((i-1)*nb+1):(i*nb),((j-1)*nb+1):(j*nb)) &
!!                + cijk(indexi)*Adgg !!Agg !!Apce((np-nb)+1:np,(np-nb)+1:np,k)!!Adif((np-nb)+1:np,(np-nb)+1:np)
!!
!!                Asir(((i-1)*(np-nb)+1):(i*(np-nb)),((j-1)*nri+1):(j*nri)) = &
!                Asir(((i-1)*(np-nb)+1):(i*(np-nb)),((j-1)*nri+1):(j*nri)) &
!!                + cijk(indexi)*Adir ! Air !!Apce(1:(np-nb),(np-nb)+1:(np-nci),k)
!!
!!                Asri(((i-1)*nri+1):(i*nri),((j-1)*(np-nb)+1):(j*(np-nb))) = &
!!                Asri(((i-1)*nri+1):(i*nri),((j-1)*(np-nb)+1):(j*(np-nb))) &
!!                + cijk(indexi)*Adri !! Ari !!Apce((np-nb)+1:(np-nci),1:(np-nb),k)
!!
!!                Asic(((i-1)*(np-nb)+1):(i*(np-nb)),((j-1)*nci+1):(j*nci)) = &
!!                Asic(((i-1)*(np-nb)+1):(i*(np-nb)),((j-1)*nci+1):(j*nci)) &
!!                + cijk(indexi)*Adic  !!Aic !!Apce(1:(np-nb),(np-nci)+1:np,k)
!!
!!                Asci(((i-1)*nci+1):(i*nci),((j-1)*(np-nb)+1):(j*(np-nb))) = &
!!                Asci(((i-1)*nci+1):(i*nci),((j-1)*(np-nb)+1):(j*(np-nb))) &
!!                + cijk(indexi)*Adci !!Aci !!Apce((np-nci)+1:np,1:(np-nb),k)
!
!!                Asrc(((i-1)*nri+1):(i*nri),((j-1)*nci+1):(j*nci)) = &
!!                Asrc(((i-1)*nri+1):(i*nri),((j-1)*nci+1):(j*nci)) &
!!                + cijk(indexi)*Adrc !!Arc !!Apce((np-nb)+1:(np-nci),(np-nci)+1:np,k)
!!
!!                Ascr(((i-1)*nci+1):(i*nci),((j-1)*nri+1):(j*nri)) = &
!!                Ascr(((i-1)*nci+1):(i*nci),((j-1)*nri+1):(j*nri)) &
!!                + cijk(indexi)*Adcr !!Acr !!Apce((np-nci)+1:np,(np-nb)+1:(np-nci),k)
!!
!!                Asrr(((i-1)*nri+1):(i*nri),((j-1)*nri+1):(j*nri)) = &
!!                Asrr(((i-1)*nri+1):(i*nri),((j-1)*nri+1):(j*nri)) &
!!                + cijk(indexi)*Adrr !!Arr !!Apce((np-nb)+1:(np-nci),(np-nb)+1:(np-nci),k)
!!
!!                Ascc(((i-1)*nci+1):(i*nci),((j-1)*nci+1):(j*nci)) = &
!!                Ascc(((i-1)*nci+1):(i*nci),((j-1)*nci+1):(j*nci)) &
!!                + cijk(indexi)*Adcc !!Acc !!Apce((np-nci)+1:np,(np-nci)+1:np,k)
!
!!                indexi = indexi+1
!!            end if
!!        end do
!!        end do
!!    end do
!!
!!!   Amat(((i-1)*(np-nb)+1):(i*(np-nb)),((j-1)*(np-nb)+1):(j*(np-nb))) = Amat(((i-1)*(np-nb)+1):(i*(np-nb)), &
!!!   ((j-1)*(np-nb)+1):(j*(np-nb))) + cijk(indexi)*Apce(1:(np-nb),1:(np-nb),k) ! Aii
!!!Amat(((i-1)*(np-nb)+1):(i*(np-nb)),((np-nb)*npceout)+((j-1)*nb+1):((np-nb)*npceout)+(j*nb)) = &
!!!   Amat(((i-1)*(np-nb)+1):(i*(np-nb)),((np-nb)*npceout)+((j-1)*nb+1):((np-nb)*npceout)+(j*nb)) &
!!!   + cijk(indexi)*Apce(1:(np-nb),(np-nb)+1:np,k) ! Aigamma
!!!Amat(((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+(i*nb),((j-1)*(np-nb)+1):(j*(np-nb))) = &
!!!   Amat(((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+(i*nb),((j-1)*(np-nb)+1):(j*(np-nb))) &
!!!   + cijk(indexi)*Apce((np-nb)+1:np,1:(np-nb),k) ! Agammai
!!!Amat(((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+(i*nb),((np-nb)*npceout)+((j-1)*nb+1): &
!!!   ((np-nb)*npceout)+(j*nb)) = Amat(((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+(i*nb), &
!!!   ((np-nb)*npceout)+((j-1)*nb+1):((np-nb)*npceout)+(j*nb)) &
!!!   + cijk(indexi)*Apce((np-nb)+1:np,(np-nb)+1:np,k) ! Agg
!!
!!!!!! First Way: Matched in Matlab
!!
!!!Amat(((i-1)*(np-nb)+1):(i*(np-nb)),((np-nb)*npceout)+((j-1)*nb+1):((np-nb)*npceout)+(j*nb)) = &
!!!   Amat(((i-1)*(np-nb)+1):(i*(np-nb)),((np-nb)*npceout)+((j-1)*nb+1):((np-nb)*npceout)+(j*nb)) &
!!!   + cijk(indexi)*Apce(1:(np-nb),(np-nb)+1:np,k) ! Aigamma
!!!Amat(((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+(i*nb),((j-1)*(np-nb)+1):(j*(np-nb))) = &
!!!   Amat(((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+(i*nb),((j-1)*(np-nb)+1):(j*(np-nb))) &
!!!   + cijk(indexi)*Apce((np-nb)+1:np,1:(np-nb),k) ! Agammai
!!!Amat(((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+(i*nb),((np-nb)*npceout)+((j-1)*nb+1): &
!!!   ((np-nb)*npceout)+(j*nb)) = Amat(((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+(i*nb), &
!!!   ((np-nb)*npceout)+((j-1)*nb+1):((np-nb)*npceout)+(j*nb)) &
!!!   + cijk(indexi)*Apce((np-nb)+1:np,(np-nb)+1:np,k) ! Agg
!!Amat( ((i-1)*(np-nb)+1):(i*(np-nb)), ((np-nb)*npceout)+((j-1)*nb+1):((np-nb)*npceout)+((j-1)*nb+nri) ) = &
!!   Amat( ((i-1)*(np-nb)+1):(i*(np-nb)), ((np-nb)*npceout)+((j-1)*nb+1):((np-nb)*npceout)+((j-1)*nb+nri) ) &
!!   + cijk(indexi)*Apce(1:(np-nb),(np-nb)+1:(np-nci),k) ! Air
!!Amat( ((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+((i-1)*nb+nri) , ((j-1)*(np-nb)+1):(j*(np-nb)) ) = &
!!   Amat( ((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+((i-1)*nb+nri) , ((j-1)*(np-nb)+1):(j*(np-nb)) ) &
!!   + cijk(indexi)*Apce((np-nb)+1:(np-nci),1:(np-nb),k) ! Ari
!!Amat( ((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+((i-1)*nb+nri) , ((np-nb)*npceout)+((j-1)*nb+1): &
!!   ((np-nb)*npceout)+((j-1)*nb+nri) ) = Amat( ((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+((i-1)*nb+nri), &
!!   ((np-nb)*npceout)+((j-1)*nb+1):((np-nb)*npceout)+((j-1)*nb+nri) ) &
!!   + cijk(indexi)*Apce((np-nb)+1:(np-nci),(np-nb)+1:(np-nci),k) ! Arr    
!!Amat( ((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+((i-1)*nb+nri) , ((np-nb)*npceout)+((j-1)*nb+nri+1): &
!!   ((np-nb)*npceout)+((j-1)*nb+nb) ) = Amat( ((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+((i-1)*nb+nri), &
!!   ((np-nb)*npceout)+((j-1)*nb+nri+1):((np-nb)*npceout)+((j-1)*nb+nb) ) &
!!   + cijk(indexi)*Apce((np-nb)+1:(np-nci),(np-nci)+1:np,k) ! Arc
!!Amat( ((i-1)*(np-nb)+1):(i*(np-nb)) , ((np-nb)*npceout)+((j-1)*nb+nri+1):((np-nb)*npceout)+((j-1)*nb+nb) ) = &
!!   Amat( ((i-1)*(np-nb)+1):(i*(np-nb)) , ((np-nb)*npceout)+((j-1)*nb+nri+1):((np-nb)*npceout)+((j-1)*nb+nb) ) &
!!   + cijk(indexi)*Apce(1:(np-nb),(np-nci)+1:np,k) ! Aic                  
!!Amat( ((np-nb)*npceout)+((i-1)*nb+nri+1):((np-nb)*npceout)+((i-1)*nb+nb) , ((j-1)*(np-nb)+1):(j*(np-nb)) ) = &
!!   Amat( ((np-nb)*npceout)+((i-1)*nb+nri+1):((np-nb)*npceout)+((i-1)*nb+nb) ,((j-1)*(np-nb)+1):(j*(np-nb))) &
!!   + cijk(indexi)*Apce((np-nci)+1:np,1:(np-nb),k) ! Aci
!!Amat( ((np-nb)*npceout)+((i-1)*nb+nri+1):((np-nb)*npceout)+((i-1)*nb+nb) , ((np-nb)*npceout)+((j-1)*nb+nri+1): &
!!   ((np-nb)*npceout)+((j-1)*nb+nb) ) = Amat( ((np-nb)*npceout)+((i-1)*nb+nri+1):((np-nb)*npceout)+((i-1)*nb+nb), &
!!   ((np-nb)*npceout)+((j-1)*nb+nri+1):((np-nb)*npceout)+((j-1)*nb+nb) ) &
!!   + cijk(indexi)*Apce((np-nci)+1:np,(np-nci)+1:np,k) ! Acc                  
!!Amat( ((np-nb)*npceout)+((i-1)*nb+nri+1):((np-nb)*npceout)+((i-1)*nb+nb) , ((np-nb)*npceout)+((j-1)*nb+1): &
!!   ((np-nb)*npceout)+((j-1)*nb+nri)) = Amat( ((np-nb)*npceout)+((i-1)*nb+nri+1):((np-nb)*npceout)+((i-1)*nb+nb), &
!!   ((np-nb)*npceout)+((j-1)*nb+1):((np-nb)*npceout)+((j-1)*nb+nri)) &
!!   + cijk(indexi)*Apce((np-nci)+1:np,(np-nb)+1:(np-nci),k) ! Acr
!!
!!
!!Amat(((i-1)*(np-nb)+1):(i*(np-nb)),((np-nb)*npceout)+((j-1)*nb+1):((np-nb)*npceout)+((j-1)*nb+nri)) = &
!!Amat(((i-1)*(np-nb)+1):(i*(np-nb)),((np-nb)*npceout)+((j-1)*nb+1):((np-nb)*npceout)+((j-1)*nb+nri)) &
!!   + cijk(indexi)*Apce(1:(np-nb),(np-nb)+1:(np-nci),k) ! Air
!!Amat(((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+((i-1)*nb+nri),((j-1)*(np-nb)+1):(j*(np-nb))) = &
!!   Amat(((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+((i-1)*nb+nri),((j-1)*(np-nb)+1):(j*(np-nb))) &
!!   + cijk(indexi)*Apce((np-nb)+1:(np-nci),1:(np-nb),k) ! Ari
!!Amat(((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+((i-1)*nb+nri),((np-nb)*npceout)+((j-1)*nb+1): &
!!   ((np-nb)*npceout)+((j-1)*nb+nri)) = Amat(((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+((i-1)*nb+nri), &
!!   ((np-nb)*npceout)+((j-1)*nb+1):((np-nb)*npceout)+((j-1)*nb+nri)) &
!!   + cijk(indexi)*Apce((np-nb)+1:(np-nci),(np-nb)+1:(np-nci),k) ! Arr    
!!Amat(((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+((i-1)*nb+nri),((np-nb)*npceout)+((j-1)*nb+nri+1): &
!!   ((np-nb)*npceout)+((j-1)*nb+nb)) = Amat(((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+((i-1)*nb+nri), &
!!   ((np-nb)*npceout)+((j-1)*nb+nri+1):((np-nb)*npceout)+((j-1)*nb+nb)) &
!!   + cijk(indexi)*Apce((np-nb)+1:(np-nci),(np-nci)+1:np,k) ! Arc
!!Amat(((i-1)*(np-nb)+1):(i*(np-nb)),((np-nb)*npceout)+((j-1)*nb+nri+1):((np-nb)*npceout)+((j-1)*nb+nb)) = &
!!   Amat(((i-1)*(np-nb)+1):(i*(np-nb)),((np-nb)*npceout)+((j-1)*nb+nri+1):((np-nb)*npceout)+((j-1)*nb+nb)) &
!!   + cijk(indexi)*Apce(1:(np-nb),(np-nci)+1:np,k) ! Aic                  
!!Amat(((np-nb)*npceout)+((i-1)*nb+nri+1):((np-nb)*npceout)+((i-1)*nb+nb),((j-1)*(np-nb)+1):(j*(np-nb))) = &
!!   Amat(((np-nb)*npceout)+((i-1)*nb+nri+1):((np-nb)*npceout)+((i-1)*nb+nb),((j-1)*(np-nb)+1):(j*(np-nb))) &
!!   + cijk(indexi)*Apce((np-nci)+1:np,1:(np-nb),k) ! Aci
!!Amat(((np-nb)*npceout)+((i-1)*nb+nri+1):((np-nb)*npceout)+((i-1)*nb+nb),((np-nb)*npceout)+((j-1)*nb+nri+1): &
!!   ((np-nb)*npceout)+((j-1)*nb+nb)) = Amat(((np-nb)*npceout)+((i-1)*nb+nri+1):((np-nb)*npceout)+((i-1)*nb+nb), &
!!   ((np-nb)*npceout)+((j-1)*nb+nri+1):((np-nb)*npceout)+((j-1)*nb+nb)) &
!!   + cijk(indexi)*Apce((np-nci)+1:np,(np-nci)+1:np,k) ! Acc                  
!!Amat(((np-nb)*npceout)+((i-1)*nb+nri+1):((np-nb)*npceout)+((i-1)*nb+nb),((np-nb)*npceout)+((j-1)*nb+1): &
!!   ((np-nb)*npceout)+((j-1)*nb+nri)) = Amat(((np-nb)*npceout)+((i-1)*nb+nri+1):((np-nb)*npceout)+((i-1)*nb+nb), &
!!   ((np-nb)*npceout)+((j-1)*nb+1):((np-nb)*npceout)+((j-1)*nb+nri)) &
!!   + cijk(indexi)*Apce((np-nci)+1:np,(np-nb)+1:(np-nci),k) ! Acr
!
!!!!!! Second Way: Don't match in Matlab: It worked finally
!!Amat(((i-1)*(np-nb)+1):(i*(np-nb)),((np-nb)*npceout)+((j-1)*nri+1):((np-nb)*npceout)+(j*nri)) = &
!!   Amat(((i-1)*(np-nb)+1):(i*(np-nb)),((np-nb)*npceout)+((j-1)*nri+1):((np-nb)*npceout)+(j*nri)) &
!!   + cijk(indexi)*Apce(1:(np-nb),(np-nb)+1:(np-nci),k) ! Air
!!Amat(((i-1)*(np-nb)+1):(i*(np-nb)),((np-nci)*npceout)+((j-1)*nci+1):((np-nci)*npceout)+(j*nci)) = &
!!   Amat(((i-1)*(np-nb)+1):(i*(np-nb)),((np-nci)*npceout)+((j-1)*nci+1):((np-nci)*npceout)+(j*nci)) &
!!   + cijk(indexi)*Apce(1:(np-nb),(np-nci)+1:np,k) ! Aic
!!Amat(((np-nb)*npceout)+((i-1)*nri+1):((np-nb)*npceout)+(i*nri),((j-1)*(np-nb)+1):(j*(np-nb))) = &
!!   Amat(((np-nb)*npceout)+((i-1)*nri+1):((np-nb)*npceout)+(i*nri),((j-1)*(np-nb)+1):(j*(np-nb))) &
!!   + cijk(indexi)*Apce((np-nb)+1:(np-nci),1:(np-nb),k) ! Ari
!!Amat(((np-nci)*npceout)+((i-1)*nci+1):((np-nci)*npceout)+(i*nci),((j-1)*(np-nb)+1):(j*(np-nb))) = &
!!   Amat(((np-nci)*npceout)+((i-1)*nci+1):((np-nci)*npceout)+(i*nci),((j-1)*(np-nb)+1):(j*(np-nb))) &
!!   + cijk(indexi)*Apce((np-nci)+1:np,1:(np-nb),k) ! Aci
!!Amat(((np-nb)*npceout)+((i-1)*nri+1):((np-nb)*npceout)+(i*nri),((np-nb)*npceout)+((j-1)*nri+1): &
!!   ((np-nb)*npceout)+(j*nri)) = Amat(((np-nb)*npceout)+((i-1)*nri+1):((np-nb)*npceout)+(i*nri), &
!!   ((np-nb)*npceout)+((j-1)*nri+1):((np-nb)*npceout)+(j*nri)) &
!!   + cijk(indexi)*Apce((np-nb)+1:(np-nci),(np-nb)+1:(np-nci),k) ! Arr
!!Amat(((np-nci)*npceout)+((i-1)*nci+1):((np-nci)*npceout)+(i*nci),((np-nci)*npceout)+((j-1)*nci+1): &
!!   ((np-nci)*npceout)+(j*nci)) = Amat(((np-nci)*npceout)+((i-1)*nci+1):((np-nci)*npceout)+(i*nci), &
!!   ((np-nci)*npceout)+((j-1)*nci+1):((np-nci)*npceout)+(j*nci)) &
!!   + cijk(indexi)*Apce((np-nci)+1:np,(np-nci)+1:np,k) ! Acc
!!Amat(((np-nb)*npceout)+((i-1)*nri+1):((np-nb)*npceout)+(i*nri),((np-nci)*npceout)+((j-1)*nci+1): &
!!   ((np-nci)*npceout)+(j*nci)) = Amat(((np-nb)*npceout)+((i-1)*nri+1):((np-nb)*npceout)+(i*nri), &
!!   ((np-nci)*npceout)+((j-1)*nci+1):((np-nci)*npceout)+(j*nci)) &
!!   + cijk(indexi)*Apce((np-nb)+1:(np-nci),(np-nci)+1:np,k) ! Arc
!!Amat(((np-nci)*npceout)+((i-1)*nci+1):((np-nci)*npceout)+(i*nci),((np-nb)*npceout)+((j-1)*nri+1): &
!!   ((np-nb)*npceout)+(j*nri)) = Amat(((np-nci)*npceout)+((i-1)*nci+1):((np-nci)*npceout)+(i*nci), &
!!   ((np-nb)*npceout)+((j-1)*nri+1):((np-nb)*npceout)+(j*nri)) &
!!   + cijk(indexi)*Apce((np-nci)+1:np,(np-nb)+1:(np-nci),k) ! Acr
!
!
!!   Aii = Amat(1:(np-nb)*npceout,1:(np-nb)*npceout)
!!   Aig = Amat(1:(np-nb)*npceout,((np-nb)*npceout+1):np*npceout)
!!   Agi = Amat(((np-nb)*npceout+1):np*npceout,1:(np-nb)*npceout)
!!   Agg = Amat(((np-nb)*npceout+1):np*npceout,((np-nb)*npceout+1):np*npceout)
!!!!!*
!!!!Amat(((i-1)*(np-nb)+1):(i*(np-nb)),((np-nb)*npceout)+((j-1)*nri+1):((np-nb)*npceout)+(j*nri)) = Air
!!!!Amat(((i-1)*(np-nb)+1):(i*(np-nb)),((np-nci)*npceout)+((j-1)*nci+1):((np-nci)*npceout)+(j*nci)) = Aic
!!   Air = Amat(1:(np-nb)*npceout,((np-nb)*npceout+1):(np-nci)*npceout)
!!   Aic = Amat(1:(np-nb)*npceout,((np-nci)*npceout+1):np*npceout)
!!   Ari = Amat(((np-nb)*npceout+1):(np-nci)*npceout,1:(np-nb)*npceout)
!!   Aci = Amat(((np-nci)*npceout+1):np*npceout,1:(np-nb)*npceout)
!!   Arr = Amat(((np-nb)*npceout+1):(np-nci)*npceout,((np-nb)*npceout+1):(np-nci)*npceout)
!!   Acc = Amat(((np-nci)*npceout+1):np*npceout,((np-nci)*npceout+1):np*npceout)
!!   Arc = Amat(((np-nb)*npceout+1):(np-nci)*npceout,((np-nci)*npceout+1):np*npceout)
!!   Acr = Amat(((np-nci)*npceout+1):np*npceout,((np-nb)*npceout+1):(np-nci)*npceout)
!
!!!Amat(((i-1)*(np-nb)+1):(i*(np-nb)),((np-nb)*npceout)+((j-1)*nri+1):((np-nb)*npceout)+(j*nri)) = Air
!!!Amat(((i-1)*(np-nb)+1):(i*(np-nb)),((np-nci)*npceout)+((j-1)*nci+1):((np-nci)*npceout)+(j*nci)) = Aic
!!Air = Amat(1:(np-nb)*npceout,(((np-nb)*npceout)+1):((np-nb)*npceout)+nri)
!!Aic = Amat(1:(np-nb)*npceout,(((np-nb)*npceout)+nri+1):((np-nb)*npceout)+nb)
!!Ari = Amat(((np-nb)*npceout)+1:((np-nb)*npceout)+nri,1:(np-nb)*npceout)
!!Aci = Amat(((np-nb)*npceout)+nri+1:((np-nb)*npceout)+nb,1:(np-nb)*npceout)
!!Arr = Amat(((np-nb)*npceout)+1:((np-nb)*npceout)+nri,((np-nb)*npceout)+1:((np-nb)*npceout)+nri)
!!Acc = Amat(((np-nb)*npceout)+nri+1:((np-nb)*npceout)+nb,((np-nb)*npceout)+nri+1:((np-nb)*npceout)+nb)
!!Arc = Amat(((np-nb)*npceout)+1:((np-nb)*npceout)+nri,((np-nb)*npceout)+nri+1:((np-nb)*npceout)+nb)
!!Acr = Amat(((np-nb)*npceout)+nri+1:((np-nb)*npceout)+nb,((np-nb)*npceout)+1:((np-nb)*npceout)+nri)
!!  
!!!-----------------------------------------------------------------------------------------------
!!! To assemble matrix & vector
!!! initializing all blocks-matrix by assinging zeros
!Asii(:,:) = 0.0d0
!!Asig(:,:) = 0.0d0
!!Asgi(:,:) = 0.0d0
!!Asgg(:,:) = 0.0d0
!!Ascc(:,:) = 0.0d0
!!Asrr(:,:) = 0.0d0
!!Asrc(:,:) = 0.0d0
!!Ascr(:,:) = 0.0d0
!!Asic(:,:) = 0.0d0
!!Asir(:,:) = 0.0d0
!!Asri(:,:) = 0.0d0
!!Asci(:,:) = 0.0d0
!! Assemble subdomian level block-matrices
!! Advection Term
!
!!!indexi = 1
!!   call blockAssembleMatrix(Adii,p,e,t,np,ne,nt,nb,nci,ndim,1,0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],&
!!                           [0.0d0,0.0d0,0.0d0,0.0d0],0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],dbounds,0,0)
!!   call blockAssembleMatrix(Adii,p,e,t,np,ne,nt,nb,nci,ndim,2,const_diff,omegas,multipliers,&
!!                            sigma,[0.0d0,0.0d0,0.0d0,0.0d0],dbounds,3,1)
!!    Adiipce(:,:,1) = Adii
!!
!!do k = 2,npcein
!!   call blockAssembleMatrix(Adii,p,e,t,np,ne,nt,nb,nci,ndim,2,const_diff,omegas,multipliers,&
!!                            sigma,[0.0d0,0.0d0,0.0d0,0.0d0],dbounds,3,1)
!!    Adiipce(:,:,k) = Adii
!!end do
!!
!!    indexi = 1
!!    do i = 1,npceout
!!        do j = 1,npceout
!!            do k = 1,npcein
!!                if ((i .eq. ijk(indexi,1)) .and. (j .eq. ijk(indexi,2)) .and. (k .eq. ijk(indexi,3))) then
!!
!!                Asii(((i-1)*(np-nb)+1):(i*(np-nb)),((j-1)*(np-nb)+1):(j*(np-nb))) =&
!!                Asii(((i-1)*(np-nb)+1):(i*(np-nb)),((j-1)*(np-nb)+1):(j*(np-nb))) &
!!                + cijk(indexi)*Adiipce(1:(np-nb),1:(np-nb),k)
!!
!!                indexi = indexi+1
!!
!!            end if
!!        end do
!!    end do
!!end do

!!!!!!
   
   

