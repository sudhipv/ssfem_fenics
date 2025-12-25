!! Main Code: By Ajit Desai, 8/2015, Mohammad Khalil, 1/2012
!! purpose: is to define inputs, to construct SSFEM and call DDM solvers 
!
!input:
!  PCE:KLE-Data:
!      ndim  : random dimension of KLE/PCE expansion used for input represenations
!      nord  : order of basis polynomials used for PCE expansion 
!      npcein: number of PCE terms : npceout
!      Apce  : stochastic-stiffness matrix 
!      Fpce  : stochastic-source vector 
!      Cijk  : multidimentional moments : tripple product 
!      sigma : σ for underlying Guassian process
!      omegas: w's of KLE expansion terms : Rogers Book 
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
!  MPI-Data: 
!      pid   : Processor ID
!      nproc : number of processors 
!      ierr  : MPI errors 
!  FEM-Data: 
!      Amat  : FEM assembly matrix 
!      Aadv  : advection matrix 
!      Adif  : diffusion matrix 
!      fvec  : FEM assembly vector
!  DDM-Data:
!      Aii   : interior-interior sub-block of Amat
!      Aig   : interior-boundary sub-block of Amat
!      Air   : interior-remainig sub-block of Amat
!      Aic   : interior-corner sub-block of Amat
!      Fi Fg : interior & boundary sub-block of Fvec
!
! output:
!      Ui Ub Ub_g  : local & global interior & boundary nodes solution vector
!!
!!---------------------------------------------------------------------------------------
PROGRAM main

   use common
   use myCommon
   use assembly
   use PETScSolvers
   implicit none

   include 'mpif.h'

!!---------------------------------------------------
!! One-Level-PCGM with PETSc Sparse Iterative Solver
!!---------------------------------------------------

!!MPI-Data:
    integer :: ierr, pid

!! Exit-Criteria:
    integer, parameter :: nvec = 3  !!For 3D
    integer, parameter :: maxiter = 1000
    double precision, parameter :: tol = 1.0d-7

!!PCE-Data:
    !! Select casep
    !! 2: for automatic
    !! 3: for 3dim, 3ord
    !! 4: for 4dim, 3ord
    !! 5: for 5dim, 3ord
    !integer, parameter :: casep = 2
    !integer, parameter :: nomga = 12
    integer            :: nord, ndim, npceout, npcein

!!Mallocs-Data:
    !! To pre-calculated "mallocsCals = 1" : First time run
    !! To use calculated "mallocsCals = 0" : Repeated runs
    integer, parameter :: mallocsCals = 1

!!KLE-Data-Inputs:
    integer            :: ncijk
    !! double precision   :: amp, const_diff, sigma, gamma
    integer, allocatable, dimension(:,:)        :: ijk
    double precision, allocatable, dimension(:) :: cijk
    !! double precision, dimension(2,2)            :: dbounds
    !! double precision, allocatable, dimension(:) :: omegas, multipliers, xi
    !! integer, allocatable, dimension(:,:)        :: mIndex, sIndex

!!MESH-Data-Inputs:
    integer :: np, ne, nt, nh, nb, npg, neg, ntg, nhg, nbg, nParts
    character(len=255) :: label, extension
    integer, allocatable, dimension(:,:)          :: hg
    double precision, allocatable, dimension(:,:) :: pg
    !! integer, allocatable, dimension(:,:)          :: e, eg
    !! integer, allocatable, dimension(:,:)          :: t, tg

!!FEM & DDM-Data-Outputs:
    double precision, allocatable, dimension(:)   :: Ui, Ub, Ub_g, U_g, U_gi

!!PCKF-Data-Inputs/Outputs:                                                       !!** PCKF ***
    integer  :: i, nmeasg !!, nmeas
    double precision, allocatable, dimension(:,:) :: Af, Afb
    double precision, allocatable, dimension(:)   :: var_g
    !! integer, allocatable, dimension(:) :: measnodes
    !! double precision, allocatable, dimension(:)   :: d

!!Time Calculation Data:
    integer  :: c1,c2,cr
    double precision :: time1
    double precision :: time2

!!-----------------------------------------------------------------------------------------------
!! MPI initiation 
   CALL MPI_INIT(ierr)
   CALL MPI_COMM_RANK(MPI_COMM_WORLD,pid,ierr)

!!-----------------------------------------------------------------------------------------------
!!!!!!!!!!> 1: Pre-Processing: Memory-Allocation / Assembly / Distribution of Amat
!!-----------------------------------------------------------------------------------------------
!!*** Pre-Processing Time T1 ***!!!
    call system_clock(count_rate=cr)
    call cpu_time(time1)
    call system_clock(c1)

!    if (pid .eq. 0) then
!      call readdbounds(dbounds)                            !!**
!    end if
!    CALL MPI_BCAST(dbounds,4,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

!!!-----------------------------------------------------------------------------------------------
!!! PCE/KLE Data
!    amp   = 45.23d0
!    sigma = 0.3d0       !!making sigma very small  0.0003d0
!    gamma = 0.05d0**2   !!0.01d0**2
!    const_diff = 0.1d0  !!0.009d0! 1e-1 for coarse mesh (meanc: in assemlby)

!!-----------------------------------------------------------------------------------------------
!! PCE data
    if (pid .eq. 0) then
        call readGlobMeshDim(nParts)
        print*, '--------------------------------------'
        print*, 'nParts    =', nParts
        call readPceData(nord, ndim, npceout, npcein)
        print*, 'Order     =', nord
        print*, 'Dimen     =', ndim
        print*, 'nPCEin    =', npcein
        print*, 'nPCEout   =', npceout
    end if

    CALL MPI_BCAST(nord,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(ndim,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(npceout,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(npcein,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)


!!!-----------------------------------------------------------------------------------------------
!!! PCE/KLE Indices
!    !!allocate(mIndex(ndim,npceout),sIndex(2,ndim),omegas(nomga),multipliers(nomga),xi(ndim))
!
!    if (pid .eq. 0) then
!        call readIndices(ndim,npceout,mIndex,sIndex)        !!**
!    end if
!
!    CALL MPI_BCAST(mIndex,ndim*npceout,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!    CALL MPI_BCAST(sIndex,2*ndim,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

!!!-----------------------------------------------------------------------------------------------
!!! Stochastic Parameters
!    if (pid .eq. 0) then
!        open(unit=1,file='../../data/meshData/params.dat')
!        write(1,*) amp, const_diff, sigma, gamma
!        close(1)
!    end if

!!-----------------------------------------------------------------------------------------------
!! Tripple Product (Cijk =  <Psi_i Psi_j Psi_k> )
    if (pid .eq. 0) then
        call readncijk(ncijk)
        print*, 'nCijk     =',ncijk!!**
        !! call readprocess(nomga,omegas,multipliers)             !!**
    end if
    !CALL MPI_BCAST(omegas,nomga,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    !CALL MPI_BCAST(multipliers,nomga,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(ncijk,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

!!-----------------------------------------------------------------------------------------------
!! Read mesh data : local & global
    if (pid .eq. 0) then
        extension = ''
        !!call readmeshdim(npg,neg,ntg,nbg,extension)         !!**
        !!call readmeshdata(pg,eg,tg,npg,neg,ntg,extension)   !!**
        !! call readnmeas(nmeasg,extension)                   !!** PCKF ***
        nmeasg = 0
        call readmeshdim3D(npg,neg,ntg,nhg,nbg,extension)         !!**
        print*, 'nNodes    =',npg
        allocate(pg(3,npg),hg(4,nhg))
        call readmeshdata3D(pg,hg,npg,nhg,extension)
    end if
    CALL MPI_BCAST(nbg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(npg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(nmeasg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)                 !!** PCKF ***

    call int2str(extension,pid+1,4)
    extension = trim(extension)
    !!call readnmeas(nmeas,extension)                         !!**                   !!** PCKF ***
    !!call readmeshdim(np,ne,nt,nb,extension)                 !!**
    call readmeshdim3D(np,ne,nt,nh,nb,extension)          !!**

!!-----------------------------------------------------------------------------------------------
!! Dynamic array allocation
    allocate(ijk(ncijk,3),cijk(ncijk)) !!,p(2,np),e(3,ne),t(3,nt),
    allocate(Ui((np-nb)*nvec*npceout),Ub(nb*nvec*npceout),U_gi(npg*nvec),Ub_g(nbg*nvec*npceout),U_g(npg*nvec))
    allocate(var_g(npg*nvec),Af(np*nvec,npceout+nmeasg),Afb(nbg*nvec,npceout+nmeasg))
    !!allocate(d(nmeas),measnodes(nmeas),uvecb(nbg),uveci(np-nb))                    !!** PCKF ***

!!-----------------------------------------------------------------------------------------------
    if (pid .eq. 0) then
        call readcijk(ncijk,ijk,cijk)                       !!**
    end if
    CALL MPI_BCAST(ijk,ncijk*3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(cijk,ncijk,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

    !! call readmeshdata(p,e,t,np,ne,nt,extension)             !!**
    !! call readmeshdim3D(np,ne,nt,nh,nb,extension)
    !! call readmeasnodes(nmeas,measnodes,extension)           !!**                   !!** PCKF ***
!!!-----------------------------------------------------------------------------------------------
!!!!!!!!!> 2: Calling Solver Only for Mallocs Calculations
!!!-----------------------------------------------------------------------------------------------
!! If you want to precalculate mallocs "mallocsCals = 1'      : First time run : Necessory
!! If have calculated mallocs then just use "mallocsCals = 0" : Repeated runs
!! mallocsCals = 1

!    if (mallocsCals .eq. 1) then
!    !call PETSc_stolpcgm(pid, p, e, t, np, ne,nt,nb,nbg,ndim,npcein,npceout,nomga,ncijk, &
!    !            ijk,cijk,casep,mIndex,sIndex,dbounds,const_diff,omegas,multipliers, &
!    !            sigma,amp,Ui,Ub,Ub_g,mallocsCals,maxiter,tol)
!    call PETSc_stolpcgm(pid,np,nb,nbg,npcein,npceout,ncijk,ijk,cijk,Ui,Ub,Ub_g, &
!                        mallocsCals,maxiter,tol)
!    end if


!!-----------------------------------------------------------------------------------------------
!!*** Pre-Processing Time T2 ***!!!
    if (pid .eq. 0) then
        call cpu_time(time2)
        call system_clock(c2)
        if (pid .eq. 0) print*, '--------------------------------------------------------------'
        print*, 'Time taken for pre-processing', dble(c2-c1)/dble(cr),'seconds'
        print*, 'CPU-Time for pre-processing', time2-time1,'seconds'
    end if

!!-----------------------------------------------------------------------------------------------
!!*** Processing Time T1 ***!!!
    call cpu_time(time1)
    call system_clock(c1)

!!-----------------------------------------------------------------------------------------------
!!!!!!!!> 3: Main Process: Calling Solvers
!!-----------------------------------------------------------------------------------------------

!    call PETSc_stolpcgm(pid, p, e, t, np, ne,nt,nb,nbg,ndim,npcein,npceout,nomga,ncijk, &
!                    ijk,cijk,casep,mIndex,sIndex,dbounds,const_diff,omegas,multipliers, &
!                    sigma,amp,Ui,Ub,Ub_g,0,maxiter,tol)

    call PETSc_stolpcgm(pid,np,nb,nbg,npcein,npceout,ncijk,ijk,cijk,Ui,Ub,Ub_g, &
                    mallocsCals,maxiter,nvec,tol)

!!-----------------------------------------------------------------------------------------------
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
!!!*** Post-Processing Time T1 ***!!!
!   call cpu_time(time1)
!   call system_clock(c1)

!!-----------------------------------------------------------------------------------------------
!!!!!!!!!!! Post-Processing: Creating Final Solution/Writing Output Files
!!-----------------------------------------------------------------------------------------------
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

!!!***** Acrual Solution *****!!!
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
!!    ELSE
!!-----------------------------------------------------------------------------------------------
!!! PCGM-Post-Processing  
!!-----------------------------------------------------------------------------------------------
!        call create_Af_pckfddm(nmeasg,npceout,np,nb,Ui,Ub,Af)
!        if (pid .eq. 0) then
!          call create_Afb_pckfddm(nmeasg,npceout,nbg,Ub_g,Afb)
!        end if

!        xi(:) = 0.0d0
!        xi(1) = -1.0d0

!!!***** Mean *****!!!
!       call construct_coln(pid,Af(1:(np-nb),1),U_gi)
!       CALL MPI_REDUCE(U_gi,U_g,npg,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!       if (pid .eq. 0) then
!          call construct_coln_boundary(Afb(:,1),U_g)
!          label='../../data/vtkOutputs/out_mean_prior_pce.vtk'                             !!**
!          call writevtk3D(label,pg,hg,U_g,npg,nhg)
!          !print*, maxval(U_g)
!       end if
!
!        if (pid .eq. 0) then
!            open(unit=2,file='../../data/vtkOutputs/mean_solutionVector.dat',status='replace')
!            write(unit=2,fmt='(ES17.8E3)') U_g
!            close(2)
!        end if


    call create_Af_pckfddm3D(nmeasg,npceout,np,nb,Ui,Ub,nvec,Af)

    if (pid .eq. 0) then
        call create_Afb_pckfddm3D(nmeasg,npceout,nbg,Ub_g,nvec,Afb)
    end if

    !xi(:) = 0.0d0
    !xi(1) = -1.0d0

    !!!***** Mean *****!!!
    call construct_coln3D(pid,np,nb,npg,Af(1:(np-nb)*nvec,1),U_gi)
    !!call construct_coln3D(pid,np,nb,npg,Ui,U_gi)

    CALL MPI_REDUCE(U_gi,U_g,npg*nvec,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    if (pid .eq. 0) then
        call construct_coln_boundary3D(nbg,npg,Afb(:,1),U_g)
        !!label='../../data/vtkOutputs/out_mean_prior_pce.vtk'
        !!call writevtk(label,pg,tg,U_g,npg,ntg)
        !!print*, maxval(U_g)
    end if

    if (pid .eq. 0) then
        open(unit=2,file='../../data/vtkOutputs/mean_solutionVector.dat',status='replace')
        write(unit=2,fmt='(ES17.8E3)') U_g
        close(2)
    end if

    !!!***** SD *****!!!
    var_g(:) = 0.0d0
    do i = 2,npceout
        call construct_coln3D(pid,np,nb,npg,Ui((i-1)*(np-nb)*nVec+1:i*(np-nb)*nVec),U_gi)
        call MPI_REDUCE(U_gi,U_g,npg*nVec,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if (pid .eq. 0) then
            call construct_coln_boundary3D(nbg,npg,Ub_g((i-1)*nbg*nVec+1:i*nbg*nVec),U_g)

            call int2str(extension,i,2)
            label= '../../data/vtkOutputs/pce_' // trim(extension) // '.dat'
            open(unit=22,file=label,status='replace')
            write(unit=22,fmt='(ES17.8E3)') U_g
            close(22)

            var_g = var_g + U_g**2
        end if
    end do
    if (pid .eq. 0) then
        var_g = sqrt(var_g)
        !!label='../../data/vtkOutputs/out_sd_prior_pce.vtk'
        !!call writevtk(label,pg,tg,var_g,npg,ntg)
        print*, maxval(var_g)

        open(unit=2,file='../../data/vtkOutputs/sd_solutionVector.dat',status='replace')
        write(unit=2,fmt='(ES17.8E3)') var_g
        close(2)

    end if

!    !!***** Higher Order PCE Coefficients *****!!!
!    do i = 2,npceout
!      !call construct_coln3D(pid,np,nb,npg,Af(1:(np-nb)*nvec,i),U_gi)
!      !CALL MPI_REDUCE(U_gi,U_g,npg*nVec,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!      call construct_coln3D(pid,np,nb,npg,Ui((i-1)*(np-nb)*nVec+1:i*(np-nb)*nVec),U_gi)
!      call MPI_REDUCE(U_gi,U_g,npg*nVec,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!
!      if (pid .eq. 0) then
!         call construct_coln_boundary3D(nbg,npg,Afb(:,i),U_g)
!         call int2str(extension,i,2)
!
!         label= '../../data/vtkOutputs/pce_' // trim(extension) // '.dat'
!         !call writevtk(label,pg,tg,U_g,npg,ntg)
!
!        open(unit=2,file=label,status='replace')
!        write(unit=2,fmt='(ES17.8E3)') U_g
!        close(2)
!
!        end if
!    end do


!!!!***** SD *****!!!
!       var_g(:) = 0.0d0
!       do i = 2,npceout
!            call construct_coln(pid,Ui((i-1)*(np-nb)+1:i*(np-nb)),U_gi)
!            call MPI_REDUCE(U_gi,U_g,npg,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!            if (pid .eq. 0) then
!                 call construct_coln_boundary(Ub_g((i-1)*nbg+1:i*nbg),U_g)
!                 var_g = var_g + U_g**2
!            end if
!       end do
!       if (pid .eq. 0) then
!            var_g = sqrt(var_g)
!            !label='../../data/vtkOutputs/out_sd_prior_pce.vtk'                                !!**
!            !call writevtk(label,pg,hg,var_g,npg,nhg)
!            !print*, maxval(var_g)
!            open(unit=2,file='../../data/vtkOutputs/sd_solutionVector.dat',status='replace')
!            write(unit=2,fmt='(ES17.8E3)') var_g
!            close(2)
!       end if

!!!!!***** Higher Order PCE Coefficients *****!!!
!       do i = 2,npceout
!          call construct_coln(pid,Af(1:(np-nb),i),U_gi)
!          CALL MPI_REDUCE(U_gi,U_g,npg,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!          if (pid .eq. 0) then
!             call construct_coln_boundary(Afb(:,i),U_g)
!             call int2str(extension,i,2)
!             label= '../../data/vtkOutputs/out_pce_prior_' // trim(extension) // '.vtk'     !!**
!             call writevtk(label,pg,hg,U_g,npg,nhg)
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
!
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
!! END IF
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
!!        print*, 'Time taken by solver', dble(c2-c1)/dble(cr),'seconds'
!!        print*, 'CPU-Time taken by olver',time2-time1,'seconds'
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
!!!!UntillHere
!
!!-----------------------------------------------------------------------------------------------
!! Dynamic array de-allocation 
    if (pid .eq. 0) then
        deallocate(pg,hg)  !!eg,tg
    end if

    deallocate(Ub,Ui,Ub_g,U_g,U_gi,ijk,cijk)  !!p,e,t,
    !! deallocate(omegas,multipliers,xi,mIndex,sIndex)
    deallocate(Af,Afb,var_g) !!d,measnodes,uveci,uvecb)                          !!** PCKF ***


   CALL MPI_FINALIZE(ierr)


END PROGRAM main

!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%** END **%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!-----------------------------------------------------------------------------------------------
!!!!!!!!!!! END !!!!!!!!!!
!!-----------------------------------------------------------------------------------------------
!!!!!!! Extra IMP Stuffs !!!!!!!
!!-----------------------------------------------------------------------------------------------
!    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
!    integer, dimension (:), allocatable :: seed
!    integer :: seed_size
!    call random_seed (seed_size)
!    allocate(seed(seed_size))
!    do i=1, seed_size
!      seed(i) = pid*100
!    end do
!    call random_seed(put = seed)
!    deallocate(seed)
!
!
!!! To assemble vector :: Method - 1
!   call assemblevector(fvec, p, e, t, np, ne, nt, amp, dbounds)
!   Fpce(:) = 0.0d0
!   Fpce(1:(np-nb)) = fvec(1:(np-nb))
!   Fpce((np-nb)*npceout+1:(np-nb)*npceout+nb) = fvec((np-nb)+1:np)
!   Fi = Fpce(1:(np-nb)*npceout)
!   Fg = Fpce((((np-nb)*npceout)+1):np*npceout)
!   !!Fr = Fpce((np-nb)*npceout+1:(np-nci)*npceout)
!   !!Fc = Fpce((np-nci)*npceout+1:np*npceout)

!!! To assemble vector : Method-4
!Fsi(:) = 0.0d0
!Fsg(:) = 0.0d0
!call SubAssembleVector(Fsi, Fsg, p, e, t, np, ne, nt, nb, npceout, amp, dbounds)
!
!!!-----------------------------------------------------------------------------------------------
!!!-----------------------------------------------------------------------------------------------
!!! To assemble matrix
!!! Stochastic Matices
!Asii(:,:) = 0.0d0
!Asig(:,:) = 0.0d0
!Asgi(:,:) = 0.0d0
!Asgg(:,:) = 0.0d0
!
!
!!!-----------------------------------------------------------------------------------------------
!!! Method : 3
!indexi = 1
!do k = 1,npcein
!    !! Deterministic Matices
!    Adii = 0.0d0
!    Adgg = 0.0d0
!    Adig = 0.0d0
!    Adgi = 0.0d0
!
!    if (k .eq. 1) then
!        !! Det-Advection Matrix
!        call SubAssembleMatrix(Adii,Adgg,Adig,Adgi,p,e,t,np,ne,nt,nb,ndim,1,0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],&
!                           [0.0d0,0.0d0,0.0d0,0.0d0],0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],dbounds,0,0)
!        !! Det-Diffusion Matrix
!        call SubAssembleMatrix(Adii,Adgg,Adig,Adgi,p,e,t,np,ne,nt,nb,ndim,2,const_diff,omegas,multipliers,&
!                            sigma,[0.0d0,0.0d0,0.0d0,0.0d0],dbounds,3,1)
!    else
!        !! Sto-Diffusion Matrices
!        call SubAssembleMatrix(Adii,Adgg,Adig,Adgi,p,e,t,np,ne,nt,nb,ndim,2,const_diff,omegas,multipliers,&
!                            sigma,[0.0d0,0.0d0,0.0d0,0.0d0],dbounds,3,k)
!    end if
!
!    do i = 1,npceout
!    do j = 1,npceout
!
!        if ((k .eq. ijk(indexi,1)) .and. (i .eq. ijk(indexi,2)) .and. (j .eq. ijk(indexi,3))) then
!            Asii(((i-1)*(np-nb)+1):(i*(np-nb)),((j-1)*(np-nb)+1):(j*(np-nb))) =&                  !! ii, remove zero
!            Asii(((i-1)*(np-nb)+1):(i*(np-nb)),((j-1)*(np-nb)+1):(j*(np-nb))) &                   !! jj, remove zero
!            + cijk(indexi)*Adii !!Aii !!Apce(1:(np-nb),1:(np-nb),k) !!Adif(1:(np-nb),1:(np-nb))   !! Val
!
!            Asig(((i-1)*(np-nb)+1):(i*(np-nb)),((j-1)*nb+1):(j*nb)) = &
!            Asig(((i-1)*(np-nb)+1):(i*(np-nb)),((j-1)*nb+1):(j*nb)) &
!            + cijk(indexi)*Adig !!Aig !!Apce(1:(np-nb),(np-nb)+1:np,k) !!Adif(1:(np-nb),(np-nb)+1:np)
!
!            Asgi(((i-1)*nb+1):(i*nb),((j-1)*(np-nb)+1):(j*(np-nb))) = &
!            Asgi(((i-1)*nb+1):(i*nb),((j-1)*(np-nb)+1):(j*(np-nb))) &
!            + cijk(indexi)*Adgi !!Agi !!Apce((np-nb)+1:np,1:(np-nb),k) !!Adif((np-nb)+1:np,1:(np-nb))
!
!            Asgg(((i-1)*nb+1):(i*nb),((j-1)*nb+1):(j*nb)) =&
!            Asgg(((i-1)*nb+1):(i*nb),((j-1)*nb+1):(j*nb)) &
!            + cijk(indexi)*Adgg !!Agg !!Apce((np-nb)+1:np,(np-nb)+1:np,k)!!Adif((np-nb)+1:np,(np-nb)+1:np)
!
!            indexi = indexi+1
!        end if
!
!    end do
!    end do
!end do
!!!-----------------------------------------------------------------------------------------------
!
!!-----------------------------------------------------------------------------------------------
!! Method : 2
!indexi = 1
!do k = 1,npcein
!
!    if (k .eq. 1) then
!    call assemblematrix(Aadv,p,e,t,np,ne,nt,ndim,1,0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],&
!                      [0.0d0,0.0d0,0.0d0,0.0d0],0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],dbounds,0,0)
!    call assemblematrix(Adif,p,e,t,np,ne,nt,ndim,2,const_diff,omegas,multipliers,&
!                       sigma,[0.0d0,0.0d0,0.0d0,0.0d0],dbounds,3,1)
!    !Apce(:,:,1) = Aadv+Adif
!    Adif = Adif+Aadv
!
!    else
!    call assemblematrix(Adif,p,e,t,np,ne,nt,ndim,2,const_diff,omegas,multipliers,&
!                       sigma,[0.0d0,0.0d0,0.0d0,0.0d0],dbounds,3,k)
!    !Apce(:,:,k) = Adif
!    end if
!
!    do i = 1,npceout
!    do j = 1,npceout
!
!        if ((k .eq. ijk(indexi,1)) .and. (i .eq. ijk(indexi,2)) .and. (j .eq. ijk(indexi,3))) then
!            Asii(((i-1)*(np-nb)+1):(i*(np-nb)),((j-1)*(np-nb)+1):(j*(np-nb))) =&
!            Asii(((i-1)*(np-nb)+1):(i*(np-nb)),((j-1)*(np-nb)+1):(j*(np-nb))) &
!            + cijk(indexi)*Adif(1:(np-nb),1:(np-nb))    !! Aii !!Apce(1:(np-nb),1:(np-nb),k)
!
!            Asig(((i-1)*(np-nb)+1):(i*(np-nb)),((j-1)*nb+1):(j*nb)) = &
!            Asig(((i-1)*(np-nb)+1):(i*(np-nb)),((j-1)*nb+1):(j*nb)) &
!            + cijk(indexi)*Adif(1:(np-nb),(np-nb)+1:np)    !! Aig !!Apce(1:(np-nb),(np-nb)+1:np,k)
!
!            Asgi(((i-1)*nb+1):(i*nb),((j-1)*(np-nb)+1):(j*(np-nb))) = &
!            Asgi(((i-1)*nb+1):(i*nb),((j-1)*(np-nb)+1):(j*(np-nb))) &
!            + cijk(indexi)*Adif((np-nb)+1:np,1:(np-nb))    !! Agi !!Apce((np-nb)+1:np,1:(np-nb),k)
!
!            Asgg(((i-1)*nb+1):(i*nb),((j-1)*nb+1):(j*nb)) =&
!            Asgg(((i-1)*nb+1):(i*nb),((j-1)*nb+1):(j*nb)) &
!            + cijk(indexi)*Adif((np-nb)+1:np,(np-nb)+1:np)    !! Agg !!Apce((np-nb)+1:np,(np-nb)+1:np,k)
!
!            indexi = indexi+1
!        end if
!    end do
!    end do
!end do
!!!-----------------------------------------------------------------------------------------------
!!-----------------------------------------------------------------------------------------------
!! First level-Data distribution for StoDDM
!! Amat(:,:) = 0.0d0
!! Method : 1
!call assemblematrix(Aadv,p,e,t,np,ne,nt,ndim,1,0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],&
!                    [0.0d0,0.0d0,0.0d0,0.0d0],0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],dbounds,0,0)
!call assemblematrix(Adif,p,e,t,np,ne,nt,ndim,2,const_diff,omegas,multipliers,&
!                    sigma,[0.0d0,0.0d0,0.0d0,0.0d0],dbounds,3,1)
!Apce(:,:,1) = Aadv+Adif
!do i = 2,npcein
!    call assemblematrix(Adif,p,e,t,np,ne,nt,ndim,2,const_diff,omegas,multipliers,&
!                        sigma,[0.0d0,0.0d0,0.0d0,0.0d0],dbounds,3,i)
!    Apce(:,:,i) = Adif
!end do
!
!
!indexi = 1
!do i = 1,npceout
!    do j = 1,npceout
!    do k = 1,npcein
!    if ((i .eq. ijk(indexi,1)) .and. (j .eq. ijk(indexi,2)) .and. (k .eq. ijk(indexi,3))) then
!    !Amat(((i-1)*(np-nb)+1):(i*(np-nb)),((j-1)*(np-nb)+1):(j*(np-nb))) =&
!    !Amat(((i-1)*(np-nb)+1):(i*(np-nb)),((j-1)*(np-nb)+1):(j*(np-nb))) &
!    !+ cijk(indexi)*Apce(1:(np-nb),1:(np-nb),k) !! Aii
!
!    Asii(((i-1)*(np-nb)+1):(i*(np-nb)),((j-1)*(np-nb)+1):(j*(np-nb))) =&
!    Asii(((i-1)*(np-nb)+1):(i*(np-nb)),((j-1)*(np-nb)+1):(j*(np-nb))) &
!    + cijk(indexi)*Apce(1:(np-nb),1:(np-nb),k) !! Aii
!
!    !Amat(((i-1)*(np-nb)+1):(i*(np-nb)),((np-nb)*npceout)+((j-1)*nb+1):((np-nb)*npceout)+(j*nb)) = &
!    !Amat(((i-1)*(np-nb)+1):(i*(np-nb)),((np-nb)*npceout)+((j-1)*nb+1):((np-nb)*npceout)+(j*nb)) &
!    !+ cijk(indexi)*Apce(1:(np-nb),(np-nb)+1:np,k)    !! Aigamma
!
!    Asig(((i-1)*(np-nb)+1):(i*(np-nb)),((j-1)*nb+1):(j*nb)) = &
!    Asig(((i-1)*(np-nb)+1):(i*(np-nb)),((j-1)*nb+1):(j*nb)) &
!    + cijk(indexi)*Apce(1:(np-nb),(np-nb)+1:np,k)    !! Aigamma
!
!    !Amat(((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+(i*nb),((j-1)*(np-nb)+1):(j*(np-nb))) = &
!    !Amat(((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+(i*nb),((j-1)*(np-nb)+1):(j*(np-nb))) &
!    !+ cijk(indexi)*Apce((np-nb)+1:np,1:(np-nb),k)    !! Agammai
!
!    Asgi(((i-1)*nb+1):(i*nb),((j-1)*(np-nb)+1):(j*(np-nb))) = &
!    Asgi(((i-1)*nb+1):(i*nb),((j-1)*(np-nb)+1):(j*(np-nb))) &
!    + cijk(indexi)*Apce((np-nb)+1:np,1:(np-nb),k)    !! Agammai
!    !
!    !Amat(((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+(i*nb),((np-nb)*npceout)+((j-1)*nb+1): &
!    !((np-nb)*npceout)+(j*nb)) = Amat(((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+(i*nb), &
!    !((np-nb)*npceout)+((j-1)*nb+1):((np-nb)*npceout)+(j*nb)) &
!    !+ cijk(indexi)*Apce((np-nb)+1:np,(np-nb)+1:np,k) !! Agammagamma
!
!    Asgg(((i-1)*nb+1):(i*nb),((j-1)*nb+1):(j*nb)) =&
!    Asgg(((i-1)*nb+1):(i*nb),((j-1)*nb+1):(j*nb)) &
!    + cijk(indexi)*Apce((np-nb)+1:np,(np-nb)+1:np,k) !! Agammagamma
!
!    indexi = indexi+1
!    end if
!    end do
!    end do
!end do
!
!Aii = Asii !Amat(1:(np-nb)*npceout,1:(np-nb)*npceout)
!Aig = Asig !Amat(1:(np-nb)*npceout,((np-nb)*npceout+1):np*npceout)
!Agi = Asgi !Amat(((np-nb)*npceout+1):np*npceout,1:(np-nb)*npceout)
!Agg = Asgg !Amat(((np-nb)*npceout+1):np*npceout,((np-nb)*npceout+1):np*npceout)
!
!!!-----------------------------------------------------------------------------------------------
!!! Second level-Data distribution for StoDDM
!   Amat(:,:) = 0.0d0        !!*
!   indexi = 1
!   do i = 1,npceout
!    do j = 1,npceout
!     do k = 1,npcein
!      if ((i .eq. ijk(indexi,1)) .and. (j .eq. ijk(indexi,2)) .and. (k .eq. ijk(indexi,3))) then
!          Amat(((i-1)*(np-nb)+1):(i*(np-nb)),((j-1)*(np-nb)+1):(j*(np-nb))) =&
!             Amat(((i-1)*(np-nb)+1):(i*(np-nb)),((j-1)*(np-nb)+1):(j*(np-nb)))&
!             + cijk(indexi)*Apce(1:(np-nb),1:(np-nb),k) ! Aii
!          Amat(((i-1)*(np-nb)+1):(i*(np-nb)),((np-nb)*npceout)+((j-1)*nri+1):((np-nb)*npceout)+(j*nri)) = &
!             Amat(((i-1)*(np-nb)+1):(i*(np-nb)),((np-nb)*npceout)+((j-1)*nri+1):((np-nb)*npceout)+(j*nri)) &
!             + cijk(indexi)*Apce(1:(np-nb),(np-nb)+1:(np-nci),k) ! Air
!          Amat(((i-1)*(np-nb)+1):(i*(np-nb)),((np-nci)*npceout)+((j-1)*nci+1):((np-nci)*npceout)+(j*nci)) = &
!             Amat(((i-1)*(np-nb)+1):(i*(np-nb)),((np-nci)*npceout)+((j-1)*nci+1):((np-nci)*npceout)+(j*nci)) &
!             + cijk(indexi)*Apce(1:(np-nb),(np-nci)+1:np,k) ! Aic
!          Amat(((np-nb)*npceout)+((i-1)*nri+1):((np-nb)*npceout)+(i*nri),((j-1)*(np-nb)+1):(j*(np-nb))) = &
!             Amat(((np-nb)*npceout)+((i-1)*nri+1):((np-nb)*npceout)+(i*nri),((j-1)*(np-nb)+1):(j*(np-nb))) &
!             + cijk(indexi)*Apce((np-nb)+1:(np-nci),1:(np-nb),k) ! Ari
!          Amat(((np-nci)*npceout)+((i-1)*nci+1):((np-nci)*npceout)+(i*nci),((j-1)*(np-nb)+1):(j*(np-nb))) = &
!             Amat(((np-nci)*npceout)+((i-1)*nci+1):((np-nci)*npceout)+(i*nci),((j-1)*(np-nb)+1):(j*(np-nb))) &
!             + cijk(indexi)*Apce((np-nci)+1:np,1:(np-nb),k) ! Aci
!          Amat(((np-nb)*npceout)+((i-1)*nri+1):((np-nb)*npceout)+(i*nri),((np-nb)*npceout)+((j-1)*nri+1): &
!             ((np-nb)*npceout)+(j*nri)) = Amat(((np-nb)*npceout)+((i-1)*nri+1):((np-nb)*npceout)+(i*nri), &
!             ((np-nb)*npceout)+((j-1)*nri+1):((np-nb)*npceout)+(j*nri)) &
!             + cijk(indexi)*Apce((np-nb)+1:(np-nci),(np-nb)+1:(np-nci),k) ! Arr
!          Amat(((np-nci)*npceout)+((i-1)*nci+1):((np-nci)*npceout)+(i*nci),((np-nci)*npceout)+((j-1)*nci+1): &
!             ((np-nci)*npceout)+(j*nci)) = Amat(((np-nci)*npceout)+((i-1)*nci+1):((np-nci)*npceout)+(i*nci), &
!             ((np-nci)*npceout)+((j-1)*nci+1):((np-nci)*npceout)+(j*nci)) &
!             + cijk(indexi)*Apce((np-nci)+1:np,(np-nci)+1:np,k) ! Acc
!          Amat(((np-nb)*npceout)+((i-1)*nri+1):((np-nb)*npceout)+(i*nri),((np-nci)*npceout)+((j-1)*nci+1): &
!             ((np-nci)*npceout)+(j*nci)) = Amat(((np-nb)*npceout)+((i-1)*nri+1):((np-nb)*npceout)+(i*nri), &
!             ((np-nci)*npceout)+((j-1)*nci+1):((np-nci)*npceout)+(j*nci)) &
!             + cijk(indexi)*Apce((np-nb)+1:(np-nci),(np-nci)+1:np,k) ! Arc
!          Amat(((np-nci)*npceout)+((i-1)*nci+1):((np-nci)*npceout)+(i*nci),((np-nb)*npceout)+((j-1)*nri+1): &
!             ((np-nb)*npceout)+(j*nri)) = Amat(((np-nci)*npceout)+((i-1)*nci+1):((np-nci)*npceout)+(i*nci), &
!             ((np-nb)*npceout)+((j-1)*nri+1):((np-nb)*npceout)+(j*nri)) &
!             + cijk(indexi)*Apce((np-nci)+1:np,(np-nb)+1:(np-nci),k) ! Acr
!         indexi = indexi+1
!      end if
!     end do
!    end do
!   end do
!
!   Air = Amat(1:(np-nb)*npceout,((np-nb)*npceout+1):(np-nci)*npceout)
!   Aic = Amat(1:(np-nb)*npceout,((np-nci)*npceout+1):np*npceout)
!   Ari = Amat(((np-nb)*npceout+1):(np-nci)*npceout,1:(np-nb)*npceout)
!   Aci = Amat(((np-nci)*npceout+1):np*npceout,1:(np-nb)*npceout)
!   Arr = Amat(((np-nb)*npceout+1):(np-nci)*npceout,((np-nb)*npceout+1):(np-nci)*npceout)
!   Acc = Amat(((np-nci)*npceout+1):np*npceout,((np-nci)*npceout+1):np*npceout)
!   Arc = Amat(((np-nb)*npceout+1):(np-nci)*npceout,((np-nci)*npceout+1):np*npceout)
!   Acr = Amat(((np-nci)*npceout+1):np*npceout,((np-nb)*npceout+1):(np-nci)*npceout)                !!*
!!!!-----------------------------------------------------------------------------------------------
!
!!   Amat(((i-1)*(np-nb)+1):(i*(np-nb)),((j-1)*(np-nb)+1):(j*(np-nb))) = Amat(((i-1)*(np-nb)+1):(i*(np-nb)), &
!!   ((j-1)*(np-nb)+1):(j*(np-nb))) + cijk(indexi)*Apce(1:(np-nb),1:(np-nb),k) ! Aii
!!Amat(((i-1)*(np-nb)+1):(i*(np-nb)),((np-nb)*npceout)+((j-1)*nb+1):((np-nb)*npceout)+(j*nb)) = &
!!   Amat(((i-1)*(np-nb)+1):(i*(np-nb)),((np-nb)*npceout)+((j-1)*nb+1):((np-nb)*npceout)+(j*nb)) &
!!   + cijk(indexi)*Apce(1:(np-nb),(np-nb)+1:np,k) ! Aigamma
!!Amat(((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+(i*nb),((j-1)*(np-nb)+1):(j*(np-nb))) = &
!!   Amat(((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+(i*nb),((j-1)*(np-nb)+1):(j*(np-nb))) &
!!   + cijk(indexi)*Apce((np-nb)+1:np,1:(np-nb),k) ! Agammai
!!Amat(((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+(i*nb),((np-nb)*npceout)+((j-1)*nb+1): &
!!   ((np-nb)*npceout)+(j*nb)) = Amat(((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+(i*nb), &
!!   ((np-nb)*npceout)+((j-1)*nb+1):((np-nb)*npceout)+(j*nb)) &
!!   + cijk(indexi)*Apce((np-nb)+1:np,(np-nb)+1:np,k) ! Agg
!
!!!!! First Way: Matched in Matlab
!
!!Amat(((i-1)*(np-nb)+1):(i*(np-nb)),((np-nb)*npceout)+((j-1)*nb+1):((np-nb)*npceout)+(j*nb)) = &
!!   Amat(((i-1)*(np-nb)+1):(i*(np-nb)),((np-nb)*npceout)+((j-1)*nb+1):((np-nb)*npceout)+(j*nb)) &
!!   + cijk(indexi)*Apce(1:(np-nb),(np-nb)+1:np,k) ! Aigamma
!!Amat(((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+(i*nb),((j-1)*(np-nb)+1):(j*(np-nb))) = &
!!   Amat(((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+(i*nb),((j-1)*(np-nb)+1):(j*(np-nb))) &
!!   + cijk(indexi)*Apce((np-nb)+1:np,1:(np-nb),k) ! Agammai
!!Amat(((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+(i*nb),((np-nb)*npceout)+((j-1)*nb+1): &
!!   ((np-nb)*npceout)+(j*nb)) = Amat(((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+(i*nb), &
!!   ((np-nb)*npceout)+((j-1)*nb+1):((np-nb)*npceout)+(j*nb)) &
!!   + cijk(indexi)*Apce((np-nb)+1:np,(np-nb)+1:np,k) ! Agg
!Amat( ((i-1)*(np-nb)+1):(i*(np-nb)), ((np-nb)*npceout)+((j-1)*nb+1):((np-nb)*npceout)+((j-1)*nb+nri) ) = &
!   Amat( ((i-1)*(np-nb)+1):(i*(np-nb)), ((np-nb)*npceout)+((j-1)*nb+1):((np-nb)*npceout)+((j-1)*nb+nri) ) &
!   + cijk(indexi)*Apce(1:(np-nb),(np-nb)+1:(np-nci),k) ! Air
!Amat( ((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+((i-1)*nb+nri) , ((j-1)*(np-nb)+1):(j*(np-nb)) ) = &
!   Amat( ((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+((i-1)*nb+nri) , ((j-1)*(np-nb)+1):(j*(np-nb)) ) &
!   + cijk(indexi)*Apce((np-nb)+1:(np-nci),1:(np-nb),k) ! Ari
!Amat( ((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+((i-1)*nb+nri) , ((np-nb)*npceout)+((j-1)*nb+1): &
!   ((np-nb)*npceout)+((j-1)*nb+nri) ) = Amat( ((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+((i-1)*nb+nri), &
!   ((np-nb)*npceout)+((j-1)*nb+1):((np-nb)*npceout)+((j-1)*nb+nri) ) &
!   + cijk(indexi)*Apce((np-nb)+1:(np-nci),(np-nb)+1:(np-nci),k) ! Arr    
!Amat( ((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+((i-1)*nb+nri) , ((np-nb)*npceout)+((j-1)*nb+nri+1): &
!   ((np-nb)*npceout)+((j-1)*nb+nb) ) = Amat( ((np-nb)*npceout)+((i-1)*nb+1):((np-nb)*npceout)+((i-1)*nb+nri), &
!   ((np-nb)*npceout)+((j-1)*nb+nri+1):((np-nb)*npceout)+((j-1)*nb+nb) ) &
!   + cijk(indexi)*Apce((np-nb)+1:(np-nci),(np-nci)+1:np,k) ! Arc
!Amat( ((i-1)*(np-nb)+1):(i*(np-nb)) , ((np-nb)*npceout)+((j-1)*nb+nri+1):((np-nb)*npceout)+((j-1)*nb+nb) ) = &
!   Amat( ((i-1)*(np-nb)+1):(i*(np-nb)) , ((np-nb)*npceout)+((j-1)*nb+nri+1):((np-nb)*npceout)+((j-1)*nb+nb) ) &
!   + cijk(indexi)*Apce(1:(np-nb),(np-nci)+1:np,k) ! Aic                  
!Amat( ((np-nb)*npceout)+((i-1)*nb+nri+1):((np-nb)*npceout)+((i-1)*nb+nb) , ((j-1)*(np-nb)+1):(j*(np-nb)) ) = &
!   Amat( ((np-nb)*npceout)+((i-1)*nb+nri+1):((np-nb)*npceout)+((i-1)*nb+nb) ,((j-1)*(np-nb)+1):(j*(np-nb))) &
!   + cijk(indexi)*Apce((np-nci)+1:np,1:(np-nb),k) ! Aci
!Amat( ((np-nb)*npceout)+((i-1)*nb+nri+1):((np-nb)*npceout)+((i-1)*nb+nb) , ((np-nb)*npceout)+((j-1)*nb+nri+1): &
!   ((np-nb)*npceout)+((j-1)*nb+nb) ) = Amat( ((np-nb)*npceout)+((i-1)*nb+nri+1):((np-nb)*npceout)+((i-1)*nb+nb), &
!   ((np-nb)*npceout)+((j-1)*nb+nri+1):((np-nb)*npceout)+((j-1)*nb+nb) ) &
!   + cijk(indexi)*Apce((np-nci)+1:np,(np-nci)+1:np,k) ! Acc                  
!Amat( ((np-nb)*npceout)+((i-1)*nb+nri+1):((np-nb)*npceout)+((i-1)*nb+nb) , ((np-nb)*npceout)+((j-1)*nb+1): &
!   ((np-nb)*npceout)+((j-1)*nb+nri)) = Amat( ((np-nb)*npceout)+((i-1)*nb+nri+1):((np-nb)*npceout)+((i-1)*nb+nb), &
!   ((np-nb)*npceout)+((j-1)*nb+1):((np-nb)*npceout)+((j-1)*nb+nri)) &
!   + cijk(indexi)*Apce((np-nci)+1:np,(np-nb)+1:(np-nci),k) ! Acr
!
!
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
   
   

