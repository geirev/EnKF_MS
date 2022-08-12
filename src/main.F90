program main
   use mod_dimensions   ! Defines state dimension
   use mod_state        ! Defines model state
   use mod_observation  ! Defines observation state
   use m_set_random_seed2
   use m_model
   use m_pseudo1D
   use m_fixsample1D
   use m_dumpsol
   use m_ensemblemean
   use m_ensemblevariance
   use m_enkf
   use m_measurements
   implicit none

! state variables
   type(state), allocatable :: mem(:)
   type(state), allocatable :: old(:)
   type(state), allocatable :: sysnoise(:)
   type(state) ana  ! analytical solution
   type(state) ave  ! ensemble average
   type(state) var  ! ensemble variance
   real, allocatable :: samples(:,:)

   type(observation), allocatable :: obs(:)

! Variables read from infile
   integer nrt                           ! Number of timesteps
   type(substate) u
   type(substate) rh
   type(substate) inivar
   type(substate) sysvar
   type(substate) obsvar


   real, parameter :: dx=1.0             ! horizontal grid spacing
   real, parameter :: dt=1.0             ! Time step of atmospheric model
   integer nrens                         ! ensemble size
   integer mode_analysis                 ! 1 standard, 2 fixed R
   logical samp_fix

   integer nro                    ! Number of ocean measurement per assimilation time
   integer nra                    ! Number of atmos measurement per assimilation time
   integer nrobs                         ! Total number of measurement per assimilation time
   logical mkobs                         ! Create or read measurements
   logical Rexact                        ! Use exact(true) or lowrank(false) R matrix
   real obsdt                            ! time between assimilation times
   logical :: lrandrot=.true.            ! random rotation in SQRT schemes
   logical :: lsymsqrt=.true.            ! Always use the symmetrical square root rather than one-sided
   real deltaobs                         ! distance between observations

! inflation
   integer inflate                       ! 0--no inflation, 1--constant inflation, 2--adaptive inflation
   real infmult                          ! constant inflation or adjustment of adaptive inflation

! local analysis
   integer local                         ! 0-no localization, 1-distance based, 2-adaptive
   real obs_radius                       ! Number of grid cells for including measurements (distance based)
   real obs_truncation                   ! Correlations for truncating measurements in adaptive scheme

! parameters for analysis scheme
   real truncation                       ! Truncation of singular values
   character(len=8) covmodel             ! Diagonal or Gaussian measurement error covariance model
   real rd                               ! Horizontal correlation of observation errors in Gaussian case

! other variables
   integer i,j,m,k
   integer iobs                          ! Counting assimilation steps (counter for records in obs.uf)
   real time
   logical leuler

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! reading input data
   open(10,file='infile.in')
      read(10,*)nrens              ; print *,'nrens=       ',nrens
      read(10,*)nrt                ; print *,'nrt=         ',nrt
      read(10,*)u%ocean            ; print *,'u%ocean=     ',u%ocean
      read(10,*)rh%ocean           ; print *,'rh%ocean=    ',rh%ocean
      read(10,*)rh%atmos           ; print *,'rh%atmos=    ',rh%atmos
      read(10,'(1x,l1)')samp_fix   ; print *,'samp_fix=    ',samp_fix
      read(10,*)inivar%ocean       ; print *,'inivar%ocean=',inivar%ocean
      read(10,*)inivar%atmos       ; print *,'inivar%atmos=',inivar%atmos
      read(10,*)sysvar%ocean       ; print *,'sysvar%ocean=',sysvar%ocean
      read(10,*)sysvar%atmos       ; print *,'sysvar%atmos=',sysvar%atmos
      read(10,*)obsvar%ocean       ; print *,'obsvar%ocean=',obsvar%ocean
      read(10,*)obsvar%atmos       ; print *,'obsvar%atmos=',obsvar%atmos
      read(10,*)nro         ; print *,'nro=  ',nro
      read(10,*)nra         ; print *,'nra=  ',nra
      read(10,*)obsdt              ; print *,'obsdt=       ',obsdt
      read(10,'(1x,l1)')mkobs      ; print *,'mkobs=       ',mkobs
      read(10,*)mode_analysis      ; print *,'mode_ana=    ',mode_analysis
      read(10,*)truncation        ; print *,'truncation=  ',truncation
      read(10,'(1x,a)')covmodel   ; print *,'covmodel=    ',trim(covmodel)
      read(10,*)rd                ; print *,'rd      =    ',rd
      read(10,'(1x,l1)')Rexact    ; print *,'Rexact=      ',Rexact
      read(10,'(1x,l1)')lrandrot  ; print *,'lrandrot=    ',lrandrot
      read(10,*)inflate,infmult   ; print *,'inflation=   ',inflate,infmult
      read(10,*)local,obs_radius,obs_truncation; print *,'localization=',local,obs_radius,obs_truncation
   close(10)

   u%atmos=1.0

   call set_random_seed2

   call system('rm -f eigenvalues.dat')

   allocate (mem(nrens))
   allocate (old(0:nrens))
   allocate (sysnoise(nrens))
   allocate (samples(nx,nrens))

   nrobs=nro+nra
   allocate (obs(nrobs))

   iobs=0
   time=0.0



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Allocate and define measurement setup

! Uniformly distributed measurements for ocean
   deltaobs=nint(real(nx)/real(nro))
   obs(1)%pos=nint(deltaobs/2.0)
   do m=2,nro
      obs(m)%pos=min(nint(obs(1)%pos+real(m-1)*deltaobs) , nx)
   enddo

! Uniformly distributed measurements for atmosphere
   deltaobs=nint(real(nx)/real(nra))
   i=nro
   obs(i+1)%pos=nint(deltaobs/2.0)
   do m=2,nra
      obs(i+m)%pos=min(nint(obs(i+1)%pos+real(m-1)*deltaobs) , nx)
   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The true solution smooth pseudo random field drawn from  N(0,1,rh).
   call pseudo1D(ana%ocean,nx,1,rh%ocean,dx,nx)
   call pseudo1D(ana%atmos,nx,1,rh%atmos,dx,nx)
   print *,'main: ana ok'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First guess solution stored in ave
   call pseudo1D(ave%ocean,nx,1,rh%ocean,dx,nx)
   call pseudo1D(ave%atmos,nx,1,rh%atmos,dx,nx)
   ave=(ave + ana)*(1.0/sqrt(2.0))
   print *,'main: fg ok'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialization of ensemble
   call pseudo1D(samples,nx,nrens,rh%ocean,dx,nx)
   if (samp_fix) call fixsample1D(samples,nx,nrens)
   do j=1,nrens
      mem(j)%ocean(1:nx)=samples(1:nx,j)
   enddo

   call pseudo1D(samples,nx,nrens,rh%atmos,dx,nx)
   if (samp_fix) call fixsample1D(samples,nx,nrens)
   do j=1,nrens
      mem(j)%atmos(1:nx)=samples(1:nx,j)
   enddo

   do j=1,nrens
      mem(j)%ocean=ave%ocean + sqrt(inivar%ocean)*mem(j)%ocean
      mem(j)%atmos=ave%atmos + sqrt(inivar%atmos)*mem(j)%atmos
   enddo

   print *,'main: ensemble ok'

   call ensemblemean(mem,ave,nrens)
   call ensemblevariance(mem,ave,var,nrens)
   call dumpsol(time,ana,ave,var,nx,dx,obs,nro,nra,mem,nrens,'I')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Time stepping
   print *,'main: start time stepping'
   do k=1,nrt
      time=time+dt
      print '(a,i6,f10.2)','timestep ',k,time

! Advection
      if (k==1 .or. mod(time-dt,obsdt).lt.0.05*dt) then
         leuler=.true.
      else
         leuler=.false.
      endif
      print *,'euler:',leuler
      do j=1,nrens
         call model(mem(j),old(j),u,dx,dt,leuler)

      enddo
      call model(ana,old(0),u,dx,dt,leuler)

! System noise
      if (sysvar%ocean > 0.0 .or. sysvar%atmos > 0.0) then
         call pseudo1D(samples,nx,nrens,rh%ocean,dx,nx)
         if (samp_fix) call fixsample1D(samples,nx,nrens)
         do j=1,nrens
            mem(j)%ocean=mem(j)%ocean+sqrt(2.0*sysvar%ocean*dt)*samples(:,j)
         enddo

         call pseudo1D(samples,nx,nrens,rh%atmos,dx,nx)
         if (samp_fix) call fixsample1D(samples,nx,nrens)
         do j=1,nrens
            mem(j)%atmos=mem(j)%atmos+sqrt(2.0*sysvar%atmos*dt)*samples(:,j)
         enddo
      endif

! Assimilation step
      if (mod(time,obsdt).lt.0.05*dt) then
         iobs=iobs+1


         call measurements(ana,obs,obsvar,nro,nra,iobs,mkobs,time)

         call ensemblemean(mem,ave,nrens)
         call ensemblevariance(mem,ave,var,nrens)
         call dumpsol(time,ana,ave,var,nx,dx,obs,nro,nra,mem,nrens,'F')

         call enkf(mem,nrens,obs,nro,nra,mode_analysis,&
                  &truncation,covmodel,dx,rh,Rexact,rd,lrandrot,lsymsqrt,&
                  &inflate,infmult,local,obs_radius,obs_truncation)

         call ensemblemean(mem,ave,nrens)
         call ensemblevariance(mem,ave,var,nrens)
         call dumpsol(time,ana,ave,var,nx,dx,obs,nro,nra,mem,nrens,'A')

      endif

   enddo


end program main

