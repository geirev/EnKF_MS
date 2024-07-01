program main
   use mod_dimensions                      ! Defines state dimension
   use mod_state                           ! Defines model state
   use mod_observation                     ! Defines observation state
   use m_readinfile                        ! Reading infile.in
   use m_set_random_seed2                  ! Sets new random seed unless 'seed.dat' exists
   use m_model                             ! Model time stepping
   use m_frobenius                         ! Frobenius norm between two matrices
   use m_pseudo1D                          ! Samples pseudo-random fields in 1-D
   use m_fixsample1D                       ! Ensures zero mean and unit variance of sampled ensembles
   use m_obsxloc                           ! Computes x position of measurements
   use m_obstloc                           ! Computes time of measurements
   use m_obscount                          ! Counts the number of measurements in a DA window
   use m_obspoints                         ! Saves the measurement locations in space and time for plotting
   use m_ensemblemean                      ! Compute ensemble mean at current time (from mem)
   use m_ensemblevariance                  ! Compute ensemble variance at current time (from mem)
   use m_ensemblecovariance                ! Compute ensemble covariance at current time (from mem)
   use m_enkfprep                          ! Setting up predicted measurements and matrices for analysis
   use m_assimilation_old                  ! Old analysis assimilation step gathered in one routine
   use m_assimilation                      ! New analysis assimilation step gathered in one routine
   use m_prepY                             ! Setting up predicted measurements for IES analysis
   use m_scaling                           ! Setting up predicted measurements for IES analysis
   use m_prepD                             ! Setting up predicted measurements for IES analysis
   use m_pertE                             ! Setting up predicted measurements for IES analysis
   use m_random                            ! Generate random normal numbers
   use m_tecfld                            ! Tecplot output (not used)
   use mod_shapiro                         ! Shapiro filter in case it is needed
   use m_windowstat                        ! Compute ensemble mean and variance in a DA window (from win)
   use m_covstat                           ! Computes space time covariances in case you have a lot of memory
   use m_gnuplot                           ! Generate space time plots for plotting in gnuplot (c.gnu)
   use m_ies                               ! The ies analysis routine
   use m_print_ies_status                  ! Prints status of ies convergence
   use m_ies_steplength                    ! Computes the new step lengths for ies
   use, intrinsic :: omp_lib
   use m_ansi_colors
   implicit none

   type(state), allocatable :: win(:,:)    ! The ensemble of realizations over an assimilation window.
   type(state), allocatable :: win0(:,:)   ! Initial ensemble of realizations over an assimilation window.
   type(state), allocatable :: mem(:)      ! The ensemble of realizations
   type(state), allocatable :: sysnoise(:) ! The ensemble of sampled system noise updated every timestep
   type(state) ref                         ! reference solution
   type(state) ave                         ! ensemble average

   real, allocatable :: samples(:,:)       ! work array used when sampling in pseudo1D

! Spacew time statistics diagnostic variables
   type(state), allocatable :: refout(:)   ! ensemble average as a function of space and time
   type(state), allocatable :: mean(:)     ! ensemble average as a function of space and time
   type(state), allocatable :: stdt(:)     ! ensemble std dev as a function of space and time
   type(state), allocatable :: rest(:)     ! ensemble std dev as a function of space and time
   type(state), allocatable :: covo(:)     ! ensemble std dev as a function of space and time
   type(state), allocatable :: cova(:)     ! ensemble std dev as a function of space and time
   type(substate), allocatable :: rmse(:)  ! time series of rmse values (mean - referece)
   type(substate), allocatable :: rmss(:)  ! time series of rms std.dev values
   type(substate) rmset(1)                  ! Total rmse value (mean - referece)
   type(substate) rmsst(1)                  ! Total rmss value rms of std.dev.

! Observation location variables
   integer, allocatable :: obsoloc(:)      ! location of ocean observations in space
   integer, allocatable :: obsaloc(:)      ! location of atmos observations in space
   integer, allocatable :: obsotimes(:)    ! location of ocean observations in time
   integer, allocatable :: obsatimes(:)    ! location of atmos observations in time

   type(observation), allocatable :: obs(:)! Stores all observation information in a DA window

! EnKF analysis variables
   real, allocatable :: Y(:,:)             ! Predicted meaurements
   real, allocatable :: S(:,:)             ! Predicted meaurements anomalies
   real, allocatable :: E(:,:)             ! Measurement perturbations
   real, allocatable :: D(:,:)             ! Perturbed measurements/Innovations
   real, allocatable :: DA(:,:)            ! Measurements
   real, allocatable :: W(:,:)             ! ies iteration matrix
   real, allocatable :: Wold(:,:)          ! ies iteration matri for reducing steplength
   real, allocatable :: Yold(:,:)          ! ies predicted measurements for reducing steplength
   real, allocatable :: X(:,:)            ! The X matrix :-)


! other variables
   integer tini                            ! start time of an assimilation window
   integer tfin                            ! end   time of an assimilation window
   integer nrobs                           ! Total number of measurement per assimilation window

   integer j,k,l,iter
   real time
   real :: dxsamp=1.0
   real :: steplength
   real, allocatable :: costf(:,:)

   call readinfile()

   call set_random_seed2

   call system('rm -f eigenvalues.dat')
   call system('rm -f obsloc?.dat')
!   call system('touch obsloca.dat obsloco.dat')

! Now that we know all dimensions, allocate the main arrays
   allocate (win(0:nrt,nrens))
   allocate (win0(0:nrt,nrens))
   allocate (refout(0:nrt))
   allocate (mean(0:nrt))
   allocate (stdt(0:nrt))
   allocate (rest(0:nrt))
   allocate (covo(0:nrt))
   allocate (cova(0:nrt))
   allocate (mem(nrens))
   allocate (sysnoise(nrens))
   allocate (samples(nx,nrens))
   allocate (rmse(0:nrt))
   allocate (rmss(0:nrt))

   allocate (W(nrens,nrens))
   allocate (X(nrens,nrens))
   allocate (Wold(nrens,nrens))
   allocate (Yold(nrens,nrens))

   allocate (obsoloc(nro))
   allocate (obsaloc(nra))
   allocate (obsotimes(nrt))
   allocate (obsatimes(nrt))

   allocate(costf(nmda,nrwindows))
   costf=0.0
   time=0.0

!   call omp_set_num_threads(10)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Uniformly distributed measurements for ocean and atmosphere in space and time
   print *,'Ocean observation locations'
   call obsxloc(nro,obsoloc)
   call obstloc(nrt,obst0o,obsdto,obsotimes)
   call obspoints(trim(outdir)//'/obsloco',obsotimes,nrt,obsoloc,nro)

   print *,'Atmos observation locations'
   call obsxloc(nra,obsaloc)
   call obstloc(nrt,obst0a,obsdta,obsatimes)
   call obspoints(trim(outdir)//'/obsloca',obsatimes,nrt,obsaloc,nra)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The true solution smooth pseudo random field drawn from  N(0,1,rh).
   call pseudo1D(ref%ocean,nx,1,rh%ocean,dxsamp,nx)
   call pseudo1D(ref%atmos,nx,1,rh%atmos,dxsamp,nx)
   print *,'main: ref ok'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First guess solution stored in ave
   call pseudo1D(ave%ocean,nx,1,rh%ocean,dxsamp,nx)
   call pseudo1D(ave%atmos,nx,1,rh%atmos,dxsamp,nx)
   ave=(ave + ref)*(1.0/sqrt(2.0))
   print *,'main: fg ok'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialization of ensemble
   call pseudo1D(samples,nx,nrens,rh%ocean,dxsamp,nx)
   if (samp_fix) call fixsample1D(samples,nx,nrens)
   do j=1,nrens
      mem(j)%ocean(1:nx)=samples(1:nx,j)
   enddo

   call pseudo1D(samples,nx,nrens,rh%atmos,dxsamp,nx)
   if (samp_fix) call fixsample1D(samples,nx,nrens)
   do j=1,nrens
      mem(j)%atmos(1:nx)=samples(1:nx,j)
   enddo

   do j=1,nrens
      mem(j)%ocean=ave%ocean + sqrt(inivar%ocean)*mem(j)%ocean
      mem(j)%atmos=ave%atmos + sqrt(inivar%atmos)*mem(j)%atmos
   enddo
   print *,'main: ensemble ok'


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Loop over assimilation windows
   print '(a)','---------------------------------------------------------------------------'
   tfin=0
   do l=1,nrwindows                                   ! Integrate model over the assimilation window
      tini=tfin                                  ! time at beginning of assimilation window
      if (l==1) then
         tfin=50
      else
         tfin=tfin+nrw                                      ! time at end of assimilation window
      endif
      print *
      print '(a)','---------------------------------------------------------------------------'
      print '(a,i3.3,a,i5.5,a,i5.5)','DA window number ',l,' running from ',tini,' to ',tfin
      print *

! Reference solution over assimilation window stored in refout
      refout(tini)=ref
      do k=tini+1,tfin
         call model(ref)
         refout(k)=ref
      enddo

      print '(tr5,a,i3,a,f10.2,a,i3)','main: prior ensemble simulation:  window=',l
!$OMP PARALLEL DO
      do j=1,nrens
         win(tini,j)=mem(j)                                 ! store ensemble at tini in win
         do k=tini+1,tfin
            call model(mem(j))
            win(k,j)=mem(j)
         enddo
      enddo
!$OMP END PARALLEL DO

      win0(tini:tfin,:)=win(tini:tfin,:) ! Store prior ensemble over window for IES algorithm as win=win0*X

! Counting the number of measurements in the current DA window
      nrobs=obscount(nrt,tini,tfin,obsotimes,obsatimes,nro,nra)
      print '(tr5,a,i5)','main: Total number of measurements in the DA window= ',nrobs
      print *

      if ((nrobs > 0).and.(tini >= 50)) then
! Allocate DA variables
         allocate(obs(nrobs))
         allocate(DA(nrobs,nrens))
         allocate(E(nrobs,nrens))
         allocate(D(nrobs,nrens))
         allocate(Y(nrobs,nrens))
         allocate(S(nrobs,nrens))

! Generate the observations in obs and DA
         print '(tr5,a)','main: -> Calling prepD to get DA'
         call prepD(obs,DA,refout,nrobs,tini,tfin,obsoloc,obsaloc,obsotimes,obsatimes,iter)

! Loop over nmda steps
         steplength=steplength0
         W=0.0

         ! For localization compute correlation matrix between state and predicted measurements
         do iter=1,nmda

            if (oldana) then
            ! Old analysis update for ES ans ESMDA
               call assimilation_old(obs,win,refout,DA,tini,tfin,nrobs,iter,obsoloc,obsaloc,obsotimes,obsatimes)
            else
            ! New IES implmentation supporting ES, ESMDA and IES
               call assimilation(obs,win,win0,DA,E,D,Y,S,X,W,Wold,Yold,costf,&
                        tini,tfin,nrobs,iter,obsoloc,obsaloc,obsotimes,obsatimes,l,steplength)
               if (steplength < 0.01) exit
            endif

! timestep loop over assimilation window
            if ((iter < nmda).or.lsim) then
               write(*,'(tr5,a,i3,a,i3,i5,a,i5)',advance='no')'main: iter=',iter,' -> post-update-RESIM: window=',l,tini,'->',tfin
!$OMP PARALLEL DO
               do j=1,nrens
                  mem(j)=win(tini,j)
                  do k=tini+1,tfin
                     call model(mem(j))
                     win(k,j)=mem(j)
                  enddo
               enddo
!$OMP END PARALLEL DO
               print '(a,i3,a)','..... done'
               print *
            else
               write(*,'(tr5,a,i3,a,i3,i5,a,i5)',advance='no')'main: iter=',iter,' -> post-update-NOSIM: window=',l,tini,'->',tfin
            endif

! update prior ensemble in MDA
            if (cmethod(1:3)=="MDA") win0(tini:tfin,:)=win(tini:tfin,:)

         enddo ! end loop over mda steps
         deallocate(obs,Y,S,E,D,DA)
      endif

! Continue integration in next assimilation window from final value of current window
      mem(:)=win(tfin,:)

   enddo ! End loop over assimilation windows





! Compute ensemble mean and variance over current DA window
   call windowstat(win,nrt,nrens,mean,stdt)

! Compute root-mean-squared errors and std. dev.
   rmset(1)=0.0
   rmsst(1)=0.0
   do k=0,nrt
      rmse(k)=sqrt(average((mean(k)-refout(k))*(mean(k)-refout(k))))
      rmss(k)=sqrt(average(stdt(k)*stdt(k)))
      if (k > 200) rmset(1)=rmset(1)+rmse(k)*rmse(k)
      if (k > 200) rmsst(1)=rmsst(1)+rmss(k)*rmss(k)
   enddo
   rmset(1)=rmset(1)/real(nrt-200)
   rmsst(1)=rmsst(1)/real(nrt-200)
   rmset(1)=sqrt(rmset(1))
   rmsst(1)=sqrt(rmsst(1))

   open(10,file=trim(outdir)//'/rmse.dat')
      do k=0,nrt
         write(10,'(i5,1024g12.4)')k,rmse(k)%Atmos,rmse(k)%Ocean,rmss(k)%Atmos,rmss(k)%Ocean
      enddo
   close(10)
   open(10,file=trim(outdir)//'/rmset.dat')
      write(10,'(10g12.4)')rmset%Atmos,rmset%Ocean,rmsst%Atmos,rmsst%Ocean
   close(10)

! Print the cost functions values for each window.
   if (cmethod == 'IES') then
      open(10,file=trim(outdir)//'/costf.dat')
         do iter=1,nmda
            write(10,'(i5,1000g12.4)')iter,costf(iter,1:nrwindows)
         enddo
      close(10)
   endif

! Dumping reference solution, mean and standard deviations for plotting
   call gnuplot('gnu_ref',refout,nrt,outdir)
   call gnuplot('gnu_ave',mean,nrt,outdir)
   call gnuplot('gnu_std',stdt,nrt,outdir)
   do k=0,nrt
      rest(k)=mean(k)-refout(k)
      rest(k)=abs(rest(k))
   enddo
   call gnuplot('gnu_res',rest,nrt,outdir)
   if (lglobstat) then
      call covstat(win,nrt,nrens,mean,stdt,covo,cova,outdir)
      call gnuplot('gnu_covo',covo,nrt,outdir)
      call gnuplot('gnu_cova',cova,nrt,outdir)
   endif

end program main

