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

   type(state), allocatable :: full(:,:)   ! The ensemble of realizations over the whole simulation
   type(state), allocatable :: win(:,:)    ! The ensemble of realizations over an assimilation window.
   type(state), allocatable :: win0(:,:)   ! Initial ensemble of realizations over an assimilation window.
   type(state), allocatable :: winref(:)   ! The reference solution over an assimilation window.
   type(state), allocatable :: mem(:)      ! The ensemble of realizations
   type(state), allocatable :: sysnoise(:) ! The ensemble of sampled system noise updated every timestep
   type(state) ref                         ! reference solution
   type(state) ave                         ! ensemble average

   real, allocatable :: samples(:,:)       ! work array used when sampling in pseudo1D

! Spacew time statistics diagnostic variables
   type(state), allocatable :: refout(:)   ! ensemble average as a function of space and time
   type(state), allocatable :: mean(:)     ! ensemble average as a function of space and time
   type(state), allocatable :: stdt(:)     ! ensemble std dev as a function of space and time
   type(state), allocatable :: covo(:)     ! ensemble std dev as a function of space and time
   type(state), allocatable :: cova(:)     ! ensemble std dev as a function of space and time
   type(substate), allocatable :: rmse(:)  ! time series of rmse values (mean - referece)
   type(substate), allocatable :: rmss(:)  ! time series of rms std values

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
   real, allocatable :: D(:,:)             ! Innovations
   real, allocatable :: DA(:,:)            ! Innovations
   real, allocatable :: R(:,:)             ! Measurement error covariance matrix
   real, allocatable :: meanS(:)           ! Mean of predicted measurements
   real, allocatable :: innov(:)           ! Mean innovation
   real, allocatable :: W(:,:)             ! ies iteration matrix
   real, allocatable :: Wold(:,:)          ! ies iteration matri for reducing steplength
   real, allocatable :: Yold(:,:)          ! ies predicted measurements for reducing steplength
   real, allocatable :: X(:,:)            ! The X matrix :-)


! other variables
   integer tini                            ! start time of an assimilation window
   integer tfin                            ! end   time of an assimilation window
   integer nrobs                           ! Total number of measurement per assimilation window

   integer j,k,l,iter,ldw,iprt
   real time
   real fac
   real :: dxsamp=1.0
   real :: steplength
   real, allocatable :: costf(:,:)

   call readinfile()

   call set_random_seed2

   call system('rm -f eigenvalues.dat')
   call system('rm -f obsloc?.dat')
!   call system('touch obsloca.dat obsloco.dat')

! Now that we know all dimensions, allocate the main arrays
   if (lglobstat) allocate (full(0:nrt,nrens))
   allocate (win(0:nrw,nrens))
   allocate (win0(0:nrw,nrens))
   allocate (winref(0:nrw))
   allocate (refout(0:nrt))
   allocate (mean(0:nrt))
   allocate (stdt(0:nrt))
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
   do l=1,nrwindows                                   ! Integrate model over the assimilation window
      tini=(l-1)*nrw                                  ! time at beginning of assimilation window
      tfin=l*nrw                                      ! time at end of assimilation window
      print *
      print '(a)','---------------------------------------------------------------------------'
      print '(a,i3.3,a,i5.5,a,i5.5)','DA window number ',l,' running from ',tini,' to ',tfin
      print *

! Reference solution over assimilation window stored in refout
      winref(0)=ref
      do k=1,nrw
         call model(ref)
         winref(k)=ref
      enddo
      refout(tini:tfin)=winref

! prior ensemble integration over assimilation window
      win(0,:)=mem(:)                                 ! store ensemble at tini in win
      print '(tr5,a,i3,a,f10.2,a,i3)','main: prior ensemble simulation:  window=',l
!$OMP PARALLEL DO
      do j=1,nrens
         do k=1,nrw
            call model(mem(j))
            win(k,j)=mem(j)
         enddo
      enddo
!$OMP END PARALLEL DO
      win0=win ! Store prior ensemble over window for IES algorithm as win=win0*X

! Counting the number of measurements in the current DA window
      nrobs=obscount(nrt,tini,tfin,obsotimes,obsatimes,nro,nra)
      print '(tr5,a,i5)','main: Total number of measurements in the DA window= ',nrobs
      print *

      if (nrobs > 0) then
! Allocate DA variables
         allocate(obs(nrobs))
         allocate(E(nrobs,nrens))
         allocate(D(nrobs,nrens))
         allocate(DA(nrobs,nrens))
         allocate(Y(nrobs,nrens))

! oldana vars
         allocate(S(nrobs,nrens))
         allocate(meanS(nrobs))
         allocate(innov(nrobs))
         allocate(R(nrobs,nrobs))

! Generate the observations in obs and DA
         winref=refout(tini:tfin)
         print '(tr5,a)','main: -> Calling prepD to get DA'
         call prepD(obs,DA,winref,nrobs,l,tini,tfin,obsoloc,obsaloc,obsotimes,obsatimes,iter)
         print *

! Loop over nmda steps
         steplength=steplength0
         W=0.0
         fac=1.0
         do iter=1,nmda
            if (oldana) then
               print '(tr5,a,i3,a)','main: iter=',iter,' -> Calling enkf preprep'
               ! Returns all matrices UNSCALED by sqrt(nrens-1) for use in analysis.F90
               call enkfprep(mem,obs,Y,S,E,D,DA,meanS,R,innov,winref,win,nrobs,l,tini,tfin,obsoloc,obsaloc,obsotimes,obsatimes)
               D=D-Y
               print '(tr5,a,2(i3,a))','main: iter=',iter,' -> Calling analysis with mode: ',mode_analysis
               call analysis(win, R, E, S, D, innov, ndim*(nrw+1), nrens, nrobs, .true., truncation, mode_analysis, &
                           lrandrot, lupdate_randrot, lsymsqrt, inflate, infmult, 1)
            else

               if ((cmethod(1:3)=='MDA') .or. ((cmethod(1:3)=='IES').and.(iter==1))) then
!                 Using ESMDA: We simullate a new E for every DA step and then scale it with sqrt(nmda),
!                              thus, D=DA+sqrt(nmda)*E is updated every MDA step.
!                 Using IES  : We simulate E only in the first iteration and compute D=DA+E which we use unchanged
!                              in all following IES iterations.  We are passing the predicted measurement Y and the
!                              perturbed measurements D to IES
                  call pertE(E,nrobs,l   ,tini,tfin,obsoloc,obsaloc,obsotimes,obsatimes)
                  if (cmethod(1:3) == 'MDA') then
                     E=sqrt(real(nmda))*E
                     W=0.0
                     steplength=1.0
                  endif
                  D=DA+E
                  call scaling(D,E,nrobs,nrens)
               endif

               if ((cmethod(1:3)=='IES').and.LM) then
               ! Levenberg–Marquardt algorithm with IES
                  fac=10.0+100.0*exp(-10.0*real(iter-1)/(real(nmda)))
               endif

               call prepY(Y,win,nrobs,l,tini,tfin,obsoloc,obsaloc,obsotimes,obsatimes)
               call scaling(Y,E,nrobs,nrens)

               if (cmethod(1:3) == 'IES') then
                  call ies_steplength(steplength,costf,nrwindows,nmda,W,Wold,D,Y,Yold,nrens,nrobs,iter,l)
                  if (steplength < 0.01) exit
               endif

               print '(tr5,a,2(i3,a))','main: iter=',iter,' -> Calling ies with mode: ',mode_analysis
               call ies(Y,D,W,nrens,nrobs,steplength,mode_analysis,fac)

               if (cmethod(1:3) == 'IES') iprt=print_ies_status(fac,steplength,costf(iter,l),Wold,W,nrens,iter)

! X = I + W/sqrt(N-1)
               X=W/sqrt(real(nrens-1))
               do j=1,nrens
                  X(j,j)=X(j,j)+1.0
               enddo

! Window update: for last iteration and in ES we are updating the whole window not just the intial condition.
! For linear dynamics, this doesn't change anything, but for strongly unstable dynamics it may give a better
! posterior estimate over the window, and better starting point for the next window.
               write(*,'(tr5,a,i3,a)',advance='no')'main: iter=',iter,' -> Ensemble update'
               if ((iter==nmda).and.(.not.lsim)) then
                  ldw=ndim*(nrw+1)
                  call dgemm('N','N',ldw,nrens,nrens,1.0,win0,ldw,X,nrens,0.0,win,ldw)
               else
                  ldw=ndim
                  call dgemm('N','N',ldw,nrens,nrens,1.0,win0(0,:),ldw,X,nrens,0.0,win(0,:),ldw)
               endif
               print '(a,i3,a)','..... done'
            endif

! timestep loop over assimilation window
            if ((iter < nmda).or.lsim) then
               write(*,'(tr5,a,i3,a,i3,i5,a,i5)',advance='no')'main: iter=',iter,' -> post-update-RESIM: window=',l,tini,'->',tfin
               mem(:)=win(0,:)
!$OMP PARALLEL DO
               do j=1,nrens
                  do k=1,nrw
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
            if (cmethod(1:3)=="MDA") win0=win

         enddo ! end loop over mda steps
         deallocate(obs,Y,S,E,D,DA,meanS,R,innov)
      endif

! Continue integration in next assimilation window from final value of current window
      mem(:)=win(nrw,:)

! Only needed to compute space-time covariances
      if (lglobstat) full(tini:tfin,:)=win(0:nrw,:)

! Compute ensemble mean and variance over current DA window
      call windowstat(win,nrw,nrens,mean(tini:tfin),stdt(tini:tfin))
   enddo ! End loop over assimilation windows

! Compute root-mean-squared errors and std. dev.
   do k=0,nrt
      rmse(k)=sqrt(average((mean(k)-refout(k))*(mean(k)-refout(k))))
      rmss(k)=sqrt(average(stdt(k)*stdt(k)))
   enddo
   open(10,file=trim(outdir)//'/rmse.dat')
      do k=0,nrt
         write(10,'(i5,1024g12.4)')k,rmse(k)%Atmos,rmse(k)%Ocean,rmss(k)%Atmos,rmss(k)%Ocean
      enddo
   close(10)

! Print the cost functions values for each window.
   if (cmethod == 'IES') then
      open(10,file=trim(outdir)//'/costf.dat')
         do iter=1,nmda
            write(10,'(i5,100g12.4)')iter,costf(iter,1:nrwindows)
         enddo
      close(10)
   endif

! Dumping reference solution, mean and standard deviations for plotting
   call gnuplot('gnu_ref',refout,nrt,outdir)
   call gnuplot('gnu_ave',mean,nrt,outdir)
   call gnuplot('gnu_std',stdt,nrt,outdir)
   if (lglobstat) then
      print *,'calling full covariance statistics'
      call covstat(full,nrt,nrens,mean,stdt,covo,cova,outdir)
      call gnuplot('gnu_covo',covo,nrt,outdir)
      call gnuplot('gnu_cova',cova,nrt,outdir)
   endif

end program main

